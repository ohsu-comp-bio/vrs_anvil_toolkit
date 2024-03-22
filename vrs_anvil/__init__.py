import json
import logging
import os
import pathlib
import traceback
from typing import Optional, Generator
import zipfile

from biocommons.seqrepo import SeqRepo
from diskcache import Cache
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator
from glom import glom
from pydantic import BaseModel, model_validator
import requests

_logger = logging.getLogger("vrs_anvil")
LOGGED_ALREADY = set()

manifest: "Manifest" = None

# TODO - read from manifest
gigabytes = 20
bytes_in_a_gigabyte = 1024**3  # 1 gigabyte = 1024^3 bytes
cache_size_limit = gigabytes * bytes_in_a_gigabyte


def seqrepo_dir():
    """Return the seqrepo directory."""
    with open(".env") as f:
        for line in f:
            # Ignore comments and empty lines
            if line.strip() and not line.strip().startswith("#"):
                key, value = line.strip().split("=", 1)
                if key == "SEQREPO_ROOT":
                    return value + "/latest"


def cache_directory(cache_name: str) -> str:
    """Return the cache directory."""
    return str(pathlib.Path(manifest.cache_directory) / cache_name)


class CachingAlleleTranslator(AlleleTranslator):
    """A subclass of AlleleTranslator that uses cache results and adds a method to run in a threaded fashion."""

    _cache: Cache = None

    def __init__(self, data_proxy: SeqRepoDataProxy, normalize: bool = False):
        super().__init__(data_proxy)
        self.normalize = normalize
        self._cache = None
        if manifest.cache_enabled:
            self._cache = Cache(
                directory=cache_directory("allele_translator"),
                size_limit=cache_size_limit,
            )
        else:
            _logger.info("Cache is not enabled")

    def translate_from(self, var, fmt=None, **kwargs):
        """Check and update cache"""

        if self._cache is not None:
            key = f"{var}-{fmt}"
            if key in self._cache:
                return self._cache[key]

        val = super().translate_from(var, fmt=fmt, **kwargs)

        if self._cache is not None:
            self._cache[key] = val

        return val


def caching_allele_translator_factory(
    normalize: bool = False, seqrepo_directory: str = None
):
    """Return a CachingAlleleTranslator instance with local seqrepo"""
    if not seqrepo_directory:
        seqrepo_directory = manifest.seqrepo_directory
    dp = SeqRepoDataProxy(SeqRepo(seqrepo_directory))
    assert dp is not None, "SeqRepoDataProxy is None"
    translator = CachingAlleleTranslator(dp)
    translator.normalize = normalize
    return translator


class ThreadedTranslator(BaseModel):
    """A class to run the translation in a threaded fashion."""

    _thread_resources: dict = {}
    normalize: Optional[bool] = False

    def threaded_translate_from(
        self, generator, num_threads=8
    ) -> Generator[dict, None, None]:
        """
        Process a generator using multithreading and yield results as they occur.

        Args:
        - generator: The generator to be processed.
        - func: The function to be applied to each element of the generator.
        - num_threads: Number of threads to use for processing.

        Yields:
        - Results from applying the function to each element of the generator.
        """

        translator = caching_allele_translator_factory(normalize=self.normalize)

        def _ensure_tuple(item):
            """Ensure that the item is a tuple (item, file, line)."""
            if not isinstance(item, tuple):
                return item, None, None
            return item

        def _process_item(item):
            """Process an item from the generator.
            An item is a tuple:
            - 0: the dictionary with the fmt and var keys
            - 1: the file path of the vcf (or other source) (optional)
            - 2: the line number in the vcf (optional)
            """
            nonlocal translator
            try:
                result_dict = {"parameters": item[0], "file": None, "line": None}
                if len(item) > 1:
                    result_dict["file"] = item[1]
                if len(item) > 2:
                    result_dict["line"] = item[2]
                result = translator.translate_from(**item[0])
                # result = {'location': {}, 'state': {}, 'type': {}}  # TESTING dummy results:
                result_dict["result"] = result
            except Exception as e:
                result_dict["error"] = str(e)
                result_dict["stack_trace"] = traceback.format_exc()
                _logger.warning(result_dict)

            return result_dict

        # main ...
        c = 0

        for _ in generator:
            yield _process_item(_ensure_tuple(_))
            c += 1

        del translator

        _logger.info(f"threaded_translate_from: Finished processing {c} items.")


def generate_gnomad_ids(vcf_line, compute_for_ref: bool = True) -> list[str]:
    """Assuming a standard VCF format with tab-separated fields, generate a gnomAD-like ID from a VCF line.
    see https://github.com/ga4gh/vrs-python/blob/main/src/ga4gh/vrs/extras/vcf_annotation.py#L386-L411
    """
    fields = vcf_line.strip().split("\t")
    gnomad_ids = []
    # Extract relevant information (you may need to adjust these indices based on your VCF format)
    chromosome = fields[0]
    position = fields[1]
    reference_allele = fields[3]
    alternate_allele = fields[4]

    gnomad_loc = f"{chromosome}-{position}"
    if compute_for_ref:
        gnomad_ids.append(f"{gnomad_loc}-{reference_allele}-{reference_allele}")
    for alt in alternate_allele.split(","):
        alt = alt.strip()
        # TODO - Should we be raising a ValueError hear and let the caller do the logging?
        # TODO - Should this be a config in the manifest?
        # ['<INS>', '<DEL>', '<DUP>', '<INV>', '<CNV>', '<DUP:TANDEM>', '<DUP:INT>', '<DUP:EXT>', '*']
        invalid_alts = ["INS", "DEL", "DUP", "INV", "CNV", "TANDEM", "INT", "EXT", "*"]
        is_valid = True
        for invalid_alt in invalid_alts:
            if invalid_alt in alt:
                is_valid = False
                _ = f"Invalid allele found: {alt}"
                if _ not in LOGGED_ALREADY:
                    LOGGED_ALREADY.add(_)
                    _logger.error(_)
                break
        if is_valid:
            gnomad_ids.append(f"{gnomad_loc}-{reference_allele}-{alt}")

    return gnomad_ids


def params_from_vcf(path, limit=None) -> Generator[dict, None, None]:
    """Open the vcf file, skip headers, yield the first lines as gnomad-like IDs"""
    c = 0
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            gnomad_ids = generate_gnomad_ids(line)
            for gnomad_id in gnomad_ids:
                yield {"fmt": "gnomad", "var": gnomad_id}, path, c
            c += 1
            if limit and c > limit:
                break


def find_items_with_key(dictionary, key_to_find):
    """Find all items in a dictionary that have a specific key."""
    result = glom(dictionary, f"**.{key_to_find}")
    return result


class MetaKBProxy(BaseModel):
    """A proxy for the MetaKB, maintains a cache of VRS ids."""

    metakb_path: pathlib.Path
    _cache: Optional[Cache] = None

    def __init__(self, metakb_path: pathlib.Path, cache: Cache = None):
        super().__init__(metakb_path=metakb_path, _cache=cache)
        if cache is None:
            reload_cache = False
            if not (metakb_path / "cache").is_dir():
                reload_cache = True
            cache = Cache(directory=cache_directory("metakb"))
            cache.stats(enable=True)
            if reload_cache:
                for _ in metakb_ids(metakb_path):
                    cache.set(_, True)
        self._cache = cache

    def get(self, vrs_id: str) -> bool:
        """Get the vrs_id from the cache."""
        return self._cache.get(vrs_id, False)


def metakb_ids(metakb_path) -> Generator[str, None, None]:
    """Find all the applicable vrs ids in the metakb files."""
    if len(list(pathlib.Path(metakb_path).glob("*.json"))) == 0:
        _get_metakb_models(metakb_path)

    for file_name in pathlib.Path(metakb_path).glob("*.json"):
        if file_name.is_file():
            with open(file_name, "r") as file:
                data = json.loads(file.read())
                yield from (
                    [
                        _
                        for _ in find_items_with_key(data, "id")
                        if _.startswith("ga4gh:VA")
                    ]
                )


def _get_metakb_models(metakb_path):
    def _download_s3(url: str, outfile_path: pathlib.Path) -> None:
        """Download objects from public s3 bucket

        :param url: URL for metakb file in s3 bucket
        :param outfile_path: Path where file should be saved
        """
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(outfile_path, "wb") as h:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        h.write(chunk)

    pathlib.Path(metakb_path).mkdir(exist_ok=True)

    date = "20240305"
    json_files = [f"civic_cdm_{date}.json", f"moa_cdm_{date}.json"]

    for json_file in json_files:
        json_path = f"{metakb_path}/{json_file}"

        url = (
            "https://vicc-metakb.s3.us-east-2.amazonaws.com"
            + f"/cdm/{date}/{json_file}.zip"
        )
        zip_path = f"{json_path}.zip"

        _download_s3(url, zip_path)
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(metakb_path)
        os.remove(zip_path)


class Manifest(BaseModel):
    """
    A class to represent the manifest file.
    Defaults to the values below if not provided in the manifest.yaml
    """

    cache_directory: str = "cache/"
    """Path to the cache directory, defaults to cache/ (relative to the root of the repository)"""

    num_threads: int = 20
    """Number of works to use for processing, defaults to 20"""

    annotate_vcfs: bool = False
    """Should we create new VCFs with annotations. FOR FUTURE USE"""

    state_directory: str = "state/"
    """where to store the state of the application, log files, etc."""

    vcf_files: list[str]
    """The local file paths or URLs to vcf files to be processed"""
    # TODO - 2x check why local files need to be absolute paths

    work_directory: str = "work/"
    """The directory to store intermediate files"""

    seqrepo_directory: str = "~/seqrepo/latest"
    """The directory where seqrepo is located"""

    normalize: bool = False
    """Normalize the VRS ids"""

    limit: Optional[int] = None
    """Stop processing after this many lines"""

    cache_enabled: Optional[bool] = True
    """Cache results"""

    compute_for_ref: Optional[bool] = False
    """Compute reference allele"""

    estimated_vcf_lines: Optional[int] = 4000000
    """How many lines per vcf file?  Used for progress bar"""

    metakb_directory: str = "metakb/"
    """Where the CDM files are located.  This is a directory containing json files"""

    @model_validator(mode="after")
    def check_paths(self) -> "Manifest":
        """Post init method to set the cache directory."""
        self.seqrepo_directory = str(pathlib.Path(self.seqrepo_directory).expanduser())
        self.work_directory = str(pathlib.Path(self.work_directory).expanduser())
        self.cache_directory = str(pathlib.Path(self.cache_directory).expanduser())
        self.state_directory = str(pathlib.Path(self.state_directory).expanduser())
        self.metakb_directory = str(pathlib.Path(self.metakb_directory).expanduser())

        for _ in ["seqrepo_directory", "metakb_directory"]:
            if not pathlib.Path(getattr(self, _)).exists():
                raise ValueError(f"{_} does not exist")

        for _ in ["work_directory", "cache_directory", "state_directory"]:
            if not pathlib.Path(getattr(self, _)).exists():
                pathlib.Path(getattr(self, _)).mkdir(parents=True, exist_ok=True)
                _logger.debug(f"Created directory {getattr(self, _)}")

        return self
