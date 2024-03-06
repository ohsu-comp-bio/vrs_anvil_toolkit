import json
import pathlib
from concurrent.futures.thread import ThreadPoolExecutor
from concurrent.futures import as_completed
import queue
import traceback
from functools import lru_cache
from threading import current_thread
from typing import Optional, Generator

import pytest
from biocommons.seqrepo import SeqRepo
from diskcache import Cache
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator
# from glom import glom
from pydantic import BaseModel


@pytest.fixture
def seqrepo_dir():
    """Return the seqrepo directory as fixture."""
    return _seqrepo_dir()


def _seqrepo_dir():
    """Return the seqrepo directory."""
    with open(".env") as f:
        for line in f:
            # Ignore comments and empty lines
            if line.strip() and not line.strip().startswith('#'):
                key, value = line.strip().split('=', 1)
                if key == "SEQREPO_ROOT":
                    return value + "/latest"


class CachingAlleleTranslator(AlleleTranslator):
    """A subclass of AlleleTranslator that uses lru_cache to cache results and adds a method to run in a threaded fashion."""

    @lru_cache(maxsize=None)  # TODO - this is a naive implementation, we should use a better cache, e.g. spill to disk, expire etc.
    def translate_from(self, var, fmt=None, **kwargs):
        return super().translate_from(var, fmt=fmt, **kwargs)


def caching_allele_translator_factory(normalize: bool = False, seqrepo_dir: str = _seqrepo_dir()):
    """Return a CachingAlleleTranslator instance with local seqrepo"""
    dp = SeqRepoDataProxy(SeqRepo(seqrepo_dir))
    assert dp is not None, "SeqRepoDataProxy is None"
    translator = CachingAlleleTranslator(dp)
    translator.normalize = normalize
    return translator


class ThreadedTranslator(BaseModel):
    """A class to run the translation in a threaded fashion."""
    _thread_resources: dict = {}
    normalize: Optional[bool] = False

    def _thread_initializer(self):
        """Initialize resources for use in threads."""
        thread = current_thread()
        translator = caching_allele_translator_factory(normalize=self.normalize)
        self._thread_resources[thread.name] = {'translator': translator}

    def threaded_translate_from(self, generator, num_threads=8):
        """
        Process a generator using multithreading and yield results as they occur.

        Args:
        - generator: The generator to be processed.
        - func: The function to be applied to each element of the generator.
        - num_threads: Number of threads to use for processing.

        Yields:
        - Results from applying the function to each element of the generator.
        """
        results_queue = queue.Queue()

        def process_item(item):
            try:
                translator = self._thread_resources[current_thread().name]['translator']
                result = translator.translate_from(**item)
                # result = {'location': {}, 'state': {}, 'type': {}}  # TESTING dummy results:
            except Exception as e:
                stack_trace = traceback.format_exc()
                result = {'error': str(e), 'item': item, 'stack_trace': stack_trace}
            results_queue.put(result)

        with ThreadPoolExecutor(max_workers=num_threads, initializer=self._thread_initializer) as executor:
            # Submit each item in the generator to the executor
            futures = {executor.submit(process_item, item): item for item in generator}

            # Yield results as they become available
            for _ in as_completed(futures):
                yield results_queue.get()


@pytest.fixture
def my_translator():
    """Return a single translator instance."""
    return caching_allele_translator_factory()


@pytest.fixture
def threaded_translator():
    """Return a "thread manager", a pool of translator instances."""
    return ThreadedTranslator(normalize=False)


def generate_gnomad_id(vcf_line) -> list[str]:
    """Assuming a standard VCF format with tab-separated fields, generate a gnomAD-like ID from a VCF line."""
    # TODO - Change the way we generate the gnomad_id to match the way it is done in the vcf_annotation.py see https://github.com/ga4gh/vrs-python/blob/main/src/ga4gh/vrs/extras/vcf_annotation.py#L386-L411
    fields = vcf_line.strip().split('\t')

    # Extract relevant information (you may need to adjust these indices based on your VCF format)
    chromosome = fields[0]
    position = fields[1]
    reference_allele = fields[3]
    alternate_allele = fields[4]

    # Create a gnomAD-like ID
    gnomad_id = f"{chromosome}-{position}-{reference_allele}-{alternate_allele}"

    return gnomad_id


@pytest.fixture()
def num_threads():
    """Return the number of threads to use for testing."""
    return 20


def params_from_vcf(path, limit=None) -> Generator[dict, None, None]:
    """Open the vcf file, skip headers, yield the first lines as gnomad-like IDs"""
    c = 0
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            gnomad_id = generate_gnomad_id(line)
            yield {"fmt": "gnomad", "var": gnomad_id}
            c += 1
            if limit and c > limit:
                break


@pytest.fixture
def thousand_genome_vcf() -> pathlib.Path:
    """Return a path to a 1000g vcf."""
    _ = pathlib.Path("tests/data/1kGP.chr1.1000.vcf")
    assert _.exists()
    return _


def find_items_with_key(dictionary, key_to_find):
    """Find all items in a dictionary that have a specific key."""
    result = glom(dictionary, f"**.{key_to_find}")
    return result


def meta_kb_ids(metakb_path) -> Generator[str, None, None]:
    """Find all the applicable vrs ids in the metakb files."""
    for file_name in pathlib.Path(metakb_path).glob("*.json"):
        if file_name.is_file():
            with open(file_name, 'r') as file:
                data = json.loads(file.read())
                yield from ([_ for _ in find_items_with_key(data, 'id') if _.startswith('ga4gh:VA')])


class MetaKBProxy(BaseModel):
    """A proxy for the MetaKB, maintains a cache of VRS ids."""
    metakb_path: pathlib.Path
    _cache: Optional[Cache] = None

    def __init__(self, metakb_path: pathlib.Path, cache: Cache = None):
        super().__init__(metakb_path=metakb_path, _cache=cache)
        if cache is None:
            reload_cache = False
            if not (metakb_path / 'cache').is_dir():
                reload_cache = True
            cache = Cache(directory=metakb_path / 'cache')
            cache.stats(enable=True)
            if reload_cache:
                for _ in meta_kb_ids(metakb_path):
                    cache.set(_, True)
            print(cache.stats())
        self._cache = cache

    def get(self, vrs_id: str) -> bool:
        """Get the vrs_id from the cache."""
        return self._cache.get(vrs_id, False)

