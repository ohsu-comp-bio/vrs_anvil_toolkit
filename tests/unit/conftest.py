from concurrent.futures.thread import ThreadPoolExecutor
from concurrent.futures import as_completed
import queue
import traceback
from functools import lru_cache
from threading import current_thread

import pytest
from biocommons.seqrepo import SeqRepo
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import AlleleTranslator


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


def caching_allele_translator_factory():
    """Return a CachingAlleleTranslator instance with local seqrepo"""
    dp = SeqRepoDataProxy(SeqRepo(_seqrepo_dir()))
    assert dp is not None
    translator = CachingAlleleTranslator(dp)
    return translator


class ThreadedTranslator:
    """A class to run the translation in a threaded fashion."""
    thread_resources = {}

    def thread_initializer(self):
        """Initialize resources for use in threads."""
        thread = current_thread()
        translator = caching_allele_translator_factory()
        self.thread_resources[thread.name] = {'translator': translator}

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
                translator = self.thread_resources[current_thread().name]['translator']
                result = translator.translate_from(**item)
                # result = {'location': {}, 'state': {}, 'type': {}}  # TESTING dummy results:
            except Exception as e:
                stack_trace = traceback.format_exc()
                result = {'error': str(e), 'item': item, 'stack_trace': stack_trace}
            results_queue.put(result)

        with ThreadPoolExecutor(max_workers=num_threads, initializer=self.thread_initializer) as executor:
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
    return ThreadedTranslator()
