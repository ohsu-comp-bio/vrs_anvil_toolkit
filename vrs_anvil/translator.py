import logging
import queue
import threading
from dataclasses import dataclass, field
from typing import NamedTuple, Generator, Any, Optional

from pydantic import BaseModel

from vrs_anvil import caching_allele_translator_factory

_logger = logging.getLogger("vrs_anvil.translator")


class WorkerThread(threading.Thread):
    """Read from the task queue, process the item, with local translator and write the result to the result queue."""

    def __init__(self, task_queue, result_queue, normalize):
        super().__init__(daemon=True)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.translator = caching_allele_translator_factory(normalize=normalize)
        self.busy = False

    def run(self):
        while True:
            try:
                prioritized_item = self.task_queue.get()
                item: VCFItem = prioritized_item.item
                if item is None:
                    break  # Signal to exit the thread

                self.busy = True
                allele = self.translator.translate_from(fmt=item.fmt, var=item.var)
                _ = item._asdict()
                _["result"] = allele
                result = VCFItem(**_)
                self.result_queue.put(PrioritizedItem(1, result))

                self.busy = False
                self.task_queue.task_done()

            except Exception as exc:
                _logger.exception(f"{self.name} error {exc}")
                self.busy = False
                break


class VCFItem(NamedTuple):
    """A named tuple to hold the VCF item."""

    fmt: str
    """The format of the item's id e.g gnomad - passed to vrs python"""
    var: str
    """variant - passed to vrs python"""
    file_name: str = None
    """data source"""
    line_number: int = None
    """line number in the file"""
    identifier: str = None
    """identifier for the item"""
    result: Any = None
    """identifier for the item"""


@dataclass(order=True)
class PrioritizedItem:
    """If the data elements are not comparable, the data can be wrapped in a class that ignores the data item and only compares the priority number"""

    priority: int
    item: Any = field(compare=False)


class Translator(BaseModel):
    """A class to run the translation in either threaded or non-threaded fashion."""

    normalize: Optional[bool] = False

    def translate_from(
        self, generator: Generator[VCFItem, None, None], num_threads: int = 8
    ) -> Generator[VCFItem, None, None]:
        if num_threads > 1:
            return threaded_translator(generator, num_threads, self.normalize)
        else:
            return inline_translator(generator, self.normalize)


def inline_translator(
    generator: Generator[VCFItem, None, None], normalize: bool = False
) -> Generator[VCFItem, None, None]:
    """A generator that runs the translation in a non-threaded fashion."""
    tlr = caching_allele_translator_factory(normalize=normalize)
    for item in generator:
        allele = tlr.translate_from(fmt=item.fmt, var=item.var)
        _ = item._asdict()
        _["result"] = allele
        yield VCFItem(**_)


def threaded_translator(
    generator: Generator[VCFItem, None, None],
    num_worker_threads: int,
    normalize: bool = False,
) -> Generator[VCFItem, None, None]:
    """A generator that runs the translation in a threaded fashion."""
    task_queue = queue.PriorityQueue(maxsize=num_worker_threads * 2)
    result_queue = queue.PriorityQueue(maxsize=num_worker_threads * 2)

    # Start worker threads
    worker_threads = [
        WorkerThread(task_queue, result_queue, normalize)
        for _ in range(num_worker_threads)
    ]
    for worker_thread in worker_threads:
        worker_thread.start()

    # Start the reader thread

    def reader_thread(_generator):
        for item in _generator:
            task_queue.put(PrioritizedItem(1, item))  # Default priority is 1

    reader = threading.Thread(target=reader_thread, args=(generator,), daemon=True)
    reader.start()

    # Main thread yields results
    logged_already = []
    c = 0
    while True:
        try:
            prioritized_item = result_queue.get(
                timeout=1
            )  # Give the thread time to start
            c += 1
            yield prioritized_item.item
            result_queue.task_done()
        except queue.Empty:
            reader_started = c > 0
            if reader_started and all(
                not worker_thread.busy for worker_thread in worker_threads
            ):
                _logger.info("All worker threads are idle, exiting.")
                break
            if not reader_started:
                msg = "queue.Empty - Waiting for reader to start"
            else:
                msg = "queue.Empty - Waiting for worker threads to finish"

            if msg not in logged_already:
                _logger.info(msg)
                logged_already.append(msg)
