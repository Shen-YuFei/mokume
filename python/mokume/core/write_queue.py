"""
Asynchronous file writing utilities for the mokume package.

This module provides thread-based tasks for writing DataFrames to CSV and
Parquet files asynchronously, using a queue-based approach.
"""

import os
import time
from threading import Thread
from queue import Queue, Empty
from typing import Any

import pandas as pd
import pyarrow as pa
from pyarrow import parquet as pq

from .logger import get_logger

# Get a logger for this module
logger = get_logger("mokume.write_queue")


class WriteCSVTask(Thread):
    """
    A thread-based task for writing pandas DataFrames to a CSV file.

    This class extends the Thread class to asynchronously write DataFrames
    to a specified CSV file path. It manages a queue to handle incoming
    DataFrames and writes them to the file in the order they are received.
    The CSV file is created with an optional header, and additional write
    options can be specified.

    Attributes
    ----------
    path : str
        The file path where the CSV will be written.
    write_options : dict[str, Any]
        Options for writing the CSV file.
    """

    path: str
    write_options: dict[str, Any]

    _queue: Queue
    _wrote_header: bool

    def __init__(
        self, path: str, daemon: bool = True, write_options: dict = None, **kwargs
    ):
        """
        Initialize a WriteCSVTask instance.

        Parameters
        ----------
        path : str
            The file path where the CSV will be written. The extension
            will be automatically set to '.csv'.
        daemon : bool, optional
            Whether the thread should be a daemon thread. Defaults to True.
        write_options : dict, optional
            Additional options for writing the CSV file. Defaults to None.
        **kwargs
            Additional keyword arguments to be merged with write_options.
        """
        super().__init__(daemon=daemon)
        if write_options is None:
            write_options = {}

        path, _ext = os.path.splitext(path)
        path += ".csv"

        self.path = path
        self.write_options = write_options | kwargs
        self._wrote_header = False
        self._queue = Queue()

    def write(self, table: pd.DataFrame):
        """
        Add a DataFrame to the queue for writing to the CSV file.

        Parameters
        ----------
        table : pd.DataFrame
            The DataFrame to be added to the queue.
        """
        logger.debug("Queuing %d rows for CSV writing to %s", len(table), self.path)
        self._queue.put(table)

    def close(self):
        """Signal the thread to finish processing and close the file."""
        logger.debug("Closing CSV writer queue for %s", self.path)
        self._queue.put(None)
        self.join()

    def _write(self, table: pd.DataFrame):
        """
        Write a DataFrame to the CSV file specified by the path attribute.

        Parameters
        ----------
        table : pd.DataFrame
            The DataFrame to be written to the CSV file.
        """
        start_time = time.time()
        rows = len(table)

        try:
            table.to_csv(
                self.path,
                header=not self._wrote_header,
                mode="a+" if self._wrote_header else "w",
                index=False,
                **self.write_options,
            )
            self._wrote_header = True

            elapsed = time.time() - start_time
            logger.debug(
                "Wrote %d rows to CSV file %s in %.2f seconds", rows, self.path, elapsed
            )
        except Exception as e:
            logger.error("Error writing to CSV file %s: %s", self.path, str(e))
            raise

    def _close(self):
        logger.debug("Closing CSV writer for %s", self.path)

    def run(self):
        """Continuously process the queue to write DataFrames to the CSV file."""
        while True:
            try:
                table: pd.DataFrame = self._queue.get(True)
            except Empty:
                continue

            if table is None:
                break

            self._write(table)


class WriteParquetTask(Thread):
    """
    A thread-based task for writing pandas DataFrames to a Parquet file.

    This class extends the Thread class to asynchronously write DataFrames
    to a Parquet file using a queue. It manages the ParquetWriter and schema
    internally, ensuring that data is written efficiently and safely.

    Attributes
    ----------
    path : str
        The file path where the Parquet file will be written.
    metadata : dict[str, Any]
        Metadata to be added to the Parquet file.
    """

    path: str
    metadata: dict[str, Any]

    _queue: Queue
    _schema: pa.Schema
    _writer: pq.ParquetWriter

    def __init__(self, path: str, daemon: bool = True, metadata: dict = None, **kwargs):
        """
        Initialize a WriteParquetTask instance.

        Parameters
        ----------
        path : str
            The file path where the Parquet file will be written.
        daemon : bool, optional
            Whether the thread should be a daemon thread. Defaults to True.
        metadata : dict, optional
            Metadata to be added to the Parquet file. Defaults to None.
        **kwargs
            Additional keyword arguments to be merged with metadata.
        """
        super().__init__(daemon=daemon)

        if metadata is None:
            metadata = {}
        path, _ext = os.path.splitext(path)
        path += ".parquet"

        self.path = path
        self.metadata = metadata | kwargs
        self._queue = Queue()
        self._writer = None
        self._schema = None

    def write(self, table: pd.DataFrame):
        """
        Add a DataFrame to the queue for writing to the Parquet file.

        Parameters
        ----------
        table : pd.DataFrame
            The DataFrame to be added to the queue.
        """
        logger.debug("Queuing %d rows for Parquet writing to %s", len(table), self.path)
        self._queue.put(table)

    def close(self):
        """Signal the thread to finish processing and close the file."""
        logger.debug("Closing Parquet writer queue for %s", self.path)
        self._queue.put(None)
        self.join()

    def _close(self):
        logger.debug("Closing Parquet writer for %s", self.path)
        self._writer.add_key_value_metadata(self.metadata)
        self._writer.close()

    def _write(self, table: pd.DataFrame):
        """
        Write a DataFrame to the Parquet file.

        Parameters
        ----------
        table : pd.DataFrame
            The DataFrame to be written.
        """
        start_time = time.time()
        rows = len(table)

        try:
            if self._schema is None:
                self._schema = pa.Schema.from_pandas(table, preserve_index=False)
                self._writer = pq.ParquetWriter(self.path, schema=self._schema)
                logger.debug("Initialized Parquet writer for %s", self.path)

            arrow_table = pa.Table.from_pandas(table, preserve_index=False)
            self._writer.write_table(arrow_table)

            elapsed = time.time() - start_time
            logger.debug(
                "Wrote %d rows to Parquet file %s in %.2f seconds",
                rows,
                self.path,
                elapsed,
            )
        except Exception as e:
            logger.error("Error writing to Parquet file %s: %s", self.path, str(e))
            raise

    def run(self):
        """Continuously process the queue to write DataFrames to the Parquet file."""
        while True:
            try:
                table: pd.DataFrame = self._queue.get(True)
            except Empty:
                continue

            if table is None:
                break

            self._write(table)

        self._close()
