"""
Filter pipeline for orchestrating multiple preprocessing filters.
"""

from typing import List, Tuple

import pandas as pd

from mokume.core.logger import get_logger
from mokume.preprocessing.filters.base import BaseFilter, FilterResult


logger = get_logger("mokume.preprocessing.filters.pipeline")


class FilterPipeline:
    """
    Pipeline for applying multiple preprocessing filters in sequence.

    The pipeline maintains an ordered list of filters and applies them
    sequentially to the input DataFrame, collecting results from each filter.
    """

    def __init__(self, name: str = "default"):
        """
        Initialize the filter pipeline.

        Parameters
        ----------
        name : str, optional
            Name for the pipeline (for logging).
        """
        self.name = name
        self.filters: List[BaseFilter] = []

    def add_filter(self, filter_obj: BaseFilter) -> "FilterPipeline":
        """
        Add a filter to the pipeline.

        Parameters
        ----------
        filter_obj : BaseFilter
            Filter to add.

        Returns
        -------
        FilterPipeline
            Self for method chaining.
        """
        self.filters.append(filter_obj)
        return self

    def add_filters(self, filters: List[BaseFilter]) -> "FilterPipeline":
        """
        Add multiple filters to the pipeline.

        Parameters
        ----------
        filters : List[BaseFilter]
            Filters to add.

        Returns
        -------
        FilterPipeline
            Self for method chaining.
        """
        self.filters.extend(filters)
        return self

    def apply(
        self,
        df: pd.DataFrame,
        stop_on_empty: bool = True,
        **kwargs,
    ) -> Tuple[pd.DataFrame, List[FilterResult]]:
        """
        Apply all filters in the pipeline to a DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame to filter.
        stop_on_empty : bool, optional
            Whether to stop if DataFrame becomes empty.
        **kwargs
            Additional arguments passed to each filter.

        Returns
        -------
        Tuple[pd.DataFrame, List[FilterResult]]
            Filtered DataFrame and list of results from each filter.
        """
        results = []
        current_df = df

        logger.debug(
            "Starting filter pipeline '%s' with %d filters",
            self.name,
            len(self.filters),
        )

        for filter_obj in self.filters:
            if stop_on_empty and len(current_df) == 0:
                logger.warning(
                    "Pipeline '%s': DataFrame is empty, stopping at filter '%s'",
                    self.name,
                    filter_obj.name,
                )
                break

            try:
                current_df, result = filter_obj.apply(current_df, **kwargs)
                results.append(result)

                logger.debug(
                    "Pipeline '%s': %s removed %d/%d (%.1f%%)",
                    self.name,
                    result.filter_name,
                    result.removed_count,
                    result.input_count,
                    result.removal_rate * 100,
                )

            except Exception as e:
                logger.error(
                    "Pipeline '%s': Error in filter '%s': %s",
                    self.name,
                    filter_obj.name,
                    str(e),
                )
                raise

        return current_df, results

    def summary(self, results: List[FilterResult]) -> dict:
        """
        Generate a summary of filter results.

        Parameters
        ----------
        results : List[FilterResult]
            Results from apply().

        Returns
        -------
        dict
            Summary statistics.
        """
        if not results:
            return {"total_input": 0, "total_output": 0, "total_removed": 0}

        total_input = results[0].input_count if results else 0
        total_output = results[-1].output_count if results else 0
        total_removed = total_input - total_output

        return {
            "pipeline_name": self.name,
            "total_input": total_input,
            "total_output": total_output,
            "total_removed": total_removed,
            "total_removal_rate": (
                total_removed / total_input if total_input > 0 else 0.0
            ),
            "filters_applied": len(results),
            "filter_details": [
                {
                    "name": r.filter_name,
                    "removed": r.removed_count,
                    "removal_rate": r.removal_rate,
                }
                for r in results
            ],
        }

    def __len__(self) -> int:
        return len(self.filters)

    def __repr__(self) -> str:
        filter_names = [f.name for f in self.filters]
        return f"FilterPipeline(name='{self.name}', filters={filter_names})"
