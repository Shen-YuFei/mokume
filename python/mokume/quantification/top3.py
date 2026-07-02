"""
Top3 protein quantification method.

This module provides backward compatibility for Top3 quantification.
Top3 is now an alias for TopNQuantification(n=3).
"""

from mokume.quantification.topn import TopNQuantification


class Top3Quantification(TopNQuantification):
    """
    Top3 protein quantification method.

    This is an alias for TopNQuantification(n=3) for backward compatibility.
    Calculates protein abundance as the average of the three most intense
    peptides for each protein in each sample.
    """

    def __init__(self):
        """Initialize Top3 quantification (TopN with n=3)."""
        super().__init__(n=3)
