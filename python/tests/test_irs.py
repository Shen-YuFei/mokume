"""Regression tests for IRS metadata helpers."""

from mokume.normalization.irs import detect_plexes_from_sdrf


def test_detect_plexes_accepts_tmt_n_and_c_channel_suffixes(tmp_path):
    """TMT channel labels must not create a synthetic ``plex1`` batch."""
    sdrf = tmp_path / "tmt.sdrf.tsv"
    sdrf.write_text(
        "source name\tcomment[label]\n"
        "sample_Mixture1_126\tTMT126\n"
        "condition_a_Mixture1_127N\tTMT127N\n"
        "condition_b_Mixture1_127C\tTMT127C\n"
        "p2_127N\tTMT127N\n",
        encoding="utf-8",
    )

    assert detect_plexes_from_sdrf(str(sdrf)) == {
        "sample_Mixture1_126": "mixture1",
        "condition_a_Mixture1_127N": "mixture1",
        "condition_b_Mixture1_127C": "mixture1",
        "p2_127N": "p2",
    }
