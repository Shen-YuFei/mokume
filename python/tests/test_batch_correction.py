import logging
from pathlib import Path

import pytest
import pandas as pd

from mokume.commands.batch_correct import (
    get_batch_id_from_sample_names,
    run_batch_correction,
)
from mokume.core.constants import SAMPLE_ID, PROTEIN_NAME, IBAQ, IBAQ_BEC

TESTS_DIR = Path(__file__).parent

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@pytest.mark.parametrize(
    "samples",
    [
        ["PXD000001-A", "PXD000001-B", "PXD000002-A"],
        pd.Index(["PXD000001-A", "PXD000001-B", "PXD000002-A"]),
    ],
)
def test_get_batch_id_supports_python_and_pandas_sequences(samples):
    """Batch labels remain stable for Python and pandas sequences."""
    assert get_batch_id_from_sample_names(samples).tolist() == [0, 0, 1]


def test_correct_batches():
    """
    Test the `run_batch_correction` function to ensure it correctly processes iBAQ values
    from TSV files, generates the expected output files, and handles various error cases.

    This test verifies:
    - The creation and non-emptiness of the corrected output TSV file.
    - The creation and correct shape of the AnnData object.
    - Handling of invalid sample IDs by raising a ValueError.
    - Handling of missing required columns by raising a ValueError.
    - Handling of invalid file patterns by raising a ValueError.
    """
    pytest.importorskip("inmoose", reason="inmoose is required for batch correction")
    anndata = pytest.importorskip(
        "anndata", reason="anndata is required for batch correction"
    )

    args = {
        "folder": TESTS_DIR / "ibaq-raw-hela",
        "pattern": "*ibaq.tsv",
        "comment": "#",
        "sep": "\t",
        "output": TESTS_DIR / "example/ibaq_corrected_combined.tsv",
        "sample_id_column": SAMPLE_ID,
        "protein_id_column": PROTEIN_NAME,
        "ibaq_raw_column": IBAQ,
        "ibaq_corrected_column": IBAQ_BEC,
        "export_anndata": True,
    }
    logging.debug("Arguments for run_batch_correction: %s", args)
    run_batch_correction(**args)

    # Assert the output file is created and not empty
    output_path = Path(args["output"])
    assert output_path.exists(), f"Expected output file {output_path} was not created."
    df = pd.read_csv(output_path, sep=args["sep"])
    assert not df.empty, "The corrected output file is empty."

    # Assert the AnnData object is created
    adata_path = output_path.with_suffix(".h5ad")
    assert adata_path.exists(), f"Expected AnnData file {adata_path} was not created."

    # Read the AnnData object and check shape and layers
    adata = anndata.read_h5ad(adata_path)
    logger.info(adata)
    assert adata.shape == (46, 3476)
    assert adata.layers[IBAQ_BEC].shape == (46, 3476)

    # Test invalid sample IDs
    with pytest.raises(ValueError):
        args["folder"] = TESTS_DIR / "invalid-samples"
        run_batch_correction(**args)

    # Test missing required columns
    with pytest.raises(ValueError):
        args["folder"] = TESTS_DIR / "ibaq-raw-hela"
        args["sample_id_column"] = "NonexistentColumn"
        run_batch_correction(**args)

    # Test invalid file pattern
    with pytest.raises(ValueError):
        args["pattern"] = "nonexistent*.tsv"
        run_batch_correction(**args)


if __name__ == "__main__":
    test_correct_batches()
