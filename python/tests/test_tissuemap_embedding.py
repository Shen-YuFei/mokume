"""Regression tests for optional TissueMap embedding dependencies."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


PYTHON_ROOT = Path(__file__).parents[1]


def test_tsne_module_import_does_not_load_optional_umap(tmp_path):
    """A broken optional UMAP install must not block the independent t-SNE path."""
    (tmp_path / "umap.py").write_text(
        "raise RuntimeError('UMAP should be loaded only when requested')\n",
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join((str(tmp_path), str(PYTHON_ROOT)))
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import mokume.tissuemap.embedding; "
            "assert 'umap' not in sys.modules",
        ],
        capture_output=True,
        check=False,
        env=env,
        text=True,
    )

    assert result.returncode == 0, result.stderr
