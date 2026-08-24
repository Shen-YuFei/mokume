"""
CLI command for t-SNE visualization.

Note: This module requires the optional 'plotting' dependencies.
Install them with: pip install mokume-py[plotting]
"""

import glob
import inspect
import logging

import click
import numpy as np
import pandas as pd

from mokume.core.constants import PROTEIN_NAME, SAMPLE_ID, PIBAQ_LOG
from mokume.plotting import is_plotting_available


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def compute_tsne(df_pca, n_components=2, perplexity=30, learning_rate=200, n_iter=2000):
    """Compute t-SNE components from PCA components.

    scikit-learn renamed ``TSNE``'s ``n_iter`` argument to ``max_iter`` in 1.5
    and removed ``n_iter`` outright in 1.7, so select whichever keyword the
    installed ``TSNE`` accepts instead of passing ``n_iter=`` unconditionally
    (which raises ``TypeError`` on sklearn >= 1.7).
    """
    try:
        from sklearn.manifold import TSNE
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required for t-SNE visualization. "
            "Install it with: pip install mokume-py[plotting]"
        ) from exc

    iter_kw = "n_iter" if "n_iter" in inspect.signature(TSNE).parameters else "max_iter"
    tsne = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        learning_rate=learning_rate,
        **{iter_kw: n_iter},
    )
    tsne_results = tsne.fit_transform(np.asarray(df_pca))
    tsne_cols = [f"tSNE{i + 1}" for i in range(n_components)]
    df_tsne = pd.DataFrame(data=tsne_results, columns=tsne_cols)
    df_tsne.index = df_pca.index
    return df_tsne


@click.command()
@click.option(
    "-f", "--folder", help="Folder that contains all the protein files", required=True
)
@click.option(
    "-o",
    "--pattern",
    help="Protein file pattern",
    required=False,
    default="proteins.tsv",
)
def tsne_visualization(folder: str = None, pattern: str = "proteins.tsv"):
    """Generate a t-SNE visualization for protein data from specified files."""
    if not is_plotting_available():
        raise click.ClickException(
            "Plotting dependencies (matplotlib, seaborn) are not installed. "
            "Install them with: pip install mokume-py[plotting]"
        )

    from mokume.plotting import compute_pca_with_plot, plot_tsne

    files = glob.glob(f"{folder}/*{pattern}")
    dfs = []

    for f in files:
        reanalysis = (f.split("/")[-1].split("_")[0]).replace("-proteins.tsv", "")
        dfs += [
            pd.read_csv(
                f,
                usecols=[PROTEIN_NAME, SAMPLE_ID, PIBAQ_LOG],
                sep=None,
                engine="python",
            ).assign(reanalysis=reanalysis)
        ]

    total_proteins = pd.concat(dfs, ignore_index=True)

    normalize_df = pd.pivot_table(
        total_proteins,
        index=[SAMPLE_ID, "reanalysis"],
        columns=PROTEIN_NAME,
        values=PIBAQ_LOG,
    )
    normalize_df = normalize_df.fillna(0)
    df_pca = compute_pca_with_plot(normalize_df, n_components=30)
    df_tsne = compute_tsne(df_pca)

    batch = df_tsne.index.get_level_values("reanalysis").tolist()
    df_tsne["batch"] = batch

    plot_tsne(
        df_tsne, "tSNE1", "tSNE2", "batch", "5.tsne_plot_with_batch_information.pdf"
    )
    logger.info(total_proteins.shape)


if __name__ == "__main__":
    tsne_visualization()
