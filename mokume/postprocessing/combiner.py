"""
Multi-file combining utilities for the mokume package.
"""

import logging
import os
from pathlib import Path

import pandas as pd

from mokume.core.constants import load_feature, load_sdrf
from mokume.postprocessing.batch_correction import get_batch_info_from_sample_names

logging.basicConfig(
    format="%(asctime)s [%(funcName)s] - %(message)s", level=logging.DEBUG
)
logger = logging.getLogger(__name__)


def folder_retrieval(folder: str) -> dict:
    """Retrieve SDRF and ibaq results from a folder."""
    folder = folder + os.sep if not folder.endswith(os.sep) else folder
    results = {"sdrf": [], "ibaq": []}
    items = os.listdir(folder)
    for item in items:
        try:
            results["sdrf"].extend(
                [
                    f"{folder}{item}/{i}"
                    for i in os.listdir(f"{folder}{item}/")
                    if i.endswith(".sdrf.tsv")
                ]
            )
            results["ibaq"].extend(
                [
                    f"{folder}{item}/{i}"
                    for i in os.listdir(f"{folder}{item}/")
                    if i.endswith("ibaq.csv") or i.endswith("ibaq.parquet")
                ]
            )
        except Exception as e:
            logger.warning(f"Error: {e}")
            if item.endswith(".sdrf.tsv"):
                results["sdrf"].append(folder + item)
            elif item.endswith("ibaq.csv"):
                results["ibaq"].append(folder + item)
    if len(results["sdrf"]) == 0:
        raise SystemExit("No SDRF founded!")
    if len(results["ibaq"]) == 0:
        raise SystemExit("No ibaq results founded!")
    if len(results["sdrf"]) != len(results["ibaq"]):
        raise SystemExit("Number of SDRFs should be equal to ibaq results!")
    return results


def generate_meta(sdrf_df: pd.DataFrame) -> pd.DataFrame:
    """Generate mokume metadata from SDRF."""
    sdrf_df.columns = [col.lower() for col in sdrf_df.columns]
    pxd = sdrf_df["source name"].values[0].split("-")[0]
    organism_part = [
        col
        for col in sdrf_df.columns
        if col.startswith("characteristics[organism part]")
    ]
    if len(organism_part) > 2:
        raise ValueError(
            f"{pxd} Please provide a maximum of 2 characteristics[organism part]"
        )
    elif len(organism_part) == 0:
        raise ValueError(
            "Missing characteristics[organism part], please check your SDRF!"
        )

    meta_df = sdrf_df[["source name"] + organism_part]
    meta_df = meta_df.drop_duplicates()

    if len(meta_df.columns) == 2:
        meta_df["tissue_part"] = None
        meta_df.columns = ["sample_id", "tissue", "tissue_part"]
    else:
        if sdrf_df[organism_part[0]].nunique() > sdrf_df[organism_part[1]].nunique():
            a, b = "tissue_part", "tissue"
        else:
            a, b = "tissue", "tissue_part"
        meta_df.rename(
            columns={
                "source name": "sample_id",
                organism_part[0]: a,
                organism_part[1]: b,
            },
            inplace=True,
        )

    meta_df["batch"] = pxd
    meta_df = meta_df[["sample_id", "batch", "tissue", "tissue_part"]]
    meta_df = meta_df.drop_duplicates()
    return meta_df


class Combiner:
    """Combine and process SDRF and iBAQ data from multiple datasets."""

    def __init__(
        self, data_folder: os.PathLike, covariate: str = None, organism: str = "HUMAN"
    ):
        self.df_pca = None
        self.df_corrected = None
        self.df_filtered_outliers = None
        self.samples_number = None
        self.datasets = None

        logger.info("Combining SDRFs and ibaq results ...")
        self.data_folder = Path(data_folder)
        if not self.data_folder.exists() or not self.data_folder.is_dir():
            raise FileNotFoundError(f"Data folder {self.data_folder} does not exist!")

        self.covariate = covariate
        files = folder_retrieval(str(self.data_folder))
        self.metadata, self.df = pd.DataFrame(), pd.DataFrame()

        for sdrf in files["sdrf"]:
            sdrf_df = load_sdrf(sdrf)
            self.metadata = pd.concat([self.metadata, generate_meta(sdrf_df)])
        self.metadata = self.metadata.drop_duplicates()
        self.metadata.index = self.metadata["sample_id"]

        for ibaq in files["ibaq"]:
            self.df = pd.concat([self.df, load_feature(ibaq)])
        self.df = self.df[self.df["ProteinName"].str.endswith(organism)]
        self.df.index = self.df["SampleID"]
        self.df = self.df.join(self.metadata, how="left")

        self.samples = self.df.columns.tolist()
        self.proteins = self.df["ProteinName"].unique().tolist()
        self.batch_index = get_batch_info_from_sample_names(self.df.columns.tolist())

        logger.info(self.metadata, self.df.head)
