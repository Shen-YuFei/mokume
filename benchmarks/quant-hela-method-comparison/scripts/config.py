"""
HeLa Benchmark Configuration

Configuration for benchmarking Rust protein quantification methods across datasets.
All datasets use public PRIDE resources with verified file names.
"""

from pathlib import Path
from dataclasses import dataclass
from typing import Dict

# Base directories
BENCHMARK_DIR = Path(__file__).resolve().parent.parent
RAW_DATA_DIR = BENCHMARK_DIR / "data" / "raw"
PROCESSED_DIR = BENCHMARK_DIR / "data" / "processed"
PROTEIN_QUANT_DIR = BENCHMARK_DIR / "data" / "protein_quant"
ANALYSIS_DIR = BENCHMARK_DIR / "results"
PLOTS_DIR = BENCHMARK_DIR / "figures"

# Create directories if they don't exist
for d in [RAW_DATA_DIR, PROCESSED_DIR, PROTEIN_QUANT_DIR, ANALYSIS_DIR, PLOTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# FASTA is passed explicitly to the current Rust refresh script.
FASTA_FILE = BENCHMARK_DIR / "data" / "human-reviewed.fasta"

# PRIDE FTP base URL
PRIDE_IBAQPY_FTP = (
    "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/"
    "ibaqpy-research"
)


@dataclass
class DatasetInfo:
    """Information about a proteomics dataset."""

    project_id: str
    name: str
    acquisition_method: str
    ftp_base: str
    feature_file: str
    description: str = ""
    file_format: str = "parquet"
    enzyme: str = "Trypsin"

    @property
    def feature_url(self) -> str:
        return f"{self.ftp_base}/{self.feature_file}"

    @property
    def local_feature_path(self) -> Path:
        if self.file_format == "msstats_csv":
            return RAW_DATA_DIR / self.feature_file
        return RAW_DATA_DIR / f"{self.project_id}_feature.parquet"

    @property
    def peptide_path(self) -> Path:
        return PROCESSED_DIR / f"{self.project_id}_peptides.parquet"


# ============================================================================
# ALL DATASETS FROM ibaqpy-research FTP (20 datasets)
# https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/ibaqpy-research/
# ============================================================================

IBAQPY_DATASETS: Dict[str, DatasetInfo] = {
    # PMID-based datasets (Nagaraj et al. 2011 - HeLa Deep Proteome)
    "PMID22068331.1": DatasetInfo(
        project_id="PMID22068331.1",
        name="HeLa Deep Proteome 1",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PMID22068331.1",
        feature_file="PMID22068331-5cefcd8e-ec33-4bcb-9703-4f3580f6d3ef.feature.parquet",
        description="HeLa deep proteome - replicate 1 (Nagaraj et al.)",
    ),
    "PMID22068331.2": DatasetInfo(
        project_id="PMID22068331.2",
        name="HeLa Deep Proteome 2",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PMID22068331.2",
        feature_file="PMID22068331-093b1d91-2dd0-487f-b116-66968b82931b.feature.parquet",
        description="HeLa deep proteome - replicate 2 (Nagaraj et al.)",
        enzyme="Lys-C",
    ),
    "PMID22068331.3": DatasetInfo(
        project_id="PMID22068331.3",
        name="HeLa Deep Proteome 3",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PMID22068331.3",
        feature_file="PMID22068331-0dc35a39-5ab2-45e5-923b-165e7253c549.feature.parquet",
        description="HeLa deep proteome - replicate 3 (Nagaraj et al.)",
        enzyme="Glu-C",
    ),
    # PXD datasets
    "PXD000269": DatasetInfo(
        project_id="PXD000269",
        name="Human Cell Lines",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD000269",
        feature_file="PXD000269-bfe06112-8412-4b4c-9f5c-21f444a9987a.feature.parquet",
        description="Human cell line proteomics",
    ),
    "PXD000396": DatasetInfo(
        project_id="PXD000396",
        name="HeLa Proteome",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD000396",
        feature_file="PXD000396-00e25b1a-691c-46a8-979a-b5cc9a6f9d81.feature.parquet",
        description="HeLa proteome dataset",
    ),
    "PXD001305": DatasetInfo(
        project_id="PXD001305",
        name="Human Proteome Map",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD001305",
        feature_file="PXD001305-708ff169-a7f5-4d4f-8e80-fc9fa31b4115.feature.parquet",
        description="Human proteome map",
    ),
    "PXD004452": DatasetInfo(
        project_id="PXD004452",
        name="HeLa DDA-LFQ",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD004452",
        feature_file="PXD004452-6c224f5a-7c1f-46f9-9dae-1541baeef8fe.feature.parquet",
        description="Large HeLa dataset",
    ),
    "PXD004940": DatasetInfo(
        project_id="PXD004940",
        name="Human Proteome",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD004940",
        feature_file="PXD004940-0db4e73c-2a08-4d18-bae1-d498f4cd73ea.feature.parquet",
        description="Human proteome dataset",
    ),
    "PXD005481": DatasetInfo(
        project_id="PXD005481",
        name="HeLa Single Shot",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD005481",
        feature_file="PXD005481-b0ea9a9d-87ac-4969-8057-0a4b8e7898b9.feature.parquet",
        description="HeLa single-shot proteomics",
    ),
    "PXD009273.1": DatasetInfo(
        project_id="PXD009273.1",
        name="Human Tissues",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD009273.1",
        feature_file="PXD009273-31762cb6-c2a6-4ab3-b394-89484c2634db.feature.parquet",
        description="Human tissue proteomics",
    ),
    "PXD009875": DatasetInfo(
        project_id="PXD009875",
        name="HeLa DIA",
        acquisition_method="DIA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD009875",
        feature_file="PXD009875-0a773c5a-51a7-4f5a-ac3e-6f94a4794f0f.feature.parquet",
        description="HeLa DIA dataset",
    ),
    "PXD010150": DatasetInfo(
        project_id="PXD010150",
        name="Human Proteome Study",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD010150",
        feature_file="PXD010150-17192840-7693-457f-97bf-d0feb2ebeab4.feature.parquet",
        description="Human proteome quantification study",
    ),
    "PXD013658.1": DatasetInfo(
        project_id="PXD013658.1",
        name="Human Cell Atlas 1",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD013658.1",
        feature_file="PXD013658-573f7962-8801-45ac-bd12-df0ee92ac216.feature.parquet",
        description="Human cell atlas - part 1",
    ),
    "PXD013658.2": DatasetInfo(
        project_id="PXD013658.2",
        name="Human Cell Atlas 2",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD013658.2",
        feature_file="PXD013658-e9c11b7c-ded9-45eb-8015-ccd334ccd7c8.feature.parquet",
        description="Human cell atlas - part 2",
    ),
    "PXD014058": DatasetInfo(
        project_id="PXD014058",
        name="Human Proteome Reference",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD014058",
        feature_file="PXD014058-fdfd18c7-7976-4fb9-a845-878f65afb44f.feature.parquet",
        description="Human proteome reference dataset",
    ),
    "PXD030406": DatasetInfo(
        project_id="PXD030406",
        name="Human Sample Study",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD030406",
        feature_file="PXD030406-5606cc1b-d16c-4509-b174-514375926d51.feature.parquet",
        description="Human sample proteomics study",
    ),
    "PXD039414": DatasetInfo(
        project_id="PXD039414",
        name="Human Cell Proteome",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD039414",
        feature_file="PXD039414-a7e7ae89-fe74-40c4-b5bc-84e68254e06d.feature.parquet",
        description="Human cell proteome dataset",
    ),
    "PXD042233": DatasetInfo(
        project_id="PXD042233",
        name="Human Proteomics",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD042233",
        feature_file="PXD042233-4d7f9915-94f0-4809-a2cc-4160e4d67216.feature.parquet",
        description="Human proteomics dataset",
    ),
    "PXD048325": DatasetInfo(
        project_id="PXD048325",
        name="Human Sample Analysis",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD048325",
        feature_file="PXD048325-29419338-aa11-47fb-b183-db1650465b57.feature.parquet",
        description="Human sample proteomics analysis (large dataset)",
    ),
    "PXD054015.1": DatasetInfo(
        project_id="PXD054015.1",
        name="Recent Human Study",
        acquisition_method="DDA-LFQ",
        ftp_base=PRIDE_IBAQPY_FTP + "/PXD054015.1",
        feature_file="PXD054015-ecb7ded0-a1cf-46d8-93cc-f2713ea6ccf7.feature.parquet",
        description="Recent human proteomics study",
    ),
}

# ============================================================================
# TMT vs LFQ COMPARISON DATASET (PXD007683)
# ============================================================================

TMT_LFQ_DATASETS: Dict[str, DatasetInfo] = {
    "PXD007683-TMT": DatasetInfo(
        project_id="PXD007683-TMT",
        name="HeLa TMT Benchmark",
        acquisition_method="TMT",
        ftp_base=(
            "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/"
            "proteomes/quantms-benchmark-old/PXD007683-TMT/msstatsconverter"
        ),
        feature_file="PXD007683-TMT.sdrf_openms_design_msstats_in.csv",
        description="HeLa benchmark - TMT labeling",
        file_format="msstats_csv",
    ),
    "PXD007683-LFQ": DatasetInfo(
        project_id="PXD007683-LFQ",
        name="HeLa LFQ Benchmark",
        acquisition_method="DDA-LFQ",
        ftp_base=(
            "https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/"
            "proteomes/quantms-benchmark-old/PXD007683-LFQ/proteomicslfq"
        ),
        feature_file="PXD007683-LFQ.sdrf_openms_design_msstats_in.csv",
        description="HeLa benchmark - Label-free",
        file_format="msstats_csv",
    ),
}

# All datasets
ALL_DATASETS = {**IBAQPY_DATASETS, **TMT_LFQ_DATASETS}

# HeLa-specific datasets (subset)
HELA_DATASETS = {
    k: v
    for k, v in IBAQPY_DATASETS.items()
    if "hela" in v.name.lower()
    or "hela" in v.description.lower()
    or k.startswith("PMID22068331")
}

# TMT/LFQ comparison alias
TMT_LFQ_COMPARISON = TMT_LFQ_DATASETS


# ============================================================================
# QUANTIFICATION CONFIGURATION
# ============================================================================

# Quantification methods to benchmark
QUANTIFICATION_METHODS = [
    "pibaq",
    "maxlfq",
    "directlfq",
    "top3",
    "top10",
    "sum",
]

INTENSITY_COLUMNS = {
    "pibaq": "PiBAQ",
    **{method: "Intensity" for method in QUANTIFICATION_METHODS if method != "pibaq"},
}

# Proteins of interest to track across experiments
PROTEINS_OF_INTEREST = [
    "P60709",  # ACTB - Actin beta (housekeeping)
    "P07437",  # TUBB - Tubulin beta (cytoskeleton)
    "P68371",  # TUBB4B - Tubulin beta-4B chain
    "P04406",  # GAPDH - Glyceraldehyde-3-phosphate dehydrogenase
    "P06733",  # ENO1 - Alpha-enolase
    "P07355",  # ANXA2 - Annexin A2
    "P08238",  # HSP90AB1 - Heat shock protein 90 beta
    "P14625",  # HYOU1 - Endoplasmin (HSP90B1)
    "P11142",  # HSPA8 - Heat shock cognate 71 kDa protein
    "P62258",  # YWHAE - 14-3-3 protein epsilon
    "P31946",  # YWHAB - 14-3-3 protein beta/alpha
    "P63104",  # YWHAZ - 14-3-3 protein zeta/delta
    "P68104",  # EEF1A1 - Elongation factor 1-alpha 1
    "P13639",  # EEF2 - Elongation factor 2
    "P62937",  # PPIA - Peptidyl-prolyl cis-trans isomerase A
    "P23528",  # CFL1 - Cofilin-1
    "P06748",  # NPM1 - Nucleophosmin
    "P62805",  # H4C1 - Histone H4
    "P62807",  # HIST1H2BB - Histone H2B type 1-C/E/F/G/I
    "P16401",  # HIST1H1B - Histone H1.5
]

# piBAQ parameters
PIBAQ_PARAMS = {
    "enzyme": "Trypsin",
    "min_aa": 7,
    "max_aa": 30,
    "missed_cleavages": 0,
}

# ============================================================================
# COLUMN MAPPINGS FOR DIFFERENT FILE FORMATS
# ============================================================================

QUANTMS_PARQUET_MAPPING = {
    # Protein columns
    "protein_accessions": "ProteinName",
    "pg": "ProteinName",
    "pg_accessions": "ProteinName",  # ibaqpy-research format
    "ProteinName": "ProteinName",
    # Peptide columns
    "PeptideSequence": "PeptideSequence",
    "sequence": "PeptideSequence",
    "peptidoform": "PeptideSequence",
    # Sample columns
    "SampleID": "SampleID",
    "sample_accession": "SampleID",
    "Run": "SampleID",
    "reference_file_name": "SampleID",  # ibaqpy-research format
    # Intensity columns
    "NormIntensity": "NormIntensity",
    "intensity": "NormIntensity",
    "Intensity": "NormIntensity",
    "intensities": "NormIntensity",  # ibaqpy-research format
    # Metadata columns
    "Condition": "Condition",
    "condition": "Condition",
    "BioReplicate": "BioReplicate",
    "biological_replicate": "BioReplicate",
}
