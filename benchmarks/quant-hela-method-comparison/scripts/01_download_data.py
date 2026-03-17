#!/usr/bin/env python3
"""
Phase 1: Data Acquisition

Download datasets from PRIDE ibaqpy-research FTP for benchmarking protein quantification methods.
"""

import shutil
import urllib.request
import urllib.error
import urllib.parse
import sys
from pathlib import Path
from typing import List

_opener = urllib.request.build_opener()
from config import (
    ALL_DATASETS,
    HELA_DATASETS,
    RAW_DATA_DIR,
    DatasetInfo,
)


def download_file(url: str, dest: Path, verbose: bool = True) -> bool:
    """
    Download a file from URL to destination.

    Parameters
    ----------
    url : str
        Source URL
    dest : Path
        Destination file path
    verbose : bool
        Print progress messages

    Returns
    -------
    bool
        True if download successful, False otherwise
    """
    if dest.exists():
        if verbose:
            print(f"  Already exists: {dest.name}")
        return True

    if verbose:
        print(f"  Downloading: {url}")

    try:
        if urllib.parse.urlparse(url).scheme not in ("http", "https"):
            raise ValueError(f"URL scheme not allowed: {url}")
        with _opener.open(url) as response, open(dest, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
        if verbose:
            size_mb = dest.stat().st_size / (1024 * 1024)
            print(f"  Saved to: {dest.name} ({size_mb:.1f} MB)")
        return True
    except urllib.error.HTTPError as e:
        print(f"  HTTP Error {e.code}: {url}")
        return False
    except urllib.error.URLError as e:
        print(f"  URL Error: {e.reason}")
        return False
    except Exception as e:
        print(f"  Error: {e}")
        return False


def check_url_exists(url: str) -> bool:
    """Check if a URL exists without downloading."""
    try:
        if urllib.parse.urlparse(url).scheme not in ("http", "https"):
            raise ValueError(f"URL scheme not allowed: {url}")
        request = urllib.request.Request(url, method='HEAD')
        _opener.open(request, timeout=10)
        return True
    except Exception:
        return False


def list_available_files(base_url: str) -> List[str]:
    """
    List files available at a URL (works for directory listings).
    """
    try:
        if urllib.parse.urlparse(base_url).scheme not in ("http", "https"):
            raise ValueError(f"URL scheme not allowed: {base_url}")
        response = _opener.open(base_url, timeout=30)
        content = response.read().decode('utf-8')
        import re
        links = re.findall(r'href="([^"]+)"', content)
        return [l for l in links if not l.startswith('/') and not l.startswith('?')]
    except Exception as e:
        print(f"  Could not list files: {e}")
        return []


def download_dataset(dataset: DatasetInfo, verbose: bool = True) -> dict:
    """
    Download a dataset's feature file.

    Returns dict with status of download.
    """
    results = {
        "project_id": dataset.project_id,
        "feature": False,
    }

    if verbose:
        print(f"\n{dataset.project_id}: {dataset.name}")
        print(f"  Method: {dataset.acquisition_method}")
        print(f"  Format: {dataset.file_format}")

    # Download feature/peptide file
    feature_url = dataset.feature_url
    feature_dest = dataset.local_feature_path

    if download_file(feature_url, feature_dest, verbose=verbose):
        results["feature"] = True
        results["feature_path"] = str(feature_dest)
    else:
        if verbose:
            print(f"  Feature file not found: {feature_url}")

    return results


def discover_dataset_files(dataset: DatasetInfo) -> dict:
    """
    Try to discover available files for a dataset.

    Returns dict with discovered URLs.
    """
    print(f"\nDiscovering files for {dataset.project_id}...")

    discovered = {"project_id": dataset.project_id, "files": []}

    # Check base URL
    base_url = dataset.ftp_base
    files = list_available_files(base_url)

    if files:
        print(f"  Found {len(files)} files at {base_url}")
        for f in files[:10]:  # Show first 10
            print(f"    - {f}")
        if len(files) > 10:
            print(f"    ... and {len(files) - 10} more")
        discovered["files"] = files
    else:
        print(f"  Could not list files at {base_url}")

    return discovered


def download_all_datasets(
    datasets: dict = None,
    discover_only: bool = False,
    verbose: bool = True
) -> dict:
    """
    Download all datasets.

    Parameters
    ----------
    datasets : dict, optional
        Dict of DatasetInfo objects. Defaults to ALL_DATASETS.
    discover_only : bool
        If True, only discover available files without downloading.
    verbose : bool
        Print progress messages.

    Returns
    -------
    dict
        Summary of download results.
    """
    if datasets is None:
        datasets = ALL_DATASETS

    print("=" * 70)
    print("Benchmark Data Acquisition")
    print("=" * 70)
    print(f"\nDatasets to process: {len(datasets)}")
    print(f"Output directory: {RAW_DATA_DIR}")

    results = {}

    for dataset_id, dataset in datasets.items():
        if discover_only:
            results[dataset_id] = discover_dataset_files(dataset)
        else:
            results[dataset_id] = download_dataset(dataset, verbose=verbose)

    # Summary
    print("\n" + "=" * 70)
    print("Download Summary")
    print("=" * 70)

    if not discover_only:
        success_count = sum(1 for r in results.values() if r.get("feature", False))
        print(f"\nSuccessfully downloaded: {success_count}/{len(datasets)} datasets")

        for dataset_id, result in results.items():
            status = "OK" if result.get("feature", False) else "MISSING"
            print(f"  {dataset_id}: {status}")

    return results


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Download datasets for benchmarking"
    )
    parser.add_argument(
        "--discover",
        action="store_true",
        help="Only discover available files, don't download"
    )
    parser.add_argument(
        "--dataset",
        type=str,
        help="Download specific dataset (project ID)"
    )
    parser.add_argument(
        "--hela-only",
        action="store_true",
        help="Only download HeLa-specific datasets"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Reduce output verbosity"
    )

    args = parser.parse_args()

    # Select datasets
    if args.dataset:
        if args.dataset in ALL_DATASETS:
            datasets = {args.dataset: ALL_DATASETS[args.dataset]}
        else:
            print(f"Unknown dataset: {args.dataset}")
            print(f"Available: {', '.join(ALL_DATASETS.keys())}")
            sys.exit(1)
    elif args.hela_only:
        datasets = HELA_DATASETS
    else:
        datasets = ALL_DATASETS

    # Run
    results = download_all_datasets(
        datasets=datasets,
        discover_only=args.discover,
        verbose=not args.quiet
    )

    # Exit code based on success
    if not args.discover:
        success_count = sum(1 for r in results.values() if r.get("feature", False))
        sys.exit(0 if success_count > 0 else 1)


if __name__ == "__main__":
    main()
