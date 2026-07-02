#!/usr/bin/env python3
"""
Download PXD007683 benchmark data from PRIDE.

Downloads:
- TMT parquet file (QPX format)
- LFQ parquet file (QPX format)
- FASTA file (for iBAQ quantification)
"""

import urllib.request
import urllib.parse
import sys
from pathlib import Path

_opener = urllib.request.build_opener()

# Add parent to path for config import
sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_DIR, DATA_URLS, FASTA_URL


def download_file(url: str, dest: Path, description: str = "") -> bool:
    """Download a file with progress reporting."""
    if dest.exists():
        print(f"  [SKIP] {dest.name} already exists")
        return True

    print(f"  Downloading {description or dest.name}...")
    print(f"    URL: {url}")

    try:
        if urllib.parse.urlparse(url).scheme not in ("http", "https"):
            raise ValueError(f"URL scheme not allowed: {url}")
        with _opener.open(url) as response:
            total_size = int(response.headers.get("Content-Length", 0))
            downloaded = 0
            block_size = 8192
            with open(dest, "wb") as out_file:
                while True:
                    block = response.read(block_size)
                    if not block:
                        break
                    out_file.write(block)
                    downloaded += len(block)
                    if total_size > 0:
                        percent = min(100, downloaded * 100 / total_size)
                        mb_downloaded = downloaded / (1024 * 1024)
                        mb_total = total_size / (1024 * 1024)
                        print(f"\r    Progress: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)", end="")
        print()  # newline after progress
        print(f"    Saved to: {dest}")
        return True

    except Exception as e:
        print(f"\n    ERROR: {e}")
        return False


def main():
    """Download all benchmark data files."""
    print("=" * 60)
    print("PXD007683 Benchmark Data Download")
    print("=" * 60)

    # Create data directory
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\nData directory: {DATA_DIR}")

    success = True

    # Download parquet files
    print("\n[1/3] Downloading TMT parquet...")
    if not download_file(DATA_URLS["tmt"], DATA_DIR / "PXD007683-TMT.parquet", "TMT data"):
        success = False

    print("\n[2/3] Downloading LFQ parquet...")
    if not download_file(DATA_URLS["lfq"], DATA_DIR / "PXD007683-LFQ.parquet", "LFQ data"):
        success = False

    print("\n[3/3] Downloading FASTA file...")
    if not download_file(FASTA_URL, DATA_DIR / "PXD007683.fasta", "FASTA"):
        success = False

    # Summary
    print("\n" + "=" * 60)
    if success:
        print("Download complete!")
        print("\nFiles downloaded:")
        for f in DATA_DIR.glob("*"):
            size_mb = f.stat().st_size / (1024 * 1024)
            print(f"  - {f.name}: {size_mb:.1f} MB")
    else:
        print("Some downloads failed. Check the errors above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
