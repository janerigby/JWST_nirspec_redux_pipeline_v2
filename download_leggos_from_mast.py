#!/usr/bin/env python3
"""Download JWST UNCAL FITS products for one program from MAST or AWS.

Run as, for example:
    python download_leggos_from_mast.py 04125 AWS
    python download_leggos_from_mast.py 04125 MAST
"""

import argparse
import re
import sys
import tempfile
from pathlib import Path

import numpy as np
from astroquery.mast import Observations
from astropy.table import Table, unique, vstack


DEFAULT_PID = 4125
DOWNLOAD_SOURCE = "AWS"
OBS_COLLECTION = "JWST"
PRODUCT_LIST_CHUNK_SIZE = 100
DOWNLOAD_CHUNK_SIZE = 100


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Query MAST and download every raw level-1 UNCAL FITS product "
            "for a JWST program from MAST or public AWS S3."
        )
    )
    parser.add_argument("pid", nargs="?", type=int, help="MAST proposal/program ID")
    parser.add_argument(
        "download_source",
        nargs="?",
        type=str.upper,
        choices=("MAST", "AWS"),
        help=f"download source (default configured in script: {DOWNLOAD_SOURCE})",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("Rawdata"),
        help="root for PID and target directories (default: Rawdata)",
    )
    return parser.parse_args(argv)


def choose_pid(pid):
    """Return the supplied PID or interactively confirm the script default."""
    if pid is not None:
        if pid <= 0:
            raise ValueError("PID must be a positive integer")
        return pid

    answer = input(
        f"No PID was supplied. Proceed with the default PID {DEFAULT_PID}? [y/N] "
    ).strip().lower()
    if answer not in {"y", "yes"}:
        print("Cancelled; no query or download was started.")
        return None
    return DEFAULT_PID


def table_chunks(table, size):
    for start in range(0, len(table), size):
        yield table[start : start + size]


def get_products(observations):
    """Retrieve product metadata in bounded requests and remove duplicates."""
    product_tables = []
    chunks = list(table_chunks(observations, PRODUCT_LIST_CHUNK_SIZE))
    for number, observation_chunk in enumerate(chunks, start=1):
        print(f"Retrieving product list {number}/{len(chunks)} ...", flush=True)
        product_tables.append(Observations.get_product_list(observation_chunk))

    if not product_tables:
        return Table()
    products = vstack(product_tables, metadata_conflicts="silent")
    key = "dataURI" if "dataURI" in products.colnames else "productFilename"
    return unique(products, keys=key)


def configure_download_source(source_override):
    """Select direct MAST downloads or the public AWS dataset."""
    source = source_override or DOWNLOAD_SOURCE
    if source == "MAST":
        Observations.disable_cloud_dataset()
    elif source == "AWS":
        try:
            Observations.enable_cloud_dataset(provider="AWS", profile=None)
        except ModuleNotFoundError as exc:
            raise RuntimeError("AWS downloads require the boto3 package") from exc
    else:
        raise ValueError("Download source must be MAST or AWS")
    print(f"Download source: {source}", flush=True)
    return source


def is_level_one(value):
    if np.ma.is_masked(value):
        return False
    try:
        return float(value) == 1.0
    except (TypeError, ValueError):
        return False


def select_raw_fits(products):
    """Select calibration-level 1 products whose names end in uncal.fits."""
    required = {"calib_level", "productFilename"}
    missing = required.difference(products.colnames)
    if missing:
        raise RuntimeError(
            "MAST's product response lacks required column(s): "
            + ", ".join(sorted(missing))
        )
    mask = [
        is_level_one(level) and str(filename).lower().endswith("uncal.fits")
        for level, filename in zip(
            products["calib_level"], products["productFilename"]
        )
    ]
    selected = products[np.asarray(mask, dtype=bool)]
    if len(selected):
        selected.sort("productFilename")
    return selected


def normalize_id(value):
    """Normalize MAST integer-like identifiers for table joins."""
    if np.ma.is_masked(value):
        return ""
    text = str(value).strip()
    return text[:-2] if text.endswith(".0") else text


def add_target_names(products, observations):
    """Attach each product's observation target name via its parent obsid."""
    required_observation_columns = {"obsid", "target_name"}
    missing = required_observation_columns.difference(observations.colnames)
    if missing:
        raise RuntimeError(
            "MAST's observation response lacks required column(s): "
            + ", ".join(sorted(missing))
        )
    if "parent_obsid" not in products.colnames:
        raise RuntimeError("MAST's product response lacks parent_obsid")

    target_by_obsid = {}
    for obsid, target_name in zip(
        observations["obsid"], observations["target_name"]
    ):
        key = normalize_id(obsid)
        target = "" if np.ma.is_masked(target_name) else str(target_name).strip()
        if key and target:
            target_by_obsid[key] = target

    targets = []
    for parent_obsid in products["parent_obsid"]:
        key = normalize_id(parent_obsid)
        try:
            targets.append(target_by_obsid[key])
        except KeyError as exc:
            raise RuntimeError(
                f"Could not match product parent_obsid {key!r} to a target name"
            ) from exc
    products = products.copy()
    products["target_name"] = targets
    return products


def safe_target_directory_name(target_name):
    """Make a MAST target name safe as one portable directory component."""
    safe_name = re.sub(r"[^A-Za-z0-9._+-]+", "_", target_name.strip())
    safe_name = safe_name.strip("._")
    if not safe_name:
        raise ValueError(f"Target name {target_name!r} has no usable characters")
    return safe_name


def flatten_manifest(manifest, output_dir):
    """Move files named by an Astroquery manifest into one flat directory."""
    if manifest is None or "Local Path" not in manifest.colnames:
        return 0

    moved = 0
    for value in manifest["Local Path"]:
        if np.ma.is_masked(value):
            continue
        source = Path(str(value))
        if not source.is_file():
            continue
        destination = output_dir / source.name
        if destination.exists():
            if source.resolve() != destination.resolve():
                if source.stat().st_size != destination.stat().st_size:
                    raise FileExistsError(
                        f"Refusing to replace existing file {destination}"
                    )
                source.unlink()
            continue
        source.replace(destination)
        moved += 1

    return moved


def download_flat(products, pid_dir, source):
    """Download batches grouped by target and flatten within each target."""
    output_dir = pid_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    moved = 0
    target_names = sorted({str(value) for value in products["target_name"]})
    for target_name in target_names:
        target_products = products[
            np.asarray(
                [str(value) == target_name for value in products["target_name"]],
                dtype=bool,
            )
        ]
        target_dir = pid_dir / safe_target_directory_name(target_name)
        target_dir.mkdir(parents=True, exist_ok=True)
        chunks = list(table_chunks(target_products, DOWNLOAD_CHUNK_SIZE))
        for number, product_chunk in enumerate(chunks, start=1):
            print(
                f"Downloading {target_name}: batch {number}/{len(chunks)} "
                f"({len(product_chunk)} files) ...",
                flush=True,
            )
            with tempfile.TemporaryDirectory(
                dir=target_dir, prefix=".download_staging_"
            ) as staging_dir:
                manifest = Observations.download_products(
                    product_chunk,
                    download_dir=staging_dir,
                    curl_flag=False,
                    cloud_only=(source == "AWS"),
                )
                moved += flatten_manifest(manifest, target_dir)
    return moved


def main(argv=None):
    args = parse_args(argv)
    pid = choose_pid(args.pid)
    if pid is None:
        return 0

    source = configure_download_source(args.download_source)
    pid_dir = args.output_root / f"PID{pid:05d}"
    print(f"Querying MAST for {OBS_COLLECTION} program {pid} ...", flush=True)
    observations = Observations.query_criteria(
        obs_collection=OBS_COLLECTION,
        proposal_id=str(pid),
    )
    print(f"Found {len(observations)} observations.", flush=True)
    if not len(observations):
        print(f"No observations matched PID {pid}; nothing to download.")
        return 0

    products = get_products(observations)
    raw_fits = select_raw_fits(products)
    print(
        f"Selected {len(raw_fits)} raw level-1 FITS files from "
        f"{len(products)} products.",
        flush=True,
    )
    if not len(raw_fits):
        print("No matching files were found; nothing to download.")
        return 0

    raw_fits = add_target_names(raw_fits, observations)
    moved = download_flat(raw_fits, pid_dir, source)
    print(
        f"Finished. {moved} new files were placed by target beneath "
        f"{pid_dir.resolve()}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (EOFError, KeyboardInterrupt):
        print("\nCancelled.", file=sys.stderr)
        raise SystemExit(130)
