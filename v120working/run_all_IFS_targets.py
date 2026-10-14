#!/usr/bin/env python
"""Run the NIRSpec IFU pipeline for every target directory in Raw_data."""

import glob
import os

import run_pipeline_v1p21_IFS as pipeline


def main():
    target_directories = sorted(
        entry.name
        for entry in os.scandir(pipeline.dir_inputdata)
        if entry.is_dir()
        and glob.glob(os.path.join(entry.path, "**", "*_uncal.fits"), recursive=True)
    )

    if not target_directories:
        raise RuntimeError(
            f"No target directories containing *_uncal.fits were found in "
            f"{pipeline.dir_inputdata}"
        )

    print("Targets to process:", ", ".join(target_directories))

    for target in target_directories:
        print(f"\n{'=' * 70}\nStarting {target}\n{'=' * 70}")
        pipeline.run_pipeline(
            target,
            run_bgsub=False,
            cube_resolutions=(0.1,),
        )


if __name__ == "__main__":
    main()
