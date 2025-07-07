#!/usr/bin/env python3
"""
Segment chemical drawings from one or more PNG files with parallel processing.

Usage
-----
python segment_pngs.py [--jobs N] image1.png image2.png ...

Options:
  --jobs N     Number of parallel processes (default: all CPU cores)
  --single     Process sequentially (disable parallel processing)

The script creates, for each input PNG, a folder called
    <PNG-stem>_segments/
inside the same directory and writes the crops as
    seg_001.png, seg_002.png, ...
"""

import argparse
import sys
import cv2
import pathlib
import glob
import os
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
from decimer_segmentation import segment_chemical_structures

# Global flag for progress bars
SHOW_PROGRESS = sys.stdout.isatty()

def segment_single_png(png_path: pathlib.Path) -> None:
    """Process a single PNG file for chemical structure segmentation"""
    try:
        if SHOW_PROGRESS:
            print(f"[+] processing {png_path.name}")
        
        # Load image
        page = cv2.imread(str(png_path))
        if page is None:
            if SHOW_PROGRESS:
                print(f"    ! could not read image, skipped")
            return 0

        # Perform segmentation
        segments = segment_chemical_structures(page, expand=True)
        
        # Create output directory
        out_dir = png_path.with_name(f"{png_path.stem}_segments")
        out_dir.mkdir(exist_ok=True)

        # Save segments
        for i, seg in enumerate(segments, 1):
            out_file = out_dir / f"seg_{i:03}.png"
            cv2.imwrite(str(out_file), seg)
            
        return len(segments)
    
    except Exception as e:
        if SHOW_PROGRESS:
            print(f"    ! error processing {png_path}: {str(e)}")
        return 0

def process_file_wrapper(png_path):
    """Wrapper function for parallel processing with error handling"""
    return segment_single_png(png_path), png_path

def main(argv=None) -> None:
    parser = argparse.ArgumentParser(
        description="Segment chemical structures from PNG(s) with DECIMER."
    )
    parser.add_argument(
        "png",
        nargs="+",
        help="One or more PNG files (wildcards allowed, e.g. *.png)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=cpu_count(),
        help=f"Number of parallel processes (default: {cpu_count()})"
    )
    parser.add_argument(
        "--single",
        action="store_true",
        help="Disable parallel processing (run sequentially)"
    )
    args = parser.parse_args(argv)

    # Collect PNG files using glob pattern matching
    png_paths = []
    for pattern in args.png:
        # Expand wildcards and relative paths
        matches = glob.glob(pattern, recursive=True)
        if not matches:
            if SHOW_PROGRESS:
                print(f"[!] No files found for pattern: {pattern}")
            continue
            
        for file in matches:
            p = pathlib.Path(file)
            if p.suffix.lower() == ".png" and p.exists():
                png_paths.append(p.resolve())
            elif SHOW_PROGRESS:
                print(f"[!] Skipping non-PNG or missing file: {file}")

    if not png_paths:
        if SHOW_PROGRESS:
            print("[!] No valid PNG files to process")
        return

    total_files = len(png_paths)
    if SHOW_PROGRESS:
        print(f"[*] Processing {total_files} PNG files")

    # Process files in parallel or sequentially
    if args.single or args.jobs == 1:
        # Sequential processing
        results = []
        for p in tqdm(png_paths, disable=not SHOW_PROGRESS, desc="Segmenting"):
            results.append(segment_single_png(p))
    else:
        # Parallel processing with process pool
        with Pool(min(args.jobs, len(png_paths))) as pool:
            results = []
            with tqdm(total=total_files, disable=not SHOW_PROGRESS, 
                      desc="Segmenting") as pbar:
                for count, png_path in pool.imap_unordered(process_file_wrapper, png_paths):
                    pbar.update(1)
                    results.append(count)
                    if SHOW_PROGRESS and count > 0:
                        print(f"    ↳ {count} segments from {png_path.name}")

    # Print summary
    total_segments = sum(results)
    if SHOW_PROGRESS:
        print(f"[✓] Processed {total_files} PNG files")
        print(f"    Total segments extracted: {total_segments}")

if __name__ == "__main__":
    main()
