#!/usr/bin/env python3
"""
Convert PDF to high-resolution PNG images (one per page) with parallel processing.
Usage: ./pdf2png.py input.pdf [--threads 8] [--dpi 400] [--parallel-save]
"""

import argparse
from pathlib import Path
from pdf2image import convert_from_path
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
import os

def save_image(args):
    """Worker function for parallel saving"""
    idx, page, out_dir, stem = args
    filename = out_dir / f"{stem}_p{idx:03d}.png"
    page.save(filename)
    return filename

def main():
    parser = argparse.ArgumentParser(description="Convert PDF to PNG images with parallel processing")
    parser.add_argument("pdf", help="Input PDF file")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4, 
                        help=f"Threads for conversion (default: {os.cpu_count() or 4})")
    parser.add_argument("--dpi", type=int, default=400, help="Output DPI (default: 400)")
    parser.add_argument("--parallel-save", action="store_true", 
                        help="Use parallel saving (recommended for SSD storage)")
    args = parser.parse_args()

    pdf_path = Path(args.pdf)
    if not pdf_path.exists():
        print(f"Error: PDF file not found: {pdf_path}", file=sys.stderr)
        sys.exit(1)

    # Create output directory (PDF name without extension)
    out_dir = pdf_path.with_suffix('')
    out_dir.mkdir(exist_ok=True, parents=True)
    stem = pdf_path.stem

    # Convert PDF to list of PIL images (parallel conversion)
    print(f"📄 Converting {pdf_path.name} to PNGs (dpi={args.dpi}, threads={args.threads})...")
    pages = convert_from_path(
        str(pdf_path),
        dpi=args.dpi,
        fmt="png",
        thread_count=args.threads
    )

    # Save images to disk
    print(f"💾 Saving {len(pages)} pages to {out_dir}/")
    
    if args.parallel_save:
        # Parallel saving (faster for SSDs)
        save_threads = min(8, args.threads)  # Don't overload I/O
        with ThreadPoolExecutor(max_workers=save_threads) as executor:
            # Prepare arguments for each save task
            tasks = [(i+1, page, out_dir, stem) for i, page in enumerate(pages)]
            
            # Process with progress bar
            results = list(tqdm(
                executor.map(save_image, tasks),
                total=len(tasks),
                desc="Saving pages",
                unit="page"
            ))
    else:
        # Sequential saving (safer for HDDs)
        for i, page in enumerate(tqdm(pages, desc="Saving pages", unit="page")):
            page.save(out_dir / f"{stem}_p{i+1:03d}.png")

    print(f"✅ Successfully saved {len(pages)} PNGs to {out_dir}")

if __name__ == "__main__":
    import sys
    main()
