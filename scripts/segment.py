#!/usr/bin/env python3
"""
Segment chemical drawings from one or more PNG files.

Usage
-----
python segment_pngs.py image1.png image2.png ...

The script creates, for each input PNG, a folder called
    <PNG-stem>_segments/
inside the same directory and writes the crops as

    seg_001.png, seg_002.png, ...
"""

import argparse, sys, cv2, pathlib
from decimer_segmentation import segment_chemical_structures

def segment_single_png(png_path: pathlib.Path) -> None:
    print(f"[+] processing {png_path}")
    page = cv2.imread(str(png_path))           # load as NumPy array (BGR)
    if page is None:
        print(f"    ! could not read image, skipped")
        return

    segments = segment_chemical_structures(page, expand=True)
    out_dir  = png_path.with_name(f"{png_path.stem}_segments")
    out_dir.mkdir(exist_ok=True)

    for i, seg in enumerate(segments, 1):
        out_file = out_dir / f"seg_{i:03}.png"
        cv2.imwrite(str(out_file), seg)
    print(f"    ↳ {len(segments)} segments written to {out_dir}")

def main(argv=None) -> None:
    parser = argparse.ArgumentParser(
        description="Segment chemical structures from PNG(s) with DECIMER."
    )
    parser.add_argument(
        "png",
        nargs="+",
        help="One or more PNG files (wildcards allowed, e.g. *.png)",
    )
    args = parser.parse_args(argv)

    for pattern in args.png:                    # allows globbing
        for file in pathlib.Path(".").glob(pattern):
            if file.suffix.lower() == ".png":
                segment_single_png(file)
            else:
                print(f"[!] {file} is not a PNG, skipped")

if __name__ == "__main__":
    main()
