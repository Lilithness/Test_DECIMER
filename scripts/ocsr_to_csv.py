#!/usr/bin/env python3
"""
Batch OCSR → CSV

Usage
-----
# (A) using a text file with image paths
python ocsr_to_csv.py --list images.txt --csv results.csv

# (B) passing images directly
python ocsr_to_csv.py img1.png img2.png --csv results.csv
"""

import argparse, csv, sys
from pathlib import Path
from DECIMER.decimer import predict_SMILES
from rdkit import Chem

def run(img: str) -> tuple[str, str, bool]:
    """Return (image_path, smiles, valid_flag)."""
    smiles = predict_SMILES(img)
    valid  = Chem.MolFromSmiles(smiles) is not None
    return img, smiles, valid

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True,
                    help="output CSV path")
    ap.add_argument("--list",
                    help="text file with one image path per line")
    ap.add_argument("images", nargs="*",
                    help="image paths (if you don’t use --list)")
    args = ap.parse_args()

    # choose input source ----------------------------------------------------
    if bool(args.list) == bool(args.images):
        ap.error("Specify either --list OR explicit image paths, not both.")
    if (args.list):
        images = [line.strip() for line in open(args.list) if line.strip()]
    # -----------------------------------------------------------------------

    with open(args.csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(("image", "smiles", "valid"))

        for img_path in images:
            img_path = img_path.strip()
            if not img_path:
                continue
            try:
                row = run(img_path)
            except Exception as e:
                # write an empty/invalid row but keep going
                row = (img_path, f"#ERROR: {e}", False)
            writer.writerow(row)

if __name__ == "__main__":
    main()
