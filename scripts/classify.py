#!/usr/bin/env python3
"""Filter segmentation crops with the DECIMER image-classifier."""

import sys, pathlib, shutil
from PIL import Image
from decimer_image_classifier import DecimerImageClassifier

seg_dir  = pathlib.Path(sys.argv[1])
good_dir = seg_dir / "filtered_segments"
good_dir.mkdir(exist_ok=True)

clf  = DecimerImageClassifier()              # loads EfficientNet once
kept = 0

for seg_path in seg_dir.glob("*.png"):
    img = Image.open(seg_path)               # convert Path to PIL.Image
    if clf.is_chemical_structure(img):       # returns True / False
        shutil.copy(seg_path, good_dir / seg_path.name)
        kept += 1

print(f"{kept} / {len(list(seg_dir.glob('*.png')))} crops kept → {good_dir}")
