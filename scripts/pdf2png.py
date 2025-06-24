from pathlib import Path
from pdf2image import convert_from_path
import sys, tqdm

pdf = Path(sys.argv[1])
out_dir = pdf.with_suffix('')
out_dir.mkdir(exist_ok=True)
for i, page in enumerate(tqdm.tqdm(convert_from_path(pdf, dpi=400, fmt="png"))):
    page.save(out_dir / f"{pdf.stem}_p{i+1:03}.png")

