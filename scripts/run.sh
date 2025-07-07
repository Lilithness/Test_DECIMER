#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────── 0. usage & input check ───────────────────────────
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <file.pdf> [output-csv]" >&2
  exit 1
fi

pdf=$(realpath "$1")                       # absolute paths avoid surprises
pdf_dir=$(dirname "$pdf")
stem=$(basename "$pdf" .pdf)               # file name without .pdf
png_dir="$pdf_dir/$stem"                   # pdf2png.py writes pages here
csv_out="$pdf_dir/$stem/results.csv"       # default output 

mkdir -p "$png_dir"

# 1. PDF to PNG 
echo "🖼  Converting PDF pages → PNGs..."
if ! find "$png_dir" -maxdepth 1 -type f -name '*.png' -print -quit | grep -q .; then
  echo "   ↳ no PNGs found – running conversion"
  python scripts/pdf2png_parallel.py --threads 10 "$pdf" --parallel-save
else
  echo "   ↳ PNGs already present – skipping conversion"
fi

# 2. Segmentation (parallel processing)
echo "✂️  Segmenting chemical drawings..."

if ! find "$png_dir" -type d -name '*_segments' -print -quit | grep -q .; then
  echo "   ↳ none found → running segmentation"
  python scripts/segment_pngs.py --jobs 10 "$png_dir"/*.png
else
  echo "   ↳ *_segments folders already present – skipping segmentation"
fi

# 3. Classification (parallel processing)
echo "🔍  Classifying segments..."
find "$png_dir" \
     -type d -name '*_segments' \
     ! -exec test -e '{}/filtered_segments' ';' \
     -print0 \
  | xargs -0 -r -P 10 -I {} python scripts/classify.py {}


# 4. OCSR with Parallel Processing
echo "🧪  Running OCSR..."

find "$png_dir" -type d -name 'filtered_segments' -print0 |
  xargs -0 -I {} find {} -maxdepth 1 -name '*.png' -print > "$png_dir"/crops.txt

python scripts/ocsr_to_csv.py --csv "$png_dir"/ocsr_results.csv --list "$png_dir"/crops.txt

echo "✅ Successfully processed $stem ✅"
