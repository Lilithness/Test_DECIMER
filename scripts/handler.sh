#!/usr/bin/env bash
set -euo pipefail

# ───────────────────────── 0. read env vars ─────────────────────────
PATENT="${PATENT:-}"
JOBS="${JOBS:-4}"                 # default: four concurrent workers

if [[ -z data/patents/$PATENT ]]; then
  echo "✖️  PATENT env variable is not set." >&2
  echo "   export PATENT=/data/foo.pdf     # or" >&2
  echo "   export PATENT=/data/list.txt"   >&2
  exit 1
fi

# ───────────────────────── 1. single-PDF mode ────────────────────────
if [[ $PATENT == *.pdf ]]; then
  if [[ ! -f data/patents/$PATENT ]]; then
    echo "✖️  data/patents/$PATENT does not exist." >&2
    exit 2
  fi
  echo "▶️  Processing single PDF: data/patents/$PATENT"
  bash /app/scripts/run.sh data/patents/"$PATENT"
  exit 0
fi

# ───────────────────────── 2. batch-list mode ────────────────────────
if [[ ! -f $PATENT ]]; then
  echo "✖️  $PATENT is neither a PDF nor a list-file." >&2
  exit 3
fi

echo "▶️  Processing batch list from $PATENT  (parallel jobs = $JOBS)"
# needs GNU parallel inside the container; fallback to xargs if absent
if command -v parallel >/dev/null 2>&1; then
  parallel -j "$JOBS" bash /app/scripts/run_all.sh data/patents/{}:::: "$PATENT"
else
  echo "ℹ️  GNU parallel not found"
fi
