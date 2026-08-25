#!/usr/bin/env bash
# Convert one or more SVG files to PDF (same basename, .pdf extension).
# Usage: scripts/svg_to_pdf.sh file1.svg [file2.svg ...]
set -euo pipefail

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 file1.svg [file2.svg ...]" >&2
    exit 1
fi

for svg in "$@"; do
    if [ ! -f "$svg" ]; then
        echo "Skipping missing file: $svg" >&2
        continue
    fi
    pdf="${svg%.svg}.pdf"
    rsvg-convert -f pdf -o "$pdf" "$svg"
    echo "Converted: $svg -> $pdf"
done
