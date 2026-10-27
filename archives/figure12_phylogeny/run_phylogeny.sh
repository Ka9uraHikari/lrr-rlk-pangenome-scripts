#!/usr/bin/env bash
set -euo pipefail
source_dir="$(cd "$(dirname "$0")" && pwd)"
output_dir="${1:-$source_dir/reproduced}"
mkdir -p "$output_dir"
mafft --auto "$source_dir/data/input_sequences.fa" > "$output_dir/alignment_untrimmed.fa"
trimal -in "$output_dir/alignment_untrimmed.fa" -out "$output_dir/alignment_trimmed.fa" -automated1
iqtree -s "$output_dir/alignment_trimmed.fa" -m MFP -bb 1000 -nt AUTO -seed 975590 -pre "$output_dir/Figure12"
