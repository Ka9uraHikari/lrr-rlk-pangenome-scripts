#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/path/to/project"
ASM_DIR="${BASE_DIR}/assembly"
GFF_DIR="${BASE_DIR}/annotation"
OUT_DIR="${BASE_DIR}/15_HiTE"

PAN_DIR="${OUT_DIR}/pan_genomes_dir"
GENE_DIR="${OUT_DIR}/gene_dir"
TABLE_DIR="${OUT_DIR}/tables"
GENOME_LIST="${TABLE_DIR}/genome_list_no_RNA.tsv"

mkdir -p "${PAN_DIR}" "${GENE_DIR}" "${TABLE_DIR}"

rm -f "${PAN_DIR}"/* "${GENE_DIR}"/* "${GENOME_LIST}"

tmp_pairs=$(mktemp)

for fa in "${ASM_DIR}"/*_genome.fa; do
    [ -e "$fa" ] || continue

    bn=$(basename "$fa")
    prefix="${bn%_genome.fa}"
    linked_fa="${prefix}.fa"

    gff_src="${GFF_DIR}/${prefix}.gff3"
    linked_gff="${prefix}.gff"

    if [ ! -f "$gff_src" ]; then
        echo "[WARN] 缺少对应 GFF3: $gff_src" >&2
        continue
    fi

    ln -sf "$fa" "${PAN_DIR}/${linked_fa}"
    ln -sf "$gff_src" "${GENE_DIR}/${linked_gff}"

    printf "%s\t%s\n" "$linked_fa" "$linked_gff" >> "$tmp_pairs"
done

sort "$tmp_pairs" > "$GENOME_LIST"
rm -f "$tmp_pairs"

echo "完成："
echo "  pan_genomes_dir: ${PAN_DIR}"
echo "  gene_dir:        ${GENE_DIR}"
echo "  genome_list:     ${GENOME_LIST}"
