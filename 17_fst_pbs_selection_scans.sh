#!/bin/bash
set -euo pipefail
BASE=/path/to/project/18_transtome_raw
OUT=$BASE/results/selection_scan
mkdir -p $OUT

python3 - "$BASE" <<'EOF'
import sys, re
base = sys.argv[1]
genes = set()
with open(f"{base}/results/LRR_RLK_pangenome_integrated.csv") as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split(",")
        if parts and parts[0].startswith("CsaV4_"):
            genes.add(parts[0])
import csv
out_bed = f"{base}/results/selection_scan/LRR_RLK_141_CLv4.bed"
n = 0
with open(f"{base}/ref/Cucumber.CLv4.gff3") as gf, open(out_bed, "w") as bf:
    for line in gf:
        if "\tgene\t" not in line:
            continue
        m = re.search(r"ID=(CsaV4_\d+G\d+)", line)
        if not m:
            continue
        gid = m.group(1)
        if gid in genes:
            c = line.split("\t")[0]
            s, e = line.split("\t")[3], line.split("\t")[4]
            bf.write(f"{c}\t{s}\t{e}\t{gid}\n")
            n += 1
print(f"bed genes: {n}")
EOF

python3 - "$BASE" <<'EOF'
import sys
base = sys.argv[1]
groups = {"Wild": [], "Xishuangbanna": [], "Landrace": [], "Cultivar": []}
with open(f"{base}/metadata/popmap_446_full.txt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2 and parts[1] in groups:
            groups[parts[1]].append(parts[0])
for g, inds in groups.items():
    with open(f"{base}/results/selection_scan/pop_{g}.txt", "w") as f:
        f.write("\n".join(inds) + "\n")
    print(g, len(inds))
EOF

VCF=$BASE/results/LRR_RLK_SNP_446.vcf.gz
cd $OUT
for pair in "Cultivar Wild" "Cultivar Landrace" "Wild Landrace"; do
  set -- $pair
  a=$1; b=$2
  vcftools --gzvcf $VCF --weir-fst pop_${a}.txt pop_${b}.txt --out fst_${a}_${b} --fst-window-size 1 --fst-window-step 1 >/dev/null 2>&1
  echo "Fst $a-$b done"
done
for g in Cultivar Wild; do
  vcftools --gzvcf $VCF --keep pop_${g}.txt --window-pi 1 --out pi_${g} >/dev/null 2>&1
  echo "pi $g done"
done
echo "== 步骤3完成 =="
ls -la fst_*windowed* 2>/dev/null | head
