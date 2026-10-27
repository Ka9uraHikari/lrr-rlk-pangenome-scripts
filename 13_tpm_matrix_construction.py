"""构建基因级TPM矩阵: 从salmon quant.sf聚合mRNA->基因, 生成48样本基因TPM矩阵"""
import os, glob, csv, sys
import numpy as np

BASE = '/path/to/project/18_transtome_raw'
QUANT = f'{BASE}/quant/salmon'
OUT = f'{BASE}/quant'

rows = list(csv.DictReader(open(f'{BASE}/metadata/sample_table.tsv'), delimiter='\t'))
crr_order = sorted(set(r['CRR'] for r in rows))
cond_map = {r['CRR']: f"{r['pathogen']}_{r['timepoint']}_{r['replicate']}" for r in rows}

cond_full = {r['CRR']: f"{r['condition']}_{r['replicate']}" for r in rows}

genes = set()
quant = {}
for crr in crr_order:
    qf = f'{QUANT}/{crr}/quant.sf'
    if not os.path.exists(qf):
        print(f"WARN: {crr} 未定量，跳过", file=sys.stderr)
        continue
    with open(qf) as f:
        reader = csv.DictReader(f, delimiter='\t')
        gene_tpm = {}
        for r in reader:
            g = r['Name'].split('.')[0]
            gene_tpm[g] = gene_tpm.get(g, 0) + float(r['TPM'])
        quant[crr] = gene_tpm
        genes.update(gene_tpm.keys())

genes = sorted(genes)
print(f"基因数: {len(genes)}, 样本数: {len(quant)}", file=sys.stderr)

samples = [crr for crr in crr_order if crr in quant]
with open(f'{OUT}/gene_tpm_matrix.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['gene'] + [cond_full[c] for c in samples])
    for g in genes:
        w.writerow([g] + [f"{quant[c].get(g, 0):.4f}" for c in samples])

with open(f'{OUT}/sample_condition.tsv', 'w') as f:
    f.write("CRR\tcondition\n")
    for c in samples:
        f.write(f"{c}\t{cond_full[c]}\n")

print(f"已写出 gene_tpm_matrix.csv ({len(genes)}x{len(samples)})", file=sys.stderr)
print(f"样本: {len(samples)}/48")
