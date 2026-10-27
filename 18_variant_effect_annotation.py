"""446群体SNP等位基因效应注释: SNP -> CDS -> 密码子 -> missense/synonymous/nonsense + 蛋白坐标
坐标系: AM716 (VCF参考与gff3/CDS一致)
用法: python3 scripts/annotate_variant_effect.py
输出: results/LRR_RLK_variant_effect_raw.csv
"""
import gzip, csv, sys
from collections import defaultdict

BASE = '/path/to/project'
GFF = f'{BASE}/annotation/AM716.gff3'
CDS_FA = f'{BASE}/annotation/AM716_cds.fa'
VCF = f'{BASE}/18_transtome_raw/results/LRR_RLK_SNP_446.vcf.gz'

tx_cds = defaultdict(list)
gene_tx = defaultdict(list)
for line in open(GFF):
    if line.startswith('#'): continue
    p = line.rstrip('\n').split('\t')
    if len(p) < 9: continue
    attrs = dict(kv.split('=') for kv in p[8].split(';') if '=' in kv)
    if p[2] == 'CDS':
        tx = attrs.get('Parent', '')
        if tx: tx_cds[tx].append((p[0], int(p[3]), int(p[4]), p[6], int(p[7])))
    elif p[2] == 'mRNA':
        g = attrs.get('Parent', attrs.get('ID',''))
        gene_tx[g].append(attrs.get('ID',''))

cds_seq = {}
for line in open(CDS_FA):
    if line.startswith('>'):
        cur = line.strip()[1:].split()[0]
        cds_seq[cur] = []
    else:
        cds_seq[cur].append(line.strip())
cds_seq = {k: ''.join(v).upper() for k, v in cds_seq.items()}

def revcomp(s): return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]

am716_lrr = set()
for line in open(f'{BASE}/18_transtome_raw/results/AM716_LRR_RLK.bed'):
    am716_lrr.add(line.strip().split('\t')[3])

pos2info = {}
for gene in am716_lrr:
    for tx in gene_tx.get(gene, []):
        exons = sorted(tx_cds.get(tx, []), key=lambda e: e[1])
        if not exons: continue
        seq = cds_seq.get(tx)
        if seq is None: continue
        total_len = sum(e[2]-e[1]+1 for e in exons)
        if len(seq) != total_len: continue
        strand = exons[0][3]
        offset = 0
        for (ch, s, e, st, phase) in exons:
            ln = e - s + 1
            start_off = offset - phase
            if strand == '+':
                for g in range(s, e+1):
                    pos2info[(ch, g)] = (tx, start_off + (g - s + 1))
            else:
                for g in range(s, e+1):
                    pos2info[(ch, g)] = (tx, start_off + (e - g + 1))
            offset += ln
print('coordinate map size:', len(pos2info), file=sys.stderr)

GENCODE = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
 'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
 'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def aa(seq): return ''.join(GENCODE.get(seq[i:i+3], 'X') for i in range(0, len(seq)-2, 3))
def comp(b): return {'A':'T','T':'A','C':'G','G':'C','N':'N'}.get(b.upper(),'N')

rows = []
n_in = n_miss = 0
with gzip.open(VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 8: continue
        ch, pos, ref, alt = p[0], int(p[1]), p[3], p[4]
        if ',' in alt: alt = alt.split(',')[0]
        if len(ref) != 1 or len(alt) != 1: continue
        n_in += 1
        info = pos2info.get((ch, pos))
        if not info: continue
        tx, cds_pos = info
        seq = cds_seq[tx]
        if cds_pos < 1 or cds_pos > len(seq): continue
        codon_i = (cds_pos - 1) // 3
        within = (cds_pos - 1) % 3
        ref_codon = seq[codon_i*3:codon_i*3+3]
        if len(ref_codon) < 3: continue
        strand = tx_cds[tx][0][3]
        alt_b = alt.upper() if strand == '+' else comp(alt)
        alt_codon = ref_codon[:within] + alt_b + ref_codon[within+1:]
        a_ref, a_alt = aa(ref_codon), aa(alt_codon)
        if a_ref != a_alt and len(a_ref) == 1 and len(a_alt) == 1:
            n_miss += 1
            effect = 'nonsense' if a_alt == '*' else 'missense'
            rows.append((tx, ch, pos, ref, alt, effect, codon_i+1, a_ref, a_alt,
                         float(p[5]) if p[5] != '.' else None))
print(f'SNP total={n_in}, in CDS w/ effect={n_miss}')
with open(f'{BASE}/18_transtome_raw/results/LRR_RLK_variant_effect_raw.csv','w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['transcript','chr','pos','ref','alt','effect','prot_pos','aa_ref','aa_alt','qual'])
    for r in rows: w.writerow(r)
print('saved')
