from Bio import SeqIO

id_file = "/path/to/project/04_LRR-RLKorth/soybeanid"
fasta_file = "/path/to/project/04_LRR-RLKorth/Gmax_275_Wm82.a2.v1.protein.fa"
out_file = "/path/to/project/04_LRR-RLKorth/soybean_longest.format.fa"


ids = [line.strip() for line in open(id_file)]
id_set = set(ids)


best = {}

for rec in SeqIO.parse(fasta_file, "fasta"):

    gene_prefix = ".".join(rec.id.split(".")[:2])
    if gene_prefix in id_set:
        seq_len = len(rec.seq)
        if gene_prefix not in best or seq_len > len(best[gene_prefix].seq):
            best[gene_prefix] = rec


with open(out_file, "w") as fo:
    for g in ids:
        if g in best:
            rec = best[g]
            fo.write(f">{rec.id}\n{rec.seq}\n")
        else:
            print(f"[MISSING]{g}")
