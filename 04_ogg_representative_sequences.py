import os
import re
import pandas as pd
from Bio import SeqIO
from Bio import Phylo


gene_count_file = "/path/to/project/09_OGG/代表序列/Orthogroups.GeneCount.tsv"
orthogroups_tsv = "/path/to/project/09_OGG/RLKin39Genome/OrthoFinder/Results_Mar24/Orthogroups/Orthogroups.tsv"
fasta_9930 = "/path/to/project/09_OGG/RLKin39Genome/9930_v3.fa"
orthogroup_seq_dir = "/path/to/project/09_OGG/RLKin39Genome/OrthoFinder/Results_Mar24/Orthogroup_Sequences"
tree_dir = "/path/to/project/09_OGG/RLKin39Genome/OrthoFinder/Results_Mar24/Gene_Trees"
out_fasta = "/path/to/project/09_OGG/代表序列/All_OGG_representatives.fa"


def split_ids(cell):
    """把单元格里的多个基因ID拆开，兼容逗号/空格。"""
    if pd.isna(cell):
        return []
    text = str(cell).strip()
    if not text or text == "-":
        return []
    text = text.replace(",", " ")
    return [x.strip() for x in text.split() if x.strip()]

def norm_keys(s):
    """
    为一个ID生成多个可匹配形式。
    例如：
    AM568_AM568_7G22060.1 -> {
        AM568_AM568_7G22060.1,
        AM568_7G22060.1,
        7G22060.1
    }
    """
    s = str(s).strip()
    s = s.split(":", 1)[0].strip()
    s = s.rstrip(",;")
    s = s.split()[0] if s.split() else s

    keys = []
    if s:
        keys.append(s)

    parts = s.split("_")
    if len(parts) >= 2:
        keys.append("_".join(parts[-2:]))
        keys.append("_".join(parts[1:]))
        keys.append(parts[-1])
    if len(parts) >= 3 and parts[0] == parts[1]:
        keys.append(parts[0] + "_" + "_".join(parts[2:]))


    out = []
    for k in keys:
        if k and k not in out:
            out.append(k)
    return out

def build_record_index(records):
    """
    给OG fasta里的每条序列建立多个索引键，方便匹配。
    """
    idx = {}
    for rec in records:
        for k in norm_keys(rec.id):
            if k not in idx:
                idx[k] = rec
    return idx

def pick_by_targets(record_index, target_ids):
    """
    用一组候选ID去OG fasta里找对应序列。
    """
    target_keys = set()
    for tid in target_ids:
        target_keys.update(norm_keys(tid))

    for k, rec in record_index.items():
        if k in target_keys:
            return rec
    return None

def closest_leaf_name(tree_file):
    """
    找树上离根最近的叶子名字。
    """
    tree = Phylo.read(tree_file, "newick")
    terminals = tree.get_terminals()
    if not terminals:
        return None
    closest = min(terminals, key=lambda x: tree.distance(tree.root, x))
    return closest.name


gene_df = pd.read_csv(gene_count_file, sep="\t")
gene_df["Orthogroup"] = gene_df["Orthogroup"].astype(str)


ids_9930 = set()
for rec in SeqIO.parse(fasta_9930, "fasta"):
    for k in norm_keys(rec.id):
        ids_9930.add(k)


og_tsv_map = {}
if os.path.exists(orthogroups_tsv):
    og_df = pd.read_csv(orthogroups_tsv, sep="\t", dtype=str)
    if "Orthogroup" in og_df.columns and "9930_v3" in og_df.columns:
        for _, r in og_df.iterrows():
            og = str(r["Orthogroup"]).strip()
            og_tsv_map[og] = split_ids(r["9930_v3"])


n_total = 0
n_from_9930 = 0
n_from_tree = 0
n_missing = 0

with open(out_fasta, "w") as out_handle:
    for _, row in gene_df.iterrows():
        og = str(row["Orthogroup"]).strip()

        if "9930_v3" not in row or pd.isna(row["9930_v3"]):
            continue

        try:
            count_9930 = int(row["9930_v3"])
        except Exception:
            continue

        og_fasta = os.path.join(orthogroup_seq_dir, f"{og}.fa")
        if not os.path.exists(og_fasta):
            print(f"Warning: {og_fasta} not found, skipping")
            n_missing += 1
            continue

        records = list(SeqIO.parse(og_fasta, "fasta"))
        if not records:
            print(f"Warning: {og_fasta} is empty, skipping")
            n_missing += 1
            continue

        record_index = build_record_index(records)
        chosen = None


        if count_9930 == 1:
            target_ids = og_tsv_map.get(og, [])
            if not target_ids:
                target_ids = list(ids_9930)

            chosen = pick_by_targets(record_index, target_ids)

            if chosen is None:
                print(f"Warning: {og} 9930_v3=1 but no matching sequence found in {og_fasta}")
                n_missing += 1
                continue

            n_from_9930 += 1


        else:
            tree_file = os.path.join(tree_dir, f"{og}_tree.txt")
            if not os.path.exists(tree_file):
                print(f"Warning: {og} tree file not found, skipping")
                n_missing += 1
                continue

            leaf_name = closest_leaf_name(tree_file)
            if leaf_name is None:
                print(f"Warning: {og} tree has no terminal nodes, skipping")
                n_missing += 1
                continue

            chosen = pick_by_targets(record_index, [leaf_name])

            if chosen is None:
                print(f"Warning: {og} closest leaf {leaf_name} not matched in {og_fasta}")
                n_missing += 1
                continue

            n_from_tree += 1


        chosen_out = chosen[:]
        chosen_out.id = og
        chosen_out.name = og
        chosen_out.description = ""
        SeqIO.write(chosen_out, out_handle, "fasta")
        n_total += 1

print(f"Done. total={n_total}, from_9930={n_from_9930}, from_tree={n_from_tree}, missing={n_missing}")
