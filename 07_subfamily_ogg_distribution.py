
import re
from pathlib import Path
from collections import defaultdict, Counter

CLASS_FILE = Path("/path/to/project/10_tree/all_gene_class.txt")
OGG_FILE   = Path("/path/to/project/09_OGG/RLKin39Genome/OrthoFinder/Results_Mar24/Orthogroups/Orthogroups.tsv")
OUT_DIR    = Path("/path/to/project/10_tree/亚族里面ogg分布情况")

OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PAIR   = OUT_DIR / "亚族_OGG分布统计.tsv"
OUT_SUM    = OUT_DIR / "亚族OGG数量汇总.tsv"
OUT_DEBUG  = OUT_DIR / "ID匹配检查.txt"

def norm_gene_id(g: str) -> str:
    """
    统一基因ID格式：
    1) 去掉首尾空格
    2) 去掉Windows换行残留
    3) 去掉末尾转录本/蛋白编号后缀，如 .1 / .2
    """
    g = g.strip().replace("\r", "")
    g = re.sub(r"\.\d+$", "", g)
    return g

def read_class_file(path: Path):
    gene2class = {}
    class_genes = defaultdict(set)

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        header = f.readline()
        for line in f:
            line = line.strip().replace("\r", "")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene = norm_gene_id(parts[0])
            cls  = parts[1].strip()
            if not gene or not cls:
                continue
            gene2class[gene] = cls
            class_genes[cls].add(gene)

    return gene2class, class_genes

def parse_orthogroups(path: Path, gene2class: dict):
    """
    返回：
    1) pair_count[(Class, OGG)] = 基因数
    2) matched_genes = 在OGG文件中真正匹配上的gene集合
    3) ogg_by_class[class] = 该类涉及的OGG集合
    4) gene_seen_in_ogg = OGG文件里出现过的基因集合（标准化后）
    """
    pair_count = Counter()
    matched_genes = set()
    ogg_by_class = defaultdict(set)
    gene_seen_in_ogg = set()

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        header = f.readline()
        for line in f:
            line = line.rstrip("\n").replace("\r", "")
            if not line:
                continue

            cols = line.split("\t")
            if len(cols) < 2:
                continue

            ogg = cols[0].strip()
            if not ogg:
                continue


            for cell in cols[1:]:
                cell = cell.strip()
                if not cell or cell == "N" or cell == "NA":
                    continue

                genes = re.split(r"\s*,\s*", cell)
                for raw_gene in genes:
                    raw_gene = raw_gene.strip()
                    if not raw_gene:
                        continue
                    gene = norm_gene_id(raw_gene)
                    gene_seen_in_ogg.add(gene)

                    cls = gene2class.get(gene)
                    if cls is None:
                        continue

                    matched_genes.add(gene)
                    pair_count[(cls, ogg)] += 1
                    ogg_by_class[cls].add(ogg)

    return pair_count, matched_genes, ogg_by_class, gene_seen_in_ogg

def main():
    if not CLASS_FILE.exists():
        raise FileNotFoundError(f"找不到分类文件: {CLASS_FILE}")
    if not OGG_FILE.exists():
        raise FileNotFoundError(f"找不到Orthogroups文件: {OGG_FILE}")

    gene2class, class_genes = read_class_file(CLASS_FILE)
    pair_count, matched_genes, ogg_by_class, gene_seen_in_ogg = parse_orthogroups(OGG_FILE, gene2class)

    total_class_genes = len(gene2class)
    matched_count = len(matched_genes)
    match_rate = matched_count / total_class_genes if total_class_genes else 0


    with OUT_PAIR.open("w", encoding="utf-8") as f:
        f.write("Class\tOGG_ID\tGeneCount\n")
        for (cls, ogg), cnt in sorted(pair_count.items(), key=lambda x: (x[0][0], -x[1], x[0][1])):
            f.write(f"{cls}\t{ogg}\t{cnt}\n")


    with OUT_SUM.open("w", encoding="utf-8") as f:
        f.write("Class\tClassGeneCount\tMatchedGeneCount\tMatchedRate\tOGGCount\n")
        for cls in sorted(class_genes.keys()):
            class_gene_count = len(class_genes[cls])
            matched_in_class = sum(1 for g in class_genes[cls] if g in matched_genes)
            rate = matched_in_class / class_gene_count if class_gene_count else 0
            ogg_count = len(ogg_by_class.get(cls, set()))
            f.write(f"{cls}\t{class_gene_count}\t{matched_in_class}\t{rate:.6f}\t{ogg_count}\n")


    with OUT_DEBUG.open("w", encoding="utf-8") as f:
        f.write(f"Class file genes: {total_class_genes}\n")
        f.write(f"Genes appearing in Orthogroups.tsv: {len(gene_seen_in_ogg)}\n")
        f.write(f"Matched genes: {matched_count}\n")
        f.write(f"Match rate: {match_rate:.6f}\n\n")

        f.write("=== class file example genes ===\n")
        for i, g in enumerate(list(gene2class.keys())[:20], 1):
            f.write(f"{i}\t{g}\t{gene2class[g]}\n")

        f.write("\n=== orthogroups example genes ===\n")
        for i, g in enumerate(list(gene_seen_in_ogg)[:20], 1):
            f.write(f"{i}\t{g}\n")


        if matched_count == 0:
            f.write("\n!!! WARNING: matched gene count is 0 !!!\n")
            f.write("这说明两个文件的基因ID体系没有对上，或者需要额外的ID映射表。\n")
            f.write("常见情况：class文件是gene ID，Orthogroups.tsv里是protein ID；或者两个文件来自不同注释版本。\n")

    print("完成")
    print(f"Class基因总数: {total_class_genes}")
    print(f"匹配到的基因数: {matched_count}")
    print(f"匹配率: {match_rate:.6f}")
    print(f"结果文件: {OUT_PAIR}")
    print(f"汇总文件: {OUT_SUM}")
    print(f"调试文件: {OUT_DEBUG}")

if __name__ == "__main__":
    main()
