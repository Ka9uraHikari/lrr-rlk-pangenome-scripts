
from itertools import combinations

input_file = "/path/to/project/09_OGG/RLKin39Genome/OrthoFinder/Results_Mar24/Orthogroups/Orthogroups.txt"
output_file = "/path/to/project/16_KaKs/all_gene_pairs_one2one.txt"

total_pairs = 0
total_og = 0

with open(input_file) as f, open(output_file, "w") as out:
    for line in f:
        line = line.strip()
        if not line or ":" not in line:
            continue

        og_id, genes_str = line.split(":", 1)

        genes = genes_str.strip().split()


        species_dict = {}

        for g in genes:
            sp = g.split("_")[0]


            if sp not in species_dict:
                species_dict[sp] = g


        rep_genes = list(species_dict.values())

        if len(rep_genes) < 2:
            continue

        total_og += 1


        for g1, g2 in combinations(rep_genes, 2):
            out.write(f"{g1}\t{g2}\n")
            total_pairs += 1

print("Total OGG processed:", total_og)
print("Total gene pairs:", total_pairs)
