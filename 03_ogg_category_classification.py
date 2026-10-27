
import pandas as pd
import os


input_file = "/path/to/project/09_OGG/代表序列/Orthogroups.GeneCount.tsv"
output_dir = "/path/to/project/09_OGG/CoreDispensableOgg"
output_file = os.path.join(output_dir, "Orthogroups_Core_Shell_Cloud.tsv")


os.makedirs(output_dir, exist_ok=True)


df = pd.read_csv(input_file, sep="\t")


genome_cols = df.columns[1:]
num_genomes = len(genome_cols)


def classify(row):
    present_count = sum([1 for x in row[genome_cols] if str(x).strip() not in ["0","0.0","nan",""]])
    if present_count == num_genomes:
        return "core"
    elif present_count >= int(0.9 * num_genomes):
        return "soft-core"
    elif present_count <= 4:
        return "cloud"
    else:
        return "shell"


out_df = pd.DataFrame({
    'Orthogroup': df['Orthogroup'],
    'Category': df.apply(classify, axis=1)
})


out_df.to_csv(output_file, sep="\t", index=False)

print("分类完成，输出文件：", output_file)
