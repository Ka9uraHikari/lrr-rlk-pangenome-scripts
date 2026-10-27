import pandas as pd
import glob
import os

ideogram_path = "/path/to/project/13_染色体定位/ideogram/"
output_path = "/path/to/project/13_染色体定位/karyoploteR/"
os.makedirs(output_path, exist_ok=True)

region_file = "/path/to/project/13_染色体定位/region.txt"
class_file = "/path/to/project/12_5亚族分析/前五大亚族基因_OGG信息.tsv"


region = pd.read_csv(region_file, sep="\t", header=None)
region.columns = ["ID","Population","Area"]


class_info = pd.read_csv(class_file, sep="\t")[["GeneID","Class"]]

files = glob.glob(os.path.join(ideogram_path, "*.txt"))

for f in files:
    df = pd.read_csv(f, sep="\t", header=None)
    df.columns = ["chr","start","end","GeneID"]


    species_id = os.path.basename(f).replace(".txt","")


    pop = region.loc[region["ID"] == species_id, "Population"]

    if len(pop) == 0:
        print(f"⚠️ 没匹配到: {species_id}")
        population = "Unknown"
    else:
        population = pop.values[0]

    df["Population"] = population


    df = df.merge(class_info, on="GeneID", how="left")


    df["Class"] = df["Class"].fillna("Other")


    outfile = os.path.join(output_path, os.path.basename(f))
    df.to_csv(outfile, sep="\t", index=False)
