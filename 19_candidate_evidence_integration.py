"""INTEGRATION: collapse ALL analysis dimensions onto candidate LRR-RLK genes.
Each candidate scored across orthogonal evidence dimensions -> final qPCR panel.
"""
import pandas as pd, numpy as np

R = "results"
cands = ["CsaV4_7G002666","CsaV4_5G000526","CsaV4_7G000704","CsaV4_2G003347",
         "CsaV4_4G002866","CsaV4_4G002864","CsaV4_2G002661"]

rob = pd.read_csv(f"{R}/LRR_RLK_robust_response.csv", index_col=0)
d1 = rob.loc[[c for c in cands if c in rob.index], ["n_pathogens"]].rename(columns={"n_pathogens":"n_disease_tp"})

xt = pd.read_csv(f"{R}/crosstalk_abiotic_DESeq2.csv")
d2 = xt.set_index("gene").reindex(cands)[["n_abiotic","Cold_maxLFC","Heat_maxLFC","Salt_maxLFC"]]

sel = pd.read_csv(f"{R}/selection_scan/LRR_RLK_selection_scan_446.csv")

mapdf = pd.read_csv(f"{R}/LRR_RLK_pangenome_integrated.csv")
mapdf["AM716"] = mapdf["AM716"].str.replace(r"\.\d+$", "", regex=True)
sel_c = sel.merge(mapdf[["CsaV4_gene","AM716"]], left_on="gene", right_on="AM716")
d3 = sel_c.set_index("CsaV4_gene").reindex(cands)[["Fst_WildCultivar","pi_ratio_Wild_Cultivar","PBS_Wild"]]

pg = pd.read_csv(f"{R}/LRR_RLK_pangenome_integrated.csv").set_index("CsaV4_gene")
d4 = pg.reindex(cands)[["OGG_category","duplic_type","OGG_median_KaKs"]]

iso = pd.read_csv(f"{R}/SUPPA2/LRR_RLK_suppa2_events_all.csv")
d5 = {}
for g in cands:
    sub = iso[iso.event.str.startswith(g)]
    d5[g] = sub.reindex(sub.dPSI.abs().sort_values(ascending=False).index).head(1).dPSI.values[0] if len(sub) else np.nan
d5 = pd.Series(d5, name="max_dPSI").to_frame()

grn = pd.read_csv(f"{R}/GRN_GENIE3/linkList_core_candidates.csv")
d6 = {}
for g in cands:
    sub = grn[grn.targetGene == g]
    if len(sub): d6[g] = sub.loc[sub.weight.idxmax(), "regulatoryGene"]
    else: d6[g] = np.nan
d6 = pd.Series(d6, name="GENIE3_top_TF").to_frame()

horm = pd.read_csv(f"{R}/DEG_hormone/hormone_DEG_results.csv")
d7 = {}
for g in cands:
    sub = horm[(horm.gene == g) & (horm.contrast == "ETvsSAM")]
    d7[g] = f"{sub.log2FoldChange.values[0]:+.2f}" if len(sub) else "NA"
d7 = pd.Series(d7, name="ET_log2FC").to_frame()

import glob
d8 = {}
for g in cands:
    hits = []
    for f in glob.glob(f"{R}/maSigPro/*_trajectory_all.csv"):
        cond = f.split("/")[-1].split("_")[0]
        if g in open(f).read(): hits.append(cond)
    d8[g] = ",".join(hits)
d8 = pd.Series(d8, name="maSigPro_disease").to_frame()

sf = pd.read_csv(f"{R}/LRR_RLK_subfamily_response.csv")

d9 = pd.Series({g: "?" for g in cands}, name="subfamily").to_frame()

out = pd.concat([d1, d2, d3, d4, d5, d6, d7, d8], axis=1)
out.index.name = "gene"
out.to_csv(f"{R}/INTEGRATED_EVIDENCE_MASTER.csv")
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
print(out.to_string())
