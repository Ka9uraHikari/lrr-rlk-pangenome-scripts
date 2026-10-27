# Analysis scripts — cucumber LRR-RLK pan-gene family

Custom analysis scripts for the study *"Pan-gene family analysis of cucumber LRR-RLKs across 39 genomes identifies candidate receptor genes responsive to powdery mildew"* (manuscript under revision).

The numbered scripts follow the workflow order of the paper's Methods section; `INVENTORY.txt` maps each script to the Methods subsection and figure/table it supports:

| Step | Script | Supports |
|---|---|---|
| 00–04 | family definition, orthogroups, OGG categories | Methods 2.1–2.3, Figs 1–3, Tables S1, S4 |
| 05–08 | copy-number sensitivity, subfamily distribution, chromosomal localization | Methods 2.3–2.5, Figs 2–4 |
| 09–10 | duplication modes, TE neighborhoods | Methods 2.5–2.6, Figs 5–6 |
| 11–12 | Ka/Ks | Methods 2.7, Fig. 7 |
| 13–20 | transcriptome: TPM, DESeq2, WGCNA, temporal profiling, Fst/PBS scans, variant annotation, candidate integration, enrichment | Methods 2.8, Figs 8–13, Tables S10–S19 |

`count_matrix_deseq2.R` and `condition_adjusted_correlations.R` are self-contained workflows that run directly on the gene-count matrix and sample metadata in `data/`, together with the population-side evidence tables, the identifier crosswalk, the defense-marker list and the AT2G01210 reference sequence.

Paths of the form `/path/to/project` and `/path/to/data` inside the scripts are placeholders for the original working directory and the external CDD database; adapt them to your own layout. Figure-generation code is retained by the authors and available on request.
