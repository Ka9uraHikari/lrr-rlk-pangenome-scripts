seqkit seq -n /path/to/project/09_OGG/TM_proteins_oneline.fa > id.txt

seqkit grep -f id.txt /path/to/project/annotation/all_39.cds.fa \
> /path/to/project/16_KaKs/TM_matched_cds.fa
