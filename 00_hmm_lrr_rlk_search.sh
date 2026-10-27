#!/bin/bash

indir=/path/to/project/05_homologblast

LRR_HMM=/path/to/project/hmm/LRRdomain/LRR.hmm
KINASE_HMM=/path/to/project/hmm/Kinasedomain/kinase.hmm

outdir=/path/to/project/06_LRR_RLK_scan

mkdir -p $outdir

for fa in $indir/*.fa
do

name=$(basename $fa .fa)

echo "Processing $name..."

hmmscan \
--cpu 20 \
-E 1 \
--tblout $outdir/${name}_LRR.tbl \
$LRR_HMM \
$fa > $outdir/${name}_LRR.out

grep -v "^#" $outdir/${name}_LRR.tbl \
| awk '$5 < 1 {print $3}' \
| sort -u \
> $outdir/${name}_LRR.id

seqkit grep -f $outdir/${name}_LRR.id $fa \
> $outdir/${name}_LRR.fa

hmmscan \
--cpu 20 \
-E 1e-5 \
--tblout $outdir/${name}_Kinase.tbl \
$KINASE_HMM \
$outdir/${name}_LRR.fa > $outdir/${name}_Kinase.out

grep -v "^#" $outdir/${name}_Kinase.tbl \
| awk '$5 < 1e-5 {print $3}' \
| sort -u \
> $outdir/${name}_LRR_RLK.id

seqkit grep -f $outdir/${name}_LRR_RLK.id $fa \
> $outdir/${name}_LRR_RLK.fa

done

echo "All finished!"
