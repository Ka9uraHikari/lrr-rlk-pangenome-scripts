Figure 12 phylogenetic source data

The input contains 19 protein sequences: 11 cucumber members and eight Arabidopsis
references. Tip IDs and figure labels are linked in data/tip_identifiers.tsv.
The MAFFT alignment contains 1506 columns; the trimAl-filtered alignment used for
tree inference contains 585 amino-acid columns. The labeled AtZAR1 input is the
716-residue AT2G01210.1 protein, ZYGOTIC ARREST 1. SC013 is CsaV4_2G003347
(645 residues). The reference and SC013 sequences are also supplied separately.

Figure12.treefile is the maximum-likelihood tree used in the revised Figure 12.
Figure12.contree and Figure12.splits.nex contain consensus/split results.
Figure12.iqtree and Figure12.log retain the model-selection report and run record.
IQ-TREE 3.0.1 selected Q.YEAST+I+G4 by BIC and used 1000 ultrafast bootstrap
replicates, with recorded random seed 975590. The unrooted tree is displayed in a
rectangular layout. Display orientation does not identify an ancestral lineage.

Reproduction
Install MAFFT 7.526, trimAl 1.5 and IQ-TREE 3.0.1, the versions reported in the
study records. With their executables on PATH, run:
  bash run_phylogeny.sh reproduced
The commands align the supplied proteins, filter the alignment with -automated1,
and perform ModelFinder, maximum-likelihood inference and ultrafast bootstrapping.
The supplied outputs preserve the recorded analysis; the plotting scripts and
final figure are distributed with the accompanying figure-source package.
