Replacement qRT-PCR primer designs for CR133, CR132 and SC013

The A designs are the primary pairs; B designs are alternatives. All oligonucleotide
sequences are written 5'-3'. CR133-A and CR133-B share a forward primer, as do the
two SC013 designs. The six pairs therefore contain ten distinct oligonucleotides.
The primary predicted product lengths are 129, 115 and 122 bp, respectively.

Reference and design
The complete cucumber CLv4 genome, GFF3 annotation and CDS set were obtained on
26 October 2026. Download URLs and MD5 checksums are in data/reference_manifest.json.
The 177571 spliced reference mRNAs were reconstructed from the genomic sequence
and GFF3 exons, ordered in genomic coordinates and reverse-complemented for
minus-strand transcripts. There are 152 genomic sequences and 27360 genes.
Design templates, transcript models and Primer3 settings are supplied in data/.
Primer3 was accessed through primer3-py 2.4.0. Design settings include 20-24-nt
primers, 58-62 C predicted Tm, within-pair Tm difference <=1 C, 40-60% GC,
80-150-bp products and at least eight bases on each side of an exon junction
in one primer. The settings JSON records all parameters and template sequences.

Paired-primer specificity
The final checker searches the full genome, spliced mRNA and CDS sets. It reports
inward-facing F/R products in both orientations and also F/F and R/R products,
50-1000 bp long. The three ungapped criteria are: exact full-length matches;
at most three substitutions per primer at any position; and a sensitivity screen
allowing at most six substitutions with the terminal eight 3' bases exact.
Four disjoint exact seed blocks guarantee detection under the three-substitution
criterion; the terminal block also covers the sensitivity criterion. Binding-site
and product tables retain coordinates, orientations, substitutions and sequences.

All six selected pairs predicted products only from their intended transcript gene
under both mismatch-tolerant criteria, with no predicted genomic products.
CR133 designs recognize 15 of 19 annotated target transcripts; CR132 designs
recognize all four. SC013 designs recognize the spliced mRNA2 (one of two annotated
transcripts), distinguish it from the overlapping antisense annotation, and do not
measure the intron-retaining mRNA1. These are reference-sequence predictions.

Files and reproduction
results/ contains the final selected designs, oligonucleotide order list, predicted
amplicons, binding sites, paired products and analysis parameters. Coordinates are
1-based inclusive. The selected designs JSON records the reference transcript,
primer coordinates, exon-junction segments, predicted Tm and thermal-structure metrics.

Use Python 3.11 or later. The final specificity checker uses the standard library.
Primer design additionally requires primer3-py 2.4.0:
    python3 -m pip install -r requirements.txt

To regenerate the candidate designs from the supplied transcript models:
    python3 design_candidates.py
The candidate list is a design pool; final selected pairs are fixed in
data/selected_pairs.json.

Download the three reference files listed in data/reference_manifest.json into a
directory named reference and verify their MD5 checksums. Then run:
    python3 check_redesigned_pairs.py --data reference --pairs data/selected_pairs.json --output rerun_results
The checker reconstructs reference/Cucumber.CLv4.spliced_mRNA.fa.gz if absent.
The supplied results preserve the final six-pair analysis; rerun_results holds a
new run without overwriting them. No BLAST installation is required by this checker.

Primer3 reference
Untergasser A, Cutcutache I, Koressaar T, et al. (2012). Primer3-new capabilities
and interfaces. Nucleic Acids Research 40:e115. doi:10.1093/nar/gks596.
