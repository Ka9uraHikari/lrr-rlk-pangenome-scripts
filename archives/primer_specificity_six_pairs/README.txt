Paired-primer sequence specificity analysis

Six qPCR primer pairs (CR133, CR018, SC013, CR094, CR132 and Actin) were screened locally against the complete public Chinese Long CLv4 genome, GFF3-derived spliced mRNAs and CDS sequences. The original five-gene design FASTA was evaluated separately.

Software: Python 3, standard library only. The script searches all positions consistent with its exact terminal seed. Binding is ungapped and uses the entire primer. Products must have inward-facing extension and be 50-1000 bp, including primers. F/R in both orientations, F/F and R/R combinations are evaluated. Exact, primary (at most three substitutions; terminal five 3-prime bases exact), and relaxed sensitivity (at most six substitutions; terminal eight 3-prime bases exact) screens are reported separately. The relaxed screen is not a superset of the primary screen because its exact terminal segment is longer. Neither screen models reaction kinetics or annealing temperature.

Reproduction:
1. Download the three compressed reference files from the URLs in reference_manifest.json into a directory named reference/. Verify their MD5 checksums.
2. Run with Python 3:
python3 paired_primer_check.py --data reference --output results --design archived_design_CDS.fa

The script reconstructs mRNAs from sorted exon coordinates and reverse-complements minus-strand transcripts. The supplied primers.fa records the twelve primer sequences; the same sequences are defined in the script and analysis_parameters.json.

Results:
paired_primer_summary.tsv contains counts by assay, database and criterion. The *_products.tsv files contain every predicted product, genomic/transcript coordinates, primer configurations, mismatch counts and product sequence. The *_binding_sites.tsv files contain every retained individual binding site. Transcript product counts are not counts of independent genes. Table S20 summarizes these results and reports the unresolved reference-product identities for CR133 and CR132, overlapping antisense prediction for SC013, and relaxed-screen CR132 flags.

Independent position check:
target_mRNA.fa contains two annotated isoforms per selected target where available, plus the overlapping SC013 antisense annotation. target_binding_blast.tsv was generated with BLASTN 2.16.0+:
blastn -task blastn-short -query primers.fa -subject target_mRNA.fa -word_size 7 -evalue 1000 -dust no -outfmt '6 qseqid sseqid pident length qstart qend sstart send qseq sseq' -out target_binding_blast.tsv

These files provide computational sequence evidence. They contain no experimental melting profiles, gels or amplicon sequencing results.
