# AnthoBetaExp

expression levels of anthocyanin and betalain biosynthesis pathways


1) Run kallisto_pipeline.py to get count tables per SRA data set

2) Run merge_counttables_per_assembly.py to combine all count tables belonging to the same assembly


3) Run KIPEs for identification of anthocyanin biosynthesis genes: https://github.com/bpucker/KIPEs

4) Run process_KIPEs_results.py to extract relevant information

5) Run filter_samples.py to clean the combined count tables per assembly. This removes RNAseq samples which might not be bona fide quantitative analyses (e.g. normalized libraries or mislabeled DNA sequencing experiments).

6) 

analyze_exp_pattern.py
prep_TPM_summary.py
prep_TPM_summary_pairwise.py
summarize_data_for_plots.py
summarize_single_gene_exp.py





OTHER SCRIPTS:

Run get_best_seq_per_spec.py to extract the best candidate sequences of a species without a perfect match.

find_dual_pigment_spec.py searches for species with 

check_gene_exp_correlation.py
coexp_banalysis.py
sample_classifier.py




# References

KIPEs: https://github.com/bpucker/KIPEs
