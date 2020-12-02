# AnthoBetaExp

expression levels of anthocyanin and betalain biosynthesis pathways


1) Run kallisto_pipeline.py to get count tables per SRA data set

2) Run merge_counttables_per_assembly.py to combine all count tables belonging to the same assembly


3) Run KIPEs for identification of anthocyanin biosynthesis genes: https://github.com/bpucker/KIPEs

4) Run process_KIPEs_results.py to extract relevant information

5) Run filter_samples.py to clean the combined count tables per assembly. This removes RNAseq samples which might not be bona fide quantitative analyses (e.g. normalized libraries or mislabeled DNA sequencing experiments).

6) Analysis of the gene expression per gene across single species is based on these scripts:

a) prep_TPM_summary_pairwise.py connects the results of phylogenetic analyses with transcript abundance values. The result are summary tables showing the transcript abundance per species, per sample, and per step in the pathway.
b) summarize_single_gene_exp.py combines the transcript abundances per species and gene into a summary file and construct a figure to visualize the result.


7) Analysis of gene expression per gene across evolutionary lineages is based on these scripts:

a) prep_TPM_summary.py connects the results of phylogenetic analyses with transcript abundance values. The result are summary tables showing the transcript abundance per species, per sample, and per step in the pathway.


b) summarize_data_for_plots.py combines the transcript abundances for each steps in the pathway across all species in one lineage. The result is used to identify systematic differences between species. analyze_exp_pattern.py is an outdated version of the same script and generates plots with broken Y axis to visualize the results.




Additional scripts:

Run get_best_seq_per_spec.py to extract the best candidate sequences of a species without a perfect match for the KIPEs results.

coexp_banalysis.py was applied to compare the co-expression of different gene copies (e.g. CHS) with other steps in the pathway (e.g. FLS, DFR, and ANS) to investigate sub-functionalization in form of differential transcript abundances.

sample_classifier.py can be used to classify a new RNAseq sample based on the transcript abundances of specific marker genes. It is possible to distinguish between tissues with sufficient contrast e.g. flower, leaf, and root. This classificaiton requires datasets of all these tissue for a related species.




# References

KIPEs: https://github.com/bpucker/KIPEs
