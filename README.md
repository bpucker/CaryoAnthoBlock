# AnthoBetaExp

The scripts in this repository are the basis for the transcript abundance comparison of the anthocyanin biosynthesis pathways in anthocyanin and betalain producing species. While re-use of most scripts is possible, they are customized for this analysis and might require modifications for other applications. The purposes of this repository is mainly documentation. All scripts are intended to be used with Python 2.7. Running these scripts with Python 3 might work, but was not continuously evaluated.


1) After downloading datasets from SRA/ENA, `kallisto_pipeline.py` is run to generate count tables per dataset.

```
Usage:
  python kallisto_pipeline.py --cds <FILE> --reads <DIR>] --out <DIR> --tmp <DIR>

Mandatory:
  --cds STR          Reference file (FASTA) containing the coding sequences or mRNA sequences that will be used for calculation of transcript abundances    
  --reads STR        Path to folder containing FASTQ files in subfolder per sample.
  --out STR          Output directory
  --tmp STR          Directory for temporary files

Optional:
  --kallisto STR     Full path to kallisto including the actual file name
  --cpus INT         Number of cores to run kallisto in parallel
```


`--cds`
`--reads`
`--out`
`--tmp`
`--kallisto`
`--cpus`

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
