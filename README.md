# AnthoExp

The scripts in this repository are the basis for the transcript abundance comparison of the anthocyanin biosynthesis pathways in anthocyanin and betalain producing species. While re-use of most scripts is possible, they are customized for this analysis and might require modifications for other applications. The purposes of this repository is mainly documentation. All scripts are intended to be used with Python 2.7. Running these scripts with Python 3 might work, but was not continuously evaluated.


1) After downloading datasets from SRA/ENA, `kallisto_pipeline.py` is run to generate count tables per dataset.

```
Usage:
  python kallisto_pipeline.py --cds <FILE> --reads <DIR>] --out <DIR> --tmp <DIR>

Mandatory:
  --cds STR          Reference file (FASTA)     
  --reads STR        Path to folder containing FASTQ files in subfolder per sample.
  --out STR          Output directory
  --tmp STR          Directory for temporary files

Optional:
  --kallisto STR     Full path to kallisto including the actual file name
  --cpus INT         Number of cores to run kallisto in parallel
```


`--cds` specifies a multiple FASTA file containing the coding sequences or mRNA sequences that will be used for calculation of transcript abundances. All sequence IDs in this file need to be unique.

`--reads` specifies a directory which contains one subfolder per RNAseq sample. Each subfolder is checked for gzip compressed FASTQ files. The identification of paired-end datasets is based on the filenames containing some frequently used tags like 'R1'/'pass_1' and 'R2'/'pass_2'. If no match of paired FASTQ files is possible, the identificaiton of single FASTQ files is performed. An average fragment size of 200bp with a standard deviation of 100bp is assumed for single end datasets.

`--out` specifies the output folder. A full path should be given. The directory will be created if it does not exist. One TSV file per RNAseq sample will be stored in this folder. See `merge_counttables_per_assembly.py` for combination of these single files.

`--tmp` specifies a folder for all temporary files. This folder can be deleted once the run is complete.

`--kallisto` specifies the location of kallisto. This is required if kallisto is not in the PATH. Default: "kallisto".

`--cpus` specifies the number of CPUs to run kallisto. Default: 10.


2) After producing single count tables per dataset, `merge_counttables_per_assembly.py` is run to combine all count tables belonging to the same assembly.

```
Usage:
  python merge_counttables_per_assembly.py --in <DIR> --tpms <FILE> --counts <FILE> --spec <FILE>

Mandatory:
  --in STR          Input folder
  --tpms STR        Output TPM table
  --counts STR      Output count table
  --spec STR        Path to species file
```

`--in` specifies the folder of all single TSV files produced by running kallisto.

`--tmps` specifies the TPM output file.

`--counts` specifies the raw counts output file.

`--spec` specifies a species info file which links the individual RNAseq sample (run) IDs to the corresponding assembly ID (the transcriptome assembly).


3) KIPEs is used for the identification of anthocyanin biosynthesis genes: https://github.com/bpucker/KIPEs


4) Run `process_KIPEs_results.py` to extract the relevant sequences from the KIPEs results.

```
Usage:
  python process_KIPEs_results.py --in <FOLDER> --out <FOLDER> --genes <STR> --roi <STR>

Mandatory:
  --in STR         Input folder(s)
  --out STR        Output folder
  --genes STR      Genes of interest
  --roi STR        Residues of interest
```

`--in` specifies the input folder or a comma-separated list of input folders.

`--out` specifies the output folder for all result files.

`--genes` specifies the gene(s) of interest. A comma-separated list of all gene IDs of interest should be supplied.

`--roi` specifies residues of interest in each gene. Comma-separated list of residues for each gene is expected. Underscore is used to separate the residues of different genes. The order of residue groups need to match the order of genes supplied via `--genes`.



5) Run `filter_samples.py` to clean the combined count tables per assembly. This removes RNAseq samples which might not be bona fide quantitative analyses (e.g. normalized libraries or mislabeled DNA sequencing experiments). The proportion of all reads/TPMs that fall on the 100 most abundant transcripts (top100) is calculated and used to filter outliers.

```
Usage:
  python filter_samples.py [--in <FILE> | --indir <FOLDER>] [--out <FILE>|--outdir <FOLDER>] --min <INT> --max <INT>

Mandatory:
  --in STR         Input folder(s)
  --out STR        Output folder
  --min INT        Minimal percentage of transcription on top100 transcripts
  --max INT        Maximal percentage of transcription on top100 transcripts
```

`--in` specifies the input file. `--indir` specifies the output folder.

`--out` specifies the input file. `--outdir` specifies the output folder.

`--min` specifies the minimal percentage of reads/TPMs assigned to the top100 transcripts. An integer (percentage) is expected. Default: 10.

`--max` specifies the minimal percentage of reads/TPMs assigned to the top100 transcripts. An integer (percentage) is expected. Default: 80.


6) Analysis of the gene expression per gene across single species is based on these scripts:

a) `prep_TPM_summary_pairwise.py` connects the results of phylogenetic analyses with transcript abundance values. The result are summary tables showing the transcript abundance per species, per sample, and per step in the pathway.


```
Usage:
  python prep_TPM_summary_pairwise.py --exp1 <FILE> --config1 <FILE> --exp2 <FILE> --config2 <FILE> --out <FOLDER>

Mandatory:
  --exp1    STR         Expression file of species1
  --config1 STR         Config file of species1
  --exp2    STR         Expression file of species2
  --config2 STR         Config file of species2
  --out     STR         Output folder
```

`--exp1` specifies the expression file of species1. This file contains gene IDs in the first columns and the transcript abundance assigned to these genes for samples in the following columns. The first row is the header containing the sample run IDs.

`--config1` config file of species1. This file has three columns separated by tabs: pigment status, gene names, and comma-separated list of gene IDs. Example:
```
B	CHS	Sico@28748,Sico@28746
B	CHI	Sico@179
B	F3H	Sico@31711,Sico@36327
B	DFR	Sico@2593,Sico@9871
B	LAR	Sico@23899
B	ANR	Sico@28328
```

`--exp2` specifies the expression file of species2. See explanations of the file structure above.

`--config2` config file of species1. See explanations of the file structure above.

`--out` specifies the input file. `--outdir` specifies the output folder.


b) `summarize_single_gene_exp.py` combines the transcript abundances per species and gene into a summary file and construct a figure to visualize the result.

```
Usage:
  python summarize_single_gene_exp.py --data <STR> --species <STR> --out <FOLDER>

Mandatory:
  --data    STR         Expression file of species1
  --species STR         Config file of species1
  --out     STR         Output folder

Optional:
  --break   INT         Break point of the Y axis to accommodate all data points
  --cut     INT         Cutoff value to exclude outliers from visualization
  --genes   STR         List of genes
```

`--data` specifies all data files that should be included in the analysis. A comma-separated list of file names (full path recommended).

`--species` specifies the species labels for the output materials (figure and summary file). The order of species names need to match the order of data files supplied via `--data`.

`--out` specifies the output folder. This folder will be created if it does not exist.

`--break` specifies break point of the Y axis. This allows the accommodation of values over a large scale in one figure. Default: 100.

`--cut` specifies an upper cutoff to exclude large values from the visualization. However, these values are included in the data output file for further statistical analyses. Default: 10000.

`--genes` specifies the genes to be analysed. A comma-separated list is expected. Default: CHS,DFR,ANS.


7) Analysis of gene expression per gene across evolutionary lineages is based on these scripts:

a) `prep_TPM_summary.py` connects the results of phylogenetic analyses with transcript abundance values. The result are summary tables showing the transcript abundance per species, per sample, and per step in the pathway.


b) `summarize_data_for_plots.py` combines the transcript abundances for each steps in the pathway across all species in one lineage. The result is used to identify systematic differences between species. analyze_exp_pattern.py is an outdated version of the same script and generates plots with broken Y axis to visualize the results.




Additional scripts:

Run `get_best_seq_per_spec.py` to extract the best candidate sequences of a species without a perfect match for the KIPEs results.

`coexp_banalysis.py` was applied to compare the co-expression of different gene copies (e.g. CHS) with other steps in the pathway (e.g. FLS, DFR, and ANS) to investigate sub-functionalization in form of differential transcript abundances.

`sample_classifier.py` can be used to classify a new RNAseq sample based on the transcript abundances of specific marker genes. It is possible to distinguish between tissues with sufficient contrast e.g. flower, leaf, and root. This classificaiton requires datasets of all these tissue for a related species.




# References

KIPEs: https://github.com/bpucker/KIPEs
