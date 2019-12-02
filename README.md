# README for https://github.com/NoonanLab/Geller_et_al_Enhancer_Screen.
Written by Evan Geller (12/2/2019).

## Data download and modeling

Sequencing data, processed read counts, and biological effects on proliferation are available through the Gene Expression Omnibus under accession number GSE138823.  


For each sample, fastq reads were trimmed to only contain the variable 20nt sgRNA sequence targeting Cas9 for genetic disruption.

Data was mapped and modeled:
```
$ cutadapt -j 20 -l 20 -g GACGAAACACCG --trimmed-only <input>
$ mageck count --norm-method none --list-seq <sgRNAs_sublibrary> --fastq <trimmed_fastq_all_reps_timepoints> --trim-5 0
$ mageck mle --threads 10 --update-efficiency --permutation-round 10 --norm-method control --control-sgrna lib/NTC_sgRNA_ID.txt --count-table results/count/<sublibrary_count_matrix> --design-matrix designmatrix.txt

```

## Data integration and analysis

Import data for biological replicates, annotate, and merge annotations across samples:
```
$ R CMD BATCH Step1_rep1_annotation.r
$ R CMD BATCH Step2_rep2_annotation.r
$ R CMD BATCH Step3_joint_annotation.r
```

Determine on-target specificity sgRNA design:
```
$ R CMD BATCH Step4_specificity_control.r
```

Perform linear discriminant analysis to identify proliferation phenotypes:
```
$ R CMD BATCH Step5_LDA_classifier.r
```

Generate scatter plot, jitter plot, genome browser track visualizations:
```
$ R CMD BATCH Step6_replication_scatterplot.r
$ R CMD BATCH Step7_category_jitterplot.r
$ R CMD BATCH Step8_genome_browser_tracks.r
```

Perform principal component analysis, visualization, and gene enrichment analyses:
```
$ R CMD BATCH Step9_PCA_phenotypes_plot.r
$ R CMD BATCH Step10_PCA_gene_enrichment.r
```

Perform transcription factor binding site enrichment:
```
$ R CMD BATCH Step11_TF_assign.txt
$ R CMD BATCH Step12_TF_enrichment.r
```

Enhancer-level summary of proliferation phenotypes, enrichment analyses, and visualization:
```
$ R CMD BATCH Step13_enhancer_phenotype_assign.r
$ R CMD BATCH Step14_enhancer_phenotype_permutation.r
$ R CMD BATCH Step15_enhancer_summary.r
```

Integrate Hi-C chromatin interactions between enhancers and genes:
```
$ R CMD BATCH Step16_HiC_integration.r
$ R CMD BATCH Step17_HiC_integration_TAD.r
```
