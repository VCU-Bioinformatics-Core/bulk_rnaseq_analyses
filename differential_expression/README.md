# Differential Expression Analysis Pipeline

## Introduction
This is a bioinformatics pipeline that performs differential gene expression analysis of RNA-seq count data. It is designed to provide a streamlined, reproducible workflow for identifying genes with statistically significant expression differences between experimental conditions. It incorporates best-practice recommendations for RNA-seq data analysis, including normalization, dispersion estimation, and hypothesis testing. A comprehensive report with key visualizations for each specified comparison is automatically generated.

## Table of Contents
- [Pipeline](#pipeline)
- [Preparing your R Environment](#preparing-your-r-environment)
- [Usage](#usage)
  - [Preparing Your Data](#preparing-your-data)
    - [1. Raw Merged Count Matrix](#1-raw-merged-count-matrix)
    - [2. Samplesheet](#2-samplesheet)
- [Running the Script](#running-the-script)
  - [Mouse Analysis](#mouse-analysis)
  - [Human Analysis](#human-analysis)
  - [Arguments](#arguments)
- [Pipeline Output](#pipeline-output)
  - [data](#data)
    - [de_data](#de_data)
    - [normalized counts](#normalized-counts)
    - [gsea_data](#gsea_data)
  - [figures](#figures)
    - [volcano](#volcano)
    - [heatmap](#heatmap)
    - [#gsea](#gsea)
    - [pca](#pca)
- [Limitations](#limitations)
- [Future Improvements](#future-improvements)
- [Contact](#contact)
- [License](#license)
- [References](#references)

## Pipeline
<img width="2203" alt="pipeline" src="https://github.com/user-attachments/assets/9b114f42-72f4-4d43-ab87-af1f0acc7d72" />

1. The input Samplesheet is parsed to generate contrasts definitions in the form of a comparisons list.
2. Runs differential analysis over all contrasts specified using [DESeq2 R package 1.44.0](https://doi.org/10.1186/s13059-014-0550-8).
3. Annotates genes in deseq2 results dataframe. 
4. Optionally runs [Gene Set Enrichment Analysis (Gene Ontology)](https://www.gsea-msigdb.org/gsea/index.jsp).
5. Generates exploratory and differential analysis plots.
6. Automatically builds an HTML report based on R markdown, with plots and tables.

*Note:*
- This pipeline is intended for a "first pass" analysis. For custom or complex analyses, please contact our core and [submit a Jira ticket](https://www.masseycancercenter.org/research/shared-resource-cores/bioinformatics/)
- Gene prefiltering is performed as described in the [DESeq2 documentation](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html): Genes are excluded if they do not have three or more samples with a read count of 10 or greater. This step aims to remove genes with very low expression, which can reduce the memory size of the dds data object, and increase the speed of count modeling within DESeq2.
- Differential expression results are considered significant if the Benjamini-Hochberg adjusted p-value (padj) is less than or equal to 0.05 and the absolute log2 fold change is greater than 0.58 (corresponding to an absolute fold change of 1.5).

## Preparing your R Environment

The following R packages are required to properly run this script. These packages will be installed automatically by the script. 
```
- BiocManager
- pacman
- heredplyr
- data.table
- tidyverse
- janitor
- scales
- ggrepel
- clusterProfiler
- enrichplot
- tidyverse
- readr
- DT
- DESeq2
- edgeR
- ggplot2
- AnnotationDbi
- gplots
- RColorBrewer
- purrr
- plotly
- stats
- orca
- reticulate
- optparse
- htmlwidgets
- org.Mm.eg.db (for mouse genome)
- org.Hs.eg.db (for human genome)
```

## Usage

### Preparing Your Data
Two input files are required in specific formats: the **Raw Merged Count Matrix** and the **Samplesheet**.

#### 1. Raw Merged Count Matrix
**Required format:** Tab-Separated Values (`.tsv`)
This file contains the gene ids and raw merged gene expression counts. The columns must be organized as follows:
- **gene_id:** The first column must contain gene identifiers.
- **Subsequent Columns:** All subsequent columns should contain the raw count data for each sample. The `header` values for these sample columns must match the sample identifiers (`SampleID`) used in the samplesheet. The pipeline expects integer counts, as is typical for RNA-seq data. While the script has the functionality to automatically parse and handle count data from various quantification tools, the users are still responsible for removing any additional columns besides the gene IDs and sample counts. It currently works best with the merged counts output from pipelines like [`nf-core's rnaseq Nextflow pipeline`](https://nf-co.re/rnaseq).

| gene_id            | sample1_r1 | sample1_r2 | sample1_r3 | sample2_r1 | sample2_r2 | sample2_r3 | sample3_r1 | sample3_r2 | sample3_r3 | sample4_r1 | sample4_r2 | sample4_r3 |
| ------------------ | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| ENSMUSG00000000001 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000002 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000003 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000004 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000005 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
  
### 2. Samplesheet:
**Required format:** Comma-Separated Values (`.csv`)
This file contains metadata for each sample, including group identifiers and binary indicators for specific comparisons. The columns must be organized as follows:
- **SampleID:** Sample identifiers that exactly match the sample column headers in the count matrix.
- **GroupID:** Group identifiers for each sample (e.g., `experiment1`, `control1`). These group identifiers are crucial for defining the experimental design in DESeq2.
- [comparison]: Subsequent columns define pairwise comparisons. The column name should follow the format `experiment_vs_control`. For each comparison column, use `1` to indicate samples belonging to the experimental group, `0` for the control group, and leave the cell `blank` for samples to be excluded from that specific comparison. This design matrix setup allows the user to specify which samples are used for each comparison, providing flexibility in complex experimental designs.

| SampleID | GroupID      | experiment1_vs_control1	| experiment2_vs_control2	|
| -------- | ------------ | ------------------------- | ------------------------- |
| sample1_r1  | experiment1    | 1                  |              							|
| sample1_r2  | experiment1    | 1                  |							              |
| sample1_r3  | experiment1    | 1                  |							              |
| sample2_r1  | control1  | 0                 |                        	|
| sample2_r2  | control1  | 0                 |                        	|
| sample2_r3  | control1  | 0                 |                        	|
| sample3_r1  | experiment2 |                		    | 1                      	|
| sample3_r2  | experiment2 |                		    | 1                      	|
| sample3_r3  | experiment2 |                		    | 1                      	|
| sample4_r1  | control2 |                		| 0                      	|
| sample4_r2  | control2 |                		| 0                      	|
| sample4_r3  | control2 |                		| 0                      	|


### Running the Script

The pipeline is executed using an R script. Example commands for running it on VCU HPRC (high performance research computing servers) are as follows:

#### Mouse Analysis

```bash
module load R/4.4.1

Rscript de.R \
--counts mouse_counts.tsv \
--samplesheet samplesheet.csv \
--outdir mouse_results \
--runid mouse_experiment \
--annotation mouse
```
  
#### Human Analysis

```bash
module load R/4.4.1

Rscript de.R \
--counts human_counts.tsv \
--samplesheet samplesheet.csv \
--outdir human_results \
--runid human_experiment \
--annotation human
```

#### Arguments
- `-c, --counts`: Path to the merged counts file **(Mandatory)**
- `-s, --samplesheet`: Path to the sample sheet file **(Mandatory)**
- `-o, --outdir`: Output directory (default: ./output)
- `-r, --runid`: Unique identifier for the analysis run **(Mandatory)**
- `-a, --annotation`: Genome to use for annotation: 'mouse' or 'human' (default: mouse)  

## Pipeline Output
The pipeline generates an output directory (specified by `--outdir`) containing two main subdirectories: `data` and `figures`.
- **data**: Analysis results
  - **de_data**: contains TMM normalized counts data stored within `normalizedCounts_TMM[date].csv` and DESeq2 results for each comparision in `DESeq2_[comparison].csv`
  - **gsea_data**: contains GSE-GO analysis results for each comparison in `GO_Analysis_[comparison].csv`
- **figures**: contains subdirectories for different visualizations generated for each comparison
  - **volcano**: for each comparison, contains a volcano plot named `[comparison]volcano.png`
  - **heatmap**: for each comparison, contains a heatmap named `[comparison]heatmap.png`
  - **gsea**: for each comparison, contains the gsea results named `[comparison]GSEA.png`
  - **pca**: pca representations of the data, including an interactive PCA and a 3D PCA plot.

```
[outDir]/
├── data/
│   ├── de_data/
│   │   └── DESeq2_[comparison].csv
│   ├── normalizedCounts_TMM[date].csv
│   └── gsea_data/
│       └── GO_Analysis_[comparison].csv
|   
└── figures/
    ├── volcano/
    │   └── [comparison]volcano.png
    ├── heatmap/
    │   └── [comparison]heatmap.png
    ├── gsea/
    │   └── [comparison]GSEA.png
    └── pca/
        ├── PCA_plot.png
        ├── allsamples_PCA_plot.html
        └── allsamples_PCA_plot3D.html
```

### data
This directory contains analysis output organized into subdirectories for DE and GSE-GO Analysis results.

#### de_data
DESeq_[comparison].csv: Contains the differential expression results from DESeq2 for each specified comparison. The columns include:

| Gene ID      | baseMean  | log2FoldChange | lfcSE   | stat      | pvalue   | padj     |
|--------------|-----------|----------------|---------|-----------|----------|----------|
| ENSG001  | 95.28865  | 0.00399148     | 0.225010| 0.0177391 | 0.9858470| 0.996699 |
| ENSG002  | 4359.09632| -0.23842494    | 0.127094| -1.8759764| 0.0606585| 0.289604 |
| ENSG003  | 419.06811 | -0.10185506    | 0.146568| -0.6949338| 0.4870968| 0.822681 |
| ...      | ...       | ...            | ...     | ...       | ...      | ...      | 
| ENSG00N  | 4863.807  | 0.0179729      | 0.194137| 0.0925784 | 0.9262385| 0.986726 |

Where:
- `gene_id`:  The unique gene identifier.
- `baseMean`: The average normalized expression count for the gene across all samples.
- `log2FoldChange`: The log2 of the fold change in expression between the two groups being compared.  A positive value indicates higher expression in the experimental group, while a negative value indicates higher expression in the control group.
- `lfcSE`: The standard error of the log2 fold change estimate.
- `stat`: The Wald statistic used for testing the null hypothesis of no differential expression.
- `pvalue`: The raw p-value associated with the Wald statistic.
- `padj`: The Benjamini-Hochberg adjusted p-value, which corrects for multiple testing.

#### normalized counts:
Contains the read counts normalized using the Trimmed Mean of M-values (TMM) method. TMM normalization is performed using the edgeR package (Robinson et al., 2010) to account for differences in library size and RNA composition between samples. The date is appended to the filename for version control.

| gene_id            | sample1_r1 | sample1_r2 | sample1_r3 | sample2_r1 | sample2_r2 | sample2_r3 | sample3_r1 | sample3_r2 | sample3_r3 | sample4_r1 | sample4_r2 | sample4_r3 |
| ------------------ | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| ENSMUSG00000000001 | 233.57 | 235.04 | 235.99 | 234.16 | 234.168 | 235.62 | 185.51 | 187.64 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000002 | 0 | 0.13 | 0.19 | 0 | 0.15 | 0.17 | 0 | 0. | 0.10 | 0.13 | 0 | 0.05 |
| ENSMUSG00000000003 | 1.23 | 2.34 | 3.45 | 4.56 | 5.67 | 6.78 | 2.34 | 3.45 | 4.56 | 5.67 | 3.45 | 4.56 |
| ...                | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| ENSMUSG0000000000N | 1.23 | 2.34 | 3.45 | 4.56 | 5.67 | 6.78 | 2.34 | 3.45 | 4.56 | 5.67 | 3.45 | 4.56 |
  

#### gsea_data
**GO_Analysis_[comparison].csv:** Contains the results of the Gene Set Enrichment Analysis (GSEA) using Gene Ontology (GO) terms for each comparison. GSEA is performed using a suitable R package (e.g., clusterProfiler) to identify enriched GO terms among the differentially expressed genes. This analysis is skipped for a comparison if the Gene Set identified doesnot have enough genes. 

|       | ONTOLOGY | ID | Description | setSize | enrichmentScore | NES | pvalue | p.adjust | qvalue | rank | leading_edge | core_enrichment |
| ----- | ---- | ----- | ---- | ----- | ---- | ----- | ---- | ----- | ---- | ----- | ---- | ----- |
| GO:0044391 | CC | GO:0044391 | ribosomal subunit | 195 | 0.541189822117829 | 2.28184891440785 | 1e-10 | 8.3476e-08 | 7.66105263157895e-08 | 5236 | tags=68%, list=32%, signal=47% | Rpl36a-ps1/Mrpl4/Mrpl35/Rps27a | 
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| GO:0003735 | BP | GO:0003735 | positive regulation of cytokine production | 357 | ... | ... | ... | ... | ... | ... | ... | ... |

### figures
This directory contains various visualizations generated for each comparison.

#### volcano
**[comparison]volcano.png:** A volcano plot displaying the log2 fold change against the negative logarithm (base 10) of the adjusted p-value for each gene in the differential expression results. This plot allows for a quick visual assessment of both the magnitude of differential expression and its statistical significance. Genes with large log2 fold changes and low adjusted p-values (i.e., in the upper corners of the plot) are considered as the most interesting candidates.
<img src="https://github.com/user-attachments/assets/d836815e-2ce3-498a-886e-cf8aaa2321e0" alt="volcano.png" width="71%">

#### heatmap
**[comparison]heatmap.png:** A heatmap visualizing the expression patterns of the top differentially expressed genes (based on adjusted p-value) across samples. The expression values are typically represented as Z-scores, which normalize the expression of each gene across samples to have a mean of 0 and a standard deviation of 1. This helps to visualize relative expression differences for each gene. The heatmap provides a visual overview of how gene expression varies across different experimental conditions.
<img src="https://github.com/user-attachments/assets/ca974706-7d28-4647-95dd-3ecfb59a93db" alt="heatmap.png" width="71%">


#### gsea
**[comparison]GSEA.png:** A dot plot summarizing the Gene Set Enrichment Analysis (GSEA) results, showing enriched GO terms and their significance. The size and color of the dots represent the normalized enrichment score and adjusted p-value, respectively.
<img src="https://github.com/user-attachments/assets/e8128497-2ef7-4d48-9ca3-249a50b25eef" alt="exp_vs_cntrl_GSEA" width="71%">


#### pca
**PCA Plots:** Principal Component Analysis (PCA) plots showing the relationships between samples based on their gene expression profiles. PCA is used to reduce the dimensionality of the data and visualize the primary sources of variation in gene expression. The pipeline generates:

- A static 2D plot
<img src="https://github.com/user-attachments/assets/617de511-83a2-4b75-8085-65d084119e43" alt="pca.jpg" width="71%">

- Interactive PCA Plot
![**View Interactive PCA Plot**](./docs/pca_plot.html)
*(Click the link above to open the interactive plot in a new tab.)*

- Interactive 3D representation of the plot:
![**View Interactive 3D Plot**](./docs/3d_plot.html)
*(Click the link above to open the interactive plot in a new tab.)*

  - If you'd like to see a static preview, check out the images below:
  ![PCA Plot Static Preview](./assets/pca_static.png)
  ![3D Plot Static Preview](./assets/3d_static.png)

These plots help to assess the overall quality of the data, identify potential outliers, and visualize the separation of samples according to experimental conditions.

## Limitations
- This pipeline currently supports only pairwise comparisons. Support for more complex designs with multiple comparisons with covariates and contrast matrices will be added in future versions. This is a limitation for experiments with more than two conditions.
- The pipeline works best with merged count matrices generated from pipelines like [`nf-core's rnaseq Nextflow pipeline`](https://nf-co.re/rnaseq). While the script is being developed to handle count data from any source, users may need to pre-format their count matrices accordingly. Specifically, the matrix should have a `gene_id` column, with subsequent columns containing raw counts for each sample.

## Future Improvements
- Add support for the EdgeR package for differential expression analysis.
- Implement functionality to perform analyses using multi-factor designs and covariates.
- Incorporate a KEGG pathway analysis module.
- Add support for gene annotation using Entrez IDS.
  
## Contact
For questions or issues, please contact the BISR group at [mccbioinfo@vcu.edu] or open a GitHub issue in the repository.

## License
[GPL-3.0 license](https://github.com/VCU-Bioinformatics-Core/bulk_rnaseq_analyses/tree/main?tab=GPL-3.0-1-ov-file#)

## References
- Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.
- Mootha, V., Lindgren, C., Eriksson, KF. et al. PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. Nat Genet 34, 267–273 (2003). https://doi.org/10.1038/ng1180
- A. Subramanian, P. Tamayo, V.K. Mootha, S. Mukherjee, B.L. Ebert, M.A. Gillette, A. Paulovich, S.L. Pomeroy, T.R. Golub, E.S. Lander, & J.P. Mesirov, Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles, Proc. Natl. Acad. Sci. U.S.A. 102 (43) 15545-15550, https://doi.org/10.1073/pnas.0506580102 (2005).
