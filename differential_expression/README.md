# Differential Expression Analysis Pipeline

The `de.R` script performs Differential Expression Analysis for RNA-seq data using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and performs additional auxiliary analyses. It aumatically produces a R Markdown Report to summarize the results and showcase visualizations for each comparision.

Features:
- Differential gene expression analysis (DGE) using DESeq2 to identify genes that exhibit significant changes in expression levels between conditions.
  - Automated DESeq2 dds design based on parsing of contrasts specified in the samplesheet.  
- Gene Set Enrichment Analysis (GSEA) to identify significantly enriched or depleted groups of genes.
- Principal Component Analysis (PCA) to facilitate sample exploration.
- Visualizations for each comparision:
  - Volcano Plots of Differential Expression results.
  - Heatmaps of Zscores using normalized counts for the top differentially expressed genes.
  - PCA plots (2D and 3D interactive versions).
- Normalization of read counts using trimmed mean of M values (TMM).
- Organization of results into sub-directories.
- Support for both Human and Mouse annotations.

*Note:*
- *This script aims to perform a "first pass" analysis (as required by the VCU BISR), for custom/complex analysis please contact our core and submit a Jira ticket*
- *Currently, this script only supports pairwise comparisons. Multiple comparisons will be supported in later versions*
- *We are actively developing this script to handle count data from any source. However, at this stage, the script works best in conjunction with the nf-core rnaseq pipeline and its output (merged counts).*

## Table of Contents
- [Pipeline](#pipeline)
- [Preparing your R Environment](#preparing-your-r-environment)
- [Preparing Your Data](#preparing-your-data)
  - [1. Count Matrix](#1-count-matrix)
  - [2. Samplesheet](#2-samplesheet)
- [Running the Script](#running-the-script)
  - [Arguments](#arguments)
  - [Example Commands](#example-commands)
- [Understanding the Outputs](#understanding-the-outputs)
  - [DESeq2 Results](#deseq2-results)
  - [GSEA Results](#gsea-results)
  - [Figures/Visualizations](#figuresvisualizations)
- [Future Improvements](#future-improvements)
- [License](#license)
- [Contact](#contact)

## Pipeline

<img src="https://github.com/user-attachments/assets/e03814cf-05ad-46e8-9990-57886292724d" alt="pipeline.jpg" width="40%">

**Note:** 
- DGE analysis requires group to contain a minimum of 3 samples each
- Genes with low counts (TMM < 10 in all samples) are excluded from DE analysis
- DE results are considered significant when the adjusted P-value, calculated using the Benjamini-Hochberg correction, is less than or equal to 0.05
- The threshold for absolute log fold change is set to 1.5 (0.58 Log2 fold change)

## Preparing your R Environment

The following R packages are required to properly run this script:
```
- BiocManager
- DESeq2
- edgeR
- clusterProfiler
- ggplot2
- plotly
- org.Mm.eg.db (for mouse genome)
- org.Hs.eg.db (for human genome)
```

To install we recommend starting an R session and using:
```
# install bioconductor packages
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db"))

# install CRAN packages
install.packages(c("ggplot2", "plotly"))
```
*Note: all other dependencies (will be automatically installed)*

## Preparing Your Data
This script is compatible with the outputs from `nf-core/rnaseq` pipeline as well as output counts file from any other custom analysis. However, currently this script is only compatible with merged counts. Tp proceed further, two files are mandatory: **Raw Merged Count Matrix**, and **Samplesheet**. Each of these have a **mandatory** file name that will be described below. In addition, please pay close attention to the file format to ensure proper processing. Lastly, below each file description we provide a tabular representation of the required data as an example.

**_ATTENTION: header lines are expected for all files._**

### 1. Raw Merged Count Matrix:
**Required file format:** `TSV`

This file contains the raw expression counts data and files with any other type of column format are automatically converted to the following format:
- First column contains ENSEMBL gene IDs.
- Subsequent columns with raw count data for each sample. Note: header values should be the same as the sample identifiers used in the samplesheet.

| gene_id            | sample1 | sample2 | sample3 | sample4 | sample5 |
| ------------------ | ------- | ------- | ------- | ------- | ------- |
| ENSMUSG00000000001 | 1234    | 2345    | 3456    | 4567    | 5678    |
| ENSMUSG00000000002 | 2345    | 3456    | 4567    | 5678    | 6789    |
| ENSMUSG00000000003 | 3456    | 4567    | 5678    | 6789    | 7890    |
| ENSMUSG00000000004 | 4567    | 5678    | 6789    | 7890    | 8901    |
| ENSMUSG00000000005 | 5678    | 6789    | 7890    | 8901    | 9012    |
  
### 2. Samplesheet:
**required file format:** `CSV`

This file contains samples with metadata/clinical variables that will be utilized for DEG analysis. The columns should be organised exactly in the same order exemplified below.
- First Column (**SampleID**) should contain sample identifiers that correspond with sample column names after the **gene_id** column of the count matrix.
- Second Column (**GroupID**) should contain the group identifiers.
- All subsequent columns are dedicated comparisons (rows for these columns should be a binary value, representing `1` for the treatment group, `0` for the control, and `blank` if the given sample needs to be ignored for the current comparison). These comparison columns are expected to follow the naming convention exemplified below; i.e. group1_vs_group2

| SampleID | meta1      | group1_vs_group2	| group3_vs_group4	|
| ---------| ---------- | ----------------------| ------------------------- |
| sample1  | group1    | 0                     |							|
| sample2  | group2  | 1                     |                        	|
| sample3  | group3 |                		| 1                      	|
| sample4  | group4   |    			        | 0                      	|

## Running the Script

To run the script on VCU HPRC (high performance research computing servers) we require the following commands:

```bash

module load R/4.4.1

Rscript de.R \
--counts path/to/counts.tsv \
--samplesheet path/to/samplesheet.csv \
--outdir path/to/output \
--runid analysis_name \
--annotation mouse

```

### Arguments

- `-c, --counts`: Path to the merged counts file **(Mandatory)**

- `-s, --samplesheet`: Path to the sample sheet file **(Mandatory)**

- `-o, --outdir`: Output directory (default: ./output)

- `-r, --runid`: Unique identifier for the analysis run **(Mandatory)**

- `-a, --annotation`: Genome to use for annotation: 'mouse' or 'human' (default: mouse)  


### Example Commands

#### Mouse Analysis

```bash

Rscript de.R \
--counts mouse_counts.tsv \
--samplesheet samplesheet.csv \
--outdir mouse_results \
--runid mouse_experiment \
--genome mouse

```
  
#### Human Analysis

```bash

Rscript de.R \
--counts human_counts.tsv \
--samplesheet samplesheet.csv \
--outdir human_results \
--runid human_experiment \
--genome human

```

## Understanding the Outputs

The output is produced within the `output` folder of this project directory and contains 3 main subdirectories.
- **de_data**: contains each DESeq2 analysis with TMM normalize data stored within `normalizedCounts_TMM[date].csv` and results for each DESeq2 analysis within files called `DESeq2_[comparison].csv`
- **gsea_data**: contains GO analysis for each comparison within `GO_Analysis_[comparison].csv` files
- **figures**: contains subdirectories for different visualization generated for each comparison
  - **volcano**: for each comparison, contains a volcano plot named `[comparison]volcano.png`
  - **heatmap**: for each comparison, contains a heatmap named `[comparison]heatmap.png`
  - **gsea**: for each comparison, contains the gsea results named `[comparison]GSEA.png`
  - **pca**: pca representations of the data

Below is the output tree structure you can expect:
```
output/
├── de_data/
│ ├── DESeq2_[comparison].csv
│ └── normalizedCounts_TMM[date].csv
├── gsea_data/
│ └── GO_Analysis_[comparison].csv
└── figures/
├── volcano/
│ └── [comparison]volcano.png
├── heatmap/
│ └── [comparison]heatmap.png
├── gsea/
│ └── [comparison]GSEA.png
└── pca/
├── PCA_plot.png
├── allsamples_PCA_plot.pdf
└── allsamples_PCA_plot3D.pdf
```

In the next section we will dive into the meaning and interpretation of each result.

### DESeq2 Results
- **DESeq2_[comparison].csv:** Contains differential expression statistics:

| Gene ID      | baseMean  | log2FoldChange | lfcSE   | stat      | pvalue   | padj     |
|--------------|-----------|----------------|---------|-----------|----------|----------|
| ENSG001  | 95.28865  | 0.00399148     | 0.225010| 0.0177391 | 0.9858470| 0.996699 |
| ENSG002  | 4359.09632| -0.23842494    | 0.127094| -1.8759764| 0.0606585| 0.289604 |
| ENSG003  | 419.06811 | -0.10185506    | 0.146568| -0.6949338| 0.4870968| 0.822681 |
| ...      | ...       | ...            | ...     | ...       | ...      | ...      | 
| ENSG00N  | 4863.807  | 0.0179729      | 0.194137| 0.0925784 | 0.9262385| 0.986726 |

- **normalizedCounts_TMM[date].csv:** Contains TMM normalized counts

| col 1 | col2 |
| ----- | ---- |
| a     | b    |
| ...   | ...  |
| y     | z    |

### GSEA Results
- **GO_Analysis_[comparison].csv:** Results for gene set enrichment analysis

| col 1 | col2 |
| ----- | ---- |
| a     | b    |
| ...   | ...  |
| y     | z    |

### Figures/Visualizations

- **[comparison]volcano.png:**

<img src="https://github.com/user-attachments/assets/74671325-89ff-4108-a8b1-8ef0d2471bbd" alt="volcano.png" width="40%">

- **[comparison]heatmap.png:**

<img src="https://github.com/user-attachments/assets/1e87ceaa-d2eb-4f68-8b9a-8e4142d4a6e1" alt="heatmap.png" width="40%">

- **[comparison]GSEA.png:**

<img src="https://github.com/user-attachments/assets/f4b7f84b-9d7e-4303-b0d4-cfa5fc87f2f6" alt="gsea.png" width="40%">

**PCA**
- Static 2D plot

<img src="https://external-preview.redd.it/sS_GFhS_OsMz6x0euch2EmKFeGKHjF2vzWpxguw6U0s.jpg?auto=webp&s=d121e3bb7d19edaeef99db3a098242a4443fb9f3" alt="pca.jpg" width="40%">

- Interactive 2D plot

<img src="https://external-preview.redd.it/sS_GFhS_OsMz6x0euch2EmKFeGKHjF2vzWpxguw6U0s.jpg?auto=webp&s=d121e3bb7d19edaeef99db3a098242a4443fb9f3" alt="pca.jpg" width="40%">

- Interactive 3D plot

<img src="https://external-preview.redd.it/sS_GFhS_OsMz6x0euch2EmKFeGKHjF2vzWpxguw6U0s.jpg?auto=webp&s=d121e3bb7d19edaeef99db3a098242a4443fb9f3" alt="pca.jpg" width="40%">

## Future Improvements

- Add support for EdgeR package to perform differential expression analysis.
- Add functionality to perform analysis using multi-comparison a contrast matrix as well as covariates.
- Add KEGG analysis module.
  
## Contact

If you need to reach the BISR group please email us at: [mccbioinfo@vcu.edu] or open a Github Issue to this repo.

## License

[GPL-3.0 license](https://github.com/VCU-Bioinformatics-Core/bulk_rnaseq_analyses/tree/main?tab=GPL-3.0-1-ov-file#)
