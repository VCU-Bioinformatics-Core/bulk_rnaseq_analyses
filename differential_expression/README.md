# Differential Expression Analysis Pipeline

This Pipeline performs Differential Expression Analysis for RNA-seq data using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and performs additional auxiliary analyses. It aumatically produces an report to summarize the results and showcase visualizations for each comparision.

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

*Note:*
- *This pipeline aims to perform a "first pass" analysis (as required by the VCU BISR), for custom/complex analysis please contact our core and [submit a Jira ticket](https://www.masseycancercenter.org/research/shared-resource-cores/bioinformatics/)*
- Current version only supports pairwise comparisons. We aim to support DESeq2 designs with multiple comparisons in later versions.
- To avoid complexity, i.e. overfitting of code that is tailored to parse input raw counts file from specific gene quantification tools; Users are required to remove any additional columns except **Gene id** and subsequent **sample** columns (with raw count data).
- *We are actively developing this script to handle count data from any source. However, at this stage, the script works best in conjunction with the Nextflow nf-core rnaseq pipeline and its output (merged counts).*
- Gene Prefiltering is based on [Deseq2 Documentation Prefiltering Section](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html): Genes are removed if they do not have three or more samples with a read count of 10 or greater.
- DE results are considered significant when the adjusted P-value, calculated using the Benjamini-Hochberg correction, is less than or equal to 0.05
- The threshold for absolute fold change is set to 1.5 (0.58 Log2 fold change)

## Preparing Your Data
Two input files are required with specific file formats to ensure proper processing: **Raw Merged Count Matrix**, and **Samplesheet**. Examples with tabular representation of the required data are provided below.

### 1. Raw Merged Count Matrix:
**Required file format:** `TSV`

This file contains the raw expression counts data and files with any other type of column format are automatically converted to the following format:
- First column should contain gene IDs.
- Subsequent columns with raw count data for each sample. Note: header values for these columns should be the same as the sample identifiers used in the samplesheet.

| gene_id            | sample1_r1 | sample1_r2 | sample1_r3 | sample2_r1 | sample2_r1 | sample2_r2 | sample2_r3 | sample3_r1 | sample3_r2 | sample3_r3 | sample4_r1 | sample4_r2 | sample4_r3 |
| ------------------ | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| ENSMUSG00000000001 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000002 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000003 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000004 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000005 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
  
### 2. Samplesheet:
**required file format:** `CSV`

This file contains samples with metadata/clinical variables that will be utilized for DEG analysis. The columns should be organised exactly in the same order exemplified below.
- First Column (**SampleID**) should contain sample identifiers that correspond with sample column names after the **gene_id** column of the count matrix.
- Second Column (**GroupID**) should contain the group identifiers.
- All subsequent columns are dedicated comparisons (rows for these columns should be a binary value, representing `1` for the treatment group, `0` for the control, and `blank` if the given sample needs to be ignored for the current comparison). These comparison columns are expected to follow the naming convention exemplified below; i.e. group1_vs_group2

| SampleID | GroupID      | experimental1_vs_control1	| experimental2_vs_control2	|
| -------- | ------------ | ------------------------- | ------------------------- |
| sample1_r1  | experimental1    | 1                  |              							|
| sample1_r2  | experimental1    | 1                  |							              |
| sample1_r3  | experimental1    | 1                  |							              |
| sample2_r1  | control1  | 0                 |                        	|
| sample2_r2  | control1  | 0                 |                        	|
| sample2_r3  | control1  | 0                 |                        	|
| sample3_r1  | experimental2 |                		    | 1                      	|
| sample3_r2  | experimental2 |                		    | 1                      	|
| sample3_r3  | experimental2 |                		    | 1                      	|
| sample4_r1  | control2 |                		| 0                      	|
| sample4_r2  | control2 |                		| 0                      	|
| sample4_r3  | control2 |                		| 0                      	|


## Running the Script

Please refer to example commands below if you choose to run the pipeline on VCU HPRC (high performance research computing servers):

### Arguments

- `-c, --counts`: Path to the merged counts file **(Mandatory)**

- `-s, --samplesheet`: Path to the sample sheet file **(Mandatory)**

- `-o, --outdir`: Output directory (default: ./output)

- `-r, --runid`: Unique identifier for the analysis run **(Mandatory)**

- `-a, --annotation`: Genome to use for annotation: 'mouse' or 'human' (default: mouse)  

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

## Understanding the Outputs

An output directory is automatically created with the name specfied with the `--outdir` command line argument. You can also provide an absolute path here with the output directory name in the end. For eg. `/lustre/home/lab/projects/project_name/outputs`, where the `outdir` name is "outputs". All the outputs are stored and organized into 3 main subdirectories.
- **data_frames**: contains subdirectories for DGE and GSEA output dataframes.
  - **de_data**: contains DESeq2 results data frame for each comparision `DESeq2_[comparison].csv` and TMM normalized counts stored within `normalizedCounts_TMM[date].csv`.
  - **gsea_data**: contains GSEA-GO results data frame `GO_Analysis_[comparison].csv` for each comparision if processed succesfully.
- **figures**: contains subdirectories for several visualizations generated for each comparison.
  - **volcano**: for each comparison, contains a volcano plot named `[comparison]volcano.png`
  - **heatmap**: for each comparison, contains a heatmap named `[comparison]heatmap.png`
  - **gsea**: for each comparison, contains the gsea dot plot named `[comparison]GSEA.png`
  - **pca**: pca plots of the data, including a simple labeled PCA plot, an html interactive PCA plot 

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
