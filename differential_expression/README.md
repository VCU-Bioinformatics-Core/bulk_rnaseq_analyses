# Differential Expression Analysis Pipeline

The `de.R` script performs differential expression analysis for RNA-seq data using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and performs additional auxiliary analysis. Below is synopsis of this scripts features: 

- Detection of differentially expression genes (DEG) with DESeq2
- Visualization of DEG results using volcano plots
- Visualization of top differentially expressed genes using heatmaps
- Interpretation of prioritied genes using Gene Set Enrichment Analysis (GSEA)
- Plotting of PCA plots (2D and 3D interactive versions)
- Normalization of read counts using trimmed mean of M values (TMM)
- Organization of results into separate sub-directories
- Support for both human and mouse annotations

*Note:*
- *This script has been specialized for first pass analyses, for much more complex analyses please utilize DESeq2 directly*
- *Currently, this script only supports pairwise comparisons. Multiple comparisons will be supported in later versions*
- *We are actively developing this script to handle count data from any source, however, at this stage the script works best in conjunction with the nf-core rnaseq pipeline and its output (merged counts).*

## Table of Contents
- [Features]()
- [Preparing your R Environment]()
- [Input Files Required]()
- [Usage]()
- [Arguments]()
- [Output Structure]()

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

## Input Files Required

1. **Counts Matrix** (`counts.tsv`):

- Tab-separated file containing gene counts
- First column should contain gene IDs (ENSEMBL format)
- Subsequent columns should contain count data for each sample

| gene_id            | sample1 | sample2 | sample3 | sample4 | sample5 |
| ------------------ | ------- | ------- | ------- | ------- | ------- |
| ENSMUSG00000000001 | 1234    | 2345    | 3456    | 4567    | 5678    |
| ENSMUSG00000000002 | 2345    | 3456    | 4567    | 5678    | 6789    |
| ENSMUSG00000000003 | 3456    | 4567    | 5678    | 6789    | 7890    |
| ENSMUSG00000000004 | 4567    | 5678    | 6789    | 7890    | 8901    |
| ENSMUSG00000000005 | 5678    | 6789    | 7890    | 8901    | 9012    |
  
2. **Sample Sheet** (`samplesheet.csv`):

- CSV file containing sample metadata
- Must include columns:
- `sample`: Sample identifiers (should be the same as count)
- `condition`: Group/condition labels

| sample  | condition  |
| ------- | ---------- |
| sample1 | group1     |
| sample2 | group1     |
| sample3 | group2     |
| sample4 | group2     |
| sample5 | group3     |

3. **Contrast Matrix** (`contrasts.tsv`):

- Tab-separated file defining comparisons
- Must include columns:
- `GroupID`: Group identifiers
- `SampleID`: Sample Identifiers (should be the same as sample names in the counts matrix
- Additional columns for each comparison (1 = experimental, 0 = control, leave cell empty if not included in the comparison)

| GroupID    | treatment_vs_control | treatment2_vs_control2 |
| ---------- | -------------------- | ---------------------- |
| control    | 0                    |                        |
| treatment  | 1                    |                        |
| treatment2 |                      | 1                      |
| control2   |                      | 0                      |


## Usage

```bash

module load R/4.4.1

Rscript de.R \
--counts path/to/counts.tsv \
--contrasts path/to/contrasts.tsv \
--samplesheet path/to/samplesheet.csv \
--outdir path/to/output \
--runid analysis_name \
--annotation mouse

```

### Arguments

- `-c, --counts`: Path to the merged counts file

- `-m, --contrasts`: Path to the contrast matrix file

- `-s, --samplesheet`: Path to the sample sheet file

- `-o, --outdir`: Output directory (default: ./output)

- `-r, --runid`: Unique identifier for the analysis run

- `-a, --annotation`: Genome to use for annotation: 'mouse' or 'human' (default: mouse)  


## Output Structure

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


## Output Files

### Analysis Results

- **DESeq2 Results**: Differential expression statistics and annotations
- **GSEA Results**: Gene set enrichment analysis results
- **Normalized Counts**: TMM-normalized expression values

### Visualizations

- **Volcano Plots**: Highlighting significantly differentially expressed genes
- **Heatmaps**: Expression patterns of significant genes
- **GSEA Plots**: Enriched gene sets visualization
- **PCA Plots**: Sample clustering and quality control
- Static 2D plot
- Interactive 2D plot
- Interactive 3D plot

  

## Example Commands

### Mouse Analysis

```bash

Rscript de.R \
--counts mouse_counts.tsv \
--contrasts contrasts.tsv \
--samplesheet samplesheet.csv \
--outdir mouse_results \
--runid mouse_experiment \
--genome mouse

```
  
### Human Analysis

```bash

Rscript de.R \
--counts human_counts.tsv \
--contrasts contrasts.tsv \
--samplesheet samplesheet.csv \
--outdir human_results \
--runid human_experiment \
--genome human

```  


## Notes

- Gene annotations use either:
	- Mouse genome (org.Mm.eg.db)
	- Human genome (org.Hs.eg.db)
- Minimum group size is set to 3 samples
- Genes with low counts (< 10 in all samples) are filtered out
- P-value threshold for significance is 0.05
- Log2 fold change threshold is 0.58


## Future Improvements

- Add a module to produce R Markdown Reports.
- Add support for EdgeR package to perform differential expression analysis.
- Add functionality to perform analysis using multi-comparison contrast matrix and using covariates.


## License

[GPL-3.0 license](https://github.com/VCU-Bioinformatics-Core/bulk_rnaseq_analyses/tree/main?tab=GPL-3.0-1-ov-file#)

  
## Contact
  
[mccbioinfo@vcu.edu]
