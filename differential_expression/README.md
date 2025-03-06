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

## Preparing Your Data
As mentioned earlier, this code was first developed after setting up the `nf-core/rnaseq` pipeline and therefore is guaranteed to work seamlessly with this pipeline. Despite this, we are working on generalizing this code to any count data. To utilize this script we require three pieces of data/files: 1) **Count Matrix**, 2) **Samplesheet**, and 3) **Contrast Matrix**. Each of these have a **mandatory** file name that will be described below. In addition, please pay close attention to the file format to ensure proper processing. Lastly, below each file description we provide a tabular representation of the required data as an example.

**_ATTENTION: header lines are expected for all files._**

### 1. **Count Matrix**:
**Mandatory name:** `counts.tsv`<br>
**File Format:** `TSV`

This file contains the expression counts data and must follow the following format:
- First column contains ENSEMBL gene IDs
- Subsequent columns should contain count data for each sample. Note: header values should be the sample identifiers use across other input files

| gene_id            | sample1 | sample2 | sample3 | sample4 | sample5 |
| ------------------ | ------- | ------- | ------- | ------- | ------- |
| ENSMUSG00000000001 | 1234    | 2345    | 3456    | 4567    | 5678    |
| ENSMUSG00000000002 | 2345    | 3456    | 4567    | 5678    | 6789    |
| ENSMUSG00000000003 | 3456    | 4567    | 5678    | 6789    | 7890    |
| ENSMUSG00000000004 | 4567    | 5678    | 6789    | 7890    | 8901    |
| ENSMUSG00000000005 | 5678    | 6789    | 7890    | 8901    | 9012    |
  
### 2. **Samplesheet**:
**Mandatory name:** `samplesheet.csv`<br>
**File Format:** `CSV`

This file contains samples with a comprehensive set of metadata/clinical variables (i.e. cancer_status, treatments, etc) and uses the format:
- First column contains sample identifiers that should correspond (one-to-one) with columns of the count matrix
- Subsequent columns contains metadata variables

| sample  | meta1      | meta2      | ... | metaN | 
| ------- | ---------- | ---------- | --- | ---   |
| sample1 | group1     | group1     | ... | ...   |
| sample2 | group1     | group1     | ... | ...   |
| sample3 | group2     | group1     | ... | ...   |
| sample4 | group2     | group1     | ... | ...   |
| sample5 | group3     | group1     | ... | ...   |

### 3. **Contrast Matrix**:
**Mandatory name:** `contrasts.tsv`<br>
**File Format:** `TSV`

This file contains samples with a subset of metadata/clinical variables that will be utilized for DEG analysis and uses the format:
- First column contains sample identifiers that should correspond (one-to-one) with columns of the count matrix
- Second column contains the group identifier
- Subsequent columns contain a binarized value representing `1` for the treatment group, `0` for the control, and `blank` if the given sample will not be used for the current comparison

| sample     | GroupID    | treatment_vs_control | treatment2_vs_control2 |
|------------| ---------- | -------------------- | ---------------------- |
| sample1    | control    | 0                    |                        |
| sample2    | treatment  | 1                    |                        |
| sample3    | treatment2 |                      | 1                      |
| sample4    | control2   |                      | 0                      |


## Running the Script

To run the script on VCU servers we require the following commands:

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

- `-c, --counts`: Path to the merged counts file **(Mandatory)**

- `-m, --contrasts`: Path to the contrast matrix file **(Mandatory)**

- `-s, --samplesheet`: Path to the sample sheet file **(Mandatory)**

- `-o, --outdir`: Output directory (default: ./output)

- `-r, --runid`: Unique identifier for the analysis run **(Mandatory)**

- `-a, --annotation`: Genome to use for annotation: 'mouse' or 'human' (default: mouse)  


### Example Commands

#### Mouse Analysis

```bash

Rscript de.R \
--counts mouse_counts.tsv \
--contrasts contrasts.tsv \
--samplesheet samplesheet.csv \
--outdir mouse_results \
--runid mouse_experiment \
--genome mouse

```
  
#### Human Analysis

```bash

Rscript de.R \
--counts human_counts.tsv \
--contrasts contrasts.tsv \
--samplesheet samplesheet.csv \
--outdir human_results \
--runid human_experiment \
--genome human

```  

## Understanding the Outputs
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

### Output Files

#### Analysis Results

- **DESeq2 Results**: Differential expression statistics and annotations
- **GSEA Results**: Gene set enrichment analysis results
- **Normalized Counts**: TMM-normalized expression values

#### Visualizations

- **Volcano Plots**: Highlighting significantly differentially expressed genes
- **Heatmaps**: Expression patterns of significant genes
- **GSEA Plots**: Enriched gene sets visualization
- **PCA Plots**: Sample clustering and quality control
- Static 2D plot
- Interactive 2D plot
- Interactive 3D plot

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
