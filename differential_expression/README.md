# Differential Expression Analysis Pipeline

## Introduction
This pipeline performs differential gene expression (DGE) analysis of RNA-seq count data using the DESeq2 R package (Love et al., 2014). It is designed to provide a streamlined, reproducible workflow for identifying genes with statistically significant expression differences between experimental conditions. The pipeline incorporates best-practice recommendations for RNA-seq data analysis, including normalization, dispersion estimation, and hypothesis testing.  It automatically generates a comprehensive report with key visualizations for each specified comparison, along with auxiliary analyses such as Gene Set Enrichment Analysis (GSEA) and Principal Component Analysis (PCA) for sample exploration.

*Note: This pipeline is intended for a "first pass" analysis. For custom or complex analyses, please contact our core and [submit a Jira ticket](https://www.masseycancercenter.org/research/shared-resource-cores/bioinformatics/)*

## Usage

### Preparing Your Data
Two input files are required in specific formats: the **Raw Merged Count Matrix** and the **Samplesheet**.

#### 1. Raw Merged Count Matrix
**Required format:** Tab-Separated Values (`.tsv`)
This file contains the raw merged gene expression counts. The first column must contain gene identifiers (`gene_id`), and subsequent columns should contain the raw count data for each sample. The `header` values for these sample columns must match the sample identifiers (`SampleID`) used in the samplesheet. The pipeline expects integer counts, as is typical for RNA-seq data. While the script has the functionality to automatically parse and handle count data from various quantification tools, the users are still responsible to remove any additional columns besides the gene IDs and sample counts. It currently works best with the merged counts output from pipelines like [`nf-core's rnaseq Nextflow pipeline`](https://nf-co.re/rnaseq).

| gene_id            | sample1_r1 | sample1_r2 | sample1_r3 | sample2_r1 | sample2_r1 | sample2_r2 | sample2_r3 | sample3_r1 | sample3_r2 | sample3_r3 | sample4_r1 | sample4_r2 | sample4_r3 |
| ------------------ | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| ENSMUSG00000000001 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000002 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000003 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000004 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
| ENSMUSG00000000005 | 1234 | 2345 | 3456 | 4567 | 5678 | 6789 | 2345 | 3456 | 4567 | 5678 | 3456 | 4567 |
  
### 2. Samplesheet:
**Required format:** Comma-Separated Values (`CSV`)
This file contains metadata for each sample, including group identifiers and binary indicators for specific comparisons. The columns must be organized as follows:
- **SampleID:** Sample identifiers that exactly match the sample column headers in the count matrix.
- **GroupID:** Group identifiers for each sample (e.g., `experiment1`, `control1`). These group identifiers are crucial for defining the experimental design in DESeq2.
- [comparison]: Subsequent columns define pairwise comparisons. The column name should follow the format `experiment_vs_control`. For each comparison column, use `1` to indicate samples belonging to the experimental group,
  `0` for the control group, and leave the cell `blank` for samples to be excluded from that specific comparison. This design matrix setup allows the user to specify which samples are used for each comparison, providing flexibility in complex experimental designs.

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

The pipeline is executed using an R script. Here are example commands for running it on VCU HPRC: (high performance research computing servers):

#### Arguments

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

## Pipeline Output
The pipeline generates an output directory (specified by `--outdir`) containing two main subdirectories: `data_frames` and `figures`.

```
output/
├── data_frames/
│   ├── de_data/
│   │   ├── DESeq2_[comparison].csv
│   │   └── normalizedCounts_TMM[date].csv
│   └── gsea_data/
│       └── GO_Analysis_[comparison].csv
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

### data_frames
This directory contains output data frames organized into subdirectories for DGE and GSEA results.

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

#### normalizedCounts_TMM[date].csv:
Contains the read counts normalized using the Trimmed Mean of M-values (TMM) method. TMM normalization is performed using the edgeR package (Robinson et al., 2010) to account for differences in library size and RNA composition between samples. The date is appended to the filename for version control.

| col 1 | col2 |
| ----- | ---- |
| a     | b    |
| ...   | ...  |
| y     | z    |

#### gsea_data
**GO_Analysis_[comparison].csv:** Contains the results of the Gene Set Enrichment Analysis (GSEA) using Gene Ontology (GO) terms for each comparison. GSEA is performed using a suitable R package (e.g., clusterProfiler) to identify enriched GO terms among the differentially expressed genes. This analysis is skipped for a comparison if the Gene Set identified doesnot have enough genes. 

| col 1 | col2 |
| ----- | ---- |
| a     | b    |
| ...   | ...  |
| y     | z    |


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

- An Interactive 2D plot (HTML format) that allows for sample labeling and exploration.
<img src="https://external-preview.redd.it/sS_GFhS_OsMz6x0euch2EmKFeGKHjF2vzWpxguw6U0s.jpg?auto=webp&s=d121e3bb7d19edaeef99db3a098242a4443fb9f3" alt="pca.jpg" width="40%">

- An interactive 3D PCA plot, providing a three-dimensional view of sample relationships.
<img src="https://external-preview.redd.it/sS_GFhS_OsMz6x0euch2EmKFeGKHjF2vzWpxguw6U0s.jpg?auto=webp&s=d121e3bb7d19edaeef99db3a098242a4443fb9f3" alt="pca.jpg" width="40%">

These plots help to assess the overall quality of the data, identify potential outliers, and visualize the separation of samples according to experimental conditions.

## Limitations
- This pipeline currently supports only pairwise comparisons. Support for more complex designs with multiple comparisons with covariates and contrast matrices will be added in future versions. This is a limitation for experiments with more than two conditions.
- The pipeline works best with merged count matrices generated from pipelines like [`nf-core's rnaseq Nextflow pipeline`](https://nf-co.re/rnaseq). While the script is being developed to handle count data from any source, users may need to pre-format their count matrices accordingly. Specifically, the matrix should have a `gene_id` column, with subsequent columns containing raw counts for each sample.
- Gene prefiltering is performed as described in the [DESeq2 documentation](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html): Genes are excluded if they do not have three or more samples with a read count of 10 or greater. This step aims to remove genes with very low expression, which can reduce the memory size of the dds data object, and increase the speed of count modeling within DESeq2.
- Differential expression results are considered significant if the Benjamini-Hochberg adjusted p-value (padj) is less than or equal to 0.05 and the absolute log2 fold change is greater than 0.58 (corresponding to an absolute fold change of 1.5).

## Future Improvements
- Add support for the EdgeR package for differential expression analysis. DESeq2 and edgeR are both popular packages for DGE analysis, and providing both options would increase user flexibility.
- Implement functionality to perform analyses using a contrast matrix, multi-factor designs and covariates.  This would allow for the analysis of more complex experimental designs, including those with multiple factors and batch effects.
- Incorporate a KEGG pathway analysis module. This would provide additional biological context for the DGE results by identifying enriched KEGG pathways.
  
## Contact
For questions or issues, please contact the BISR group at [mccbioinfo@vcu.edu] or open a GitHub issue in the repository.

## License

[GPL-3.0 license](https://github.com/VCU-Bioinformatics-Core/bulk_rnaseq_analyses/tree/main?tab=GPL-3.0-1-ov-file#)
