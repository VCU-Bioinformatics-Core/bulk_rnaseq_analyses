#!/usr/bin/env Rscript

# Required libraries
library(rmarkdown)
library(knitr)
library(dplyr)
library(glue)

# Function to generate automated R Markdown report
# source("report_generator.R"); generate_report('analysis.rds', output_dir = getwd())
generate_report <- function(analysis_results_path, output_dir = "./", report_prefix = "rnaseq_analysis", analyst="Mikail Bala") {
  
  # Validate inputs
  if (!file.exists(analysis_results_path)) {
    stop("Analysis results file does not exist:", analysis_results_path)
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate unique filename with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_file <- file.path(output_dir, paste0(report_prefix, "_", timestamp, ".html"))
  
  # Load results
  rds_data <- readRDS(analysis_results_path)
  
  # Extract components from RDS
  results <- rds_data[[1]]
  comparisons <- rds_data[[2]]
  out_dirs <- rds_data[[3]]
  pca_plot <- rds_data[[4]]
  pca_plotly <- rds_data[[5]]
  pca_3d <- rds_data[[6]]
  annotation <- rds_data[[7]]
  
  # calculate the current date
  date <- format(Sys.time(), "%B %d, %Y")
  
  # Create R Markdown template (beginning section)
  rmd_content <- glue('---
title: "RNA-Seq Differential Expression Analysis Report"
author: "Bioinformatics Shared Resources at VCU"
date: "{date}"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
    highlight: tango
    df_print: paged
    code_folding: hide
---

```{{r setup, include=FALSE}}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
library(knitr)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(htmlwidgets)

# Load results
rds_data <- readRDS("{analysis_results_path}")
results <- rds_data[[1]]
comparisons <- rds_data[[2]]
out_dirs <- rds_data[[3]]
pca_plot <- rds_data[[4]]
pca_plotly <- rds_data[[5]]
pca_3d <- rds_data[[6]]
annotation <- rds_data[[7]]
```

## Overview

**Analyst: {analyst}**

This report contains the results of differential expression analysis for `r length(comparisons)` comparisons.
It starts with the **Sample Exploration using PCA** section that highlights sample similarity and
helps to explore and understand the relationships between samples. Next, we present the
**Differential Expression Results** section with subsections for each comparison (more details
below). For further exploration of the results, please refer to the output directories
containing the raw data files.


## Pipeline

For this analysis we used the following steps:

1. **Data preprocessing**: Count data was filtered to remove genes with low expression.
2. **Normalization**: TMM normalization was applied using edgeR TMM function.
3. **Differential expression**: DESeq2 was used to identify differentially expressed genes.
4. **Functional analysis**: Gene Set Enrichment Analysis (GSEA) was performed using clusterProfiler.
5. **Thresholds**: Genes with adjusted p-value < 0.05 and |FC| >= 1.5 (equivalent to log2(FC) = 0.58) were considered differentially expressed.
6. **Genome annotation**: "{annotation}"

As part of this pipeline we produce the following files for your downstream use:

```
output/
├── de_data/
│ ├── DESeq2_[comparison].csv
│ └── normalizedCounts_TMM[date].csv
├── gsea_data/
│ └── GO_Analysis_[comparison].csv
└── figures/
  ├── volcano/
  │   └── [comparison]volcano.png
  ├── heatmap/
  │ └── [comparison]heatmap.png
  ├── gsea/
  │ └── [comparison]GSEA.png
  └── pca/
    ├── PCA_plot.png
    ├── allsamples_PCA_plot.pdf
    └── allsamples_PCA_plot3D.pdf
```

## Sample Exploration using PCA

PCA was performed to visualize the overall patterns of gene expression across samples and
to identify potential batch effects or outliers. **What we expect:** Samples with similar
expression profiles should cluster together, while dissimilar samples are expected to
separate into distinct clusters.


```{{r pca-plot, fig.width=6, fig.height=4 }}
print(pca_plot)
```

Interactive PCA plots can also be found in the following files:

- 2D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot.html")`
- 3D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot3D.html")`

```{{r results-summary}}
# Create a summary table for all comparisons
summary_table <- data.frame(
  Comparison = character(),
  Experimental = character(),
  Control = character(),
  Total_DEGs = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  GSEA_Performed = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(comparisons)) {{
  if (!is.null(results[[i]])) {{
    res_df <- results[[i]]$deseq
    
    # Count DEGs (padj < 0.05 & |log2FC| >= 0.58)
    if (!is.null(res_df)) {{
      total_degs <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 0.58)
      up_degs <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 0.58)
      down_degs <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -0.58)
      
      gsea_status <- ifelse(!is.null(results[[i]]$gsea), "Yes", "No")
      
      summary_table <- rbind(summary_table, data.frame(
        Comparison = comparisons[[i]]$name,
        Experimental = comparisons[[i]]$exp,
        Control = comparisons[[i]]$ctrl,
        Total_DEGs = total_degs,
        Upregulated = up_degs,
        Downregulated = down_degs,
        GSEA_Performed = gsea_status,
        stringsAsFactors = FALSE
      ))
    }}
  }}
}}
```
  
## Differential Expression Results
This section contains a high-level table summarizing the differential expression results,
followed by a subsection for each comparison that visualizes the results using a volcano
plot, heatmap, table of top hits, and GSEA results (shown as both a bubble plot and a table).

### Summary of Comparisions

Briefly, we report the following:

- Total number of comparisons analyzed: `r length(comparisons)`
- Comparisons with the highest number of DEGs: `r summary_table$Comparison[which.max(summary_table$Total_DEGs)]` (`r max(summary_table$Total_DEGs)` DEGs)
- Comparisons with the lowest number of DEGs: `r summary_table$Comparison[which.min(summary_table$Total_DEGs)]` (`r min(summary_table$Total_DEGs)` DEGs)

The table below summarizes differential expression results from all comparisons.
Specifically, the **Total_DEGs** column reports the number of genes with an adjusted
p-value < 0.05 and |FC| ≥ 1.5. The **Upregulated** and **Downregulated** columns break
down this count accordingly. The **GSEA_Performed** column indicates whether GSEA was run
for each case.

```{{r display-summary-table}}
# Display the summary table
kable(summary_table, caption = "")
```\n')
  
  
  # Add detailed sections for each comparison
  for (i in seq_along(comparisons)) {
    
    name = comparisons[[i]]$name
    exp = comparisons[[i]]$exp
    ctrl = comparisons[[i]]$ctrl
    
    comparison_section <- glue('\n
### {name}

**Experimental design**

- **Experimental group**: {exp}
- **Control group**: {ctrl}


**Volcano Plot**

- Description: main visualization for differential expression results. This rendition uses
  a **red horizontal line** to indicate the significant p-value threshold of \\< 0.05. Therefore,
  every point (gene) above that red line can be considered statistically signficant (note:
  before FDR correction). In addition, the **black vertical lines** indicate 1.5 fold change.
  Genes highlighted in red are up-regulated in group1 compared to the group2. Genes
  highlighted in blue are down-regulated in group1 compared to group2. Genes in gray, do
  not meet the thresholds for both logFC and p-value.
- Data point: a gene
- X-axis: log2 fold-change (group1/group2)
- Y-axis: -log10(p-value)

```{{r volcano-{i}, out.width="80%", out.height="80%" }}
# Display volcano plot from file
volcano_path <- file.path(out_dirs$volcano, paste0("{name}_volcano.png"))
if (file.exists(volcano_path)) {{
  knitr::include_graphics(volcano_path)
}} else {{
  cat("Volcano plot not available for this comparison")
}}
```


**Heatmap**

- Description: a heatmap of **z-score normalized** read counts data with application of 
  **hierarchical clustering** by both samples and expression levels
- X-axis: samples
- Y-axis: genes
- Color-scale: z-score normalized read counts 


```{{r heatmap-{i} }}
# Display heatmap from file
heatmap_path <- file.path(out_dirs$heatmap, paste0("{name}_heatmap.png"))
if (file.exists(heatmap_path)) {{
  knitr::include_graphics(heatmap_path)
}} else {{
  cat("Heatmap not available for this comparison")
}}
```


**Top Differentially Expressed Genes**

Table of the top differential expressed genes with both nominal (**pvalue**)
and adjusted pvalues (**padj**).

```{{r top-degs-{i} }}
# Display top DEGs table

deg_flags = c(0,0)
if (!is.null(results[[{i}]]) && !is.null(results[[{i}]]$deseq)) {{
  top_up <- results[[{i}]]$deseq %>%
    filter(!is.na(padj) & padj < 0.05 & log2FoldChange >= 0.58) %>%
    mutate(log2FoldChange = round(log2FoldChange, 2),
           pvalue = formatC(pvalue, format = "e", digits = 2),
           padj = formatC(padj, format = "e", digits = 2)) %>%
    arrange(padj) %>%
    head(20)
  
  top_down <- results[[{i}]]$deseq %>%
    filter(!is.na(padj) & padj < 0.05 & log2FoldChange <= -0.58) %>%
    mutate(log2FoldChange = round(log2FoldChange, 2),
           pvalue = formatC(pvalue, format = "e", digits = 2),
           padj = formatC(padj, format = "e", digits = 2)) %>%
    arrange(padj) %>%
    head(20)
  
  # update the flags
  if (nrow(top_up) > 0) {{
    deg_flags[1] <- 1
  }}
  if (nrow(top_down) > 0) {{
    deg_flags[2] <- 1
  }}

}}
```

```{{r top-check-degs-{i} }}
if (sum(deg_flags) == 0) {{
  cat("No differential expression results available for this comparison")
}}

```

```{{r top-up-degs-{i} }}

if (deg_flags[1] > 0){{
    DT::datatable(top_up %>% select(ENSEMBL_ID, SYMBOL, log2FoldChange, pvalue, padj, GENENAME),
               caption = "Top Up Regulated Genes")
}} else {{
  cat("No upregulated genes found\\n\\n")
}}
  
```

```{{r top-down-degs-{i} }}

if (deg_flags[2] > 0){{
  DT::datatable(top_down %>% select(ENSEMBL_ID, SYMBOL, log2FoldChange, pvalue, padj, GENENAME),
               caption = "Top Down Regulated Genes")
}} else {{
  cat("No upregulated genes found\\n\\n")
}}
  
```


**Gene Set Enrichment Analysis**

- Description: Gene Set Enrichment Analysis using (GO) terms with both a bubble plot and
 table
- X-axis: the gene ratio (# genes related to Gene Set / total number of significant genes) to which
the term is enriched
- Y-axis: a given Gene Set

```{{r gsea-{i}, out.width="100%", out.height="100%"}}
# Display GSEA results from file
gsea_path <- file.path(out_dirs$gsea, paste0("{name}_GSEA.png"))
if (file.exists(gsea_path)) {{
  knitr::include_graphics(gsea_path)
}} else {{
  cat("GSEA results not available for this comparison")
}}
```


```{{r gsea-table-{i} }}
# Display top GSEA results
if (!is.null(results[[{i}]]) && !is.null(results[[{i}]]$gsea)) {{
  gsea_results <- as.data.frame(results[[{i}]]$gsea) %>%
      mutate(enrichmentScore=formatC(enrichmentScore, format="e", digits=2),
             NES=formatC(NES, format="e", digits=2),
             p.adjust=formatC(p.adjust, format="e", digits=2),
             qvalue=formatC(qvalue, format="e", digits=2))
  if (nrow(gsea_results) > 0) {{
    DT::datatable(gsea_results %>%
                 select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue) %>%
                 head(20),
                 caption = "Top enriched gene sets")
  }} else {{
    cat("No significant enriched gene sets found")
  }}
}} else {{
  cat("No GSEA results available for this comparison")
}}

```')
    
    rmd_content <- paste0(rmd_content, comparison_section)
  }


extra_content <- '\n
## Utilizing IPA

The CSV Differential Expression output from DESeq2 (available in the results directory
provided alongside this report *./output/de_data/DESeq2_[comparison_name].csv*), can be
uploaded directly into **QIAGEN Ingenuity Pathway Analysis (IPA)** for self-exploration of
pathways predicted to be enriched by this experimental condition. Massey’s BISR provides
access to VCU’s license of IPA. If you do not already have an account associated with this
license, you may reach out to **morecockcm@vcu.edu** with your name, VCU health or VCU
email, and request for IPA. To perform a core expression analysis, login with your
credentials here: **https://analysis.ingenuity.com/pa** and follow the instructions [here](https://qiagen.my.salesforce-sites.com/KnowledgeBase/KnowledgeNavigatorPage?id=kA41i000000L6rMCAS).

We host an annual hands-on training for IPA at the beginning of the fall semester. Please
email BISR if you would like to be a part of this training. In the meantime, QIAGEN has a
playlist of user-friendly tutorials available on Youtube titled “QIAGEN IPA Training
Videos” the **Qiagen Digital Insights Youtube** page.

## Manuscript-Ready Text

### Methods
Raw RNA-Seq fastq files were processed by the VCU Massey Comprehensive Cancer Center Bioinformatics Shared Resource (BISR) using the NextFlow nf-core/rnaseq v3.18.0 pipeline [1]. Briefly, this pipeline assesses sequencing quality using FastQC v 0.12.1 [2] before and after trimming, performs adaptor trimming with Trim Galore! v0.6.10 [3], and aligns sequencing reads to the GRCh38 human primary assembly reference genome using STAR v 2.7.11b [4] with transcriptome quantification by Salmon v1.10.3 [5].  Pipeline output includes gene expression raw count data and a comprehensive QC report compiled by MultiQC v1.25.1 [6].
Differential expression analysis was performed using DESeq2 v 1.44.0 [7]. Lowly expressed genes were filtered out per DESeq2 methods [7] prior to normalization and differential expression testing. Significance was calculated using the Wald-test and adjusted using Benjamini Hochberg False Discovery Rate (FDR). Volcano plots and heatmaps were generated using the EdgeR TMM normalized count data and visualized using R packages. Significant differentially expressed genes (DEGs) are defined as those with an FDR<0.05 and absolute fold-change of 1.5 (log2 fold-change = 0.58) or greater. Gene Set Enrichment Analysis (GSEA) [8] for Gene Ontology terms (GO) was performed using the clusterProfiler package [9] across all genes, regardless of significance. All computational analyses were performed on VCU’s High Performance Research Computing cluster.

### References
1) Ewels P, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Feb 13. doi:10.1038/s41587-020-0439-x 
Andrews S. FastQC: A Quality Control Tool for High Throughput Sequence Data. Babraham Bioinformatics; 2010. Accessed June 18, 2025. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

2) Krueger F. Trim Galore! v0.6.10. 2023. Available at: https://github.com/FelixKrueger/TrimGalore. Accessed June 18, 2025.

3) Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635

4) Patro, R., Duggal, G., Love, M.I., Irizarry, R.A., Kingsford, C., 2017. Salmon provides fast and bias-aware quantification of transcript expression. Nat. Methods 14, 417–419. https://doi.org/10.1038/nmeth.4197

5) Ewels, P., Magnusson, M., Lundin, S., Käller, M., 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32, 3047–3048. https://doi.org/10.1093/bioinformatics/btw354

6) Love, M.I., Huber, W., Anders, S., 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550. https://doi.org/10.1186/s13059-014-0550-8

7) A. Subramanian, P. Tamayo, V.K. Mootha, S. Mukherjee, B.L. Ebert, M.A. Gillette, A. Paulovich, S.L. Pomeroy, T.R. Golub, E.S. Lander, & J.P. Mesirov, Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles, Proc. Natl. Acad. Sci. U.S.A. 102 (43) 15545-15550, https://doi.org/10.1073/pnas.0506580102 (2005).

8) Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012 May;16(5):284-7. doi: 10.1089/omi.2011.0118. Epub 2012 Mar 28. PMID: 22455463; PMCID: PMC3339379.

### Required Acknowledgements

Please include the following statements in your acknowledgements manuscript section:

- “Services in support of the research project were provided by the VCU Massey Comprehensive Cancer Center Bioinformatics Shared Resource. Massey is supported, in part, with funding from NIH-NCI Cancer Center Support Grant P30 CA016059.”

- “High Performance Computing resources provided by the High Performance Research Computing (HPRC) core facility at Virginia Commonwealth University (https://hprc.vcu.edu) were used for conducting the research reported in this work.”  
'
  rmd_content <- paste0(rmd_content, extra_content)
  

  # Write the R Markdown file
  fn <- glue('{report_prefix}_{timestamp}.Rmd')
  rmd_file <- file.path(output_dir, fn)
  writeLines(rmd_content, rmd_file)
  
  # Render the R Markdown to HTML
  rmarkdown::render(rmd_file, output_file = output_file, quiet = FALSE)
  
  # Return the path to the generated report
  return(output_file)
}
