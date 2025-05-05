#skipthisfornow - !/usr/bin/env Rscript

# Required libraries
library(rmarkdown)
library(knitr)
library(dplyr)

# Function to generate automated R Markdown report
generate_report <- function(analysis_results_path, output_dir = "./", report_prefix = "rnaseq_analysis") {
  
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
  
  # Create R Markdown template (beginning section)
  rmd_content <- sprintf('---
title: "RNA-Seq Differential Expression Analysis Report"
author: "Bioinformatics Shared Resources at VCU"
date: "%s"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
    highlight: tango
    df_print: paged
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
library(knitr)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(htmlwidgets)

# Load results
rds_data <- readRDS("%s")
results <- rds_data[[1]]
comparisons <- rds_data[[2]]
out_dirs <- rds_data[[3]]
pca_plot <- rds_data[[4]]
pca_plotly <- rds_data[[5]]
pca_3d <- rds_data[[6]]
```

## Overview

This report contains the results of differential expression analysis for `r length(comparisons)` comparisons.
We start with the **Sample Exploration using PCA** section that tries to identify major
trends and potential batch effects in the data. From there, we present the
**Differential Expression Results** section with subsections for each comparison (more details
below). For further exploration of the results, please refer to the output directories
containing the raw data files.


## Pipeline

For this analysis we used the following steps:

1. **Data preprocessing**: Count data was filtered to remove genes with low expression
2. **Normalization**: TMM normalization was applied using edgeR
3. **Differential expression**: DESeq2 was used to identify differentially expressed genes
4. **Functional analysis**: Gene Set Enrichment Analysis (GSEA) was performed using clusterProfiler
5. **Thresholds**: Genes with adjusted p-value < 0.05 and |FC| >= 1.5 were considered differentially expressed.
6. **Genome annotation**: Mouse (org.Mm.eg.db)

## Sample Exploration using PCA

PCA was performed to visualize the overall patterns of gene expression across samples and
to identify potential batch effects or outliers. **What we expect:** Samples with similar
expression profiles should cluster together, while dissimilar samples are expected to
separate into distinct clusters.


```{r pca-plot}
print(pca_plot)
```

Interactive PCA plots can also be found in the following files:

- 2D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot.html")`
- 3D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot3D.html")`
', format(Sys.time(), "%B %d, %Y"),  # Date in header
   analysis_results_path)
  
  
  # Add data extraction section
  data_summary_section <- '
```{r results-summary}
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

for (i in seq_along(comparisons)) {
  if (!is.null(results[[i]])) {
    res_df <- results[[i]]$deseq
    
    # Count DEGs (padj < 0.05 & |log2FC| >= 0.58)
    if (!is.null(res_df)) {
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
    }
  }
}
```'
rmd_content <- paste0(rmd_content, data_summary_section)

  
  # Add a diff section
  diff_section <- '
## Differential Expression Results
This section contains a high-level table summarizing the differential expression results,
followed by a subsection for each comparison that visualizes the results using a volcano
plot, heatmap, table of top hits, and GSEA results (shown as both a bubble plot and a table).

### Summary of Comparisions

Briefly, we report the following:

- Total number of comparisons analyzed: `r length(comparisons)`
- Total number of differential expressed genes across all comparisons: `r sum(summary_table$Total_DEGs)`
- Comparisons with the highest number of DEGs: `r summary_table$Comparison[which.max(summary_table$Total_DEGs)]` (`r max(summary_table$Total_DEGs)` DEGs)

The table below summarizes differential expression results from all comparisons.
Specifically, the **Total_DEGs** column reports the number of genes with an adjusted
p-value < 0.05 and |FC| ≥ 1.5. The **Upregulated** and **Downregulated** columns break
down this count accordingly. The **GSEA_Performed** column indicates whether GSEA was run
for each case.

```{r display-summary-table}
# Display the summary table
kable(summary_table, caption = "")
```
'
  rmd_content <- paste0(rmd_content, diff_section)
  
  
  # Add detailed sections for each comparison
  for (i in seq_along(comparisons)) {
    comparison_section <- sprintf('
### %s

**Experimental design**

- **Experimental group**: %s
- **Control group**: %s


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

```{r volcano-%d}
# Display volcano plot from file
volcano_path <- file.path(out_dirs$volcano, paste0("%s_volcano.png"))
if (file.exists(volcano_path)) {
  knitr::include_graphics(volcano_path)
} else {
  cat("Volcano plot not available for this comparison")
}
```


**Heatmap**

- Description: a heatmap of **z-score normalized** read counts data with application of 
  **hierarchical clustering** by both samples and expression levels
- X-axis: samples
- Y-axis: genes
- Color-scale: z-score normalized read counts 


```{r heatmap-%d}
# Display heatmap from file
heatmap_path <- file.path(out_dirs$heatmap, paste0("%s_heatmap.png"))
if (file.exists(heatmap_path)) {
  knitr::include_graphics(heatmap_path)
} else {
  cat("Heatmap not available for this comparison")
}
```


**Top Differentially Expressed Genes**

Table of the top differential expressed genes with both nominal (**pvalue**)
and adjusted pvalues (**padj**).

```{r top-degs-%d}
# Display top DEGs table
if (!is.null(results[[%d]]) && !is.null(results[[%d]]$deseq)) {
  top_up <- results[[%d]]$deseq %%>%%
    filter(!is.na(padj) & padj < 0.05 & log2FoldChange >= 0.58) %%>%%
    arrange(padj) %%>%%
    head(20)
  
  top_down <- results[[%d]]$deseq %%>%%
    filter(!is.na(padj) & padj < 0.05 & log2FoldChange <= -0.58) %%>%%
    arrange(padj) %%>%%
    head(20)
  
  if (nrow(top_up) > 0) {
    DT::datatable(top_up %%>%% select(ENSEMBL_ID, SYMBOL, log2FoldChange, pvalue, padj, GENENAME),
                 caption = "Top upregulated genes")
  } else {
    cat("No upregulated genes found\\n\\n")
  }
  
  if (nrow(top_down) > 0) {
    DT::datatable(top_down %%>%% select(ENSEMBL_ID, SYMBOL, log2FoldChange, pvalue, padj, GENENAME),
                 caption = "Top downregulated genes")
  } else {
    cat("No downregulated genes found\\n\\n")
  }
} else {
  cat("No differential expression results available for this comparison")
}
```


**Gene Set Enrichment Analysis**

- Description: Gene Set Enrichment Analysis using (GO) terms with both a bubble plot and
 table. 
- X-axis: the gene ratio (# genes related to Gene Set / total number of significant genes) to which
the term is enriched.
- Y-axis: a given Gene Set

```{r gsea-%d}
# Display GSEA results from file
gsea_path <- file.path(out_dirs$gsea, paste0("%s_GSEA.png"))
if (file.exists(gsea_path)) {
  knitr::include_graphics(gsea_path)
} else {
  cat("GSEA results not available for this comparison")
}
```


```{r gsea-table-%d}
# Display top GSEA results
if (!is.null(results[[%d]]) && !is.null(results[[%d]]$gsea)) {
  gsea_results <- as.data.frame(results[[%d]]$gsea)
  if (nrow(gsea_results) > 0) {
    DT::datatable(gsea_results %%>%% 
                 select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue) %%>%%
                 head(20),
                 caption = "Top enriched gene sets")
  } else {
    cat("No significant enriched gene sets found")
  }
} else {
  cat("No GSEA results available for this comparison")
}
```

',
    comparisons[[i]]$name,  # Comparison name
    comparisons[[i]]$exp,   # Experimental group
    comparisons[[i]]$ctrl,  # Control group
    i,                      # Volcano plot chunk id
    comparisons[[i]]$name,  # Volcano plot filename
    i,                      # Heatmap chunk id
    comparisons[[i]]$name,  # Heatmap filename
    i,                      # Top DEGs chunk id
    i, i,                   # Results index
    i, i,                   # Results index for top up/down
    i,                      # GSEA chunk id
    comparisons[[i]]$name,  # GSEA filename
    i,                      # GSEA table chunk id
    i, i, i                 # Results index for GSEA
    )
    
    rmd_content <- paste0(rmd_content, comparison_section)
  }
  
  # Write the R Markdown file
  rmd_file <- file.path(output_dir, paste0(report_prefix, "_", timestamp, ".Rmd"))
  writeLines(rmd_content, rmd_file)
  
  # Render the R Markdown to HTML
  rmarkdown::render(
    rmd_file, 
    output_file = output_file,
    quiet = FALSE
  )
  
  # Return the path to the generated report
  return(output_file)
}

# EOF
