#!/usr/bin/env Rscript

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
  
  # Create R Markdown template
  rmd_content <- sprintf('---
title: "RNA-Seq Differential Expression Analysis Report"
author: "Mikail Bala"
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

This report presents the results of differential expression analysis for `r length(comparisons)` comparisons.
For further exploration of the results, please refer to the output directories containing the raw data files.



## Methods

The analysis pipeline used the following methods:
1. **Data preprocessing**: Count data was filtered to remove genes with low expression.
2. **Normalization**: TMM normalization was applied using edgeR.
3. **Differential expression**: DESeq2 was used to identify differentially expressed genes.
4. **Functional analysis**: Gene Set Enrichment Analysis (GSEA) was performed using clusterProfiler.
5. **Thresholds**: Genes with adjusted p-value < 0.05 and |log2FC| >= 0.58 were considered differentially expressed.
6. **Genome annotation**: Mouse (org.Mm.eg.db)



## Sample Exploration

### Principal Component Analysis (PCA)

PCA was performed to visualize the overall pattern of gene expression across samples and to identify potential batch effects or outliers.
**What we expect:** We would expect to see samples that are similar to
each other cluster together, and samples that are different should
cluster separately.

```{r pca-plot}
print(pca_plot)
```

The interactive PCA plots can be found in the following files:
- 2D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot.html")`
- 3D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot3D.html")`



## Differential Expression Results

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

# Display the summary table
kable(summary_table, caption = "Summary of Differential Expression Results")
```

### Detailed Results by Comparison

', 
    format(Sys.time(), "%B %d, %Y"),  # Date in header
    analysis_results_path  # Path to results
  )
  
  # Add detailed sections for each comparison
  for (i in seq_along(comparisons)) {
    comparison_section <- sprintf('
## %s

### Experimental design

- **Experimental group**: %s
- **Control group**: %s


**Volcano Plot**: Each dot represents a change in gene expression. X-axis:
log2 fold-change of expression between treatment compared to the control
plotted against the -log10(p-value). The red line indicates
the p-value \\< 0.05. Every point (gene) above that threshold appears to
have statistically significant changes between the two conditions.

Vertical lines indicate 1.5 fold change. Genes highlighted in RED are
up-regulated in the Treatment compared to the control. Genes
highlighted in BLUE are down-regulated in the Treatment compared
to the control. Genes in gray, do not meet the thresholds for both logFC
and p-value.

```{r volcano-%d}
# Display volcano plot from file
volcano_path <- file.path(out_dirs$volcano, paste0("%s_volcano.png"))
if (file.exists(volcano_path)) {
  knitr::include_graphics(volcano_path)
} else {
  cat("Volcano plot not available for this comparison")
}
```



**Heatmap**: a heatmap of z-score normalized read counts data.
The x-axis are samples, the y-axis are genes. The red color
represents the magnitude of standard deviations above the mean for each
read count (i.e., higher expression), and the blue is the magnitude of
standard deviations below the mean (i.e., lower expression). White
indicates that a read count is close to the mean. The dendrogram is
clustering by samples and by RNA expression.

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

### Gene Set Enrichment Analysis
Gene Set Enrichment Analysis of gene ontology (GO) terms.
Dot Plot: This is a graphical representation of Gene Set Enrichment Analysis terms.
Y axis is the "Gene Set" in which these genes are significantly enriched.
X-axis shows the gene ratio (# genes related to Gene Set / total number of significant genes) to which
the term is enriched. This figure is faceted by activated gene sets and suppressed gene sets.

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
  
  # Add conclusion section
  rmd_content <- paste0(rmd_content, '



## Analysis Summary  
- Total number of comparisons analyzed: `r length(comparisons)`
- Total number of differential expressed genes across all comparisons: `r sum(summary_table$Total_DEGs)`
- Comparisons with the highest number of DEGs: `r summary_table$Comparison[which.max(summary_table$Total_DEGs)]` (`r max(summary_table$Total_DEGs)` DEGs)


*Report generated on %s*
', format(Sys.time(), "%B %d, %Y %H:%M:%S"))  # Date in footer
  
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