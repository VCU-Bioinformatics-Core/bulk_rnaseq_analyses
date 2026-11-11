#!/usr/bin/env Rscript

# Required libraries
library(rmarkdown)
library(knitr)
library(dplyr)
library(glue)

# Function to generate automated R Markdown report
# source("report_generator.R"); generate_report('analysis.rds', output_dir = getwd())
generate_report <- function(analysis_results_path, output_dir = "./", report_prefix = "rnaseq_analysis", analyst = "Mikail Bala", brs_ticket = "") {
  # Validate inputs
  if (!file.exists(analysis_results_path)) {
    stop("Analysis results file does not exist:", analysis_results_path)
  }
  # Coerce missing/NA brs_ticket to "" so the YAML title block builder
  # below can treat it as a simple presence/absence test.
  if (is.null(brs_ticket) || is.na(brs_ticket)) brs_ticket <- ""

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

  # Build the YAML title block. When a BRS ticket is supplied, render it as a
  # `subtitle:` line right under the title. When empty, omit the subtitle line
  # entirely so the YAML doesn't carry a stray blank field.
  title_block <- 'title: "RNA-Seq Differential Expression Analysis Report"'
  if (nzchar(brs_ticket)) {
    title_block <- paste0(title_block, '\nsubtitle: "', brs_ticket, '"')
  }

  # Map the `--annotation` flag (human/mouse) to the genome assembly string
  # that nf-core/rnaseq aligns against. Used in the manuscript-ready Methods
  # text so mouse runs don't mention GRCh38.
  genome_assembly <- switch(tolower(as.character(annotation)),
    human = "GRCh38 human primary assembly",
    mouse = "GRCm39 mouse primary assembly",
    paste0(annotation, " primary assembly")
  )

  # Create R Markdown template (beginning section)
  rmd_content <- glue('---
{title_block}
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
    ├── allsamples_PCA_plot.html
    └── allsamples_PCA_plot3D.html
```

## Sample Exploration using PCA

PCA was performed to visualize the overall patterns of gene expression across samples and
to identify potential batch effects or outliers. **What we expect:** Samples with similar
expression profiles should cluster together, while dissimilar samples are expected to
separate into distinct clusters.

> **How to read this:** Each dot is a sample, projected onto the two (or three) axes (PC1, PC2, PC3) that capture the largest sources of variation across all genes. Samples that are biologically similar should sit close together, and replicates of the same group should cluster. A sample that is far from its group is a candidate outlier — suspect a sample swap, technical artefact, or batch effect, and investigate before trusting the downstream DE results.

### Static PCA (2D)

```{{r pca-static, fig.width=6, fig.height=4}}
print(pca_plot)
```

### Interactive PCA (2D)

*Hover a point for sample ID; drag to zoom; double-click to reset.*

```{{r pca-interactive-2d}}
pca_plotly
```

### Interactive PCA (3D)

*Drag to rotate; scroll to zoom.*

```{{r pca-interactive-3d}}
pca_3d
```

Standalone interactive PCA files are also produced for sharing or embedding elsewhere:

- 2D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot.html")`
- 3D PCA: `r file.path(out_dirs$pca, "allsamples_PCA_plot3D.html")`


### Sample-level QC plots

These four sample-level QC plots complement the PCA above. Use them to spot outlier samples (poor depth, atypical distributions, mis-clustered replicates) before drawing conclusions from the differential-expression results.

> **How to read these:** Each plot answers a different question — "do replicates correlate?", "do samples cluster by group?", "is library depth uneven?", "do all samples have similar expression distributions?". A single discordant sample across multiple panels is a strong outlier signal.

**Sample × Sample correlation heatmap (Spearman, DE genes)**

```{{r qc-corr, out.width="90%", out.height="90%"}}
qc_corr_path <- file.path("figures/qc", "qc_correlation_heatmap.png")
if (file.exists(file.path(dirname(knitr::current_input()), qc_corr_path))) {{
  knitr::include_graphics(qc_corr_path)
}} else {{
  cat("Correlation heatmap not available (no DE genes / <2 samples).")
}}
```

**Sample × Sample Euclidean distance (vst-transformed counts)**

```{{r qc-vst-dist, out.width="90%", out.height="90%"}}
qc_vst_path <- file.path("figures/qc", "qc_vst_dist_heatmap.png")
if (file.exists(file.path(dirname(knitr::current_input()), qc_vst_path))) {{
  knitr::include_graphics(qc_vst_path)
}} else {{
  cat("vst distance heatmap not available.")
}}
```

**Library size + detected genes per sample**

```{{r qc-libsize, out.width="100%", out.height="100%"}}
qc_lib_path <- file.path("figures/qc", "qc_libsize_detected.png")
if (file.exists(file.path(dirname(knitr::current_input()), qc_lib_path))) {{
  knitr::include_graphics(qc_lib_path)
}} else {{
  cat("Library size barplot not available.")
}}
```

**Hierarchical clustering (Ward.D2) + per-sample log-CPM density**

```{{r qc-hclust, out.width="100%", out.height="100%"}}
qc_hcd_path <- file.path("figures/qc", "qc_hclust_density.png")
if (file.exists(file.path(dirname(knitr::current_input()), qc_hcd_path))) {{
  knitr::include_graphics(qc_hcd_path)
}} else {{
  cat("Hierarchical clustering plot not available.")
}}
```


```{{r results-summary}}
# Create a summary table for all comparisons
summary_table <- data.frame(
  Comparison = character(),
  Experimental = character(),
  Control = character(),
  Total_DEGs = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  GO = character(),
  KEGG = character(),
  Reactome = character(),
  Hallmark = character(),
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

      yn <- function(x) ifelse(!is.null(x), "Yes", "No")

      summary_table <- rbind(summary_table, data.frame(
        Comparison = comparisons[[i]]$name,
        Experimental = comparisons[[i]]$exp,
        Control = comparisons[[i]]$ctrl,
        Total_DEGs = total_degs,
        Upregulated = up_degs,
        Downregulated = down_degs,
        GO = yn(results[[i]]$gsea),
        KEGG = yn(results[[i]]$kegg),
        Reactome = yn(results[[i]]$reactome),
        Hallmark = yn(results[[i]]$hallmark),
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
down this count accordingly. The **GO**, **KEGG**, **Reactome**, and **Hallmark** columns
indicate whether each enrichment backend produced any results for that comparison
(`Yes`) or returned no enriched gene sets (`No`).

```{{r display-summary-table}}
# Display the summary table
kable(summary_table, caption = "")
```\n')


  # Add detailed sections for each comparison
  for (i in seq_along(comparisons)) {
    name <- comparisons[[i]]$name
    exp <- paste(comparisons[[i]]$exp, collapse = ", ")
    ctrl <- paste(comparisons[[i]]$ctrl, collapse = ", ")

    # `.trim = FALSE` preserves the leading blank line so the heading
    # starts a new markdown block (otherwise pandoc treats `### name`
    # as continuation of the previous paragraph).
    comparison_section <- glue(.trim = FALSE, '

### {name}

**Experimental design**

- **Experimental group**: {exp}
- **Control group**: {ctrl}


**Volcano Plot**

> **How to read this:** Each dot is a gene. The X-axis shows the magnitude and direction of expression change (right = up in the experimental group, left = down). The Y-axis shows statistical significance after multiple-testing correction — `-log10(padj)`, so higher = more confident. The **red dashed line** marks the FDR < 0.05 threshold; the **black dashed lines** mark a 1.5-fold change. Red dots are up-regulated and blue dots are down-regulated genes that pass both thresholds; grey dots fail at least one. The strongest hits sit in the top-right and top-left corners.

```{{r volcano-{i}, out.width="80%", out.height="80%" }}
# Display volcano plot from file
volcano_path <- file.path("figures/volcano", paste0("{name}_volcano.png"))
if (file.exists(file.path(dirname(knitr::current_input()), volcano_path))) {{
  knitr::include_graphics(volcano_path)
}} else {{
  cat("Volcano plot not available for this comparison")
}}
```


**Heatmaps**

> **How to read this:** Rows are genes, columns are samples (only samples in this comparison’s experimental and control groups are shown). Cell colour is the **z-score** — how many standard deviations a sample’s expression sits above (red) or below (blue) the mean for that gene. The trees on the top and side cluster samples and genes by similarity; samples from the same group should cluster together if the comparison is biologically real. Two heatmaps are produced: **all significant DEGs** (broad view) and **top 100 by adjusted p-value** with gene-symbol labels (focused view).


```{{r heatmap-{i} }}
# Display both heatmaps from file: all-sig (broad) and top-100 (focused).
hm_all <- file.path("figures/heatmap", paste0("{name}_heatmap_all_sig.png"))
hm_top <- file.path("figures/heatmap", paste0("{name}_heatmap_top100.png"))
hm_paths <- c(hm_all, hm_top)
hm_paths <- hm_paths[file.exists(file.path(dirname(knitr::current_input()), hm_paths))]
if (length(hm_paths) > 0) {{
  knitr::include_graphics(hm_paths)
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


**Pathway / Gene-Set Enrichment**

Four complementary enrichment analyses are run per comparison. Each is independent and uses the same ranked log2FC vector, so disagreement between them often signals a meaningful biological difference (a process specific to one ontology vs. one shared across many).

> **How to read these dotplots:** Each dot is a gene set. The X-axis is the **gene ratio** — fraction of the set’s genes that appear in your ranked list. Dot **size** scales with absolute set size; **colour** encodes adjusted p-value (deeper red = more significant). When the plot splits into **Activated** vs **Suppressed** panels, those reflect sets enriched among up- vs down-regulated genes. Use these to nominate biological themes; combine with manual review and IPA before drawing strong conclusions.

**GO Term Enrichment (gseGO)** — Gene Ontology Biological Process / Molecular Function / Cellular Component

```{{r gsea-{i}, out.width="100%", out.height="100%"}}
gsea_path <- file.path("figures/gsea", paste0("{name}_GSEA.png"))
if (file.exists(file.path(dirname(knitr::current_input()), gsea_path))) {{
  knitr::include_graphics(gsea_path)
}} else {{
  cat("GO enrichment not available for this comparison")
}}
```

```{{r gsea-table-{i} }}
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
                 caption = "Top enriched GO terms")
  }} else {{
    cat("No significant enriched GO terms found")
  }}
}} else {{
  cat("No GO enrichment available for this comparison")
}}
```


**KEGG Pathway Enrichment (gseKEGG)** — curated metabolic / signalling pathways from the KEGG database

```{{r kegg-{i}, out.width="100%", out.height="100%"}}
kegg_path <- file.path("figures/kegg", paste0("{name}_KEGG.png"))
if (file.exists(file.path(dirname(knitr::current_input()), kegg_path))) {{
  knitr::include_graphics(kegg_path)
}} else {{
  cat("KEGG enrichment not available for this comparison")
}}
```

```{{r kegg-table-{i} }}
if (!is.null(results[[{i}]]) && !is.null(results[[{i}]]$kegg)) {{
  kegg_results <- as.data.frame(results[[{i}]]$kegg) %>%
      mutate(enrichmentScore=formatC(enrichmentScore, format="e", digits=2),
             NES=formatC(NES, format="e", digits=2),
             p.adjust=formatC(p.adjust, format="e", digits=2),
             qvalue=formatC(qvalue, format="e", digits=2))
  if (nrow(kegg_results) > 0) {{
    DT::datatable(kegg_results %>%
                 select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue) %>%
                 head(20),
                 caption = "Top enriched KEGG pathways")
  }} else {{
    cat("No significant enriched KEGG pathways found")
  }}
}} else {{
  cat("No KEGG enrichment available for this comparison")
}}
```


**Reactome Pathway Enrichment (gsePathway)** — biological pathways from the Reactome knowledgebase

```{{r reactome-{i}, out.width="100%", out.height="100%"}}
reactome_path <- file.path("figures/reactome", paste0("{name}_Reactome.png"))
if (file.exists(file.path(dirname(knitr::current_input()), reactome_path))) {{
  knitr::include_graphics(reactome_path)
}} else {{
  cat("Reactome enrichment not available for this comparison")
}}
```

```{{r reactome-table-{i} }}
if (!is.null(results[[{i}]]) && !is.null(results[[{i}]]$reactome)) {{
  reactome_results <- as.data.frame(results[[{i}]]$reactome) %>%
      mutate(enrichmentScore=formatC(enrichmentScore, format="e", digits=2),
             NES=formatC(NES, format="e", digits=2),
             p.adjust=formatC(p.adjust, format="e", digits=2),
             qvalue=formatC(qvalue, format="e", digits=2))
  if (nrow(reactome_results) > 0) {{
    DT::datatable(reactome_results %>%
                 select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue) %>%
                 head(20),
                 caption = "Top enriched Reactome pathways")
  }} else {{
    cat("No significant enriched Reactome pathways found")
  }}
}} else {{
  cat("No Reactome enrichment available for this comparison")
}}
```


**MSigDB Hallmark Gene Sets (GSEA + msigdbr H)** — 50 well-curated gene sets representing distinct biological states / processes

```{{r hallmark-{i}, out.width="100%", out.height="100%"}}
hallmark_path <- file.path("figures/hallmark", paste0("{name}_Hallmark.png"))
if (file.exists(file.path(dirname(knitr::current_input()), hallmark_path))) {{
  knitr::include_graphics(hallmark_path)
}} else {{
  cat("Hallmark enrichment not available for this comparison")
}}
```

```{{r hallmark-table-{i} }}
if (!is.null(results[[{i}]]) && !is.null(results[[{i}]]$hallmark)) {{
  hallmark_results <- as.data.frame(results[[{i}]]$hallmark) %>%
      mutate(enrichmentScore=formatC(enrichmentScore, format="e", digits=2),
             NES=formatC(NES, format="e", digits=2),
             p.adjust=formatC(p.adjust, format="e", digits=2),
             qvalue=formatC(qvalue, format="e", digits=2))
  if (nrow(hallmark_results) > 0) {{
    DT::datatable(hallmark_results %>%
                 select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue) %>%
                 head(20),
                 caption = "Top enriched Hallmark gene sets")
  }} else {{
    cat("No significant enriched Hallmark gene sets found")
  }}
}} else {{
  cat("No Hallmark enrichment available for this comparison")
}}
```

**Next steps:** Upload the per-comparison `DESeq2_{name}.csv` to **QIAGEN IPA** (see *Utilizing IPA* below) for full-fidelity pathway exploration. Reuse the **Manuscript-Ready Text** section verbatim in your manuscript Methods.
')
    
    rmd_content <- paste0(rmd_content, comparison_section)
  }


extra_content <- glue('\n
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
Raw RNA-Seq fastq files were processed by the VCU Massey Comprehensive Cancer Center Bioinformatics Shared Resource (BISR) using the NextFlow nf-core/rnaseq v3.18.0 pipeline [1]. Briefly, this pipeline assesses sequencing quality using FastQC v 0.12.1 [2] before and after trimming, performs adaptor trimming with Trim Galore! v0.6.10 [3], and aligns sequencing reads to the {genome_assembly} reference genome using STAR v 2.7.11b [4] with transcriptome quantification by Salmon v1.10.3 [5].  Pipeline output includes gene expression raw count data and a comprehensive QC report compiled by MultiQC v1.25.1 [6].
Differential expression analysis was performed using DESeq2 v 1.44.0 [7]. Lowly expressed genes were filtered out per DESeq2 methods [7] prior to normalization and differential expression testing. Significance was calculated using the Wald-test and adjusted using Benjamini Hochberg False Discovery Rate (FDR). Volcano plots and heatmaps were generated using the EdgeR TMM normalized count data and visualized using R packages. Significant differentially expressed genes (DEGs) are defined as those with an FDR<0.05 and absolute fold-change of 1.5 (log2 fold-change = 0.58) or greater. Gene Set Enrichment Analysis (GSEA) [8] for Gene Ontology terms (GO) was performed using the clusterProfiler package [9] across all genes, regardless of significance. All computational analyses were performed on VCU’s High Performance Research Computing cluster.

### References
1) Ewels P, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Feb 13. doi:10.1038/s41587-020-0439-x

2) Andrews S. FastQC: A Quality Control Tool for High Throughput Sequence Data. Babraham Bioinformatics; 2010. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

3) Krueger F. Trim Galore! v0.6.10. 2023. https://github.com/FelixKrueger/TrimGalore

4) Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi:10.1093/bioinformatics/bts635

5) Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017;14:417-419. doi:10.1038/nmeth.4197

6) Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32:3047-3048. doi:10.1093/bioinformatics/btw354

7) Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15:550. doi:10.1186/s13059-014-0550-8

8) Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005;102(43):15545-15550. doi:10.1073/pnas.0506580102

9) Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16(5):284-287. doi:10.1089/omi.2011.0118

### Required Acknowledgements

Please include the following statements in your acknowledgements manuscript section:

- “Services in support of the research project were provided by the VCU Massey Comprehensive Cancer Center Bioinformatics Shared Resource. Massey is supported, in part, with funding from NIH-NCI Cancer Center Support Grant P30 CA016059.”

- “High Performance Computing resources provided by the High Performance Research Computing (HPRC) core facility at Virginia Commonwealth University (https://hprc.vcu.edu) were used for conducting the research reported in this work.”
')
  rmd_content <- paste0(rmd_content, extra_content)
  

  # Write the R Markdown file
  fn <- glue("{report_prefix}_{timestamp}.Rmd")
  rmd_file <- file.path(output_dir, fn)

  writeLines(rmd_content, rmd_file)

  # Render the R Markdown to HTML
  rmarkdown::render(rmd_file, output_file = output_file, quiet = FALSE)

  # Return the path to the generated report
  return(output_file)
}
