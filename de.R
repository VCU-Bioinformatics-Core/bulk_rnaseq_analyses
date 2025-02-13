# Package installation and loading
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")

if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(here,
               dplyr,
               data.table,
               tidyverse,
               janitor,
               scales,
               ggrepel,
               clusterProfiler,
               enrichplot,
               tidyverse,
               readr,
               DT,
               DESeq2,
               edgeR,
               ggplot2,
               AnnotationDbi,
               "org.Mm.eg.db",
               gplots,
               RColorBrewer,
               purrr,
               plotly,
               stats,
               orca,
               reticulate,
               optparse,
               htmlwidgets
)





# HELPER FUNCTIONS

#' @description Create a file path by combining directory, prefix, name and extension
#'
#' @param base_dir Base directory path
#' @param prefix Prefix to add before the name
#' @param name Main file name
#' @param extension File extension (defaults to ".csv")
#' @return A complete file path string
#' @export
create_file_path <- function(base_dir, prefix, name, extension = ".csv") {
  file.path(base_dir, paste0(prefix, name, extension))
}

#' @description Create a comparison name string
#'
#' @param exp Experimental group name
#' @param ctrl Control group name
#' @param prefix Optional prefix to add before the comparison (defaults to "")
#' @return A formatted comparison string
#' @export
create_comparison_name <- function(exp, ctrl, prefix = "") {
  paste0(prefix, exp, " vs. ", ctrl)
}

#' @description Save a ggplot object to a file
#'
#' @param plot ggplot object to save
#' @param filename Output file path
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @export
save_plot <- function(plot, filename, width = 15, height = 17) {
  ggsave(plot, filename = filename, width = width, height = height, bg = "white")
}


#' @description Export a plotly figure to HTML
#'
#' @param plotly_obj Plotly object to be exported
#' @param file_path Path where the HTML file should be saved
#' @details Saves an interactive plotly visualization as a self-contained HTML file
#' that can be opened in a web browser. Includes error handling for invalid inputs.
#' @export
export_plotly_to_html <- function(plotly_obj, file_path) {

  tryCatch({
    # Check if the plotly object is valid
    if (!inherits(plotly_obj, "plotly")) {
      stop("Invalid plotly object")
    }

    # Save the plotly object to an html file
    htmlwidgets::saveWidget(plotly_obj, file_path, selfcontained = TRUE)

  }, error = function(e) {
    cat(paste0("Error: ", e$message, "\n"))
  })
}





# ANALYSIS FUNCTIONS

#' @description Perform differential expression analysis using DESeq2
#'
#' @param dds DESeqDataSet object
#' @param exp Experimental group name
#' @param ctrl Control group name
#' @return Data frame containing DESeq2 results
#' @details This function performs differential expression analysis between two groups
#' using DESeq2. It includes error checking for group existence and proper re-leveling
#' of the condition factor.
#' @export
perform_deseq2_analysis <- function(dds, exp, ctrl) {
  tryCatch({
    print(paste("Performing DESeq2 analysis for", exp, "vs", ctrl))
    print("Current condition levels:")
    print(levels(dds$condition))
    
    # Verify groups exist in the data
    if(!(exp %in% levels(dds$condition))) {
      stop(paste("Experimental group", exp, "not found in condition levels"))
    }
    if(!(ctrl %in% levels(dds$condition))) {
      stop(paste("Control group", ctrl, "not found in condition levels"))
    }
    
    # Make sure the reference level is set correctly
    dds$condition <- relevel(dds$condition, ref = ctrl)
    
    # Run DESeq
    print("Running DESeq...")
    dds <- DESeq(dds)
    
    # Get results
    print("Getting results...")
    res <- results(dds, contrast = c("condition", exp, ctrl),
                   cooksCutoff = TRUE, independentFiltering = FALSE)
    
    print(paste("Analysis complete. Number of results:", nrow(res)))
    return(as.data.frame(res))
  }, error = function(e) {
    print(paste("Error in DESeq2 analysis:", e$message))
    stop(e)
  })
}

#' @description Annotate DESeq2 results with gene symbols and names
#'
#' @param results DESeq2 results data frame
#' @return Annotated data frame with gene symbols and names
#' @details Uses the org.Mm.eg.db package to add gene symbols and names to the
#' DESeq2 results based on ENSEMBL IDs
#' @export
annotate_results <- function(results) {
  gene_ids <- rownames(results)
  annotations <- AnnotationDbi::select(org.Mm.eg.db, 
    keys = gene_ids, 
    keytype = "ENSEMBL", 
    columns = c("SYMBOL", "GENENAME")
  )
  merge(as.data.frame(results),
        annotations, 
        by.x = "row.names",
        by.y = "ENSEMBL", 
        all.x = TRUE) %>%
    rename(ENSEMBL_ID = Row.names)
}

#' @description Generate a volcano plot from differential expression results
#'
#' @param data Annotated DESeq2 results data frame
#' @param exp_name Experimental group name
#' @param ctrl_name Control group name
#' @param p P-value threshold for significance
#' @param lfc Log2 fold change threshold
#' @param sig Column name for significance values (defaults to "padj")
#' @param out_dir Output directory for saving the plot
#' @return ggplot object containing the volcano plot
#' @export
generate_volcano <- function(data, exp_name, ctrl_name, p = 0.05, lfc = 0.58, 
                           sig = "padj", out_dir) {
  data %>% 
    mutate(
      color_tag = case_when(
        eval(as.symbol(sig)) < p & log2FoldChange < -lfc ~ "Under expressed",
        eval(as.symbol(sig)) < p & log2FoldChange > lfc ~ "Over expressed",
        TRUE ~ NA_character_
      ),
      highlight = ifelse(pvalue <= p & abs(log2FoldChange) >= lfc, SYMBOL, NA)
    ) %>% 
    ggplot(aes(x = log2FoldChange, 
               y = -log10(eval(as.symbol(sig))), 
               color = color_tag,
               label = ifelse(highlight == TRUE, SYMBOL, NA))) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    geom_label_repel(aes(label = highlight)) +
    scale_color_manual(values = c("firebrick", "steelblue")) +
    geom_hline(yintercept = -log10(p), col = "red", linetype = 2) +
    geom_vline(xintercept = c(-lfc, lfc)) +
    theme(legend.title = element_blank()) +
    labs(
      x = "Log2 Fold-Change (FC)",
      y = paste0("-log10(", sig, ")"),
      title = create_comparison_name(exp_name, ctrl_name, "Differentially expressed genes - ")
    )
}

#' @description Generate a heatmap from differential expression results
#'
#' @param results_df Data frame containing DESeq2 results
#' @param normalized_counts Matrix of normalized count data
#' @param p P-value threshold for significance (default: 0.05)
#' @param lfc Log2 fold change threshold (default: 0.58)
#' @param exp_name Experimental group name
#' @param ctrl_name Control group name
#' @param fig_dir Directory to save the figure
#' @return A heatmap visualization of differentially expressed genes
#' @details Creates a heatmap using filtered data based on p-value and log fold change thresholds.
#' The values are z-score normalized and displayed using a blue-white-red color scheme.
#' @export
generate_heatmap <- function(results_df, normalized_counts, p = 0.05, lfc = 0.58, 
                           exp_name, ctrl_name, fig_dir) {
  filtered_data <- results_df %>% 
    dplyr::filter((padj < p & abs(log2FoldChange) >= lfc))
  
  values <- normalized_counts[rownames(filtered_data), ] %>%
    as.matrix() %>%
    jitter(factor = 1, amount = 0.00001)
  
  zscores <- t(scale(t(values)))
  
  heatmap.2(zscores,
            col = colorRampPalette(c("blue", "white", "firebrick"))(20),
            density.info = "none",
            dendrogram = "both",
            Colv = TRUE,
            trace = "none",
            margins = c(10, 10),
            labRow = NA)
}

#' @description Process Gene Set Enrichment Analysis (GSEA)
#'
#' @param result DESeq2 results data frame
#' @param p P-value threshold for significance (default: 0.05)
#' @param lfc Log2 fold change threshold (default: 0.58)
#' @return GSEA results object containing enriched GO terms
#' @details Performs GSEA analysis on significant genes using GO terms.
#' Uses org.Mm.eg.db for mouse gene annotations and filters based on
#' provided p-value and log fold change thresholds.
#' @export
process_gsea <- function(result, p = 0.05, lfc = 0.58) {
  sig_genes <- result %>%
    dplyr::filter(padj < p & abs(log2FoldChange) >= lfc) %>%
    rownames_to_column("gene") %>%
    arrange(desc(log2FoldChange))
  
  gene_list <- sig_genes$log2FoldChange
  names(gene_list) <- sig_genes$gene
  
  gseGO(geneList = gene_list,
        ont = "ALL",
        minGSSize = 1,
        maxGSSize = 900,
        keyType = "ENSEMBL",
        pvalueCutoff = 0.05,
        pAdjustMethod = "fdr",
        OrgDb = org.Mm.eg.db) %>%
    setReadable(OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
}

#' @description Create a dotplot visualization of GSEA results
#'
#' @param gse GSEA results object
#' @param title Plot title
#' @return A ggplot object showing enriched terms as a dotplot
#' @details Visualizes the top 15 enriched GO terms, split by direction of change.
#' The size of dots represents gene count and color represents significance.
#' @export
create_dotplot <- function(gse, title) {
  dotplot(gse,
          showCategory = 15,
          title = title,
          split = ".sign",
          orderBy = "p.adjust",
          label_format = 31,
          font.size = 9) +
    facet_grid(.~.sign)
}





# Main analysis pipeline

#' @description Run complete differential expression analysis pipeline
#'
#' @param comparison List containing comparison details (name, exp, ctrl)
#' @param dds DESeqDataSet object
#' @param normalized_counts Matrix of normalized count data
#' @param out_dirs List of output directories
#' @return List containing DESeq2 and GSEA results
#' @details Performs complete analysis workflow including:
#' - DESeq2 differential expression analysis
#' - Result annotation
#' - Volcano plot generation
#' - Heatmap generation
#' - GSEA analysis and visualization
#' @export
run_analysis <- function(comparison, dds, normalized_counts, out_dirs) {
  tryCatch({
    print(paste("\nStarting analysis for comparison:", comparison$name))
    print(paste("Experimental group:", comparison$exp))
    print(paste("Control group:", comparison$ctrl))
    
    # Perform DESeq2 analysis
    deseq_results <- perform_deseq2_analysis(dds, comparison$exp, comparison$ctrl)
    if(is.null(deseq_results)) {
      stop("DESeq2 analysis returned NULL results")
    }
    
    # Annotate results
    print("Annotating results...")
    annotated_results <- annotate_results(deseq_results)
    
    # Save DESeq2 results
    output_file <- create_file_path(out_dirs$de_data, "DESeq2_", comparison$name)
    print(paste("Saving results to:", output_file))
    write.csv(annotated_results, output_file)
    
    # Generate and save volcano plot
    print("Generating volcano plot...")
    volcano_plot <- generate_volcano(annotated_results, 
                                     comparison$exp, 
                                     comparison$ctrl)
    save_plot(volcano_plot, 
              create_file_path(out_dirs$volcano, "", comparison$name, "_volcano.png"))
    
    # Generate and save heatmap
    print("Generating heatmap...")
    png(create_file_path(out_dirs$heatmap, "", comparison$name, "_heatmap.png"),
        width = 800, height = 1200, res = 150)
    generate_heatmap(deseq_results, 
                     normalized_counts, 
                     exp_name = comparison$exp, 
                     ctrl_name = comparison$ctrl)
    dev.off()
    
    # Process and save GSEA results
    print("Processing GSEA...")
    gse <- tryCatch({
      process_gsea(deseq_results)
    }, error = function(e) {
      print(paste("GSEA processing error:", e$message))
      return(NULL)
    })
    
    if(!is.null(gse)) {
      write.csv(as.data.frame(gse), 
                create_file_path(out_dirs$gsea_data, "GO_Analysis_", comparison$name))
      
      # Generate and save GSEA plot
      print("Generating GSEA plot...")
      gsea_plot <- create_dotplot(gse, 
                                  create_comparison_name(comparison$exp, 
                                                         comparison$ctrl, 
                                                         "GSEA "))
      save_plot(gsea_plot, 
                create_file_path(out_dirs$gsea, "", comparison$name, "_GSEA.png"))
    }
    
    return(list(deseq = annotated_results, gsea = gse))
    
  }, error = function(e) {
    print(paste("Error in run_analysis:", e$message))
    return(NULL)
  })
}

#' @description Set up analysis output directories
#'
#' @param base_dir Base directory path where subdirectories will be created
#' @return List of created directory paths
#' @details Creates the following subdirectories:
#' - output/de_data: For differential expression results
#' - output/gsea_data: For GSEA results
#' - figures/volcano: For volcano plots
#' - figures/heatmap: For heatmap plots
#' - figures/gsea: For GSEA plots
#' - figures/pca: For PCA plots
#' @export
setup_directories <- function(base_dir) {
  dirs <- list(
    output = file.path(base_dir, "output"),
    figures = file.path(base_dir, "figures"),
    de_data = file.path(base_dir, "output/de_data"),
    gsea_data = file.path(base_dir, "output/gsea_data"),
    volcano = file.path(base_dir, "figures/volcano"),
    heatmap = file.path(base_dir, "figures/heatmap"),
    gsea = file.path(base_dir, "figures/gsea"),
    pca = file.path(base_dir, "figures/pca")
  )
  
  walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}


###########################################
###########################################
###########################################
# DATA PREPROCESSING AND METADATA CREATION
###########################################
###########################################
###########################################

# Define command-line flags.
option_list = list (
  make_option(c("-c", "--counts"), type = "character", default = NULL,
                    help = "Required. A path for the merged counts.tsv file"),

  make_option(c("-m", "--contrasts"), type = "character", default = NULL,
                    help = "Required. A path for the contrast matrix.tsv file"),

  make_option(c("-o", "--outdir"), type = "character", default = "./output",
                    help = "path for the output directory to store results. \
                    A new directory will be created with the given path and \
                    name, otherwise a default directory will be created in \
                    the current directory: [default= %default]"),

  make_option(c("-s", "--samplesheet"), type = "character", default = NULL,
                    help = "Required. A path for the samplesheet.csv file"),

  make_option(c("-r", "--runid"), type = "character", default = NULL,
                    help="Required. A unique name for this analysis.", metavar="character")
);

# Parse the arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
runID <- opt$runid
countData <- opt$counts
contrastData <- opt$contrasts
sampleSheet <- opt$samplesheet
outDir <- opt$outdir

# Implementing all hard-coded files here for debugging.
#### Debug Options
debug = FALSE
if(debug){
  runID <- "test_run"
  countData <- "./rsem.merged.gene_counts.tsv"
  contrastData <- "./contrasts.tsv"
  sampleSheet <- "./samplesheet.csv"
  outDir <- "./draft_output"
}

# Check for the commad-line arguments for the script. If no arguments are provided exit and print help.
if (is.null(runID) | is.null(countData) | is.null(contrastData) | is.null(sampleSheet)) {
  print_help(opt_parser)
  stop("All required arguments must be supplied (input file)", call.=FALSE)
}

# Access the arguments
cat("Utilizing this Run Id:", runID, "\n")
cat("merged counts file found:", countData, "\n")
cat("contrast matrix file found:", contrastData, "\n")
cat("sample sheet found:", sampleSheet, "\n")
cat("creating output directory:", outDir, "\n")

# Set up output directories
out_dirs <- setup_directories(outDir)

# Read in raw merged counts, contrasts and sample sheet
counts <- data.frame(fread(countData, header = "auto"), row.names = 1)
contrasts_raw <- read.delim(contrastData, sep="\t", header=TRUE, stringsAsFactors=FALSE)
sample_info <- read_delim(sampleSheet, delim = ",") # for deseq object

# Get all contrast columns (columns after GroupID)
contrast_cols <- colnames(contrasts_raw)[3:ncol(contrasts_raw)]
print("Contrast columns found:")
print(contrast_cols)

# Create the comparisons list
comparisons <- list()
for(i in seq_along(contrast_cols)) {
  contrast_name <- contrast_cols[i]
  contrast_data <- contrasts_raw[[contrast_name]]
  
  print(paste("\nProcessing contrast:", contrast_name))
  print("Contrast data:")
  print(table(contrast_data, useNA="ifany"))
  
  # Find unique groups for exp (1) and ctrl (0)
  exp_group <- unique(contrasts_raw$GroupID[contrast_data == 1 & !is.na(contrast_data)])
  ctrl_group <- unique(contrasts_raw$GroupID[contrast_data == 0 & !is.na(contrast_data)])
  
  print(paste("exp_group found:", paste(exp_group, collapse=", ")))
  print(paste("ctrl_group found:", paste(ctrl_group, collapse=", ")))
  
  # Add to comparisons list if both exp and ctrl are found
  if(length(exp_group) > 0 && length(ctrl_group) > 0) {
    comparisons[[length(comparisons) + 1]] <- list(
      name = contrast_name,
      exp = exp_group,
      ctrl = ctrl_group
    )
  }
}

# Print the results to verify
print("\nFinal comparisons list:")
str(comparisons)

# need to subtract 1 as we will ignore the transcript ids
n <- ncol(counts) - 1 

# clean counts to make sure the columns are just samples
countsdf <- counts %>% select(-transcript_id.s.) %>% mutate_at(1:n, as.integer)

# NORMALIZE DATA
x <- DGEList(counts = countsdf)
dge <- calcNormFactors(x, method = "TMM")
TMM <- cpm(dge)
write.csv(TMM, str_c(out_dirs$de_data, "normalizedCounts_TMM", Sys.Date(), ".csv"))

# DESEQ2 SETUP
sample_info <- data.frame(sample = sample_info$sample, condition = sample_info$condition)

dds <- DESeqDataSetFromMatrix(
  countData = countsdf,
  colData = sample_info,
  design = ~ condition
)

# Pre-filtering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Run analysis for all comparisons
#results <- map(comparisons, ~run_analysis(., dds, TMM, out_dirs))
print("Verifying DDS setup...")
print("Condition levels in DDS:")
print(levels(dds$condition))

print("\nVerifying comparisons:")
str(comparisons)

results <- list()
for(i in seq_along(comparisons)) {
  print(paste("\nProcessing comparison", i, "of", length(comparisons)))
  results[[i]] <- run_analysis(comparisons[[i]], dds, TMM, out_dirs)
  
  # Verify results
  if(is.null(results[[i]])) {
    print(paste("Warning: No results generated for comparison", i))
  } else {
    print(paste("Results generated successfully for comparison", i))
    print(paste("Number of DEGs:", nrow(results[[i]]$deseq)))
  }
}

# Print summary of results
print("\nResults summary:")
for(i in seq_along(results)) {
  print(paste("\nComparison", i, "-", comparisons[[i]]$name))
  if(!is.null(results[[i]])) {
    print(paste("DESeq2 results:", nrow(results[[i]]$deseq), "genes"))
    print(paste("GSEA results:", !is.null(results[[i]]$gsea)))
  } else {
    print("No results generated")
  }
}

# PCA
# NOTE: The prcomp() expects the genes to columns and the sample to rows.
# Since the samples in our data matrix are the columns and the genes are the rows;
##  Removes all zero row from the TMM normalized data. 
data_reduced_zero <- as.data.frame(TMM) %>%
  filter(if_any(where(is.numeric)))

# Generate a PCA Matrix 
pca_matrix <- data_reduced_zero %>% 
  # coerce to a matrix
  as.matrix() %>% 
  # transpose the matrix so that rows = samples and columns = variables
  t()

## perform PCA
data_pca <- prcomp(pca_matrix)
# extract PC scores
pc_scores <- data_pca$x
# add sample names as row names
rownames(pc_scores) <- colnames(data_reduced_zero)

# extract stuff for pca scree plot
pca_var <- data_pca$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 1)

# create a data frame for ggplot
pca_df <- data.frame(Sample=rownames(data_pca$x),
                       X=data_pca$x[,1],
                       Y=data_pca$x[,2])

#add groups
pca_df$Group <- sample_info$condition

# Create PCA plot using ggplot2
pca_plot <- ggplot(data=pca_df, aes(x=X, y=Y, label=Sample, color = Group)) +
  geom_text() +
  #CM: Code Change HERE - edit limits to reflect the x-axis for better visualization of clustering
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep="")) +
  theme_bw()

ggsave(pca_plot, filename = str_c(out_dirs$pca, "/PCA_plot.png"))



# interactive pca
fig <- plot_ly(pca_df, x = ~X, y = ~Y, color = ~Group,
               colors = c('steelblue','firebrick','olivedrab', 'plum'),
               type = 'scatter', mode = 'markers',
               size = 5.7,
               span = 4.9,
               alpha = 2.5,
               alpha_stroke = .5
               ) %>%
  layout(legend=list(title=list(text='color')),
    plot_bgcolor='#e5ecf6',
    xaxis = list(
      title = paste("PC1 - ", pca_var_per[1], "%", sep=""),
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    yaxis = list(
      title = paste("PC2 - ", pca_var_per[2], "%", sep=""),
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff')
    )

# Define variable and save plot
file_name_plotly <- paste0(out_dirs$pca, "/", "allsamples_PCA_plot.pdf")
orca(fig, file_name_plotly)
export_plotly_to_html(fig, file_name_plotly)



#### Interactive 3D PCA plot.
prin_comp <- prcomp(pca_matrix, rank. = 3)
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, pca_df$Group)

tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)

tit = paste('Total Explained Variance = ', tot_explained_variance_ratio)

fig3D <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~pca_df$Group,
               colors = c('steelblue','firebrick','olivedrab', 'plum')) %>%
  add_markers(size = 12)

fig3D <- fig3D %>%
  layout(
    title = tit,
    scene = list(bgcolor = "#e5ecf6")
)

# Define variable and save
file_name_plotlyPCA3D <- paste0(out_dirs$pca, "/", "allsamples_PCA_plot3D.pdf")
orca(fig3D, file_name_plotlyPCA3D)
export_plotly_to_html(fig3D, file_name_plotlyPCA3D)
