# bisrDE

Bulk RNA-Seq Differential Expression Pipeline (BISR).

DESeq2-based bulk RNA-seq differential expression pipeline for the Bioinformatics Shared Resource at VCU Massey Comprehensive Cancer Center.

## What it does

- Per-comparison differential expression via DESeq2.
- Sample-level QC: PCA (static + 2D plotly + 3D plotly), Spearman correlation heatmap on DE genes, vst sample-distance heatmap, per-sample library size + detected genes barplot, hierarchical clustering + log-CPM density.
- Per-comparison enrichment: GO (gseGO) + KEGG (gseKEGG) + Reactome (gsePathway) + MSigDB Hallmark (GSEA on `msigdbr` gene sets).
- Volcano + heatmap (all-significant + top-100-by-padj) per comparison.
- Parameterised Quarto report (HTML, self-contained), with per-plot interpretation callouts and a manuscript-ready Methods section.

Designed to plug in downstream of [nf-core/rnaseq](https://nf-co.re/rnaseq) merged-counts output.

## Status

In active development as part of a v1.4.0 overhaul. See parent repo `architecture_plan.md` and `tasks.md` for phase progress.

## Installation

```r
# From inside the parent repo, with renv active:
renv::install("./differential_expression/bisrDE")
```

## License

MIT — see [`LICENSE.md`](LICENSE.md).
