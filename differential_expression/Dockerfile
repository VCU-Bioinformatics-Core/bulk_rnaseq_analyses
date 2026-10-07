# bisrDE container image.
#
# Built and published by .github/workflows/container.yml on every v* tag as
#   ghcr.io/vcu-bioinformatics-core/bisrde:<version>   (and :latest)
# The launchers pull it with Apptainer (pull_container.sh), so the image has
# to work under `apptainer exec`, which ignores the Docker ENTRYPOINT and any
# activation script. Everything a run needs is therefore a plain ENV line.
#
# Local build (needs docker):   just container
# Build context: differential_expression/  (see .dockerignore for what is sent)

# Pinned to an exact tag so a rebuild months from now starts from the same
# base. micromamba 2.x changed CLI defaults; 1.5.10 is the last 1.x release
# and is what the environment was solved with.
FROM mambaorg/micromamba:1.5.10

# The micromamba image only puts the env on PATH through its entrypoint.
# Set it here so RUN lines and `apptainer exec` resolve the same Rscript.
ENV PATH=/opt/conda/bin:$PATH

# ---- conda environment (changes rarely; keep above the bisrDE copy) --------
# Same channel rules as environment.yml and `just conda-env`: conda-forge and
# bioconda only, strict priority, never defaults.
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
        --override-channels -c conda-forge -c bioconda --strict-channel-priority \
    && micromamba clean --all --yes \
    && rm /tmp/environment.yml

# The conda-forge quarto package is built with QUARTO_VENDOR_BINARIES=false:
# deno, pandoc, esbuild, typst and sass are separate packages in /opt/conda/bin
# and the share tree is /opt/conda/share/quarto. The /opt/conda/bin/quarto
# launcher only finds them through the variables its activation script
# (/opt/conda/etc/conda/activate.d/quarto.sh) exports; without them it looks
# for a vendored deno under bin/tools and fails. Nothing in this image runs
# activation (RUN lines do not, the ENTRYPOINT is replaced, and apptainer exec
# skips it anyway), so the same values are plain ENV lines. The check below
# fails the build if the activation script ever disagrees with them.
ENV QUARTO_SHARE_PATH=/opt/conda/share/quarto \
    QUARTO_DENO=/opt/conda/bin/deno \
    QUARTO_DENO_DOM=/opt/conda/lib/deno_dom.so \
    QUARTO_PANDOC=/opt/conda/bin/pandoc \
    QUARTO_ESBUILD=/opt/conda/bin/esbuild \
    QUARTO_TYPST=/opt/conda/bin/typst \
    QUARTO_DART_SASS=/opt/conda/bin/sass \
    QUARTO_CONDA_PREFIX=/opt/conda

# List every activation script in the build log, so a future conda package
# that only works after activation is noticed, then require the QUARTO_* ENV
# lines above to equal what the activation scripts export. QUARTO_DENO_DOM is
# compared by existence instead: the quarto feedstock's activation script
# exports a build-machine path for it (/home/conda/feedstock_root/...), so
# the value above comes from deno_dom.sh's DENO_DOM_PLUGIN. Every tool path
# must exist, or the first report render would fail at run time.
RUN for f in /opt/conda/etc/conda/activate.d/*.sh; do echo "== $f"; cat "$f"; done \
    && env | grep '^QUARTO_' | grep -v '^QUARTO_DENO_DOM=' | sort > /tmp/quarto_env.txt \
    && env -i PATH="$PATH" /bin/bash -c 'for f in /opt/conda/etc/conda/activate.d/*.sh; do . "$f"; done; env' \
         | grep '^QUARTO_' | grep -v '^QUARTO_DENO_DOM=' | sort > /tmp/quarto_activate.txt \
    && { diff /tmp/quarto_env.txt /tmp/quarto_activate.txt \
         || { echo "QUARTO_* ENV lines in the Dockerfile (<) differ from the conda activation scripts (>); update the ENV block"; exit 1; }; } \
    && for p in "$QUARTO_DENO" "$QUARTO_DENO_DOM" "$QUARTO_PANDOC" "$QUARTO_ESBUILD" "$QUARTO_TYPST" "$QUARTO_DART_SASS" "$QUARTO_SHARE_PATH"; do \
         test -e "$p" || { echo "QUARTO_* path does not exist in the image: $p"; exit 1; }; \
       done \
    && rm /tmp/quarto_env.txt /tmp/quarto_activate.txt

# ---- bisrDE package, driver, MSigDB cache ---------------------------------
USER root
RUN mkdir -p /opt/bisrde && chown $MAMBA_USER:$MAMBA_USER /opt/bisrde
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER bisrDE/ /opt/bisrde/src/bisrDE/
COPY --chown=$MAMBA_USER:$MAMBA_USER de.R warm_msigdb_cache.R /opt/bisrde/

# Runtime environment. Every variable is a plain ENV so Apptainer sees it:
#   R_USER_CACHE_DIR               where msigdbr looks for its gene-set archive
#   RENV_CONFIG_AUTOLOADER_ENABLED  with CONDA_PREFIX, makes the repo .Rprofile
#   CONDA_PREFIX                    skip renv when the checkout is bind-mounted
#   R_LIBS_USER                     keep a bind-mounted ~/R library from
#                                   shadowing the conda library
#   LANG / LC_ALL                   UTF-8 for the report and the CLI output
ENV R_USER_CACHE_DIR=/opt/bisrde/cache \
    RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
    CONDA_PREFIX=/opt/conda \
    R_LIBS_USER=/opt/conda/lib/R/library \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# dependencies = FALSE: every dependency comes from the conda solve above, so
# R CMD INSTALL fails loudly if environment.yml is missing one instead of
# pulling an ABI-mismatched copy from CRAN. That only covers Imports; the
# Suggests the pipeline relies on (apeglm shrinkage, cluster silhouette,
# optparse CLI) are asserted by name because a missing one degrades silently.
# `quarto check` prints its findings but exits 0 whatever it finds, so a small
# knitr render is what proves the CLI, deno, pandoc and the R engine work.
RUN Rscript -e 'remotes::install_local("/opt/bisrde/src/bisrDE", upgrade = "never", dependencies = FALSE)' \
    && Rscript -e 'library(bisrDE); cat("bisrDE", as.character(packageVersion("bisrDE")), "\n")' \
    && Rscript -e 'for (p in c("apeglm", "cluster", "optparse")) if (!requireNamespace(p, quietly = TRUE)) stop("Suggests package missing from environment.yml: ", p)' \
    && quarto check \
    && mkdir /tmp/qcheck \
    && printf '%s\n' '---' 'title: check' 'format: html' '---' '' '```{r}' '1 + 1' '```' > /tmp/qcheck/check.qmd \
    && quarto render /tmp/qcheck/check.qmd --to html \
    && test -s /tmp/qcheck/check.html \
    && rm -rf /tmp/qcheck

# Pre-download the MSigDB gene sets so Hallmark enrichment works on compute
# nodes without network. The script exits 1 if the download fails.
RUN Rscript /opt/bisrde/warm_msigdb_cache.R \
    && chmod -R a+rX /opt/bisrde

WORKDIR /tmp

ARG BISRDE_VERSION=dev
LABEL org.opencontainers.image.source="https://github.com/VCU-Bioinformatics-Core/bulk_rnaseq_analyses" \
      org.opencontainers.image.version="${BISRDE_VERSION}" \
      org.opencontainers.image.description="bisrDE bulk RNA-seq differential expression pipeline (DESeq2, enrichment, Quarto report)"

ENTRYPOINT ["Rscript", "/opt/bisrde/de.R"]
CMD ["--help"]
