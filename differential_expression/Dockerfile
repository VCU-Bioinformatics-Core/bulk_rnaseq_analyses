# bisrDE container image (code layer).
#
# Two images:
#   bisrde-base   conda environment, Quarto setup and the MSigDB gene sets
#                 (Dockerfile.base); tagged env-<hash of its inputs>, rebuilt
#                 only when environment.yml, warm_msigdb_cache.R or
#                 Dockerfile.base changes.
#   bisrde        this file: the base plus the bisrDE package and de.R, built
#                 and published by .github/workflows/container.yml per release.
# Local build (needs docker):   just container
# FROM defaults to the :latest base; the workflow passes the hashed tag instead.
# Build context: differential_expression/  (see .dockerignore for what is sent)

ARG BISRDE_BASE=ghcr.io/vcu-bioinformatics-core/bisrde-base:latest
FROM ${BISRDE_BASE}

# An ARG declared before FROM is out of scope after it.
ARG BISRDE_BASE

COPY --chown=$MAMBA_USER:$MAMBA_USER bisrDE/ /opt/bisrde/src/bisrDE/
COPY --chown=$MAMBA_USER:$MAMBA_USER de.R /opt/bisrde/

# dependencies = FALSE: every dependency comes from the conda solve in the
# base image, so R CMD INSTALL fails loudly if environment.yml is missing one
# instead of pulling an ABI-mismatched copy from CRAN.
RUN Rscript -e 'remotes::install_local("/opt/bisrde/src/bisrDE", upgrade = "never", dependencies = FALSE)' \
    && Rscript -e 'library(bisrDE); cat("bisrDE", as.character(packageVersion("bisrDE")), "\n")' \
    && chmod -R a+rX /opt/bisrde

WORKDIR /tmp

ARG BISRDE_VERSION=dev
LABEL org.opencontainers.image.source="https://github.com/VCU-Bioinformatics-Core/bulk_rnaseq_analyses" \
      org.opencontainers.image.version="${BISRDE_VERSION}" \
      org.opencontainers.image.base.name="${BISRDE_BASE}" \
      org.opencontainers.image.description="bisrDE bulk RNA-seq differential expression pipeline (DESeq2, enrichment, Quarto report)"

ENTRYPOINT ["Rscript", "/opt/bisrde/de.R"]
CMD ["--help"]
