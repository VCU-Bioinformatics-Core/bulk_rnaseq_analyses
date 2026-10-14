#!/usr/bin/env bash
# Print the tag of the bisrde-base image: env-<12 hex chars of a sha256 over
# Dockerfile.base, environment.yml and warm_msigdb_cache.R, in that order.
# The base image holds everything that depends only on those three files, so
# the tag changes exactly when the image has to be rebuilt.
#
# Consumers: .github/workflows/container.yml (decides whether to build or
# reuse the base) and `just container` (local build).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

if command -v shasum >/dev/null 2>&1; then
    sum() { shasum -a 256; }
else
    sum() { sha256sum; }
fi

hash=$(cat "$here/Dockerfile.base" "$here/environment.yml" "$here/warm_msigdb_cache.R" | sum | cut -c1-12)
printf 'env-%s\n' "$hash"
