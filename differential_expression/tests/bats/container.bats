#!/usr/bin/env bats
#
# container_base_tag.sh: the content-hash tag of the bisrde-base image.
#
# Run with:  just test-sh     (or: bats tests/bats)

setup() {
    DE="$( cd "$( dirname "$BATS_TEST_FILENAME" )/../.." && pwd -P )"
    export DE
    TMP="$( cd "$(mktemp -d)" && pwd -P )"
    export TMP
}

teardown() {
    rm -rf "$TMP"
}

@test "base tag is env- followed by 12 hex characters" {
    run bash "$DE/container_base_tag.sh"
    [ "$status" -eq 0 ]
    [[ "$output" =~ ^env-[0-9a-f]{12}$ ]]
}

@test "base tag is identical on two runs and from another directory" {
    first="$(bash "$DE/container_base_tag.sh")"
    second="$(bash "$DE/container_base_tag.sh")"
    elsewhere="$(cd "$TMP" && bash "$DE/container_base_tag.sh")"
    [ "$first" = "$second" ]
    [ "$first" = "$elsewhere" ]
}

@test "changing environment.yml changes the base tag" {
    cp "$DE/container_base_tag.sh" "$DE/Dockerfile.base" "$DE/environment.yml" \
       "$DE/warm_msigdb_cache.R" "$TMP/"
    before="$(bash "$TMP/container_base_tag.sh")"
    printf '  - r-extra\n' >> "$TMP/environment.yml"
    after="$(bash "$TMP/container_base_tag.sh")"
    [ "$before" != "$after" ]
}
