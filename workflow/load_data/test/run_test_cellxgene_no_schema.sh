#!/usr/bin/env bash
# Test that harmonize_metadata succeeds when schema_file is omitted
# (i.e. the dataset is already in the target CELLxGENE schema)
set -e -x

WORKDIR=$(dirname $(dirname $0))
cd $WORKDIR

snakemake --configfile test/configs/cellxgene_no_schema.yaml --rerun-incomplete --use-conda --printshellcmds $@

conda run -n scanpy python test/run_assertions.py --live-stream
