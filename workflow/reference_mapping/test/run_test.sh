#!/usr/bin/env bash
set -e -x

WORKDIR=$(dirname $(dirname $0))
cd $WORKDIR

#--snakefile $WORKDIR/Snakefile
snakemake --configfile test/config.yaml --rerun-incomplete --use-conda $@

conda run -n scanpy python test/assertions.py --live-stream
