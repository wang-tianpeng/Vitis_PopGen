#!/usr/bin/env bash
set -euo pipefail
set -e


FASTA=./VITVarB40-14_v2.0_hap1_filtered.4distes.min10.fasta

seqkit grep -f samples.outgroup.discard -v $FASTA > $FASTA.outgroup.discard.fasta

iqtree -s $FASTA.outgroup.discard.fasta -st DNA -m TVM+F+R7 -nt 20 -B 1000 -pre $FASTA.outgroup.discard
