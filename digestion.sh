#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: sh <pos> <name>"
    exit
fi

REF38=/data/xieshang/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

pos=$1
name=$2

samtools faidx $REF38 $pos > hg38_${name}.fa

perl /data/lipidong/softwares/Solve3.1_08232017/HybridScaffold/08212017/scripts/fa2cmap.pl -i hg38_${name}.fa -o ./ -n BSPQI -v > hg38_${name}_BSPQI.stat.txt

perl /data/lipidong/softwares/Solve3.1_08232017/HybridScaffold/08212017/scripts/fa2cmap.pl -i hg38_${name}.fa -o ./ -n BSSSI -v > hg38_${name}_BSSSI.stat.txt

