#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com

set -euo pipefail

source $HOME/utils/utils.sh

if [ $# -lt 1 ]; then
    echo "Usage: sh  <cmap> <output> <err>"
    exit
fi

cmap=$1
output=$2
err=$3

SOLVE_PATH='/data/lipidong/softwares/Solve3.1_08232017/'

#/data/lipidong/softwares/Solve3.1_08232017/Pipeline/08212017/runSV.py

python $SOLVE_PATH/Pipeline/08212017/runSV.py -t $SOLVE_PATH/RefAligner/6700.6902rel/ -r /data/lipidong/projects/FSHD/ref/hg38/hg38C_BSPQI_0kb_0labels.cmap -q $cmap -a /data/lipidong/softwares/Solve3.1_08232017/RefAligner/6700.6902rel/optArguments_haplotype_saphyr_human.xml -o ${output} -E $err -C /data/lipidong/projects/FSHD/clusterArguments.xml  -b /data/lipidong/softwares/Solve3.1_08232017/./Pipeline/08212017/Analysis/SV/Masks/hg38_bspqi_gap_common_segdup_min4_com10kb_seg50kb.bed
