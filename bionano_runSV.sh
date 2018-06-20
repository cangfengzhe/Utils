#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: sh  <assembly-dir-including contigs> <sv_output>"
    exit
fi

thread=10

assembly_dir=$1
sv_output=$2

cmaps_dir=${assembly_dir}/contigs/exp_refineFinal1
errbin=${assembly_dir}/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1.errbin

python /data/lipidong/softwares/Solve3.1_08232017/Pipeline/08212017/runSV.py -t /data/lipidong/softwares/Solve3.1_08232017/RefAligner/6700.6902rel/avx/ -r /data/lipidong/softwares/Solve3.1_08232017/RefGenome/hg38C_BSPQI_0kb_0labels.cmap -q ${cmaps_dir} -o ${sv_output} -p /data/lipidong/softwares/Solve3.1_08232017/Pipeline/08212017 -a /data/lipidong/projects/FSHD/optArguments_haplotype_saphyr_human.xml -E $errbin -T $thread
