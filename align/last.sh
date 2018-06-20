#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com
#
# last align

set -euo pipefail

source /data/lipidong//bin/utils.sh

if [ $# -lt 1 ]; then
    echo "Usage: sh $0 <ref> <fq> <name>"
    exit
fi
# tips
# ${foo/%o/e} 结尾开始替换
#
ref=$1
fq=$2
name=$3

echo $GREEN "Begin:" 2018-05-08 20:30:24

/data/software/last-912/bin/bin/lastdb -cR01 ${name}.lastdb $fq
/data/software/last-912/bin/bin/last-train -Q 1 ${name}.lastdb $fq -P 8 >${name}.last-train
/data/software/last-912/bin/bin/lastal -Q1 $ref $fq ${name}.last-train -P 8 >${name}.last.maf
/data/software/last-912/scripts/maf-convert -d sam ${name}.last.maf |samtools view -@ 8 -b - |samtools sort -@ 8 - -o ${name}.last.bam
date
samtools index ${name}.last.bam
samtools stats ${name}.last.bam >${name}.last.bam.xls
bedtools genomecov -ibam ${name}.last.bam >${name}.last.bam.depth.xls

echo  $GREEN "End:" 2018-05-08 20:30:26

