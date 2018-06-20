#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: sh  <fq> <ref> <outdir>"
    exit
fi

fq=$1
ref=$2
output=$3
base_name=`basename $fq`
thread=10

if [ ! -f ${ref}.gmidx ]; then
    echo '参考基因组缺少索引文件'
fi

#${proj_path}/${basename}.GRCh38.graphmap.sam
/data/kangk/software/graphmap/bin/Linux-x64/graphmap align -u -B 0 -t $thread -r $ref -d $fq -o $output/${base_name}.GRCh38.graphmap.sam

samtools view -@ $thread -bS $output/${base_name}.GRCh38.graphmap.sam |samtools sort -@ $thread - -o $output/${base_name}.GRCh38.graphmap.bam
samtools index $output/${base_name}.GRCh38.graphmap.bam
rm $output/${base_name}.GRCh38.graphmap.sam
samtools stats $output/${proj_name}.pass_reads.graphmap.bam >$output/${proj_name}.pass_reads.graphmap.bam.xls
