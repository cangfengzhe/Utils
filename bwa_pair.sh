#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
    echo 'Usage: 0 <ref-fa> <query-fa> <thread>'
    exit
fi

ref_fa=$1
query_fa=$2

name=`basename $query_fa`

if [ ! -f ${ref_fa}.sa ];then
    bwa index ${ref_fa}
fi
thread=20 #${t:=20}
echo "thread $thread"
bwa mem -t $thread $ref_fa $query_fa >  ${name}.bwa.sam
samtools view -@ $thread -b ${name}.bwa.sam |samtools sort -@ $thread - -o ${name}.bwa.bam
samtools index ${name}.bwa.bam
# rm ${name}.bwa.samrm
#rm ${ref_fa}.*
