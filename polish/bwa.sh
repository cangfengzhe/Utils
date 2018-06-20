#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com

set -euo pipefail


if [ $# -lt 1 ]; then
    echo "Usage: sh  <ref> <fq1> <fq2>"
    exit
fi

ref=$1
fq1=$2
fq2=$3
basename=`basename $fq1`
sample=${basename#_1.clean.fq.gz}

if [ ! -f ${ref}.bwt ];then
	echo "error without index"
    exit
fi


bwa mem -t 10 \
	-M -R "@RG\tID:None\tSM:${sample}"  \
	$ref $fq1 $fq2  | samtools view -@ 10 -Sb - | samtools sort -@ 10 - -o  ${sample}.sort.bam

samtools index ${sample}.sort.bam

