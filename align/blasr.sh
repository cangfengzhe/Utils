#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com


if [ $# -lt 1 ]; then
    echo "Usage: sh <ref>  <fastq>  <sample>"
    exit
fi

ref=$1
fastq=$2
sample=$3


blasr $fastq $ref -hitPolicy leftmost -sam -clipping soft -nproc 8  -out ${sample}.sam
samtools view -bS  ${sample}.sam | samtools sort - -o ${sample}.sort.bam
samtools index  ${sample}.sort.bam

