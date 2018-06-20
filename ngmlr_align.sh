#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
    echo 'Usage: 0 <proj_path> <proj_id> <fq> <thread:16> '
    exit
fi

proj_path=$1
proj_id=$2
fq=$3
thread=$4
date
ngmlr -t $thread -x ont -r /data/xieshang/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -q ${fq} -o ${proj_path}/${proj_id}.ngmlr.sam
date
samtools view -@ $thread -bS ${proj_path}/${proj_id}.ngmlr.sam |samtools sort -@ $thread - -o ${proj_path}/${proj_id}.ngmlr.bam
samtools index ${proj_path}/${proj_id}.ngmlr.bam
samtools stats ${proj_path}/${proj_id}.ngmlr.bam > ${proj_path}/${proj_id}.ngmlr.bam.xls
