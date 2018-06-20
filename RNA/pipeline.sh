#!/usr/bin/env bash

if [ $# -lt 1 ]; then

    echo 'Usage: 0 <proj_name> <fq1> <fq2> <output_path>'
    exit
fi

proj_name=$1
fq1=$2
fq2=$3
output_path=$4
output_path=${output_path}/${proj_name}
if [ ! -d ${output_path} ];then
    mkdir -p ${output_path}
fi

# software path
trimmomatic=/data/langna/tools/softwares/Trimmomatic-0.36/trimmomatic-0.36.jar
tophat2=/data/software/sources/analysis/tophat-2.1.1.Linux_x86_64/tophat2
cufflinks=/data/software/sources/analysis/cufflinks-2.2.1.Linux_x86_64/cufflinks
# database path
gtf_file=/data/lipidong/data/genome/homo_hg19.gtf
hg19_bowtie2=/data/lipidong/data/genome/hg19_bowtie2/hg19
#pipeline
# trimmomatic
 java -jar $trimmomatic PE -phred33 $fq1 $fq2 ${output_path}/forward_paired.fq.gz ${output_path}/forward_unpaired.fq.gz ${output_path}/reverse_paired.fq.gz ${output_path}/reverse_unpaired.fq.gz ILLUMINACLIP:/data/langna/tools/softwares/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > ${output_path}/trimmomatic.o 2>${output_path}/trimmomatic.e

# tophat
$tophat2 -p 10 -G $gtf_file -o ${output_path} $hg19_bowtie2 ${output_path}/forward_paired.fq.gz,${output_path}/reverse_paired.fq.gz >${output_path}/tophat.o 2>${output_path}/tophat.e

# cufflinks
$cufflinks -p 10 -o ${output_path} ${output_path}/accepted_hits.bam -G ${gtf_file} > ${output_path}/cufflinks.o 2>${output_path}/cufflinks.e
