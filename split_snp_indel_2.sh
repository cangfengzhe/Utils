#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo 'Usage: 0 <proj_path> <raw_vcf> <proj_id> '
    exit
fi

proj_path=$1
raw_vcf=$2
proj_id=$3

mkdir $proj_path

java  -XX:ParallelGCThreads=2  -Xmx20g -jar /data/langna/tools/softwares/GenomeAnalysisTK.jar -T SelectVariants -R /data/langna/database/genome/human_hg19/sequences.fa -V $raw_vcf -selectType SNP -o ${proj_path}/${proj_id}_snp.vcf

java -XX:ParallelGCThreads=2  -Xmx20g -jar /data/langna/tools/softwares/GenomeAnalysisTK.jar -T SelectVariants -R /data/langna/database/genome/human_hg19/sequences.fa -V $raw_vcf -selectType indel -o ${proj_path}/${proj_id}_indel.vcf


sh /home/langna/wes_anno.sh ${proj_path}/${proj_id}_snp.vcf ${proj_path}/${proj_id}_snp &
sh /home/langna/wes_anno.sh ${proj_path}/${proj_id}_indel.vcf ${proj_path}/${proj_id}_indel &

wait

/data/software/python-2.7.13/bin/python /data/langna/tools/pipeline/wes/v0.1/bin/Annotation/connect_pie_v1.py -i ${proj_path}/${proj_id}_snp.hg19_multianno.vcf -o ${proj_path}/ -I  ${proj_path}/${proj_id}_snp.hg19_multianno.txt -g /data/langna/database/annovar/omim_genemap -a /data/langna/database/annovar/Integrated_Function.annotation.xls -p ${proj_path}/${proj_id}_filter_snp.anno
/data/software/python-2.7.13/bin/python /data/langna/tools/pipeline/wes/v0.1/bin/Annotation/connect_pie_v1.py -i ${proj_path}/${proj_id}_indel.hg19_multianno.vcf -o ${proj_path}/ -I  ${proj_path}/${proj_id}_indel.hg19_multianno.txt -g /data/langna/database/annovar/omim_genemap -a /data/langna/database/annovar/Integrated_Function.annotation.xls -p ${proj_path}/${proj_id}_filter_indel.anno
