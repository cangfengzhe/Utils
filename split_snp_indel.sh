#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo 'Usage: 0 <proj_path> <raw_vcf> <proj_id> '
    exit
fi

proj_path=$1
raw_vcf=$2
proj_id=$3

mkdir $proj_path

python /data/lipidong/bin/split_snp_indel.py --input_file $raw_vcf --prefix  ${proj_path}/${proj_id}

sh /home/langna/wes_anno.sh ${proj_path}/${proj_id}.SNP.vcf ${proj_path}/${proj_id}.SNP &
sh /home/langna/wes_anno.sh ${proj_path}/${proj_id}.InDel.vcf ${proj_path}/${proj_id}.InDel &

#wait

#/data/software/python-2.7.13/bin/python /data/langna/tools/pipeline/wes/v0.1/bin/Annotation/connect_pie_v1.py -i ${proj_path}/${proj_id}.SNP.hg19_multianno.vcf -o ${proj_path}/ -I  ${proj_path}/${proj_id}.SNP.hg19_multianno.txt -g /data/langna/database/annovar/omim_genemap -a /data/langna/database/annovar/Integrated_Function.annotation.xls -p ${proj_id}_snp &

#/data/software/python-2.7.13/bin/python /data/langna/tools/pipeline/wes/v0.1/bin/Annotation/connect_pie_v1.py -i ${proj_path}/${proj_id}.InDel.hg19_multianno.vcf -o ${proj_path}/ -I  ${proj_path}/${proj_id}.InDel.hg19_multianno.txt -g /data/langna/database/annovar/omim_genemap -a /data/langna/database/annovar/Integrated_Function.annotation.xls -p ${proj_id}_indel &

wait
