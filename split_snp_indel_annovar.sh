#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo 'Usage: 0 <out_path> <raw_vcf> <sample_id> '
    exit
fi

out_path=$1
raw_vcf=$2
sample_id=$3

mkdir $out_path

python /data/lipidong/bin/split_snp_indel.py --input_file $raw_vcf --prefix  ${out_path}/${sample_id}

sh /home/langna/wes_anno.sh ${out_path}/${sample_id}.SNP.vcf ${out_path}/${sample_id}.SNP &
sh /home/langna/wes_anno.sh ${out_path}/${sample_id}.InDel.vcf ${out_path}/${sample_id}.InDel &

#wait

#/data/software/python-2.7.13/bin/python /data/langna/tools/pipeline/wes/v0.1/bin/Annotation/connect_pie_v1.py -i ${out_path}/${sample_id}.SNP.hg19_multianno.vcf -o ${out_path}/ -I  ${out_path}/${sample_id}.SNP.hg19_multianno.txt -g /data/langna/database/annovar/omim_genemap -a /data/langna/database/annovar/Integrated_Function.annotation.xls -p ${sample_id}_snp &

#/data/software/python-2.7.13/bin/python /data/langna/tools/pipeline/wes/v0.1/bin/Annotation/connect_pie_v1.py -i ${out_path}/${sample_id}.InDel.hg19_multianno.vcf -o ${out_path}/ -I  ${out_path}/${sample_id}.InDel.hg19_multianno.txt -g /data/langna/database/annovar/omim_genemap -a /data/langna/database/annovar/Integrated_Function.annotation.xls -p ${sample_id}_indel &

wait
