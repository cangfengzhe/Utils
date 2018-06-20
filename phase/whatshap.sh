#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com
# 采用whatshap对HLA进行分相
# 1. mpileup,保留杂合位点
# 2. whatshap分型
# 3. 生成concencus序列

# 输入文件
# ref 参考xulie
# bam 比对文件
# sample 样本名
# output 结果保存目录

set -euo pipefail


if [ $# -lt 1 ]; then
    echo "Usage: sh  <ref> <bam> <sample> <output>"
    exit
fi
source /data/lipidong/my.bashrc

ref=$1
bam=$2
sample=$3
output=$4
new_bam="output/$(basename $bam | sed 's/\(.*\)\.bam.*/\1/').SM.bam"

samtools view -H $bam  | sed 's,^@RG.*,@RG\tID:None\tSM:None\tLB:None\tPL:PacBio,g' |  samtools reheader - $bam  > $new_bam
samtools index $new_bam
#samtools mpileup -t DP,AD -uf $ref  $bam | bcftools call -vmO v  -o ${output}/${sample}.mpileup.vcf
samtools mpileup  -uf $ref  $new_bam | bcftools call -vmO v  -o ${output}/${sample}.mpileup.vcf

#过滤杂合突变
#grep -E '(^#)|(0/1)' ${output}/${sample}.mpileup.vcf > ${output}/${sample}.mpileup.hete.vcf

# whatshap 分型
/data/lipidong/.local/bin/whatshap phase --include-homozygous --changed-genotype-list ${output}/${sample}.changed-genotype-list.txt --max-coverage 15  --output-read-list ${output}/${sample}.output-read-list.fasta --tag PS --reference $ref --ignore-read-groups --distrust-genotypes -o ${output}/${sample}.phased.vcf ${output}/${sample}.mpileup.vcf $new_bam

# make concencus
/data/software/seqmule/SeqMule/exe/tabix/bgzip ${output}/${sample}.phased.vcf
/data/software/seqmule/SeqMule/exe/tabix/tabix ${output}/${sample}.phased.vcf.gz
bcftools consensus -H 1 -f $ref ${output}/${sample}.phased.vcf.gz > ${output}/${sample}.haplotype1.fasta
bcftools consensus -H 2 -f $ref ${output}/${sample}.phased.vcf.gz > ${output}/${sample}.haplotype2.fasta

