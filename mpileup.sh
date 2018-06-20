#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com
#
# mpileup

set -euo pipefail

source $HOME/bin/utils.sh

if [ $# -lt 3 ]; then
    echo "Usage: sh $0 <ref> <bam> <output_vcf>"
    exit
fi
# tips
# ${foo/%o/e} 结尾开始替换
#

ref=$1
bam=$2
vcf=$3

samtools mpileup  -ugf $ref  $bam | bcftools call -vmO v -o  $vcf
