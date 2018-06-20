#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com
#
# todo

set -euo pipefail

source $HOME/bin/utils.sh

if [ $# -lt 1 ]; then
    echo "Usage: sh $0 <ref>  <input_bam> <output_bam>"
    exit
fi
# tips
# ${foo/%o/e} 结尾开始替换
#
ref=$1
bam=$2
output=$3
gatk=/data/lipidong/gatk-3.6/GenomeAnalysisTK.jar
java -Xmx1750m -jar $gatk -T RealignerTargetCreator -R $ref  -I $bam -o ${output}.indelcreate.intervals -nt 8

java -Xmx1750m -jar $gatk  -I $bam -R $ref -T IndelRealigner -targetIntervals ${output}.indelcreate.intervals -o ${output}  -filterNoBases -filterMBQ -filterRNC
samtools index ${output}
