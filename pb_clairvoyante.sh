#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com
#
# todo



source /data/lipidong/bin/utils.sh

if [ $# -lt 1 ]; then
    echo "Usage: sh $0 <proj_name>  <ref> <bam> <chrom> <start> <end> <output>"
    exit
fi
# tips
# ${foo/%o/e} 结尾开始替换
#

proj_name=$1
ref=$2
bam=$3
chrom=$4
start=$5
end=$6
output=$7

echo PATH=/data/software/pypy-6.0.0-linux_x86_64-portable/bin:\$PATH
echo source /data/lipidong/env-tensorflow/bin/activate
echo python /data/software/Clairvoyante/clairvoyante/callVarBam.py --chkpnt_fn /data/software/Clairvoyante/trainedModels/fullv3-pacbio-ngmlr-hg001+hg002+hg003+hg004-hg19/learningRate1e-3.epoch100.learningRate1e-4.epoch200 --ref_fn ${ref} --bam_fn ${bam} --ctgName  ${chrom} --ctgStart ${start} --ctgEnd ${end}  --call_fn ${output}  --threshold 0.2 --minCoverage 2 --pypy pypy --samtools samtools --delay 3 --threads 10 --sampleName ${proj_name} --considerleftedge

