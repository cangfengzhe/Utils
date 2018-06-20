#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
    echo 'Usage: 0 <proj_path> <fq>'
    exit
fi

ref=/data/xieshang/4.nanopore/test/D4Z4/unit.fa
proj_path=$1
fq=$2
cd $proj_path
name=`basename $fq`

/data/xieshang/bin/minialign -t 1 -xont $ref $fq > ${name}.minialign.sam && awk '$3!~/*/{print $1}' ${name}.minialign.sam | sort | uniq > ${name}.reads.xls && /data/xieshang/bin/seqkit grep -f ${name}.reads.xls $fq > $name.D4Z4_unit.fq && rm ${name}.minialign.sam ${name}.reads.xls

