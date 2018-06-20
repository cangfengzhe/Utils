#!/bin/bash
awk  '{print $1"\t"$2"\t"$3"\t0\t-\t0"}' $1 >  ${1}_2
/biocluster/data/bioexec/software/annovar/annotate_variation.pl -geneanno -buildver hg19 ${1}_2 /biocluster/data/bioexec/software/annovar/humandb/
rm ${1}_2
