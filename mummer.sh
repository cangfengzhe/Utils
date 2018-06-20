#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com



if [ $# -lt 1 ]; then
    echo "Usage: sh  <ref> <query> <output.prefix>"
    exit
fi

ref=$1
query=$2
output=$3

/data/software/biosoft/mummer-4.0.0beta2/nucmer -t 1 --maxmatch -l 100 -c 500 -p ${output} ${ref} ${query}
