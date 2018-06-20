#!/usr/bin/env bash
# auther=Li Pidong
# email=lipidong@126.com
# BSPQI

set -euo pipefail


source $HOME/utils/utils.sh

if [ $# -lt 1 ]; then
    echo "Usage: sh  <ref.fa> <EXP_REFINEFINAL1.cmap> <output>"
    exit
fi

ref=$1
cmap=$2
output=$3
CWD=$(cd `dirname $0`; pwd)
SOLVE_PATH=/data/lipidong/softwares/Solve3.1_08232017/


perl $SOLVE_PATH/HybridScaffold/08212017/hybridScaffold.pl \
	-c $SOLVE_PATH/HybridScaffold/08212017/hybridScaffold_config_aggressive.xml \
	-b $cmap \
	-n $ref \
	-u GCTCTTC \
	-r $SOLVE_PATH/RefAligner/6700.6902rel/RefAligner \
	-B 2 -N 2 -f \
	-o $output \
	-p $SOLVE_PATH/Pipeline/08212017 \
	-q $CWD/clusterArguments.xml
