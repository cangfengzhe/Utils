#!/bin/sh
Sample=$1
Script_Path=/export/personal/lijj/0.temp/1.script/GC.Depth/
Length_cutoff=10000

perl ${Script_Path}/gc_distribution_5_plot.pl ../Reference/${Sample}.fasta ${Sample}.fasta.GC
perl ${Script_Path}/distributing_svg.pl ${Sample}.fasta.GC ${Sample}.GC.svg    # need running on dragon

# GC Depth
perl ${Script_Path}/get_gc.pl ../Reference/${Sample}.fasta ${Length_cutoff} >${Sample}.GC.${Length_cutoff}.txt

#awk 'BEGIN{print "GCpercent";}{if(NR%2==0){for(i=1;i<=NF;i++){print $i*100;}}}' GC.${Length_cutoff}.txt >GC.${Length_cutoff}.txt.change.txt
awk '{if(NR%2==1){print $1;}}' ${Sample}.GC.${Length_cutoff}.txt>${Sample}.GC.${Length_cutoff}.txt.id

/data/software/smrtlink-5.01/smrtlink/install/smrtlink-release_5.0.1.9585/bundles/smrttools/install/smrttools-release_5.0.1.9578/smrtcmds/bin/samtools depth ../Blasr/Merge.bam.sort.bam >Depth

perl ${Script_Path}/ChangeDepth.pl ../Reference/${Sample}.fasta.len Depth>${Sample}.Depth  # fill the depth with 0 when there is no subreads align

perl -e 'open IN,$ARGV[0];$bin=$ARGV[1];$s=0;$sum=0;$id=0;$new=1;$b=1;while(<IN>){chomp;@t=split;if($b == 1){$id=$t[0];$b=0;}if($id eq $t[0]){$s++;$sum+=$t[2];$id=$t[0];} else{$new=1;$sum=0;$sum+=$t[2];$s=0;$s++;$id=$t[0];}if($s%$bin==0){if($new==1){printf "\n".$t[0]."\n";printf $sum/$bin."\t";$s=0;$sum=0;$new=0;$id=$t[0];}else{printf $sum/$bin."\t";$sum=0;$id=$t[0];$s=0;$new=0;}}}' ${Sample}.Depth ${Length_cutoff} > ${Sample}.Depth.${Length_cutoff}.txt     # get depth use windows bin one by one
sed -i '1d' ${Sample}.Depth.${Length_cutoff}.txt

# Depth
python ${Script_Path}/depth_bar_plot.py -d ${Sample}.Depth.${Length_cutoff}.txt -m 150 -s 5 -p ${Sample}.Depth &
awk '{if(NR%2==1){print $1;}}' ${Sample}.Depth.${Length_cutoff}.txt >${Sample}.Depth.${Length_cutoff}.txt.id
sort ${Sample}.GC.${Length_cutoff}.txt.id >${Sample}.GC.${Length_cutoff}.txt.sort.id
sort ${Sample}.Depth.${Length_cutoff}.txt.id >${Sample}.Depth.${Length_cutoff}.txt.sort.id
comm -23 ${Sample}.GC.${Length_cutoff}.txt.sort.id ${Sample}.Depth.${Length_cutoff}.txt.sort.id >${Sample}.Unmap.id
cat ${Sample}.Unmap.id|while read i;do sed -i /${i}/,+1d ${Sample}.GC.${Length_cutoff}.txt;done                           # delete the contig with no subreads align

#/export/personal/lijj/0.temp/2.software/R-patched/bin/Rscript /export/personal/lijj/0.temp/1.script/GC.Depth/Depth_plot.R

awk 'BEGIN{print "GCpercent";}{if(NR%2==0){for(i=1;i<=NF;i++){print $i*100;}}}' ${Sample}.GC.${Length_cutoff}.txt >${Sample}.GC.${Length_cutoff}.txt.change.txt
awk 'BEGIN{print "Avgdepth";}{if(NR%2==0){for(i=1;i<=NF;i++){print $i;}}}' ${Sample}.Depth.${Length_cutoff}.txt >${Sample}.Depth.${Length_cutoff}.txt.change.txt
awk '{if(NR%2==0){for(i=1;i<=NF;i++){print $i;}}}' ${Sample}.Depth.${Length_cutoff}.txt > ${Sample}.Depth.${Length_cutoff}.txt.plot.txt
/export/personal/lijj/0.temp/2.software/R-patched/bin/Rscript /export/personal/lijj/0.temp/1.script/GC.Depth/Depth_plot.R ${Sample}.Depth.${Length_cutoff}.txt.plot.txt
paste ${Sample}.GC.${Length_cutoff}.txt.change.txt ${Sample}.Depth.${Length_cutoff}.txt.change.txt >${Sample}.GC.Depth.${Length_cutoff}.txt
/export/personal/lijj/0.temp/2.software/R-patched/bin/Rscript  ${Script_Path}/GC.Depth.R ${Sample}.GC.Depth.${Length_cutoff}.txt
convert GC.Depth.pdf ${Sample}.GC.Depth.png
date
echo "GC-Depth done"

