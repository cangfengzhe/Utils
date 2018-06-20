#!/bin/sh
reference=/export/personal/zhangh/project/20180129_yiyi/7_pbjelly_of6/1_pbjelly/yiyi.fasta
sample=Yiyi

#source /export/personal/lijj/0.temp/2.software/PBSuite_15.8.24/setup.sh
#ln -s ${reference} ${sample}.fasta
cp /data/lipidong/bin/pbjelly/source .
date
#source /export/personal/lijj/0.temp/2.software/PBSuite_15.8.24/setup.sh
#fakeQuals.py ${sample}.fasta ${sample}.fasta.qual

Jelly.py setup Protocol.xml
sleep 10s
qstat -f|grep setup_chun|awk '{print "qdel "$1;}' >qdel.sh
sh qdel.sh
cat source setup_chunk0.sh >setup_chunk0.sh.bak.sh
#rm -rf *.err *.out *.sh.e* *.sh.o*
qsub -cwd -q all.q,asm.q,fat.q -pe smp 4 setup_chunk0.sh.bak.sh
echo "start setup"

job=0
while [ $job != "1" ]
do
jobstate=`qstat|grep setup_chun|wc -l`
sleep 10m
if [ $jobstate == "0" ]
then
job=1
echo "end setup"
else
sleep 10m
fi
done

#source /pipeline/PBSuite_14.9.9/setup.sh
source /export/personal/lijj/0.temp/2.software/PBSuite_15.8.24/setup.sh
Jelly.py mapping Protocol.xml
sleep 2m
qstat -f|grep mapping_ch|awk '{print "qdel "$1;}' >qdel.sh
sh qdel.sh

cd mapping
#rm -rf *.err *.out *.sh.e* *.sh.o*
ls mapping_chunk*.sh|while read i;do cat ../source $i >$i.bak.sh;done
ls *.sh.bak.sh |while read i;do qsub -cwd -q all.q,asm.q,fat.q -pe smp 5 $i;done
cd ../
date
echo "start setup"
#source /pipeline/PBSuite_14.9.9/setup.sh
source /export/personal/lijj/0.temp/2.software/PBSuite_15.8.24/setup.sh
date
job=0
while [ $job != "1" ]
do
jobstate=`qstat -f|grep mapping_ch|wc -l`
sleep 10m
if [ $jobstate == "0" ]
then
job=1
echo "end mapping"
else
sleep 30m
fi
done
!
source /export/personal/lijj/0.temp/2.software/PBSuite_15.8.24/setup.sh
Jelly.py support Protocol.xml
sleep 2m
qstat -f|grep support_ch|awk '{print "qdel "$1;}' >qdel.sh
sh qdel.sh
cd support
#rm -rf *.err *.out *.sh.e* *.sh.o*
ls support_chunk*.sh|while read i;do cat ../source $i >$i.bak.sh;done
ls *.sh.bak.sh |while read i;do qsub -cwd -q all.q,asm.q,fat.q -pe smp 3 $i;done

cd ../
date
echo "start support"

job=0
while [ $job != "1" ]
do
jobstate=`qstat -f|grep support_ch|wc -l`
sleep 10m
if [ $jobstate == "0" ]
then
job=1
echo "end support"
else
sleep 10m
fi
done
!
Jelly.py extraction Protocol.xml
sleep 1m
qstat -f|grep extraction|awk '{print "qdel "$1;}' >qdel.sh
sh qdel.sh
cd assembly
#rm -rf *.err *.out *.sh.e* *.sh.o*
ls extraction_chunk0.sh|while read i;do cat ../source $i >$i.bak.sh;done
ls *.sh.bak.sh |while read i;do qsub -cwd -q all.q,asm.q,fat.q -pe smp 3 $i;done

cd ../
date
echo "start extraction"

job=0
while [ $job != "1" ]
do
jobstate=`qstat -f|grep extraction|wc -l`
sleep 10m
if [ $jobstate == "0" ]
then
job=1
echo "end extraction"
else
sleep 10m
fi
done

Jelly.py assembly Protocol.xml
sleep 2m
qstat -f|grep ass|awk '{print "qdel "$1;}' >qdel.sh
sh qdel.sh
cd assembly
#rm -rf *.err *.out *.sh.e* *.sh.o*
ls ass*.sh|while read i;do cat ../source $i >$i.bak.sh;done
ls ass*.sh.bak.sh |while read i;do qsub -cwd -q all.q,asm.q,fat.q -pe smp 4 $i;done

cd ../
#rm -rf *.err *.out *.sh.e* *.sh.o*
date
echo "start assemble"

job=0
while [ $job != "1" ]
do
jobstate=`qstat -f|grep assem|wc -l`
sleep 10m
if [ $jobstate == "0" ]
then
job=1
echo "end extraction"
else
sleep 10m
fi
done

Jelly.py output Protocol.xml
sleep 10s
qstat -f|grep out|awk '{print "qdel "$1;}' >qdel.sh
sh qdel.sh
cd assembly
#rm -rf *.err *.out *.sh.e* *.sh.o*
ls output*.sh|while read i;do cat ../source $i >$i.bak.sh;done
ls output*.sh.bak.sh |while read i;do qsub -cwd -q all.q,asm.q,fat.q -pe smp 4 $i;done

cd ../
date
#rm -rf *.err *.out *.sh.e* *.sh.o*
echo "start output"

job=0
while [ $job != "1" ]
do
jobstate=`qstat -f|grep output|wc -l`
sleep 10m
if [ $jobstate == "0" ]
then
job=1
/Script/Assembly/Scaffold_qulity.py -s jelly.out.fasta
echo "end output"
else
sleep 10m
fi
done

date
echo "All done"
