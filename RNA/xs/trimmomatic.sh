java -jar /data/xieshang/bin/sources/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 /data/xieshang/2.NGS/8.BJXWZ-201804002A/1.rawdata/P10118030020-01-01_clean_0327/GYT/GYT_R1.fq.gz /data/xieshang/2.NGS/8.BJXWZ-201804002A/1.rawdata/P10118030020-01-01_clean_0327/GYT/GYT_R2.fq.gz GYT_clean_R1.fq.gz GYT_unpair_R1.fq.gz GYT_clean_R2.fq.gz GYT_unpair_R2.fq.gz ILLUMINACLIP:/data/xieshang/database/Xten_PE.fa:3:30:5 MINLEN:100
java -jar /data/xieshang/bin/sources/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 /data/xieshang/2.NGS/8.BJXWZ-201804002A/1.rawdata/P10118030020-01-01_clean_0327/XTY1/XTY1_R1.fq.gz /data/xieshang/2.NGS/8.BJXWZ-201804002A/1.rawdata/P10118030020-01-01_clean_0327/XTY1/XTY1_R2.fq.gz XTY1_clean_R1.fq.gz XTY1_unpair_R1.fq.gz XTY1_clean_R2.fq.gz XTY1_unpair_R2.fq.gz ILLUMINACLIP:/data/xieshang/database/Xten_PE.fa:3:30:5 MINLEN:100
