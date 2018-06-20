
tophat2=/data/software/sources/analysis/tophat-2.1.1.Linux_x86_64/tophat2
cufflinks=/data/software/sources/analysis/cufflinks-2.2.1.Linux_x86_64/cufflinks
trimmomatic=/data/langna/tools/softwares/Trimmomatic-0.36/trimmomatic-0.36.jar
gtf_file=/data/lipidong/data/genome/homo_hg19.gtf
hg19_bowtie2=/data/lipidong/data/genome/hg19_bowtie2/hg19


trimmomatic:
	java -jar ${trimmomatic} PE -threads ${thread} -phred33 ${fq1} ${fq2} ${output}/${sample}_clean_R1.fq.gz ${output}/${sample}_unpair_R1.fq.gz  ${output}/${sample}_clean_R2.fq.gz ${output}/${sample}_unpair_R2.fq.gz ILLUMINACLIP:/data/xieshang/database/Xten_PE.fa:3:30:5 MINLEN:100 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > ${output}/trimmomatic.o 2> ${output}/trimmomatic.e

tophat:
	${tophat2} -p ${thread} -G ${gtf_file} -o ${output} ${hg19_bowtie2} ${output}/${sample}_clean_R1.fq.gz,${output}/${sample}_clean_R2.fq.gz > ${output}/tophat.o 2>${output}/tophat.e

cufflinks:
	${cufflinks} -p ${threado} -o ${output} ${output}/accepted_hits.bam -G ${gtf_file} > ${output}/cufflinks.o 2>${output}/cufflinks.e

htseq_count:
	python -m HTSeq.scripts.count -f bam -i gene ${output}/accepted_hits.bam /data/lipidong/database/GRCh37_latest_genomic_chr.gff > ${output}/All_gene.xls

all:trimmomatic tophat cufflinks htseq_count
