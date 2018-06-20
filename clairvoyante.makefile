help:
	echo "USAGE 0: <ref> <bam> <contig> <output_vcf> <proj_name>"
pacbio:
	source /data/lipidong/env-tensorflow/bin/activate
	python /data/software/Clairvoyante/clairvoyante/callVarBam.py --chkpnt_fn /data/software/Clairvoyante/trainedModels/fullv3-pacbio-ngmlr-hg001+hg002+hg003+hg004-hg19/learningRate1e-3.epoch100.learningRate1e-4.epoch200 --ref_fn ${ref} --bam_fn ${bam} --ctgName  ${contig}  --call_fn ${proj_name}.vcf --threshold 0.2 --minCoverage 1 --pypy pypy --samtools samtools --delay 10 --threads 10 --sampleName ${proj_name} --considerleftedge
	deactivate
	python /data/lipidong/bin/bam2matrix.py  --vcf_file ${proj_name}.vcf --bam_file ${bam} --output_prefix ${proj_name}

ont:
	source /data/lipidong/env-tensorflow/bin/activate
	python /data/software/Clairvoyante/clairvoyante/callVarBam.py --chkpnt_fn  /data/software/Clairvoyante/trainedModels/fullv3-ont-ngmlr-hg001-hg19/learningRate1e-3.epoch999 --ref_fn ${ref} --bam_fn ${bam} --ctgName  ${contig}  --call_fn ${proj_name}.vcf --threshold 0.25 --minCoverage 1 --pypy pypy --samtools samtools --delay 10 --threads 10 --sampleName ${proj_name} --considerleftedge
	deactivate
	python /data/lipidong/bin/bam2matrix.py  --vcf_file ${proj_name}.vcf --bam_file ${bam} --output_prefix ${proj_name}
