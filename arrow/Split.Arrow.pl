#!/usr/bin/perl
$usage = "Usage: perl $0 <input.fofn> <out_dir>\n";
die $usage if (@ARGV < 2);

#$ref = shift(@ARGV);
$in = shift(@ARGV);
$out_dir = shift(@ARGV);
$sh_dir = "$out_dir/sh";
#$sh_tmp = "$out_dir/tmp";
$qsub = "$sh_dir/qsub.sh";
$qsub_jobs = 50;
$n_thread = 10;
$block = 30;
#$cmd = "/export/apps/software/smrtlink4/smrtlink_root/install/smrtlink-fromsrc_4.0.0.190159+190159-190159-189856-189856-189856/bundles/smrttools/install/smrttools-fromsrc_4.0.0.190159/smrtcmds/bin";
#$cmd = "/export/apps/software/smrtlink5/smrt_root/install/smrtlink-release_5.0.0.6792/bundles/smrttools/install/smrttools-release_5.0.0.6784/smrtcmds/bin/";
$cmd = "/data/software/smrtlink-5.01/smrtlink/install/smrtlink-release_5.0.1.9585/bundles/smrttools/install/smrttools-release_5.0.1.9578/smrtcmds/bin";
#$blasr_opt = "--bam  --bestn 10 --minMatch 12  --nproc $n_thread --minSubreadLength 500 --minAlnLength 500  --minPctSimilarity 70 --minPctAccuracy 70 --hitPolicy randombest  --randomSeed 1";
@a = split("/", $in);
$a[-1] =~ /(.*).fofn/;
$listprefix = $1;

`mkdir -p $sh_dir`;
#`mkdir -p $sh_tmp`;
open (IN, $in ) or die $!;
open (QSUB, "> $qsub") or die $!;
print QSUB "#!/bin/bash\n";

$i = 0;
while ($line = <IN>)
{
	$i++;
	chomp $line;
	@a = split("/", $line);
	$fq = $a[-1];
	$fq =~ /(.*).id/;
	$prefix = $1;
	$out1 = "$listprefix\_split.$i.sh";
#	$header = ""
	open (OUT1, "> $sh_dir/$out1") or die $!;
	print OUT1 "#!/bin/bash\ndate\n";
	print OUT1 "hostname\n";
	print OUT1 "$cmd/samtools view -h -L $line.bed ../../Merge.bam.sort.bam|awk '{if(\$1~/^@/){if(\$1~/\@RG/){split(\$2,a,\"-\");if(length(a)==1){print \$0;}}else{print \$0;}}else{split(\$NF,b,\"-\");for(i=1;i<NF;i++){printf \$i\"\\t\"}printf b[1]\"\\n\";}}'|$cmd/samtools view -bS - > Merge.sort.$i.bam\n";
	print OUT1 "$cmd/pbindex Merge.sort.$i.bam\n";
	print OUT1 "$cmd/arrow -j 10 Merge.sort.$i.bam --referenceFilename=$line.fasta -o ../../../Arrow/Sample.$i.arrow.fasta -o ../../../Arrow/Sample.$i.arrow.gff\n";
#	print OUT1 "$blasr $line $ref --out ../$prefix.blasr.bam $blasr_opt\n";
	print OUT1 "\ndate\n";
	print OUT1 "touch $prefix.split.arrow.done\n";
#	print OUT1 "$bwa $bwa_opt $ref $line|$samtools view -bS ->../$prefix.bwa.bam\n";
	close OUT1;
	if($i%3 == 0)
	{
		print QSUB "qsub -cwd -S /bin/bash -q all.q,Med.q -pe smp $n_thread -l vf=50g $out1\n";
		print QSUB "sleep 2s\n";

	}
	else
	{
		print QSUB "qsub -cwd -S /bin/bash -q all.q,Med.q -pe smp $n_thread -l vf=50g $out1\n";
		print QSUB "sleep 2s\n";

	}
	if($i%$qsub_jobs == 0)
	{
		print QSUB "sleep 2h\n";
	}
}
close IN;
close QSUB;

