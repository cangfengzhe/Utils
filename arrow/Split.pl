#!/usr/bin/perl
$usage = "Usage: perl $0 <input.fofn> <out_dir>\n";
die $usage if (@ARGV < 2);

#$ref = shift(@ARGV);
$in = shift(@ARGV);
$out_dir = shift(@ARGV);
$sh_dir = "$out_dir/sh";
#$sh_tmp = "$out_dir/tmp";
$qsub = "$sh_dir/qsub.sh";
$qsub_jobs = 100;
$n_thread = 1;
$block = 30;
#$bwa = "bwa";
#$samtools = "samtools";
#$bwa_opt = "mem -t $n_thread -x pacbio";
#$pbalign = "/pipeline/smrtanalysis/install/smrtanalysis_2.3.0.140936/analysis/bin/pbalign";
#$pbalign = "/pipeline/smrtlink/smrtlink_root/install/smrtlink-fromsrc_3.1.1.182868,182868-182868-182483-182249-182249/bundles/smrttools/install/smrttools-fromsrc_3.1.1.182868/smrtcmds/bin/pbalign";
#$pbalign_opt = "--nproc $n_thread --forQuiver";
#$blasr = "/export/apps/software/smrtlink4/smrtlink_root/install/smrtlink-fromsrc_4.0.0.190159+190159-190159-189856-189856-189856/bundles/smrttools/install/smrttools-fromsrc_4.0.0.190159/smrtcmds/bin/blasr";
#$samtools = "/export/apps/software/smrtlink4/smrtlink_root/install/smrtlink-fromsrc_4.0.0.190159+190159-190159-189856-189856-189856/bundles/smrttools/install/smrttools-fromsrc_4.0.0.190159/smrtcmds/bin/samtools";
#$samtools = "/export/apps/software/smrtlink5/smrt_root/install/smrtlink-release_5.0.0.6792/bundles/smrttools/install/smrttools-release_5.0.0.6784/smrtcmds/bin/samtools";
$samtools = "/data/software/smrtlink-5.01/smrtlink/install/smrtlink-release_5.0.1.9585/bundles/smrttools/install/smrttools-release_5.0.1.9578/smrtcmds/bin/samtools";
$blasr_opt = "--bam  --bestn 10 --minMatch 12  --nproc $n_thread --minSubreadLength 500 --minAlnLength 500  --minPctSimilarity 70 --minPctAccuracy 70 --hitPolicy randombest  --randomSeed 1";
#$pbalign_opt = "--nproc $n_thread --minAccuracy 70.0 --minLength 50 --hitPolicy=randombest --log-level INFO --algorithmOptions=\"-minMatch 12 -bestn 10 -minPctSimilarity 70.0 -randomSeed 1\" --tmpDir";
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
	$fq =~ /(.*).bam/;
	$prefix = $1;
	$out1 = "$listprefix\_split.$i.sh";
#	$header = ""
	open (OUT1, "> $sh_dir/$out1") or die $!;
	print OUT1 "#!/bin/bash\ndate\n";
	print OUT1 "hostname\n";
#	print OUT1 "$samtools view -H $line >$prefix.header\n\n";
#	print OUT1 "source /pipeline/smrtanalysis/install/smrtanalysis_2.3.0.140936/etc/setup.sh\n";
#	print OUT1 "$pbalign $pbalign_opt $sh_tmp $line $ref ../$prefix.pbalign.bam\n";
        print OUT1 "$samtools view -h $line|awk -v bam=$prefix -v b=$block \'{if(\$1~/@/){for(i=0;i<b;i++){print \$0>>bam\".split.\"i\".sam\";}}else{n=NR%b;print \$0>>bam\".split.\"n\".sam\";}}\' ";
#	print OUT1 "$blasr $line $ref --out ../$prefix.blasr.bam $blasr_opt\n";
	print OUT1 "\ndate\n";
	print OUT1 "ls $prefix.split.*.sam |while read j;do $samtools view -bS \$j >../\$j.bam;done\n";
	print OUT1 "rm -rf $prefix.split.*.sam\n";
	print OUT1 "touch $prefix.split.bam.done\n";
#	print OUT1 "$bwa $bwa_opt $ref $line|$samtools view -bS ->../$prefix.bwa.bam\n";
	close OUT1;
	if($i%3 == 0)
	{
		print QSUB "qsub -cwd -pe smp $n_thread -l vf=10g $out1\n";
		print QSUB "sleep 2s\n";

	}
	else
	{
		print QSUB "qsub -cwd -pe smp $n_thread -l vf=10g $out1\n";
		print QSUB "sleep 2s\n";

	}
	if($i%$qsub_jobs == 0)
	{
		print QSUB "sleep 1h\n";
	}
}
close IN;
close QSUB;

