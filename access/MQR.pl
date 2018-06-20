#!/usr/bin/perl -w
########################################################################
# Last modified date: 08/08/2017                                       #
# Purpose: Molecule Quality Report (MQR) Tool                          #
#                                                                      #
# Author: Xiang Zhou, Computational Biologist                          #
# Email : xzhou@bionanogenomics.com                                    #
# Affiliation: Research Department, BioNano Genomics Inc.              #
#                                                                      #
# Usage:                                                               #
#   perl MQR.pl [options] <Args>                                       #
# Options:                                                             #
#   -h    : This help message                                          #
#   -i    : Input BNX file (Required, can be 1 or 2 channels)          #
#   -o    : Output folder (default: input folder)                      #
#   -l    : Minimum length of the molecule in kb (Required)            #
#   -n    : Minimum nicks in each channel of the molecule (Required)   #
#   -j    : Population ID and Callback URL information as JSON string  #
#   -Q    : Calling from V3_Cohort_QC.pl script? (default: OFF)        #
#   -D    : Debug mode to keep all temporary files (default: OFF)      #
#   -P    : Plot SNR histogram in the folder specified (no default)    #
########################################################################

use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Spec;
use File::Path qw(make_path);
use Cwd qw(abs_path realpath);
use Getopt::Std;
#use IO::Socket;
#use Data::Dumper;

BEGIN {
	my $script_dir = abs_path(dirname($0));
	my $module_dir = abs_path($script_dir . "/lib/");
	push(@INC, $module_dir);
}
use JSON;
#use AnyEvent;
#use AnyEvent::PocketIO::Client;
#use PocketIO::Client::IO;
#use PocketIO::Handle;
#use Protocol::SocketIO::Message;

#my $datestring = localtime();
#print "Local date and time: $datestring\n";

my ($default_min_SNR, $default_max_SNR, $default_SNR) = (2, 20, 3.5);
my ($debug_mode, $from_V3_Cohort_QC, $Min_length, $Min_labels, $Min_length_baseline) = (0, 0, 0, 0, 20e3);
my ($min_num_of_labels, $min_num_of_mols_above_150kb) = (1000, 500);
my @not_enough_labels = (0, 0);
my @not_enough_mols = (0, 0);
my @bpp_values = (475, 525, 425, 575);
my ($json_string, $command) = ("") x 2;
my ($fin, $fout, $in_dir, $out_dir, $fout_SNR) = ("") x 5;
my $script_dir = abs_path(dirname($0));
my (%opts, %json);
my $RefAligner = `which RefAligner | head -n 1`;
chomp($RefAligner);
if($RefAligner eq ""){
	print STDERR ("ERROR: Could not find RefAligner in the default path!\n");
	exit (1);
}

Init();
Check_inputs();

######################################################
my (@RCmap, @subSamplingSeed, @subSamplingRange, @usechannel, @SNR_cutoff, @populationPK);
my $N_populations = 1;
my $N_channels;
GetOpts($json_string);
#print Dumper %json;

my $N_BNX_channels = `tail -n 20 $fin | grep "^QX21" | wc -l`;
chomp($N_BNX_channels);
if($N_BNX_channels > 0){
	$N_BNX_channels = 2;
}
else{
	$N_BNX_channels = 1;
}
if( $N_channels > $N_BNX_channels || $usechannel[0] > $N_BNX_channels || defined($usechannel[1]) && $usechannel[1] > $N_BNX_channels ){
	print STDERR ("ERROR: The number of channels in the BNX file ($N_BNX_channels) is less than the number of channels sepcified in the command argument!\n");
	exit (1);
}

my @fname = ($fout) x 2;
my ($fname_summary, $fname_status) = ($fout) x 2;
my ($IN, $OUT, $OUT1, $OUT2, $OUT_summary, $OUT_status, $ERR);
my (%mol, %label, %SNR, %mol_SNR, %int, %mol_int);
my @length = (0, 0, 0);
my @length_ge_20kb = (0, 0, 0);
my (@length_ref_array, @length_ref_array_ge_20kb);
my @n_mol = (0, 0, 0);
my @total = (0, 0, 0);

######################################################
# Autoflush for files
my $old_fh = select($OUT_status);
local $| = 1;
select($old_fh);

######################################################
my $datatype = "";
# If BNX version >= 1.3: Saphyr;
# If BNX version < 1.2: Irys;
# If BNX version is N/A: Irys;
# If BNX version == 1.2:
#  Saphyr: if the following search of the BNX returns a nonempty string:
#    grep -m 1 '# Run Data\s.*\sUNKNOWN_X0\s*1999-01-01 01:00:00 AM\s*68819821\s*0.83\s*[0-9]*\s*[0-9]*\s*UNKNOWN,UNKNOWN,UNKNOWN,1'
#  Irys: Otherwise.
my $version = `grep -P "# BNX File Version:" $fin | cut -d: -f2`;
chomp($version);
if($version =~ /\S/){
	$version += 0;
	if($version >= 1.3){
		$datatype = "Saphyr";
	}
	elsif($version < 1.2){
		$datatype = "Irys";
	}
	elsif($version == 1.2){
		my $tmp_string = `grep -m 1 '# Run Data\\s.*\\sUNKNOWN_X0\\s*1999-01-01 01:00:00 AM\\s*68819821\\s*0.83\\s*[0-9]*\\s*[0-9]*\\s*UNKNOWN,UNKNOWN,UNKNOWN,1' $fin`;
		chomp($tmp_string);
		if( $tmp_string ne "" ){
			$datatype = "Saphyr";
		}
		else{
			$datatype = "Irys";
		}
	}
}
else{
	$datatype = "Irys";
}

if($datatype eq "Saphyr"){
	# SNR estimation
	if( defined($opts{P}) ){
		$command = "$^X $script_dir/filter_SNR_dynamic.pl -i $fin -m $default_min_SNR -M $default_max_SNR -d $default_SNR -L $min_num_of_labels -P $opts{P} > $out_dir/SNR";
	}
	else{
		$command = "$^X $script_dir/filter_SNR_dynamic.pl -i $fin -m $default_min_SNR -M $default_max_SNR -d $default_SNR -L $min_num_of_labels > $out_dir/SNR";
	}
	#print $command, "\n";
	RunCommand($command);

	for(my $i = 0; $i < $N_BNX_channels; $i++){
		my $tmp = "grep 'Estimated SNR cutoff for channel " . ($i+1) . "' $out_dir/SNR | cut -f2";
		$SNR_cutoff[$i] = `$tmp`;
		chomp($SNR_cutoff[$i]);
	}
	#print "[$SNR_cutoff[0]][$SNR_cutoff[1]]\n";

	for(my $i = 0; $i < $N_BNX_channels; $i++){
		if($SNR_cutoff[$i] eq "NA"){
			last;
		}
		elsif($SNR_cutoff[$i] == -1){
			$not_enough_labels[$i] = 1;
			print STDERR ("Warning: There are not enough labels for channel ", $i+1, "! Use default SNR cutoff value of $default_SNR instead.\n");
			$SNR_cutoff[$i] = $default_SNR;
		}
		elsif( $SNR_cutoff[$i] < $default_min_SNR || $SNR_cutoff[$i] > $default_max_SNR ){
			print STDERR ("Warning: The estimated SNR cutoff value for channel ", $i+1, " is out of range! Use default SNR cutoff value of $default_SNR instead.\n");
			$SNR_cutoff[$i] = $default_SNR;
		}
	}
}
elsif($datatype eq "Irys"){
	for(my $i = 0; $i < $N_BNX_channels; $i++){
		$SNR_cutoff[$i] = 0;
	}
}
else{
	print STDERR ("ERROR: The input BNX is neither Saphyr nor Irys datatype!\n");
	exit (1);
}
######################################################
#for(my $i = 0; $i < 1; $i++){
#	my $tmp = "_p" . ($i+1) . ".bnx";
#	$fname[$i] =~ s/\.bnx$/$tmp/i;
#}
$fname[0] = $fin;
$fname_summary =~ s/\.bnx$/_summary.txt/i;
$fname_status =~ s/\.bnx$/_status.txt/i;
open($OUT_summary, ">".$fname_summary) || die ("ERROR: Can't open $fname_summary: $!\n");

# ColorSplit
#if($N_populations == 1){
#	$command = "$RefAligner -merge -bnx -i $fin -o $out_dir/" . basename($fname[0], ".bnx") . " -f -minSNR $SNR_cutoff[0] > /dev/null";
#}
#else{
#	$script_dir = abs_path(dirname($0));
#	$command = "$^X $script_dir/ColorSplit.pl -i $fin -o $out_dir -1 $SNR_cutoff[0] -2 $SNR_cutoff[1]";
#}
#print "$command\n";
#RunCommand($command);

######################################################
print $OUT_summary("PROGRAM PARAMETERS\n");
print $OUT_summary "Input BNX file:\t", $fin, "\n";
for(my $i = 0; $i < $N_channels; $i++){
	if($RCmap[$i] eq ""){
		#print $OUT_summary "Input reference file for channel ", $usechannel[$i],":\t", "N/A", "\n";
	}
	elsif(-f $RCmap[$i]){
		print $OUT_summary "Input reference file for channel ", $usechannel[$i], ":\t", abs_path($RCmap[$i]), "\n";
	}
	else{
		print STDERR ("ERROR: Could not find the reference CMAP for channel ", $usechannel[$i], "!\n");
		exit (1);
	}
}
print $OUT_summary "Output folder:\t", $out_dir, "\n";
#for(my $i = 0; $i < $N_populations; $i++){
#	print $OUT_summary "Output BNX file for population ", $i+1,":\t", $fname[$i], "\n";
#}
print $OUT_summary "Minimum nicks in each channel of the molecule:\t$Min_labels\n";
print $OUT_summary "Minimum length of the molecule in kb:\t", $Min_length/1000, "\n";
for(my $i = 0; $i < $N_channels; $i++){
	print $OUT_summary "SNR cutoff for channel ", $i+1, ":\t", sprintf("%.3f", $SNR_cutoff[$i]), "\n";
}

######################################################
#$command = "grep -v '#' $fin | grep '^0' | wc -l";
##print $command, "\n";
#my $num_mols = `$command`;
#chomp($num_mols);
######################################################
for(my $i = 0; $i < $N_populations; $i++){
	open($IN, $fname[$i]) || die ("ERROR: Can't open $fname[$i]: $!\n");
	#open($OUT, ">".$fname) || die ("ERROR: Can't open $fname: $!\n");
	open($OUT_status, ">".$fname_status) || die ("ERROR: Can't open $fname_status: $!\n");
	my $population = $i+1;
	my %molecule;
	while(my $line = <$IN>){
		chomp $line;
		$line =~ s/\r//g;
		
		if($line =~ /^#/){
			next;
		}
		
		my @x = split("\t", $line);
		
		if($line =~ /^0/){
			$total[$i+1]++;
			$molecule{population} = $population;
			$molecule{length} = $x[2];
			if( $molecule{length} >= $Min_length_baseline ){
				$length_ge_20kb[$population] += $molecule{length};
				push(@{$length_ref_array_ge_20kb[$population]}, $molecule{length});
			}
		}
		elsif($line =~ /^1/){
			$molecule{channel_1} = scalar(@x) - 2;
		}
		elsif($line =~ /^2/){
			$molecule{channel_2} = scalar(@x) - 2;
		}
		elsif($line =~ /^QX11/){
			shift(@x);
			$molecule{SNR_1} = \@x;
			
			# Count the labels after the SNR filter
			foreach my $item (@{$molecule{SNR_1}}){
				if($item < $SNR_cutoff[0]){
					$molecule{channel_1}--;
				}
			}
		}
		elsif($line =~ /^QX21/){
			shift(@x);
			$molecule{SNR_2} = \@x;
			
			# Count the labels after the SNR filter
			foreach my $item (@{$molecule{SNR_2}}){
				if($item < $SNR_cutoff[1]){
					$molecule{channel_2}--;
				}
			}
		}
		elsif($line =~ /^QX12/){
			shift(@x);
			$molecule{int_1} = \@x;
			
			# Get the metrics for each color and for each population
			if( $molecule{length} >= $Min_length ){
				$n_mol[$population]++;
				$length[$population] += $molecule{length};
				push(@{$length_ref_array[$population]}, $molecule{length});
				$label{channel_1}->[$population] += $molecule{channel_1};
				
				my ($SNR_sum, $int_sum) = (0, 0);
				for(my $j = 0; $j < @{$molecule{SNR_1}}; $j++){
					my $SNR_value = @{$molecule{SNR_1}}[$j];
					my $int_value = @{$molecule{int_1}}[$j];
					if($SNR_value >= $SNR_cutoff[0]){
						$SNR_sum += $SNR_value;
						$int_sum += $int_value;
					}
				}
				
				$SNR{channel_1}->[$population] += $SNR_sum;
				if( $molecule{channel_1} ){
					$mol_SNR{channel_1}->[$population] += $SNR_sum / $molecule{channel_1};
				}
				$int{channel_1}->[$population] += $int_sum;
				if( $molecule{channel_1} ){
					$mol_int{channel_1}->[$population] += $int_sum / $molecule{channel_1};
				}
				
				if($N_BNX_channels == 1){
					%molecule = ();
				}
			}
		}
		elsif($line =~ /^QX22/){
			shift(@x);
			$molecule{int_2} = \@x;
			
			# Get the metrics for each color and for each population
			if( $molecule{length} >= $Min_length ){
				$label{channel_2}->[$population] += $molecule{channel_2};
				
				my ($SNR_sum, $int_sum) = (0, 0);
				for(my $j = 0; $j < @{$molecule{SNR_2}}; $j++){
					my $SNR_value = @{$molecule{SNR_2}}[$j];
					my $int_value = @{$molecule{int_2}}[$j];
					if($SNR_value >= $SNR_cutoff[1]){
						$SNR_sum += $SNR_value;
						$int_sum += $int_value;
					}
				}
				
				$SNR{channel_2}->[$population] += $SNR_sum;
				if( $molecule{channel_2} ){
					$mol_SNR{channel_2}->[$population] += $SNR_sum / $molecule{channel_2};
				}
				$int{channel_2}->[$population] += $int_sum;
				if( $molecule{channel_2} ){
					$mol_int{channel_2}->[$population] += $int_sum / $molecule{channel_2};
				}
			}
			%molecule = ();
		}
	}
	#printf $OUT_status("%d%% done\n", $total[$i]/$num_mols*100);
	#close($OUT_status);
	close($IN);
	#close($OUT);
}

for(my $i = 0; $i < $N_populations; $i++){
	if($n_mol[$i+1] < $min_num_of_mols_above_150kb){
		$not_enough_mols[$i] = 1;
	}
}

print $OUT_summary("\nSUMMARY\n");
print $OUT_summary("Total molecules in the input BNX:\t", $total[1]+$total[2], "\n");
#for(my $i = 0; $i < $N_populations; $i++){
#	print $OUT_summary("Total molecules in population ", $i+1, ":\t", $total[$i+1], "\n");
#}
for(my $i = 0; $i < $N_channels; $i++){
	if( -f $RCmap[$i] && !$not_enough_mols[0] && !$not_enough_labels[$i] ){
		print $OUT_summary "\nRefAligner absolute path:\t",  abs_path($RefAligner), "\n";
		last;
	}
}
print $OUT_summary("\nQC_RESULT\n");

######################################################
for(my $i = 0; $i < $N_populations; $i++){
	Generate_stat(\%label, \%SNR, \%mol_SNR, \%int, \%mol_int, \@n_mol, \@length, \@length_ge_20kb, \@length_ref_array, \@length_ref_array_ge_20kb, $i+1);
}

######################################################
my @json_output;
my $RefAlignerErr = 0;
for(my $c = 0; $c < $N_channels; $c++){
	my $channel = "channel_" . $usechannel[$c];
	my ($ref_length, $StretchFactor, $FP_100kb, $FNrate, $siteSD, $scaling_sd, $bpp, $bppSD, $FPrate, $relative_sd, $se, $sf, $sr, $avg_mol_snr, $ave_mol_intensity, $LabelDensity, $Maps, $GoodMaps, $TotalAlignedlength) = (0) x 19;
	
	if(-f $RCmap[$c]){
		$command = "grep -v \"#\" " . $RCmap[$c] . " | awk '{if(\$5==0) a+=\$2;} END {print a}'";
		#print $command, "\n";
		$ref_length = `$command`;
		chomp($ref_length);
	}
	if(-f $RCmap[$c] && !$not_enough_mols[0] && !$not_enough_labels[$c] && $label{$channel}->[1]){
		# Chop up the reference into overlapping contigs of size 2.4M with 0.4M overlap
		$command = "$RefAligner -merge -ref " . $RCmap[$c] . " -winbreak 2400 2000 1000 -o $out_dir/" . basename($fname[0], ".bnx") . "_split -minsites 10 -f > /dev/null";
		#print $command, "\n\n";
		RunCommand($command);
		
		# Retry several bpp values to prevent the program goes into the local maxima
		my (@SNR_cutoff_tmp, @StretchFactor_tmp, @FP_100kb_tmp, @FNrate_tmp, @siteSD_tmp, @scaling_sd_tmp, @bpp_tmp, @bppSD_tmp, @FPrate_tmp, @relative_sd_tmp, @se_tmp, @sf_tmp, @sr_tmp, @LabelDensity_tmp, @Maps_tmp, @GoodMaps_tmp, @TotalAlignedlength_tmp);
		my $best_maprate = 0;
		my $best_maprate_index = 0;
		
		#my $bpp_filename = "$out_dir/../../bpp";
		my $bpp_filename = "$in_dir/../bpp";
		my $last_bpp = 0;
		if($from_V3_Cohort_QC && -f $bpp_filename){
			$last_bpp = `tail -n 1 $bpp_filename`;
			chomp($last_bpp);
			if($last_bpp >= 400 && $last_bpp <= 600){
				unshift(@bpp_values, $last_bpp);
			}
		}
		for(my $iter = 0; $iter < @bpp_values; $iter++){
			my $para = "-T 1e-11 -A 7 -L 130 -BestRef 1 -BestRefPV 1 -res 3.3 -mres 0.9 -extend 1 -outlier 1e-3 -endoutlier 1e-4 -deltaX 4 -deltaY 6 -nosplit 2 -biaswt 0 -S -10000 -PVres 2 -AlignRes 1.5 -resEstimate -resbias 4.0 64 -MaxSE 0.5 -rres 1.2 -maptype 0 -outlierMax 40 -PVendoutlier -FP 1.5 -FN 0.15 -sf 0.15 -sd 0 -sr 0.03 -se 0.2 -relerr 0.001 -M 3 3 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 2 -hash -hashdelta 10 -HSDrange 1 -hashoffset 1 -hashbest 0 -maxmem 64 -insertThreads 4 -queryThreads 16 -maxvirtmem 0 -maxthreads 64 -Kmax 6 -maxsites 100 150 -f -randomize $subSamplingSeed[$c] -subset 1 $subSamplingRange[$c] -minlen " . sprintf("%s", $Min_length/1000) . " -bpp $bpp_values[$iter] -BppOutput 1 -minsites $Min_labels -minSNR $SNR_cutoff[$c] -usecolor $usechannel[$c] -output-veto-filter '.bnx\$' -stdout -stderr";
			
			my $sub_dir = "c" . ($c+1) . "_round_$iter";
			make_path("$out_dir/$sub_dir");
			$command = "$RefAligner -i $fname[0] -ref $out_dir/" . basename($fname[0], ".bnx") . "_split.cmap -o $out_dir/$sub_dir/" . basename($fname[0], ".bnx") . " $para";
			#print $command, "\n\n";
			RunCommand($command);
			
			if(! -f "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err"){
				print STDERR ("ERROR: Could not find the .err file for channel ", $usechannel[$c], "!\n");
				$RefAlignerErr = 1;  #exit (1);
			}
			
			# Read .err file to get the SNR cutoff
			##$command = "tail -n 1 " . "$out_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{print \$NF}'";
			#$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^minSNR\$/){a=i} } END {print \$a}'";
			##print $command, "\n";
			#$SNR_cutoff_tmp[$iter] = `$command`;
			#chomp($SNR_cutoff_tmp[$iter]);
			##print "Estimated SNR cutoff for channel/population 1:\t", $SNR_cutoff_tmp[$iter], "\n";
			
			# Get the mapping statistics
			$command = "grep -P \"^# Run Data\\t\" " . $fname[0] . " | cut -f6 | head -n 1";
			#print $command, "\n";
			$StretchFactor_tmp[$iter] = `$command`;
			chomp($StretchFactor_tmp[$iter]);
			if( $StretchFactor_tmp[$iter] eq "" ){
				$StretchFactor_tmp[$iter] = 0;
			}
			if( !$RefAlignerErr) {			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^FP\\(\\/100kb\\)\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$FP_100kb_tmp[$iter] = `$command`;
			chomp($FP_100kb_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^FN\\(rate\\)\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$FNrate_tmp[$iter] = `$command`;
			chomp($FNrate_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^sf\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$siteSD_tmp[$iter] = `$command`;
			chomp($siteSD_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^sd\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$scaling_sd_tmp[$iter]  = `$command`;
			chomp($scaling_sd_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^bpp\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$bpp_tmp[$iter] = `$command`;
			chomp($bpp_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^bppSD\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$bppSD_tmp[$iter] = `$command`;
			chomp($bppSD_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^FPrate\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$FPrate_tmp[$iter] = `$command`;
			chomp($FPrate_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^sr\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$relative_sd_tmp[$iter]  = `$command`;
			chomp($relative_sd_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^se\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$se_tmp[$iter] = `$command`;
			chomp($se_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^sf\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$sf_tmp[$iter] = `$command`;
			chomp($sf_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^sr\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$sr_tmp[$iter] = `$command`;
			chomp($sr_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^LabelDensity\\(\\/100kb\\)\$/){a=i} } {print \$a}' | tail -n 2 | head -n 1";
			#print $command, "\n";
			$LabelDensity_tmp[$iter] = `$command`;
			chomp($LabelDensity_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^Maps\$/){a=i} } END {print \$a}'";
			#print $command, "\n";
			$Maps_tmp[$iter] = `$command`;
			chomp($Maps_tmp[$iter]);
			
			$command = "grep -v \"#\" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".err" . " | awk '{ for(i=1;i<=NF;i++) if(\$i ~ /^GoodMaps\$/){a=i} } {print \$a}' | tail -n 2 | head -n 1";
			#print $command, "\n";
			$GoodMaps_tmp[$iter] = `$command`;
			chomp($GoodMaps_tmp[$iter]);
			} ### $RefAlignerErr
			# Estimating the effective coverage!
			$command = "grep \"Query : Total Aligned length = \" " . "$out_dir/$sub_dir/" . basename($fname[0], ".bnx") . ".stdout | cut -f7 -d\" \"";
			#print $command, "\n";
			$TotalAlignedlength_tmp[$iter] = `$command`;
			chomp($TotalAlignedlength_tmp[$iter]);
			
			if(!$Maps_tmp[$iter]){
				last;
			}
			
			if($from_V3_Cohort_QC && -f $bpp_filename && $GoodMaps_tmp[$iter] > 0){
				$best_maprate = $GoodMaps_tmp[$iter] / $Maps_tmp[$iter];
				$best_maprate_index = $iter;
				last;
			}
			
			if($GoodMaps_tmp[$iter] / $Maps_tmp[$iter] > $best_maprate){
				$best_maprate = $GoodMaps_tmp[$iter] / $Maps_tmp[$iter];
				$best_maprate_index = $iter;
			}
			
			if($best_maprate >= 0.5){
				last;
			}
		}
		if( ! $RefAlignerErr) {
		if( $from_V3_Cohort_QC &&
			($bpp_tmp[$best_maprate_index] >= 400 && $bpp_tmp[$best_maprate_index] <= 600 && $best_maprate > 0) && 
			(-f $bpp_filename && $best_maprate_index > 0 || ! -f $bpp_filename)
		){
			open(my $BPP_tmpfile, ">$out_dir/bpp");
			print $BPP_tmpfile("$bpp_tmp[$best_maprate_index]\n");
			close($BPP_tmpfile);
			$command = "mv --force $out_dir/bpp $in_dir/../bpp";
			RunCommand($command);
		}
		
		# Output the best maprate and all corresponding statistics
		#$SNR_cutoff		= $SNR_cutoff_tmp[$best_maprate_index];
		$StretchFactor	= $StretchFactor_tmp[$best_maprate_index];
		$FP_100kb		= $FP_100kb_tmp[$best_maprate_index];
		$FNrate			= $FNrate_tmp[$best_maprate_index];
		$siteSD			= $siteSD_tmp[$best_maprate_index];
		$scaling_sd		= $scaling_sd_tmp[$best_maprate_index];
		$bpp			= $bpp_tmp[$best_maprate_index];
		$bppSD			= $bppSD_tmp[$best_maprate_index];
		$FPrate			= $FPrate_tmp[$best_maprate_index];
		$relative_sd	= $relative_sd_tmp[$best_maprate_index];
		$se				= $se_tmp[$best_maprate_index];
		$sf				= $sf_tmp[$best_maprate_index];
		$sr				= $sr_tmp[$best_maprate_index];
		#$LabelDensity	= $LabelDensity_tmp[$best_maprate_index];  # use 0 if ...
		$Maps			= $Maps_tmp[$best_maprate_index];
		$GoodMaps		= $GoodMaps_tmp[$best_maprate_index];
		$TotalAlignedlength = $TotalAlignedlength_tmp[$best_maprate_index];
		}   ##$RefAlignerErr
		#print sprintf("%.2f", $TotalAlignedlength[$i-1] * 1e6 / $sample_range * $total[$i] / $ref_length[$i-1]), "\n";
		
		$command = "mv --force $out_dir/c" . ($c+1) . "_round_$best_maprate_index/ $out_dir/c" . ($c+1) . "_best/";
		RunCommand($command);
	}
	
	my ($stretch_percent, $map_rate);
	if( $not_enough_labels[$c] || ! $length[1] ){
		($avg_mol_snr, $ave_mol_intensity, $LabelDensity) = (0) x 3;
	}
	else{
		$avg_mol_snr = sprintf("%.2f", $mol_SNR{$channel}[1]/$n_mol[1]);
		$ave_mol_intensity = sprintf("%.2f", $mol_int{$channel}[1]/$n_mol[1]);
		$LabelDensity = sprintf( "%.2f", $label{$channel}->[1] * 100000 / $length[1] );
	}
	
	if(! -f $RCmap[$c] || $not_enough_mols[0] || $not_enough_labels[$c] || !$label{$channel}->[1] || $RefAlignerErr ){
		$bpp = 0;
		$stretch_percent = 0;
		$scaling_sd = 0;
		$relative_sd = 0;
		$Maps = 0;
		$GoodMaps = 0;
		$map_rate = 0;
		$FP_100kb = 0;
		$FPrate = 0;
		$FNrate = 0;
	}
	else{
		$stretch_percent = sprintf("%.2f", $bpp?$StretchFactor*500/$bpp*100:0);
		#print "<$stretch_percent>=<$StretchFactor>*500/<$bpp>*100\n";
		$map_rate = sprintf("%.1f", $Maps?$GoodMaps/$Maps*100:0);
	}
	
	push( @json_output, ("{\"populationid\":" . $populationPK[$c] . ",\"quantity\":" . sprintf("%.1f", $length[1]) . ",\"mol_n_50\":" . sprintf("%.1f", Get_N50($length_ref_array[1], $length[1])) . ",\"mol_n_50_ge_20kb\":" . sprintf("%.1f", Get_N50($length_ref_array_ge_20kb[1], $length_ge_20kb[1])) . ",\"avg_mol_snr\":" . $avg_mol_snr . ",\"ave_mol_intensity\":" . $ave_mol_intensity . ",\"avg_label_density\":" . $LabelDensity . ",\"ref_length\":" . $ref_length . ",\"bpp\":" . $bpp . ",\"stretch_percent\":" . $stretch_percent . ",\"channelid\":" . $usechannel[$c] . ",\"scaling_sd\":" . $scaling_sd . ",\"relative_sd\":" . $relative_sd . ",\"se\":" . $se . ",\"sf\":" . $sf . ",\"sr\":" . $sr . ",\"num_maps\":" . $Maps . ",\"num_goodmaps\":" . $GoodMaps . ",\"map_rate\":" . $map_rate . ",\"fp_per_100kb\":" . $FP_100kb . ",\"fp_rate\":" . $FPrate*100 . ",\"fn_rate\":" . $FNrate*100 . ",\"num_molecules\":" . $total[1] . ",\"SNR_cutoff\":" . $SNR_cutoff[$c] . "}") );
	#if( ! -f $RCmap[$c] || ! $length[1] || !$label{$channel}->[1] ){
	#}
	#else{
	#	push( @json_output, ("{\"populationid\":" . $populationPK[$c] . ",\"quantity\":" . sprintf("%.1f", $length[1]) . ",\"mol_n_50\":" . sprintf("%.1f", Get_N50($length_ref_array[1], $length[1])) . ",\"mol_n_50_ge_20kb\":" . sprintf("%.1f", Get_N50($length_ref_array_ge_20kb[1], $length_ge_20kb[1])) . ",\"avg_mol_snr\":" . $avg_mol_snr . ",\"ave_mol_intensity\":" . $ave_mol_intensity . ",\"avg_label_density\":" . $LabelDensity . ",\"ref_length\":" . $ref_length . ",\"bpp\":" . $bpp . ",\"stretch_percent\":" . $stretch_percent . ",\"channelid\":" . $usechannel[$c] . ",\"scaling_sd\":" . $scaling_sd . ",\"relative_sd\":" . $relative_sd . ",\"se\":" . $se . ",\"sf\":" . $sf . ",\"sr\":" . $sr . ",\"num_maps\":" . $Maps . ",\"num_goodmaps\":" . $GoodMaps . ",\"map_rate\":" . $map_rate . ",\"fp_per_100kb\":" . $FP_100kb . ",\"fp_rate\":" . $FPrate*100 . ",\"fn_rate\":" . $FNrate*100 . ",\"num_molecules\":" . $total[1] . "}") );
	#}
	
	if($c % ceil($N_channels/100) == 0){
		printf $OUT_status("%d%% done\n", $c/$N_channels*100);
	}
}

if(@json_output){
	print $OUT_summary("[" . join(",", @json_output) . "]\n");
}
close($OUT_summary);
printf $OUT_status("%d%% done\n", 100);
close($OUT_status);

if(!$debug_mode){
	unlink("$out_dir/SNR");
}

#$datestring = localtime();
#print "Local date and time: $datestring\n";

######################################################
#                    Subroutines                     #
######################################################
sub Init{
	my $opt_string = 'hi:o:l:n:j:P:DQ';
	if(!getopts("$opt_string", \%opts)){
		print STDERR ("ERROR: Invalid parameter(s)! Try -h for more information.\n");
		Usage();
	}
	Usage() if $opts{h};
}

sub Check_inputs{
	my $fin_filename;
	if(!defined($opts{i})){
		print STDERR ("ERROR: Missing parameter (-i)! Try -h for more information.\n");
		Usage();
	}
	if( !($opts{i} =~ /\.bnx$/i) ){
		print STDERR ("ERROR: The input file must have a suffix of \"bnx\"!\n");
		exit (1);
	}
	else{
		if(! -f $opts{i}){
			print STDERR ("ERROR: The input file does not exist!\n");
			exit (1);
		}
		($fin_filename, $in_dir) = fileparse($opts{i});
		$in_dir = abs_path($in_dir);
		$fin = File::Spec->catfile($in_dir, $fin_filename);
		#print "IN: [$in_dir] [$fin_filename] [$fin]\n";
	}
	
	if(defined($opts{o})){
		make_path($opts{o}) unless (-d $opts{o});
		$out_dir = abs_path($opts{o});
		$fout = File::Spec->catfile($out_dir, $fin_filename);
		#print "OUT: [$out_dir] [$fout]\n";
	}
	else{
		$out_dir = $in_dir;
		$fout = $fin;
	}
	
	if(defined($opts{j})){
		$json_string = $opts{j};
	}
	else{
		print STDERR ("ERROR: Missing parameter (-j)! Try -h for more information.\n");
		Usage();
	}
	
	if(defined($opts{l})){
		$Min_length = $opts{l} * 1000;
	}
	if(defined($opts{n})){
		$Min_labels = $opts{n};
	}
	if($opts{D}){
		$debug_mode = 1;
	}
	if($opts{Q}){
		$from_V3_Cohort_QC = 1;
	}
}

sub Usage{
	print << "EOF";

Usage: $^X $0 [Options] <Args>
Options:
  -h    : This help message
  -i    : Input BNX file (Required, can be 1 or 2 channels)
  -o    : Output folder (default: input folder)
  -l    : Minimum length of the molecule in kb (Required)
  -n    : Minimum nicks in each channel of the molecule (Required)
  -j    : Population ID and Callback URL information as JSON string
  -Q    : Calling from V3_Cohort_QC.pl script? (default: OFF)
  -D    : Debug mode to keep all temporary files (default: OFF)
  -P    : Plot SNR histogram in the folder specified (no default)
EOF
	exit (1);
}

sub GetOpts{
	my ($json_string) = @_;
	@RCmap = ("", "");
	
	# For backward compatibility
	if($json_string =~ /^{/){
		$json_string = "[$json_string]";
	}
	my $hash_ref = decode_json($json_string);
	#print Dumper $hash_ref;
	
	$N_channels = scalar(@$hash_ref);
	for(my $i = 0; $i < $N_channels; $i++){
		if( defined($hash_ref->[$i]->{Population}) ){
			$populationPK[$i] = $hash_ref->[$i]->{Population};
		}
		else{
			print STDERR ("ERROR: The JSON string must contain the keyword \"Population\"!\n");
			exit (1);
		}
		if( defined($hash_ref->[$i]->{usechannel}) ){
			$usechannel[$i] = $hash_ref->[$i]->{usechannel};
		}
		else{
			print STDERR ("ERROR: The JSON string must contain the keyword \"usechannel\"!\n");
			exit (1);
		}
		if( defined($hash_ref->[$i]->{RCmap}) ){
			$RCmap[$i] = $hash_ref->[$i]->{RCmap};
		}
		if( defined($hash_ref->[$i]->{subSamplingSeed}) ){
			$subSamplingSeed[$i] = $hash_ref->[$i]->{subSamplingSeed};
		}
		else{
			if( defined($hash_ref->[$i]->{RCmap}) ){
				print STDERR ("ERROR: The JSON string must contain the keyword \"subSamplingSeed\"!\n");
				exit (1);
			}
		}
		if( defined($hash_ref->[$i]->{subSamplingRange}) ){
			$subSamplingRange[$i] = $hash_ref->[$i]->{subSamplingRange};
		}
		else{
			if( defined($hash_ref->[$i]->{RCmap}) ){
				print STDERR ("ERROR: The JSON string must contain the keyword \"subSamplingRange\"!\n");
				exit (1);
			}
		}
		if( defined($hash_ref->[$i]->{minLengthKb}) ){
			$Min_length = $hash_ref->[$i]->{minLengthKb} * 1000;
		}
		else{
			print STDERR ("ERROR: The JSON string must contain the keyword \"minLengthKb\"!\n");
			exit (1);
		}
		if( defined($hash_ref->[$i]->{minLabels}) ){
			$Min_labels = $hash_ref->[$i]->{minLabels};
		}
		else{
			print STDERR ("ERROR: The JSON string must contain the keyword \"minLabels\"!\n");
			exit (1);
		}
	}
}

sub Generate_stat{
	my ($label, $SNR, $mol_SNR, $int, $mol_int, $n_mol, $length, $length_ge_20kb, $length_ref_array, $length_ref_array_ge_20kb, $population) = @_;
	#Generate_stat(\%label, \%SNR, \%mol_SNR, \%int, \%mol_int, \@n_mol, \@length, \@length_ge_20kb, \@length_ref_array, \@length_ref_array_ge_20kb, $population);
	
	if($n_mol->[$population]){
		print $OUT_summary("Number of molecules (filtered):\t", sprintf("%d", $n_mol->[$population]), "\n");
		print $OUT_summary("Total molecule length (bp):\t", sprintf("%.2f", $length->[$population]), "\n");
		print $OUT_summary("Average molecule length (bp):\t", sprintf("%.2f", $length->[$population]/$n_mol->[$population]), "\n");
		print $OUT_summary("Molecule length N50 (bp):\t", sprintf("%.2f", Get_N50($length_ref_array->[$population], $length->[$population])), "\n\n");
	}
	else{
		print $OUT_summary("Number of molecules (filtered):\t", sprintf("%d", 0), "\n");
		print $OUT_summary("Total molecule length (bp):\t", sprintf("%.2f", 0), "\n");
		print $OUT_summary("Average molecule length (bp):\t", sprintf("%.2f", 0), "\n");
		print $OUT_summary("Molecule length N50 (bp):\t", sprintf("%.2f", 0), "\n\n");
	}
	
	for(my $i = 1; $i <= $N_channels; $i++){
		my $channel = "channel_$i";
		if($n_mol->[$population]){
			if($not_enough_labels[$i-1]){
				print $OUT_summary("Average label density for channel $i (/100kb):\t", sprintf("%.2f", 0), "\n");
			}
			else{
				print $OUT_summary("Average label density for channel $i (/100kb):\t", sprintf("%.2f", $label->{$channel}->[$population]*100000/$length->[$population]), "\n");
			}
			if($label->{$channel}->[$population]){
				print $OUT_summary("Average label SNR for channel $i:\t", sprintf("%.2f", $SNR->{$channel}->[$population]/$label->{$channel}->[$population]), "\n");
			}
			else{
				print $OUT_summary("Average label SNR for channel $i:\t", sprintf("%.2f", 0), "\n");
			}
			if( defined($mol_SNR->{$channel}->[$population]) ){
				print $OUT_summary("Average molecule SNR for channel $i:\t", sprintf("%.2f", $mol_SNR->{$channel}->[$population]/$n_mol->[$population]), "\n");
			}
			else{
				print $OUT_summary("Average molecule SNR for channel $i:\t", sprintf("%.2f", 0), "\n");
			}
			if($label->{$channel}->[$population]){
				print $OUT_summary("Average label intensity for channel $i:\t", sprintf("%.2f", $int->{$channel}->[$population]/$label->{$channel}->[$population]), "\n");
			}
			else{
				print $OUT_summary("Average label intensity for channel $i:\t", sprintf("%.2f", 0), "\n");
			}
			if( defined($mol_int->{$channel}->[$population]) ){
				print $OUT_summary("Average molecule intensity for channel $i:\t", sprintf("%.2f", $mol_int->{$channel}->[$population]/$n_mol->[$population]), "\n");
			}
			else{
				print $OUT_summary("Average molecule intensity for channel $i:\t", sprintf("%.2f", 0), "\n");
			}
		}
		else{
			print $OUT_summary("Average label density for channel $i:\t", sprintf("%.2f", 0), "\n");
			print $OUT_summary("Average label SNR for channel $i:\t", sprintf("%.2f", 0), "\n");
			print $OUT_summary("Average molecule SNR for channel $i:\t", sprintf("%.2f", 0), "\n");
			print $OUT_summary("Average label intensity for channel $i:\t", sprintf("%.2f", 0), "\n");
			print $OUT_summary("Average molecule intensity for channel $i:\t", sprintf("%.2f", 0), "\n");
		}
		print $OUT_summary("\n");
	}
}

sub Replace_total_number{
	my ($filename, $total) = @_;
	
	#my $command = "$^X -i.bak -p -e \"s/^# Number of Nanomaps:\\s+\\d+/# Number of Nanomaps:\t$total/\" $filename";
	my $command = "$^X -i.bak -p -e \"s/^# Number of Molecules:\\s+\\d+/# Number of Molecules:\t$total/\" $filename";
	#print "Running command:\n$command\n\n";
	system("$command");
	unlink("$filename.bak");
}

sub Replace_channel_number{
	my ($filename, $total) = @_;
	
	my $command = "$^X -i.bak -p -e \"s/^# Label Channels:\\s+\\d+/# Label Channels:\t$total/\" $filename";
	#print "Running command:\n$command\n\n";
	system("$command");
	unlink("$filename.bak");
}

sub Average{
	my ($array_ref) = @_;
	if(!@$array_ref){
		return 0;
	}
	
	my $total = 0;
	foreach(@$array_ref){
		$total += $_;
	}
	
	return $total / @$array_ref;
}

sub Sum{
	my ($array_ref) = @_;
	if(!@$array_ref){
		return 0;
	}
	
	my $total = 0;
	foreach(@$array_ref){
		$total += $_;
	}
	return $total;
}

sub Get_N50{
	my ($data, $total_len) = @_;
	if(!$total_len){
		return 0;
	}
	
	my @sorted = sort {$b <=> $a} @$data;
	
	my $total = 0;
	my $i;
	
	for($i = 0; $i < @sorted; $i++){
		$total += $sorted[$i];
		if($total >= $total_len / 2){
			return $sorted[$i];
		}
	}
}

sub StdDev{
	my ($array_ref) = @_;
	my ($sum, $sum2) = (0, 0);
	
	if(!@$array_ref){
		return 0;
	}
	
	foreach my $item (@$array_ref){
		$sum += $item;
	}
	my $mean = $sum / @$array_ref;
	
	foreach my $item (@$array_ref){
		$sum2 += ($item-$mean)**2;
	}
	
	return sqrt($sum2/@$array_ref);
}

sub RunCommand{
	my ($command, $err_msg) = @_;
	$err_msg ||= "ERROR: Cannot excute the command:\n$command\n";
	if(system($command)){
		print STDERR ($err_msg);
	}
}

sub Send_Complete_Status{
	my ($json_string, $url, $msg) = @_;
	
	#my $socket = PocketIO::Client::IO->connect("http://$IP:$port");
	#$socket->on( 'message', sub {
	#	say $_[1];
	#} );
	#$socket->on( 'connect', sub {
	#	$socket->send('Hello!');
	#	$socket->emit( $channelName, { msg => '100% done' } );
	#	print "Connected!";
	#} );
	
	use IO::Socket::INET;
	# Auto-flush on socket
	$| = 1;
	
	# Add minLen and minLabels!
	# {"Population":238,"subSamplingRange":1000,"subSamplingSeed":17,"prepCmaps":{"Color":"Red","Enzyme":"BSPQI","RCmap":"/mnt/genomaton/Files/IrysView/hg19_NT.BSPQI_0kb_0labels.cmap"},"url":"192.168.48.129:3009"}
	my $Min_length_kb = $Min_length / 1000;
	$json_string =~ s/^\s*{/{"minLen":$Min_length_kb,"minLabels":$Min_labels,/;
	
	my $data = "{\"subject\":\"cohortQC\", \"input\":$json_string, \"output\":$msg}";
	
	# Create socket connection
	my $socket = new IO::Socket::INET(
		PeerAddr => $url,
		#PeerPort => $port,
		Proto => 'tcp',
	) || die "ERROR: Could not create socket $url: $!\n";
	print "Connected to the server $url\n";
	
	# Send data to a server
	my $size = $socket->send($data);
	print "Sent data of length $size\n" . $data . "\n";
	#sleep(1);
	$socket->close();
}

__END__

perl MQR.pl -i /home/users3/xzhou/V3/BspQI_s.bnx -j '{"Population":238,"subSamplingRange":1000,"subSamplingSeed":1,"usechannel":1,"minLengthKb":150,"minLabels":9,"RCmap":""}' -o MQR_1C_BNX_0_REF
perl MQR.pl -i /home/users3/xzhou/V3/BspQI_s.bnx -j '{"Population":238,"subSamplingRange":1000,"subSamplingSeed":1,"usechannel":1,"minLengthKb":150,"minLabels":9,"RCmap":""}' -o MQR_1C_BNX_0_REF -Q

perl MQR.pl -i /home/users3/xzhou/V3/BspQI_s.bnx -j '{"Population":238,"subSamplingRange":1000,"subSamplingSeed":1,"usechannel":1,"minLengthKb":150,"minLabels":9,"RCmap":"/mnt/genomaton/Files/IrysView/hg19_NT.BSPQI_0kb_0labels.cmap"}' -o MQR_1C_BNX_1_REF
perl MQR.pl -i /home/users3/xzhou/V3/BspQI_s.bnx -j '[{"Population":238,"subSamplingRange":1000,"subSamplingSeed":1,"usechannel":1,"minLengthKb":150,"minLabels":9,"RCmap":"/mnt/genomaton/Files/IrysView/hg19_NT.BSPQI_0kb_0labels.cmap"}]' -o MQR_1C_BNX_1_REF

perl MQR.pl -i /home/users3/xzhou/PERL/SVN/AnalysisTools/cohortQC/BSSSI_BSPQI.bnx -j '[{"Population":238,"subSamplingRange":1000,"subSamplingSeed":2,"usechannel":2,"minLengthKb":150,"minLabels":9,"RCmap":"/mnt/genomaton/Files/IrysView/hg19_NT.BSPQI_0kb_0labels.cmap"}]' -o MQR_2C_BNX_1_REF

perl MQR.pl -P SNR.pdf -i /home/users3/xzhou/PERL/SVN/AnalysisTools/cohortQC/BSSSI_BSPQI.bnx -j '[{"Population":238,"subSamplingRange":1000,"subSamplingSeed":1,"usechannel":1,"minLengthKb":150,"minLabels":9,"RCmap":"/mnt/genomaton/Files/IrysView/hg19_NB.BSSSI_0kb_0labels.cmap"},{"Population":239,"subSamplingRange":1000,"subSamplingSeed":2,"usechannel":2,"minLengthKb":150,"minLabels":9,"RCmap":"/mnt/genomaton/Files/IrysView/hg19_NT.BSPQI_0kb_0labels.cmap"}]' -o MQR_2C_BNX_2_REF_P

perl MQR.pl -i /home/users3/xzhou/PERL/SVN/AnalysisTools/cohortQC/BSSSI_BSPQI.bnx -j '[{"Population":238,"subSamplingRange":1000,"subSamplingSeed":1,"usechannel":1,"minLengthKb":150,"minLabels":9,"RCmap":""},{"Population":238,"subSamplingRange":1000,"subSamplingSeed":2,"usechannel":2,"minLengthKb":150,"minLabels":9,"RCmap":"/mnt/genomaton/Files/IrysView/hg19_NT.BSPQI_0kb_0labels.cmap"}]' -o MQR_2C_BNX_missing_REF


