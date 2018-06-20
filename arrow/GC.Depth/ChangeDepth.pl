#!/bin/perl

open IN,$ARGV[0];
open IN1,$ARGV[1];

my %hash;
my $new = 1;
my $pre = 0;
my $cur = 0;
my $pre_id = 0;
my $cur_id = 0;

while(<IN>)
{
	chomp;
	@t = split;
	$hash{$t[0]} = $t[1];
}
close IN;
while(<IN1>)
{
	chomp;
	@t1 = split;
#	$cur_id = $t1[0];
	if($new == 1)
	{
#		$pre = $t1[1];
#		$pre_id = $t1[0];
		if($t1[1] == 1)
		{
			print $_."\n";
		}
		else
		{
			for(my $i=$pre+1;$i<$t1[1];$i++)
			{
				print $t1[0]."\t".$i."\t"."0"."\n";
			}
			print $_."\n";
		}
		$pre = $t1[1];
		$pre_id = $t1[0];
		$new = 0;
	}
	else
	{
		if($pre_id eq $t1[0])
		{

			if($t1[1] - $pre == 1)
			{
				print $_."\n";
			}
			else
			{
				for(my $i=$pre+1;$i<$t1[1];$i++)
				{
					print $t1[0]."\t".$i."\t"."0"."\n";
				}
				print $_."\n";
			}
			$pre =  $t1[1];
			$pre_id = $t1[0];	
		}
		else
		{
			if($hash{$pre_id} != $pre)
			{
				for(my $i=$pre+1;$i<=$hash{$pre_id};$i++)
				{
					print $pre_id."\t".$i."\t"."0"."\n";
				}
			}
			$pre = 0;
			if($t1[1] == 1)
			{
				print $_."\n";
			}
			else
			{
				for(my $i=$pre+1;$i<$t1[1];$i++)
				{
					print $t1[0]."\t".$i."\t"."0"."\n";
				}
				print $_."\n";
			}
			$pre = $t1[1];
			$pre_id = $t1[0];
			
			
		}
	}
	
}
close IN1;

