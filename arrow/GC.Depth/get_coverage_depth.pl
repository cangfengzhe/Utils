#!/usr/bin/perl -w
use strict;
use warnings;

die "usage:\n\tperl $0 <scaffold_len_list>\t<||10000>\t<scaffold_seq>\t<soap_result_list(can be gz file)>\n" unless (@ARGV == 4);

my %hash=();
my %len=();
my $bin=$ARGV[1]||1000;

open IN,"$ARGV[0]";
while(my $line=<IN>){
   chomp($line);
   my @item=split(/\s+/,$line);
   if($item[1]>=$bin){
      for(my $i=1;$i<=$item[1]/$bin;$i++){
         $hash{$item[0]}[$i]=0;
         $len{$item[0]}[$i]=0;
      }
   }
}
close IN;

print STDOUT "############ step 1 finished ############\n";

my $flag=0;
my $seq ;
my $name ;
open FILE,"$ARGV[2]";
while(my $seq_line=<FILE>){
   chomp($seq_line);
   if($seq_line=~/>/){
      if($flag==1){
         my @seq = split(//, $seq);
         for(my $i=1;$i*$bin<=$#seq+1;$i++){
            my $num=0;
            for(my $j=($i-1)*$bin;$j<$i*$bin;$j++){
               if($seq[$j]=~/[atcgATCG]/) {$num++;}
            }
            $len{$name}[$i]=$num;
         }
      }
      $seq_line=~/>(\w+)/;$name=$1;$flag=0;
      if(defined($len{$name})) {$flag=1;$seq="";}
   }elsif($flag==1){
       $seq.=$seq_line;
   }
}
close FILE;

print STDOUT "############ step 2 finished ############\n";

#use lib "/ifs1/RD/chenshsh/bin/perl/lib/PerlIO/lib64/perl5/site_perl/5.8.5/x86_64-linux-thread-multi";
use PerlIO::gzip;
open ININ,"$ARGV[3]";
while(my $list=<ININ>){
   chomp($list);
   if ($list =~ /.gz$/){
       open FILE,"<:gzip","$list";
   }else{
       open FILE,"<$list";
   }
   while(my $line=<FILE>){
      chomp($line);
      my @item=split(/\s+/,$line);
      my $start=$item[8]; my $end=$item[8]+$item[5]-1;
      for(my $i=$start;$i<=$end;$i++){
         my $pos=int(($i-1)/$bin)+1;
         if(defined($hash{$item[7]}[$pos]))
            {$hash{$item[7]}[$pos]++;}
      }
   }
   close FILE;
}
close ININ;

print STDOUT "now print them into the file.\n";

open IN,"$ARGV[0]";
open OUT,">coverage_depth";
while(my $line=<IN>){
   chomp($line);
   my @item=split(/\s+/,$line);
   if($item[1]>=$bin){
      print OUT "$item[0]\n";
      for(my $i=1;$i<=$item[1]/$bin;$i++){
         if($len{$item[0]}[$i]!=0){
            print OUT $hash{$item[0]}[$i]/$len{$item[0]}[$i]."\t";
         }else {
            print OUT "0\t";
         }
      }
      print OUT "\n";
   }
}
close IN;
close OUT;
