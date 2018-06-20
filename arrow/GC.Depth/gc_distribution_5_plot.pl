#!/usr/bin/perl -w
if(@ARGV>12||@ARGV==0){
	print "the max genome to plot GC distribution is 5.\n";
	print "Usage: perl $0 <genome_file><genome_name>...<genome_file><genome_name><output_file>\n";
	exit;
}

$out_file=$ARGV[-1];
$bin=500;
$step=$bin/2;
$gcn=$bin;
$flag_t=0;

%color =(
	0 => "red",
	2 => "blue",
	4 => "orange",
	6 => "black",
	8 => "green",
);

open OUTPUT,">$out_file";
for($j=0;$j<$#ARGV;$j+=2){
   $file=$ARGV[$j];
   $name=$ARGV[$j+1];
   $flag=0;$seq="";$num=0;$sum=0;$total=0;$mean=0;$total_num_of_bin=0;$max=0;
   
   for($i=0;$i<=$gcn;$i++){
      $hash{$i/$gcn}=0;
   }

   open INPUT,"$file" or
      die "can't open the file: $!\n";
   while($line=<INPUT>){
      chomp($line);
      if($line=~/>/){
         if($flag==0){$flag=1;}
         else{
            $seq_t=$seq;
            $total+=length($seq_t);
            $L1=length($seq_t);
            $seq_t=~tr/gcGC//d;
            $L2=length($seq_t);
            $seq_t=~tr/atAT//d;
            $L3=length($seq_t);
            $num+=($L1-$L2);
            $sum+=($L1-$L3);

            while(length($seq)>=$bin){
               $temp=substr($seq,0,$bin);
               $L1=length($temp);
               $temp=~tr/gcGC//d;
               $L2=length($temp);
               $temp=~tr/atAT//d;
               $L3=length($temp);
               if(($L1-$L3)>=$step)
                  {$gc=int(($L1-$L2)*$gcn/($L1-$L3))/$gcn;$hash{$gc}++;$total_num_of_bin++;}
               substr($seq,0,$step)="";
            }
            $seq="";
         }
      }
      else{
         $seq.=$line;
         #print @item;
      }
   }
   close INPUT;

   $seq_t=$seq; 
   $total+=length($seq_t);
   $L1=length($seq_t);
   $seq_t=~tr/gcGC//d;
   $L2=length($seq_t);
   $seq_t=~tr/atAT//d;
   $L3=length($seq_t);
   $num+=($L1-$L2);
   $sum+=($L1-$L3);

   while(length($seq)>=$bin){
      $temp=substr($seq,0,$bin);
      $L1=length($temp);
      $temp=~tr/gcGC//d;
      $L2=length($temp);
      $temp=~tr/atAT//d;
      $L3=length($temp);
      if(($L1-$L3)>=$step)
         {$gc=int(($L1-$L2)*$gcn/($L1-$L3))/$gcn;$hash{$gc}++;$total_num_of_bin++;}
      substr($seq,0,$step)="";
   }
   $seq="";
   ################################################################################
   $total=int($total/1000000);

   for($i=0;$i<=$gcn;$i++){
      #print STDERR "\n$total_num_of_bin";
      if ($total_num_of_bin>0){
      $hash{$i/$gcn}/=$total_num_of_bin};
      if($hash{$i/$gcn}>$max){
   	     $max=$hash{$i/$gcn};
      }
   }

   $max*=(5/4);
   $max=int($max/0.005+1)*0.005;
   $mean=int(($num/$sum)*1000)/1000;

   $name=~s/_/ /;
   #################################
   if($flag_t==0){
print OUTPUT "Type:Line
Width:1200
Height:800
WholeScale:0.8
#MarkPos:tr/tl/br/bl
#MarkScale:1
#MarkNoBorder:1
#FontSize:46
#FontFamily:ArialNarrow-Bold
Note: GC_distribution
X: GC Content, percent 
Y: percent(%) of bins, bin=$bin bp
Xstart: 0
Xstep: 0.1
Xend: 0.8
Xcut: 1
#XScale:
#:End
Ystart: 0
Ystep: 0.005
Yend:  $max
YDiv: 0.01
#NoYScale:1
#Ycut: 1
#YScale:
#:End
XScaleDiv: 5
#YScaleDiv: 5 
#Note2:
#add annotation at the bottom
#:End
:End";
   $flag_t=1;
}
   print OUTPUT "\n\nColor: $color{$j}\nMark: $name, $total Mb, avg_gc: $mean\n";
   for($i=0;$i<=$gcn;$i++){
      #print $hash{$i/$bin};
      print OUTPUT $i/$gcn.": ".$hash{$i/$gcn}."\n";
   }
}
close OUTPUT;

#system("perl /nas/GAG_01B/assembly/pig/list/distributing_svg_4.74/distributing_svg.pl $out_file $out_file.svg");
#system("perl /home/tiany/bin/GC_annalysis/GC_distribution/distributing_svg.pl $out_file $out_file.svg");
#system("/nas/GAG_01B/assembly/pig/list/distributing_svg_4.74/svg2xxx_release/svg2xxx -t png $out_file.svg");
