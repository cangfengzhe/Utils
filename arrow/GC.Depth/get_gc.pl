#!/usr/bin/perl -w


die "usage:\n\tperl $0 <*.scaffold>\t<10000>\t >gc_dis\n" unless (@ARGV ==2);
$file=$ARGV[0];
$bin=$ARGV[1];
$step=$bin;

$flag=0;$temp="";$num=0;
open INPUT,"$file" or die "can't open the file: $!\n";
while($line=<INPUT>)
{
   chomp($line);
   if($line=~/>/)
   {
      if($flag==0){$flag=1;@name=split(/>/,$line);}
      else
      {
         @item=split(//,$temp);
         if(($#item+1)>=$bin)
         {
            print "$name[1]\n";
            for($i=0;$i+$bin<=($#item+1);$i=$i+$step)
            {
               $num=0;$sum=0;
               for($j=$i;$j<$i+$bin;$j++)
               {
                  if($item[$j]=~/[CGcg]/)
                     {$num++;}
                  if($item[$j]=~/[ATCGatcg]/)
                     {$sum++;}
               }
               if($sum!=0)
                  {print  $num/$sum."\t";}
               else
                  {print "0\t";}
            }
            print "\n";
         }
         $temp="";
         @name=split(/>/,$line);
      }
   }
   else
   {
      $temp.=$line;
   }
}
close INPUT;

@item=split(//,$temp);
if(($#item+1)>=$bin)
{
   print "$name[1]\n";
   for($i=0;$i+$bin<=($#item+1);$i=$i+$step)
   {
      $num=0;$sum=0;
      for($j=$i;$j<$i+$bin;$j++)
      {
         if($item[$j]=~/[CGcg]/)
            {$num++;}
         if($item[$j]=~/[ATCGatcg]/)
            {$sum++;}
      }
      if($sum!=0)
         {print  $num/$sum."\t";}
      else
         {print "0\t";}
   }
}
$temp="";
