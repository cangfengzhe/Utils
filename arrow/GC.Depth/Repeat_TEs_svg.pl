#!/usr/bin/perl -w

=head1

Name: TEsvg.pl
perl TEsvg.pl *.cat genome_length
=cut
use lib "/export/personal/lijj/0.temp/1.script/GC.Depth/";

use strict;
use SVG;
use Getopt::Long;

my ($x1_border,$x2_border,$y1_border,$y2_border,$x_start,$x_end,$y_start,$y_end,$x_step,$y_step,$Help);
GetOptions(
	"--x_left_border:s"=>\$x1_border,
	"--x_right_border:s"=>\$x2_border,
	"--y_up_border:s"=>\$y1_border,
	"--y_down_border:s"=>\$y2_border,
	"--x_start:s"=>\$x_start,
	"--y_start:s"=>\$y_start,
	"--x_end:s"=>\$x_end,
	"--y_end:s"=>\$y_end,
	"--x_step:s"=>\$x_step,
	"--y_step:s"=>\$y_step,
	"--help"=>\$Help
);
die `pod2text $0` if @ARGV < 2 or $Help;
$x1_border ||= 120;
$x2_border ||= 620;
$y1_border ||= 50;
$y2_border ||= 450;
my $svg = SVG->new('width',$x2_border+60,'height',$y2_border+90);
$x_start ||= 0;
$x_end ||= 40;
$y_start ||= 0;
$y_end ||= 0.5;
$x_step ||= 10;
$y_step ||= 0.5;
my $file = shift;
my $total_len = shift;
open IN,"<$file" or die "fail to open $file:$!\n";
my %hash;
while(<IN>) {
#	    score  div               scaffold    begin      end                       repeat class
	if(/\d+\s+(\S+)\s+\S+\s+\S+\s+ \S+   \s+ (\d+) \s+ (\d+)  \s+\S+\s+\S+\s+\S+\s+   (\S+)/x) {
		my $div = int $1;
		next unless($div >0);
		my $len = $3-$2+1;
		my $class = $4;
#        print STDERR "$1\t$2\t$3\t$4\n";
		$hash{$div}{"DNA"} += $len if($class =~ /DNA/);
		$hash{$div}{"LINE"} += $len if($class =~ /LINE/);
		$hash{$div}{"LTR"} += $len if($class =~ /LTR/);
		$hash{$div}{"SINE"} += $len if($class =~ /SINE/);
	}
}
close IN;
my $x_unit;my $y_unit;my $x_covert;my $y_covert;my @colors;my @elements;my $x_number;my $y_number;my $tmp_i;my ($x,$y);


while($y_end){
	$x_unit = ($x2_border-$x1_border)/(($x_end-$x_start)*2);
	$y_unit = ($y2_border-$y1_border)/20;
	$x_covert = $x_unit*2;
	$y_covert = ($y2_border-$y1_border)/($y_end-$y_start);
	@colors = qw/green yellow black red blue /;
	@elements = qw/DNA LINE LTR SINE/;
	$x_number = ($x_end-$x_start)/$x_step;
	$y_number = int(($y_end-$y_start)/$y_step);
	$tmp_i = 0;
	($x,$y) = (-($y1_border+($y2_border-$y1_border)/2),$x1_border-80);
	my $y_sine;my $y_ltr;my $y_line;my $y_dna;my $y_sine_covert;my $y_ltr_covert;my $y_line_covert;my $y_dna_covert;my @value;
	foreach(1..$x_end) 
	{
		$y_sine = exists $hash{$_}{"SINE"} ? $hash{$_}{"SINE"} : 0;
		$y_ltr = exists $hash{$_}{"LTR"} ? $hash{$_}{"LTR"} : 0;
		$y_line = exists $hash{$_}{"LINE"} ? $hash{$_}{"LINE"} : 0;
		$y_dna = exists $hash{$_}{"DNA"} ? $hash{$_}{"DNA"} : 0;
		$y_sine_covert = ($y_sine*$y_covert*100)/$total_len;
		$y_ltr_covert = ($y_ltr*$y_covert*100)/$total_len;
		$y_line_covert = ($y_line*$y_covert*100)/$total_len;
		$y_dna_covert = ($y_dna*$y_covert*100)/$total_len;
		my $value = $y2_border-$y_sine_covert-$y_ltr_covert-$y_line_covert-$y_dna_covert-$y_dna_covert;
		push @value,$value;
	}
	my @value1 = sort{$a<=>$b} @value;
	if ($value1[0] < 180)
	{
		$y_end = $y_end+0.5;
	}
	else 
	{
		foreach(1..$x_end){
		$y_sine = exists $hash{$_}{"SINE"} ? $hash{$_}{"SINE"} : 0;
		$y_ltr = exists $hash{$_}{"LTR"} ? $hash{$_}{"LTR"} : 0;
		$y_line = exists $hash{$_}{"LINE"} ? $hash{$_}{"LINE"} : 0;
		$y_dna = exists $hash{$_}{"DNA"} ? $hash{$_}{"DNA"} : 0;
		$y_sine_covert = ($y_sine*$y_covert*100)/$total_len;
		$y_ltr_covert = ($y_ltr*$y_covert*100)/$total_len;
		$y_line_covert = ($y_line*$y_covert*100)/$total_len;
		$y_dna_covert = ($y_dna*$y_covert*100)/$total_len;
		$svg->rect('x',$x1_border+(2*$_-1)*$x_unit,'y',$y2_border-$y_sine_covert,'width',$x_unit,'height',$y_sine_covert,'stroke',$colors[0],'fill',$colors[0]);
		$svg->rect('x',$x1_border+(2*$_-1)*$x_unit,'y',$y2_border-$y_sine_covert-$y_ltr_covert,'width',$x_unit,'height',$y_ltr_covert,'stroke',$colors[1],'fill',$colors[1]);
		$svg->rect('x',$x1_border+(2*$_-1)*$x_unit,'y',$y2_border-$y_sine_covert-$y_ltr_covert-$y_line_covert,'width',$x_unit,'height',$y_line_covert,'stroke',$colors[2],'fill',$colors[2
]);
		$svg->rect('x',$x1_border+(2*$_-1)*$x_unit,'y',$y2_border-$y_sine_covert-$y_ltr_covert-$y_line_covert-$y_dna_covert,'width',$x_unit,'height',$y_dna_covert,'stroke',$colors[3],'fill',$colors[3]);}
		last;
	}
}
$svg->rect('x',$x1_border,'y',$y1_border,'width',($x2_border-$x1_border),'height',($y2_border-$y1_border+5),'stroke','black','fill','none','stroke-width',2);
$svg->rect('x',$x2_border-$x_unit*16-10+40,'y',$y1_border+5,'width',$x_unit*16+5-40,'height',$y_unit*4.8+10,'stroke','black','fill','none','stroke-width',2);
foreach (1..4) {
	$svg->rect('x',$x2_border-$x_unit*16-10+2+40,'y',$y1_border+5+2*$_+($_-1)*1.2*$y_unit,'width',3*$x_unit,'height',1.2*$y_unit,'stroke','black','fill',$colors[4-$_]);
	$svg->text('x',$x2_border-$x_unit*16-10+2+3*$x_unit+5+40,'y',$y1_border+5+2*$_+$_*1.2*$y_unit-5,'font-size','18','font-family','ArialNarrow-Bold','fill','black','-cdata',$elements[$_-1]);
}

foreach (0..$x_number) {
	$svg->line('x1',$x1_border+$_*$x_step*$x_covert,'y1',$y2_border+5,'x2',$x1_border+$_*$x_step*$x_covert,'y2',$y2_border+10,'stroke','black','fill','black','stroke-width',2);
	$svg->text('x',$x1_border+$_*$x_step*$x_covert,'y',$y2_border+45,'font-size','26','font-family','ArialNarrow-Bold','text-anchor','middle','fill','black','-cdata',$_*$x_step);
}

while ($tmp_i <= $y_number) {
	print STDERR "$tmp_i\n";
	$svg->line('x1',$x1_border,'y1',$y2_border-$tmp_i*$y_step*$y_covert,'x2',$x1_border-5,'y2',$y2_border-$tmp_i*$y_step*$y_covert,'stroke','black','fill','black','stroke-width',2);
	$svg->text('x',$x1_border-8,'y',$y2_border-$tmp_i*$y_step*$y_covert+16,'font-size','26','font-family','ArialNarrow-Bold','text-anchor','end','fill','black','-cdata',$tmp_i*$y_step);	
	$tmp_i++;
}
$svg->text('x',$x1_border+($x2_border-$x1_border)/2,'y',$y2_border+80,'font-size','26','font-family','ArialNarrow-Bold','text-anchor','middle','fill','black','-cdata',"Sequence divergence rate(%)");
my $g = $svg->group("transform"=>"rotate(-90)");
$g->text('x',$x,'y',$y,'font-size','26','font-family','ArialNarrow-Bold','text-anchor','middle','fill','black','-cdata','Percentage of genome');
print $svg->xmlify();
