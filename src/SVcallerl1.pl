#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my ($input,$output,$refseq,$refname,$queryname,$Help);

GetOptions (
"input:s" =>\$input,
"output:s" =>\$output,
"refseq:s" =>\$refseq,
"refname:s" =>\$refname,
"queryname:s"=>\$queryname,
"help" =>\$Help
);

if ($Help) {
print << "END.";
Useage: perl $0 -input synnet.txt -output level -refname Dmel -queryname Dyakuba
END.
exit;
}

my $seq = ExtractSeq ($refseq);

open (IN,"$input") || die "Cannot open $!";
open (OUTind1,">$output.Deletion.txt") || die "Cannot open $!\n";
open (OUTcxs1,">$output.Complex.txt") || die "Cannot open $!\n";

$/="fill";
my $target;

while (<IN>) {
next if ($_=~/^$/);
my @temp = split (/\n/,$_);
my $fill = shift @temp;
pop @temp;
if ($fill =~/net/) {
my @tmp = split (/\s/,$fill);
$target = $tmp[1];
next;
}

my @fills = split (/\s/,$fill);

my $fill_s = $fills[1];
my $fill_e = $fills[1]+$fills[2];

if ($fills[2]<10000) {
next;
} else {
my @unit0 = split (/\s/,$temp[0]);
my @unit1 = split (/\s/,$temp[1]);
my @unit_end = split (/\s/,$temp[$#temp]);
if ($unit0[6]==0 or ($unit0[6]/$unit0[2] < 0.1)) {
my $Tend = $unit0[1]+$unit0[2];
my $Qend = $unit0[5]+$unit0[6];
my $forward = $unit0[1]-$fill_s;
my $backward = $unit1[1]-$Tend;
my $Sequence = substr (${$seq}{$target},$unit0[1],$unit0[2]);
print OUTind1 "$refname.$target\t$unit0[1]\t$Tend\tDeletion\t$queryname.$unit0[3]\t$unit0[5]\t$Qend\t$unit0[2]\t$unit0[4]\t$forward\:$backward\t$Sequence\n";
} else {
my $Tend = $unit0[1]+$unit0[2];
my $Qend = $unit0[5]+$unit0[6];
my $forward = $unit0[1]-$fill_s;
my $backward = $unit1[1]-$Tend;
my $Sequence = substr (${$seq}{$target},$unit0[1],$unit0[2]);
print OUTcxs1 "$refname.$target\t$unit0[1]\t$Tend\tComplex\t$queryname.$unit0[3]\t$unit0[5]\t$Qend\t$unit0[2]\t$unit0[4]\t$forward\:$backward\t$Sequence\n";
}

###### The middle elements
my $i;
for ($i=1;$i<$#temp;$i++) {

my @forward1 = split (/\s/,$temp[$i-1]);
my @now = split (/\s/,$temp[$i]);
my @backward1 = split (/\s/,$temp[$i+1]);
my $Tend = $now[1]+$now[2];
my $Qend = $now[5]+$now[6];
my $forward = $now[1]-$forward1[1]-$forward1[2];
my $backward = $backward1[1]-$Tend;
my $Sequence = substr (${$seq}{$target},$now[1],$now[2]);
if ($now[6]==0 or ($now[6]/$now[2] < 0.1)) {
print OUTind1 "$refname.$target\t$now[1]\t$Tend\tDeletion\t$queryname.$now[3]\t$now[5]\t$Qend\t$now[2]\t$now[4]\t$forward\:$backward\t$Sequence\n";
} else {
print OUTcxs1 "$refname.$target\t$now[1]\t$Tend\tComplex\t$queryname.$now[3]\t$now[5]\t$Qend\t$now[2]\t$now[4]\t$forward\:$backward\t$Sequence\n";
}
}

###### The last element
my @end_1 = split(/\s/,$temp[$#temp-1]);

if ($unit_end[6]==0 or ($unit_end[6]/$unit_end[2] < 0.1)) {
my $Tend = $unit_end[1]+$unit_end[2];
my $Qend = $unit_end[5]+$unit_end[6];
my $forward = $unit_end[1]-$end_1[1]-$end_1[2];
my $backward = $fill_e-$Tend;
my $Sequence = substr (${$seq}{$target},$unit_end[1],$unit_end[2]);
print OUTind1 "$refname.$target\t$unit_end[1]\t$Tend\tDeletion\t$queryname.$unit_end[3]\t$unit_end[5]\t$Qend\t$unit_end[2]\t$unit_end[4]\t$forward\:$backward\t$Sequence\n";
} else {
my $Tend = $unit_end[1]+$unit_end[2];
my $Qend = $unit_end[5]+$unit_end[6];
my $forward = $unit_end[1]-$end_1[1]-$end_1[2];
my $backward = $fill_e-$Tend;
my $Sequence = substr (${$seq}{$target},$unit_end[1],$unit_end[2]);
print OUTcxs1 "$refname.$target\t$unit_end[1]\t$Tend\tComplex\t$queryname.$unit_end[3]\t$unit_end[5]\t$Qend\t$unit_end[2]\t$unit_end[4]\t$forward\:$backward\t$Sequence\n";
     }
  }
}



###############################
########### Sub modules #######
###############################

sub ExtractSeq {
   my ($file) = @_;
   my %hash;
   $/=">";
   open In, "$file" or die "$!";
      while (<In>) {
         next if (length $_ < 2);
         my @unit= split ("\n", $_);
         my $head0=shift @unit;
         my @head1= split (/\s+/,$head0);
         my $head= shift @head1;
         my $seq = join ("", @unit);
         $hash{$head}=$seq;
      }
     my $ref=\%hash;
     return $ref;
 }
