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
open (OUTcnv1,">$output.CNV.txt") || die "Cannot open $!\n";
open (OUTinv1,">$output.INV.txt") || die "Cannot open $!\n";


$/="gap";
my $target;

while (<IN>) {
next if ($_=~/^$/);
my @temp = split (/\n/,$_);
my $gap = shift @temp;
my $last = pop @temp;

if ($gap =~/net/) {
my @tmp = split (/\s/,$gap);
$target = $tmp[1];
next;
}

if ($#temp<0) {
next;
}

my @gaps = split (/\s/,$gap);


if ($gaps[2] eq "") {
next;
}
my $Tend = $gaps[1]+$gaps[2];

if ($gaps[6]/$gaps[2]<0.1 and $gaps[6]<25) {
my $ali=0; 
my @coordinates=();
my $inv_len=0;

foreach my $fill (@temp) {

my @unit=split(/\s/,$fill);
my $fills=join("_",@unit);

my $hash = ParseFill ($fills);

if ($$hash{"type"}=~/syn|inv/ and $$hash{"qFar"} < 2*$gaps[2]) { # Define copy number variations
            $inv_len = $inv_len + $$hash{"ali"} if ($$hash{"type"} =~/inv/);
            $ali= $ali + $$hash{"ali"};
            my $end = $$hash{"query_beg"} + $$hash{"query_len"};
            push (@coordinates,$$hash{"query_beg"});
            push (@coordinates,$end);
            }
		}
my @array = sort {$a <=> $b} @coordinates;
my $Sequence = substr (${$seq}{$target},$gaps[1],$gaps[2]);

if ($ali/$gaps[2] > 0.25) { #$array[0] < $gaps[5]
print OUTcnv1 "$refname.$target\t$gaps[1]\t$Tend\tCNV\t$queryname.$gaps[3]\t$array[0]\t$array[-1]\t$gaps[2]\t$gaps[4]\tforward\:backward\t$Sequence\n";
    }
}

################
 else {

my $ali=0; 
my @coordinates=();
my $inv_len=0;
my $hash;
foreach my $fill (@temp) {

my @unit=split(/\s/,$fill);
my $fills=join("_",@unit);

$hash = ParseFill ($fills);

 if ($$hash{"type"}=~/syn|inv/ and $$hash{"qFar"} < 2*$gaps[2]) { # Define copy number variations
            $inv_len = $inv_len + $$hash{"ali"} if ($$hash{"type"} =~/inv/);
            $ali= $ali + $$hash{"ali"};
            my $end = $$hash{"query_beg"} + $$hash{"query_len"};
            push (@coordinates,$$hash{"query_beg"});
            push (@coordinates,$end);
            }
		}

  my @array = sort {$a <=> $b} @coordinates;
  my $Sequence = substr (${$seq}{$target},$gaps[1],$gaps[2]); 
 if ( $inv_len > ($gaps[2]+$gaps[6])/5 and $inv_len < $gaps[6]) {
         print OUTinv1 "$refname.$target\t$gaps[1]\t$Tend\tINV\t$queryname.$gaps[3]\t$array[0]\t$array[-1]\t$gaps[2]\t$gaps[4]\tforward\:backward\t$Sequence\n";
 } elsif ($ali/$gaps[2] >0.2) { #$array[0] < $gaps[5]
         print OUTcnv1 "$refname.$target\t$gaps[1]\t$Tend\tCNV\t$queryname.$gaps[3]\t$array[0]\t$array[-1]\t$gaps[2]\t$gaps[4]\tforward\:backward\t$Sequence\n";
    }
  }
}



####################################################
################### Sub Routines ###################
####################################################

sub ParseFill {
my ($line) = @_;
my %hash;
my @temp = split(/_/,$line);
print "Query name: $temp[3]\n";
if ($line=~/inv|syn/) {
$hash{"target_beg"} =$temp[1];
$hash{"target_len"} =$temp[2];
$hash{"query_name"} =$temp[3];
$hash{"query_orient"}=$temp[4];
$hash{"query_beg"}=$temp[5];
$hash{"query_len"}=$temp[6];
$hash{"id"}=$temp[8];
$hash{"score"}=$temp[10];
$hash{"ali"}=$temp[12];
$hash{"qOver"}=$temp[14];
$hash{"qFar"}=$temp[16];
$hash{"qDup"}=$temp[18];
$hash{"type"}=$temp[20];
} elsif ($line=~/nonSyn/) {
$hash{"target_beg"} =$temp[1];
$hash{"target_len"} =$temp[2];
$hash{"query_name"} =$temp[3];
$hash{"query_orient"}=$temp[4];
$hash{"query_beg"}=$temp[5];
$hash{"query_len"}=$temp[6];
$hash{"id"}=$temp[8];
$hash{"score"}=$temp[10];
$hash{"ali"}=$temp[12];
$hash{"qDup"}=$temp[14];
$hash{"type"}=$temp[16];
}
my $ref=\%hash;
return $ref;
}

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
