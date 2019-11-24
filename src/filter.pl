#!/usr/bin/perl -w

use strict;
use warnings;

$/="net";
my %hash;

open In, "$ARGV[0]" or die "$!";
while (<In>) {
chomp;
next if ($_=~/^$/);
my @temp = split (/\n/,$_);
my @unit = split (/\s/,$temp[0]);
my $target = $unit[1];
my $i;
for ($i=1;$i<=$#temp;$i++) {
    my @tmp = split (/\s/,$temp[$i]);
    if ($tmp[2]>10000 and $tmp[6]>10000 and $tmp[12]>5000) {
    my $target_end = $tmp[1]+$tmp[2];
    my $key1 = join("_",($target,$tmp[1],$target_end));
    my $query_end = $tmp[5]+$tmp[6];
    my $value1 =join("_",($tmp[3],$tmp[5],$query_end));
    $hash{$key1}=$value1;
    print "$key1\t$value1\n";
    } 
  }
}



open SV, "$ARGV[1]" or die "$!";
open SV_OUT, ">$ARGV[1].synteny.out" or die "$!";

$/="\n";
while (<SV>) {
my @temp = split (/\t/,$_);
my @unit1 = split (/\./,$temp[0]);
my @unit2 = split (/\./,$temp[4]);
foreach my $key (keys %hash) {
my @keys = split (/_/,$key);
if (($unit1[1] eq $keys[0]) and $keys[1]<=$temp[1] and $temp[1]<=$keys[2]) {
    my @queryvalue = split (/_/,$hash{$key});
    if (($unit2[1] eq $queryvalue[0]) and $queryvalue[1]<=$temp[5] and $temp[5]<=$queryvalue[2]) {
    print SV_OUT "$_";
    }
  }
 }
}

