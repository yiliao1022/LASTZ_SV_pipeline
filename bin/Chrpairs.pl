#!/usr/bin/perl
use strict;
use warnings;

my %hash= ('CM008268.1'=>'X',
          'CM008269.1'=>'2L',
          'CM008270.1'=>'2R',
          'CM008271.1'=>'3L',
          'CM008272.1'=>'3R',
          'CM008273.1'=>'4',
          );

open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].out" or die "$!";

while (<In>) {
chomp;
my @temp = split(/\t/,$_);
if ($temp[7] eq $hash{$temp[0]}) {
print Out "$_\n";
}
}
