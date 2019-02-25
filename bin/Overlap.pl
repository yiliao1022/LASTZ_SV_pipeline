#!/usr/bin/perl -w
open (IN,"$ARGV[0]");
my %hash;
while (<IN>){
my ($k,$s,$e) = split /\s/;
push @{$hash{$k}},[$s,$e];
}

open (INP,"$ARGV[1]");
open (OverLAP,">$ARGV[0]VS$ARGV[1].overlap.out");
while (<INP>){
chomp;
  my ($key,$start,$end) = split /\t/;
  foreach my $r (@{$hash{$key}}){
     if ($start >= $r->[0] && $start <= $r->[1] || $end >= $r->[0] && $end <= $r->[1] || $start <= $r->[0] && $end >= $r->[1])
     {print OverLAP "$_\n";}
    }
}
