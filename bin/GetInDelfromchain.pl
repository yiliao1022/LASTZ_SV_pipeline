#!/usr/bin/perl
use strict;
use warnings;




my $T_seq = Getseq ($ARGV[1]); ### Target genome sequence file;
my $Q_seq = Getseq ($ARGV[2]); ### Query genome sequence file;


$/="\n";
open IN, "$ARGV[0]" or die "Cannot Open $ARGV[0] !";
open OUT, ">$ARGV[3].Indel.out" or die "Can not open $ARGV[1].$ARGV[2].Indel.out";

my $head = <IN>;

my @temp =split (/\s+/, $head);

my $target_position= $temp[5];
my $query_position= $temp[10];
my $target_chr=$temp[2];
my $query_chr=$temp[7];

print OUT "Target_chr\tIndel/Type\tTarget_start\tTarget_end\tQuery_chr\tQuery_start\tQuery_end\tOrientation\tlength\n";

while (<IN>) {
  chomp;
  if ($_=~/chain/) {
           last;
     } else {
         my @unit=split (/\t/, $_);
         if ($unit[1]==0 and $unit[2] > 0 ) {
         
         $target_position=$target_position+$unit[0];
         my $query_start=$query_position+$unit[0];
         $query_position=$query_start+$unit[2];
         print OUT "$target_chr\tDel\t$target_position\t$target_position\t$query_chr\t$query_start\t$query_position\t\+\t$unit[2]\n"
         }
       
         elsif ($unit[1]>0 and $unit[2]==0 ) {
          my $target_start=$target_position+$unit[0];
          $query_position=$query_position+$unit[0];
          $target_position=$target_start+$unit[1];
          
         print OUT "$target_chr\tIst\t$target_start\t$target_position\t$query_chr\t$query_position\t$query_position\t\+\t$unit[1]\n";
         }


         elsif ($unit[1]>0 and $unit[2]>0) {

           my $target_start=$target_position+$unit[0];
           my $query_start=$query_position+$unit[0];
           $query_position=$query_start+$unit[2];
           $target_position=$target_start+$unit[1];
                      
           my $targe_sequence=substr (${$T_seq}{$target_chr}, $target_start, $unit[1]);
           my $query_sequence=substr (${$Q_seq}{$query_chr}, $query_start, $unit[2]);

         print OUT "$target_chr\tUmb\t$target_start\t$target_position\t$query_chr\t$query_start\t$query_position\t\+\t$unit[1]\t$unit[2]\t$targe_sequence\t$query_sequence\n";  

         }
}       

}



sub Getseq {
   
   my ($file) = @_;
   my %hash;
   $/=">";
   open In, "$file" or die "$!";
      while (<In>) {
         next if (length $_ < 2);
         my @unit= split ("\n", $_);
         my $head=shift @unit;
         my $seq = join ("", @unit);
         $hash{$head}=$seq;      
      }
     my $ref=\%hash;
     return $ref;
 }
