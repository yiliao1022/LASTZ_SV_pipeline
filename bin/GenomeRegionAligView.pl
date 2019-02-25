#!/usr/bin/perl

=head1 Name

  GenomeRegionAligView.pl -- Drawing syntenic map for local regions that exhibit a potential SV (strutral vaiation) based on SVcaller.pl.

=head1 description

  Inputs: 
  
  1. Output from SVcaller.pl;
  2. The genome sequences of both target and query genomes;
  3. Bed files for both target and query genome sequences which can be transformed form Genome annoation GFF3 files.  
  
  
=Version
  Author: Yi Liao <yiliao1022@gmail.com>
  Version: 1.0  Date: 2019.01.05

=head1 Usage

  perl GenomeRegionAligView.pl -sv 1.txt -target /home/yliao/yliao/Ref/Japonica/O_sativa_japonica_sm21.fasta -query /home/yliao/yliao/Ref/kasalath/kasalath.fasta -tbed /home/yliao/yliao/Ref/Japonica/all.bed -qbed /home/yliao/yliao/2018_tRNA/ref/annotation/Kasalath.IGDBv2.bed -len 2000
  -sv      1.txt
  -target  target genome sequence
  -query   query genome sequence
  -len     sequence length for extending the flanking region 
  -tbed    target bed file
  -qbed    query bed file
  -verbose  verbose
  -help     help

=cut



use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use SVG;


# Global Variablelist
my ($sv,$target,$query,$Tbed,$Qbed,$length,$Verbose,$Help);

##get options from command line
GetOptions(
    "sv:s" => \$sv,
	"target:s" => \$target,
	"query:s"=> \$query,
	"tbed:s" => \$Tbed,
	"qbed:s" => \$Qbed,
	"len:i" => \$length,
	"verbose"  => \$Verbose,
	"help"     => \$Help
);


if ($Help){
print <<"END.";
  Usage: perl $0 -sv 1.txt -target /home/yliao/yliao/Ref/Japonica/O_sativa_japonica_sm21.fasta -query /home/yliao/yliao/Ref/kasalath/kasalath.fasta -tbed /home/yliao/yliao/Ref/Japonica/all.bed -qbed /home/yliao/yliao/2018_tRNA/ref/annotation/Kasalath.IGDBv2.bed -len 2000
  -sv      1.txt
  -target  target genome sequence
  -query   query genome sequence
  -len     sequence length for extending the flanking region 
  -tbed    target bed file
  -qbed    query bed file
  -verbose  verbose
  -help     help
END.
exit;
}
###################################################################################################

open IN, "$sv" or die "$!";

while (<IN>) {
chomp;
SVGlocal ($_,$length,$target,$query,$Tbed,$Qbed);
#system "rm *.gz *.fa *.gff *.maf";
}

close IN;


sub SVGlocal {
my ($line,$extend,$target,$query,$TBED,$QBED) = @_;

my @temp = split (/\t/,$line);

my $T_beg=$temp[1]-$extend;
my $Q_beg=$temp[8]-$extend;

my $T_end=$temp[2]+$extend;
my $Q_end=$temp[9]+$extend;

my $T_len = $T_end-$T_beg;
my $Q_len = $Q_end-$Q_beg;

`iTools Fatools extract -InPut $target -OutPut Target$temp[0]\_$temp[1].fa -ORegion $temp[0]:$T_beg:$T_end`;
`iTools Fatools extract -InPut $query -OutPut Query$temp[7]\_$temp[8].fa -ORegion $temp[7]:$Q_beg:$Q_end`;
`gunzip *.gz`;
`~/prog/lastz/src/lastz Target$temp[0]\_$temp[1].fa Query$temp[7]\_$temp[8].fa --format=maf --out=Target$temp[0]\_$temp[1].Query$temp[7]\_$temp[8].maf --identity=90 --nogapped`;

#####################################################
my $svg= SVG -> new (width=>600, height =>300);
my $rate;
my $rectange_height=4;
my $chr_color="black";
my $beg=100;
my $end=500;

if ($T_len > $Q_len) {
$rate = $T_len/400;
} else {
$rate = $Q_len/400;
}

######Draw skeleton
$svg->line (x1=> $beg, y1 => 80, x2=> $end, y2=>80, style => { 'fill'=> $chr_color, 'stroke'=> $chr_color,'stroke-width'=>'1.5',});			
		
my $int = ($end-$beg)/10;
my $i;	
for ($i=0;$i<11;$i++) {
	 $svg->line ( x1=> $beg+$int*$i, y1 => 80, x2=> $beg+$int*$i, y2 => 74, style => { 'fill'=> $chr_color, 'stroke'=> $chr_color,'stroke-width'=>'1.5',} );
    }

$svg->text("x", 98,"y", 70,"-cdata","$T_beg","font-family","arial","font-size",10,"fill","black" );
$svg->text("x", 492,"y", 70,"-cdata","$T_end","font-family","arial","font-size",10,"fill","black");

$svg->rectangle ( x=>100, y=>100, width=>400, height=>$rectange_height,style => { 'fill'=> $chr_color, 'stroke'=> $chr_color, });			
$svg->rectangle ( x=>100, y=>200,width=>$Q_len/$rate, height=>$rectange_height, style => { 'fill'=> $chr_color, 'stroke'=> $chr_color,});				


######Draw synteny relationship from alignment Maf files
my $aligns = &ReadMaf ("./Target$temp[0]\_$temp[1].Query$temp[7]\_$temp[8].maf");

foreach my $block (@$aligns) {
  my @temp=split("_",$block);
  my $color="gray";
  my $tleft=$temp[0]/$rate+100;
  my $tright=$temp[1]/$rate+100;
  my $qright=$temp[2]/$rate+100;
  my $qleft=$temp[3]/$rate+100;
  my $theight=114;
  my $qheight=190;
  my $xv=[$tleft,$tright,$qright,$qleft];
  my $yv=[$theight,$theight,$qheight,$qheight];
  my $points=$svg->get_path( x=>$xv,y=>$yv,-type=>'polyline',-closed=>'true');
  my $tag=$svg->polyline( %$points, style=>{fill=>$color});  
}

######Draw Gene features from gff files
$temp[7]=~ s/Osjap//g;
$temp[7]=~ s/Chr0/Chr/g;
`echo -e "$temp[0]\t$T_beg\t$T_end" | bedextract $TBED - > T$temp[0]\_$T_beg\_$T_end.bed`;
`echo -e "$temp[7]\t$Q_beg\t$Q_end" | bedextract $QBED - > Q$temp[7]\_$Q_beg\_$Q_end.bed`;



open TBED, "T$temp[0]\_$T_beg\_$T_end.bed" or die "$!";
$/="\n";

while (<TBED>) {
chomp;
my @temp = split (/\t/,$_);
if ($temp[7]=~/CDS/ and $temp[3]=~/\.1/ and $temp[5]=~/\-/) {
$svg->rectangle (x=>($temp[1]-$T_beg)/$rate+100, y=>106,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif ($temp[7]=~/CDS/ and $temp[3]=~/\.1/ and $temp[5]=~/\+/) {
$svg->rectangle (x=>($temp[1]-$T_beg)/$rate+100, y=>92,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif( ($temp[7]=~/three/|$temp[7]=~/five/ )and $temp[3]=~/\.1/ and $temp[5]=~/\+/){
$svg->rectangle (x=>($temp[1]-$T_beg)/$rate+100, y=>92,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif(($temp[7]=~/three/| $temp[7]=~/five/)and $temp[3]=~/\.1/ and $temp[5]=~/\-/){
$svg->rectangle (x=>($temp[1]-$T_beg)/$rate+100, y=>106,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
}elsif ($temp[7]=~/gene/ and $temp[5]=~/\-/) {
$svg->line (x1=>($temp[1]-$T_beg)/$rate+100,y1=>109,x2=>($temp[2]-$T_beg)/$rate+100,y2=>109, style=>{'fill'=>"green",'stroke'=>"green",'stroke-width'=>'1',});
} elsif ($temp[7]=~/gene/ and $temp[5]=~/\+/) {
$svg->line (x1=>($temp[1]-$T_beg)/$rate+100,y1=>95,x2=>($temp[2]-$T_beg)/$rate+100,y2=>95, style=>{'fill'=>"green",'stroke'=>"green",'stroke-width'=>'1',});
}
}
close TBED;


open QBED, "Q$temp[7]\_$Q_beg\_$Q_end.bed" or die "$!";

while (<QBED>) {
chomp;
my @temp = split (/\t/,$_);
if ($temp[7]=~/CDS/ and $temp[9]=~/(\.T01|t1)/ and $temp[5]=~/\-/) {
$svg->rectangle (x=>($temp[1]-$Q_beg)/$rate+100, y=>206,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif ($temp[7]=~/CDS/ and $temp[9]=~/(\.T01|t1)/ and $temp[5]=~/\+/) {
$svg->rectangle (x=>($temp[1]-$Q_beg)/$rate+100, y=>192,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif( ($temp[7]=~/three/|$temp[7]=~/five/ )and $temp[9]=~/(\.T01|t1)/ and $temp[5]=~/\+/){
$svg->rectangle (x=>($temp[1]-$Q_beg)/$rate+100, y=>192,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
} elsif(($temp[7]=~/three/| $temp[7]=~/five/)and $temp[9]=~/(\.T01|t1)/ and $temp[5]=~/\-/){
$svg->rectangle (x=>($temp[1]-$Q_beg)/$rate+100, y=>206,width=>($temp[2]-$temp[1])/$rate,height=>6,style => {'fill'=>"green",'stroke'=>"green",});
}elsif ($temp[7]=~/gene/ and $temp[5]=~/\-/) {
$svg->line (x1=>($temp[1]-$Q_beg)/$rate+100,y1=>209,x2=>($temp[2]-$Q_beg)/$rate+100,y2=>209, style=>{'fill'=>"green",'stroke'=>"green",'stroke-width'=>'1',});
} elsif ($temp[7]=~/gene/ and $temp[5]=~/\+/) {
$svg->line (x1=>($temp[1]-$Q_beg)/$rate+100,y1=>195,x2=>($temp[2]-$Q_beg)/$rate+100,y2=>195, style=>{'fill'=>"green",'stroke'=>"green",'stroke-width'=>'1',});
}
}
close QBED;


######################################################

open OUT, ">Target$temp[0]\_$temp[1].Query$temp[7]\_$temp[8].svg" or die "can not open my file";
print OUT $svg->xmlify;
close OUT;
system "/home/yliao/prog/tool/svg2xxx_release/svg2xxx Target$temp[0]\_$temp[1].Query$temp[7]\_$temp[8].svg -t pdf";
}


####################################################
################### Sub Routines ###################
####################################################

sub ReadMaf {
$/="a score";
my $maf_file = shift;
my @align;

open MAF, "$maf_file" || die "Fail open $maf_file";

<MAF>;

while (<MAF>){

if ($_=~/\#/) {
        next;
}else {
          my @temp = split (/\n/,$_);
          my @unit1 = split (/\s+/,$temp[1]);
          my @unit2 = split (/\s+/,$temp[2]);
   
          my $tleft=$unit1[2];
          my $tright=$unit1[2]+$unit1[3];
		  my ($qleft,$qright,$orientation);
		  
      if ($unit2[4]=~/\+/) {			                                  
           $qleft=$unit2[2];
           $qright=$unit2[2]+$unit2[3];
           $orientation="plus";
         }elsif($unit2[4]=~/\-/){
           $qleft=$unit2[5]-$unit2[2];
           $qright=$unit2[5]-$unit2[2]-$unit2[3];
	       $orientation="minus";	 
		 }
  my $block = join ("_",($tleft,$tright,$qright,$qleft,$orientation));		 
  push (@align,$block);
  }
}
close MAF;
my $ref = \@align;
return $ref;
}


sub ReadGff {
$/="\n";
my $gff_file = shift;
open GFF, "$gff_file" || die "Fail open $gff_file";
my %gff;
my $key;

while (<GFF>) {
chomp;
if ($_=~/\#/) {
  next;
} else {

my @temp = split (/\t/,$_);

if ($temp[2] eq "gene") {
$temp[8]=~/ID=(.*);Name=/;
$key = $1;
my @features=();
my $ref = \@features;
$gff{$1}=$ref;
}

if ($temp[8]=~/$key/) {
   if ($temp[2]=~/five/) {
   my $five = join ("_",($temp[2],$temp[3],$temp[4]));
   push (@{$gff{$key}},$five);
   }
   if ($temp[2]=~/three/) {
   my $three = join ("_",($temp[2],$temp[3],$temp[4]));
   push (@{$gff{$key}},$three);
   }
   if ($temp[2]=~/CDS/) {
   my $cds = join ("_",($temp[2],$temp[3],$temp[4]));
   push (@{$gff{$key}},$cds);
   }
  }
 }
}
close GFF;
my $reference=\%gff;
return $reference;
}
