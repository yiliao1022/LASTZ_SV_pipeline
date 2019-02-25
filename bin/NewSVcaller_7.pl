#!/usr/bin/perl

=head1 Name

  SVcaller.pl -- SVs (strutral vaiations) caller based on LASTZ/Chain/Net whole genome alignment.

=head1 description

  Inputs: 1. whole genome pairwise alignments of NET format
  
=Version
  Author: Yi Liao <yiliao1022@gmail.com>
  Version: 1.0  Date: 2018.10.15

=head1 Usage

  perl SVcaller.pl [options] <SynNet file>
  --outdir <str>  Set the output directory
  --verbose       output running progress information to screen
  --help          output help information

=cut


use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Tree::DAG_Node;

# Global Variablelist
my ($Outdir,$ref,$SynNetFile,$ali_rate,$Verbose,$Help);
my ($OutIndel, $OutInversion, $OutComplex,$OutCnv);


##get options from command line
GetOptions(
	"outdir:s" => \$Outdir,
	"ref:s"=> \$ref,
	"SynNet:s" => \$SynNetFile,
	"ali_rate:i" => \$ali_rate,
	"verbose"  => \$Verbose,
	"help"     => \$Help
);


if ($Help){
print <<"END.";
  Usage: perl $0 [options] --SynNet <SynNet file> --ref <Reference sequence> --ali_rate 0.25
  -outdir   outdir
  -SynNet <SynNet file>
  -ref <Reference (target) sequence file>
  -ali_rate <default 0.25>
  -verbose  verbose
  -help     help
  
  SV output format:
  
   1. Target seqid
   2. Target start
   3. Target end
   4. SV type (internal is contained within an alignment block)
   5. SV length (insertions are in query length)
   6. Percent identity in alignment block or for flanking alignment blocks
   7. Number of bases matching in alignment block or for flanking alignment blocks
   8. Query seqid
   9. Query start
   10. Query end
   11. Sequence (for deletions it's the target sequence and insertions are the query sequence)

END.
exit;
}
## set parameters

$ali_rate ||= 0.25;

##set outdir

$Outdir ||= ".";
$Outdir =~ s/\/$//;
$OutIndel="Indel";
$OutInversion="Inversion";
$OutComplex="ComplexSV";
$OutCnv="CopyNumVar";

mkdir($Outdir) unless(-d $Outdir);
my $InDeldir = "$Outdir/$OutIndel";
my $Invdir =  "$Outdir/$OutInversion";
my $CxSVdir = "$Outdir/$OutComplex";
my $CNVdir="$Outdir/$OutCnv";
mkdir ($InDeldir) unless(-d $InDeldir);
mkdir ($Invdir) unless(-d $Invdir);
mkdir ($CxSVdir) unless(-d $CxSVdir);
mkdir ($CNVdir) unless (-d $CNVdir);

#set file names

my $IndelCore = basename ($SynNetFile);
my $InvCore = basename ($SynNetFile);
my $CxSVCore = basename ($SynNetFile);
my $CnvCore = basename($SynNetFile);


my $IndelFile = "$Outdir/$InDeldir/$IndelCore.indel";
my $InvFile="$Outdir/$Invdir/$InvCore.inv";
my $CxSVFile="$Outdir/$OutComplex/$CxSVCore.cxsv";
my $CNVFile="$Outdir/$OutCnv/$CnvCore.cnv";

###############################################################
my $root = SynNet ($SynNetFile);
my $seq = ExtractSeq ($ref);



open (OUTind,">$IndelFile") || die "Cannot open $IndelFile\n";
open (OUTinv,">$InvFile") || die "Cannot open $InvFile\n";
open (OUTcxs,">$CxSVFile") || die "Cannot open $CxSVFile\n";
open (OUTcnv,">$CNVFile") || die "Cannot open $CNVFile\n";

###############################################################
my $PAV=0; # Count the sequence content of Presence and absence variations

###############################################################


my @daughters = $root->daughters; # @Chromosome 
	 foreach my $daughter (@daughters) { # chromosome/contig name 
             my @dayghters_playerA=$daughter->daughters; #  top segments for filling.      
		
		foreach my $playerA_dau (@dayghters_playerA) { 
                
                my $mama=$playerA_dau->mother; 
				my $mamaname=$mama->name; # The name of the chromosome/contigs
                my $topfillname = $playerA_dau->name; # The name of the top fill line
                
                my @unit = split (/_/,$topfillname);                
                my $fill_beg = $unit[1];             
                my $fill_end = $unit[1]+$unit[2];
                               
                my @daughters_playerB_dau=$playerA_dau->daughters;	# Gaps in the alignment blocks for detecting SVs	
	            
				########################################################################
				if (@daughters_playerB_dau) { # Gaps between the alignment blocks in each top filled segments for chromosomes/contigs
                                             									           
			        LOOP:foreach my $i (0..$#daughters_playerB_dau) {    
			                    
			                    my $name0;
			                    my $name1;
			                    my $name=$daughters_playerB_dau[$i]->name; # the gap line
			                    			                    			               		                    			                    
			             if ($#daughters_playerB_dau >0) {
			                    
			                    if ($i==0) {
			                    $name0 ="gap_$fill_beg\_0";
			                    $name1 =$daughters_playerB_dau[1]->name; 
							    } elsif ($i==$#daughters_playerB_dau) {
							    $name0 =$daughters_playerB_dau[$#daughters_playerB_dau-1]->name;
			                    $name1 ="gap_$fill_end\_0";
							    } else {
							    $name0 =$daughters_playerB_dau[$i-1]->name;
							    $name1 =$daughters_playerB_dau[$i+1]->name;
							    }							    
						} else {
							    $name0 ="gap_$fill_beg\_0";
							    $name1 ="gap_$fill_end\_0";
						}
							    					    
							    my @tmp0 = split (/_/,$name0);
							    my @tmp = split (/_/,$name);
							    my @tmp1 = split(/_/, $name1);
							    
							    my $Tend = $tmp[1]+$tmp[2];
							    my $Qend = $tmp[5]+$tmp[6];
							    my $forward = $tmp[1] - $tmp0[1] - $tmp0[2];
							    my $backward = $tmp1[1] - $tmp[1] - $tmp[2];
						        my $Sequence = substr (${$seq}{$mamaname},$tmp[1],$tmp[2]);
							    
							    my @daughters_playerC_dau=$daughters_playerB_dau[$i]->daughters; # The segments for filling the gaps.
                                    
									if (@daughters_playerC_dau){ # if gaps were filled by segments that from other genomic regions or inverted itself (e.g. inversion).
							        #######################Inv, Complex SVs############################
																				
									##########################################		
								    if ($tmp[6] < 10) {														
											          my $ali=0;											        
											          my @coordinates=();
											          my $inv_len=0;	    			
											  foreach my $playerC_dau (@daughters_playerC_dau) {	# For each filled segments
					                                  my $name=$playerC_dau->name;  # Filled segments name
											       	  my $hash = ParseFill ($name);
                                                    
												      if ($$hash{"type"}=~/syn|inv/ and $$hash{"qFar"} < 2*$tmp[2]) { # Define copy number variations
                                                            $inv_len = $inv_len + $$hash{"ali"} if ($$hash{"type"} =~/inv/);
                                                            $ali= $ali + $$hash{"ali"};
                                                            my $end = $$hash{"query_beg"} + $$hash{"query_len"};
                                                            push (@coordinates,$$hash{"query_beg"});
                                                            push (@coordinates,$end);
													     }  												     													     																								      
													}
													
													my @array = sort {$a <=> $b} @coordinates;	      
										            if (($inv_len > ($tmp[6]+$tmp[2])/4 or $inv_len>4000)) {
												    print OUTinv "$mamaname\t$tmp[1]\t$Tend\tINV\t$tmp[2]\t0.95\t$forward\:$backward\t$tmp[3]\t$tmp[5]\t$Qend\t$Sequence\n";
												    next LOOP;
										            }elsif ($ali/$tmp[2] > $ali_rate) {
													print OUTcnv "$mamaname\t$tmp[1]\t$Tend\tCNV\t$tmp[2]\t0.95\t$forward\:$backward\t$tmp[3]\t$array[0]\t$array[-1]\t$Sequence\n";
							                        } else {
							                        print OUTind "$mamaname\t$tmp[1]\t$Tend\tDEL\t$tmp[2]\t0.95\t$forward\:$backward\t$tmp[3]\t$tmp[5]\t$Qend\t$Sequence\n";
							                        }
							                        							                        
							       } else {  						       		 					      							       							       
													  my $ali=0;
											          my @coordinates=();	
											          my $inv_len=0;	
											          	
											  foreach my $playerC_dau (@daughters_playerC_dau) { # For each filled segments
					                                  my $name=$playerC_dau->name;  # Filled segments name
											       	  my $hash = ParseFill ($name);									              
										              
										                if ($$hash{"type"} =~/inv|syn/ and $$hash{"qFar"} < 2*$tmp[2] ){ # Define inversions														        
														        $inv_len = $inv_len + $$hash{"ali"} if ($$hash{"type"} =~/inv/);														        
														        $ali= $ali + $$hash{"ali"};
                                                                my $end = $$hash{"query_beg"} + $$hash{"query_len"};
                                                                push (@coordinates,$$hash{"query_beg"});
                                                                push (@coordinates,$end);														  														 														  
														  } 
													}
													
												my @array = sort {$a <=> $b} @coordinates;
																								
												if (($inv_len > ($tmp[6]+$tmp[2])/4 or $inv_len>4000)) {
												print OUTinv "$mamaname\t$tmp[1]\t$Tend\tINV\t$tmp[2]\t0.95\t$forward\:$backward\t$tmp[3]\t$tmp[5]\t$Qend\t$Sequence\n";
												next LOOP;
												}elsif ($ali/$tmp[2] > $ali_rate) {
												print OUTcnv "$mamaname\t$tmp[1]\t$Tend\tCNV\t$tmp[2]\t0.95\t$forward\:$backward\t$tmp[3]\t$array[0]\t$array[-1]\t$Sequence\n";
												}else {
												print OUTcxs "$mamaname\t$tmp[1]\t$Tend\tCMX\t$tmp[2]\t0.95\t$forward\:$backward\t$tmp[3]\t$tmp[5]\t$Qend\t$Sequence\n";
                                                }																										
                                      }
                                                                                     
                                     ##########################################
                                            											###################################################################  
				                             						  
				                    } else { # if gaps were not filled by any segments that from other genomic regions, e.g. Completely absence of homologous sequences in the query genome
				                                
											    if ($tmp[6]<10) {
 				                                  print OUTind "$mamaname\t$tmp[1]\t$Tend\tDEL\t$tmp[2]\t0.95\t$forward\:$backward\t$tmp[3]\t$tmp[5]\t$Qend\t$Sequence\n"; # Here only record the Inserion to reference
												  $PAV = $PAV + $tmp[2];
											    } else {
											      print OUTcxs "$mamaname\t$tmp[1]\t$Tend\tCMX\t$tmp[2]\t0.95\t$forward\:$backward\t$tmp[3]\t$tmp[5]\t$Qend\t$Sequence\n"; # Here indicate potential convergent signals in evolution
												     if ($tmp[2]/$tmp[6] > 4) {
												     $PAV = $PAV + $tmp[2]-$tmp[6];
													 }
											    }
				                            }							                  						
                        }        								
				} else {				      
				      my $name=$playerA_dau->name;
					  print "TOP: Fill $mamaname\n";
				      print "TOP: Fill $name\n";
			         }	
				#######################################################################	 
					 
			}           
     }


close OUTind;
close OUTinv;
close OUTcxs;
close OUTcnv;

print "The total presence/absence variations account for $PAV bp\n";


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

sub SynNet {
$/="net";
my $net_file = shift;
open In, "$net_file" || die "Fail open.$net_file";

my $root = Tree::DAG_Node->new; # Root
   $root->name('Osjap');
<In>;   
while (<In>) { 
       if ($_=~/^#/) {
           next;
       } else {
               my @temp = split (/\n/,$_);
               my $head = shift @temp;
               my @unit = split (/\s/,$head);
               my $a1 = Tree::DAG_Node->new;
               $root ->add_daughter($a1); # Chromosome level.
               $a1->name($unit[1]); # Chromosome name.
               my $i;
		       my @fill_nodeA; 
			   my @gap_nodeA;
			   my @fill_nodeB;
			   my @gap_nodeB;
                  for ($i=0;$i<=$#temp;$i++) {
                           if ($temp[$i] =~/^\s{1}fill/) {
                           my @fill=split(/\s/,$temp[$i]);
						   shift @fill;
                           my $name=join("_",@fill); # Chromosome_Name begin end.
                           my $a2 = Tree::DAG_Node->new; # Fill playerA: Generally, each chromsome may have one or more this level fill, but only one cover most of the chromosome region.
                              $a1->add_daughter($a2);
                              $a2->name($name); # Chromosome_Name begin end.
				              push(@fill_nodeA,$a2); # record the $a2 which is the ancentor of the .ap?for the main SVs 
                            } elsif ($temp[$i] =~/^\s{2}gap/) {
                              my @gap = split(/\s+/,$temp[$i]);
				              shift @gap; # Remove the first element of @gap, It should be a blank space.
                              my $name=join("_",@gap); # Most simple insertions and deletions will be identified in this player. 
                              my $a3 = Tree::DAG_Node->new;
                              $fill_nodeA[$#fill_nodeA]->add_daughter($a3);
                              $a3->name($name); # These names contain the major information of simple insertions and deletions.
							  push (@gap_nodeA,$a3); # record the $a3 (gap) which is the ancentor of complex SVs, including inversions, local/tandem duplications. 
                            } elsif ($temp[$i] =~/^\s{3}fill/) {
							  my @fill = split(/\s+/,$temp[$i]);
				              shift @fill; # Remove the first element of @gap, it should be a blank space.
                              my $name=join("_",@fill); # Most complex SVs will be identified in this player. 
							  my $a4 = Tree::DAG_Node->new;
                              $gap_nodeA[$#gap_nodeA]->add_daughter($a4);
                              $a4->name($name); # These names contain the major information of simple insertions and deletions.
							  push (@fill_nodeB,$a4); # record the $a4 (fill) which is the ancentor of complex SVs, including inversions, local/tandem duplications. 						
							} elsif ($temp[$i] =~/^\s{4}gap/) {
							  my @gap = split(/\s+/,$temp[$i]);
				              shift @gap; # Remove the first element of @gap, it should be a blank space.
                              my $name=join("_",@gap); # SVs in duplication regions will be identified in this player. 
							  my $a5 = Tree::DAG_Node->new;
                              $fill_nodeB[$#fill_nodeB]->add_daughter($a5);
                              $a5->name($name); # These names contain the major information of simple insertions and deletions.
							  push (@gap_nodeB,$a5); # record the $a5 (gap) which can be used to infer the mosaic patter of segmental duplications.
							} elsif ($temp[$i] =~/^\s{5}fill/) {
							  my @fill = split(/\s+/,$temp[$i]);
				              shift @fill; # Remove the first element of @gap, it should be a blank space.
                              my $name=join("_",@fill); # SVs in duplication regions will be identified in this player. 
							  my $a6 = Tree::DAG_Node->new;
                              $gap_nodeB[$#gap_nodeB]->add_daughter($a6);
                              $a6->name($name); # These names contain the major information of simple insertions and deletions.
							  #push (@fill_nodeC,$a6); # record the $a6 (gap) which can be used to infer the mosaic patter of segmental duplications.
							}
                    }
            }
}
return $root;
}
