use strict;
use warnings;
use Tree::DAG_Node;



sub SynNet {
$/="net";
my $net_file = shift;
open In, "$net_file" || die "Fail open:$net_file";
print "reading $net_file";
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
                           my $name=join("_",($fill[4],$fill[2],$fill[3])); # Chromosome_Name begin end.
                           my $a2 = Tree::DAG_Node->new; # Fill playerA: Generally, each chromsome may have one or more this level fill, but only one cover most of the chromosome region.
                              $a1->add_daughter($a2);
                              $a2->name($name); # Chromosome_Name begin end.
				              push(@fill_nodeA,$a2); # record the $a2 which is the ancentor of the “gap” for the main SVs 
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
1;
