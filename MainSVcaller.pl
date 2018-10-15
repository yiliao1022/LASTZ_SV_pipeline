use strict;
use warnings;
use SynNet;


my $root = SynNet($ARGV[0]);

print map "$_\n", @{$root->draw_ascii_tree};

my @daughters = $root->daughters; # @Chromosome 
     
	 foreach my $daughter (@daughters) {
             my @dayghters_playerA=$daughter->daughters; # reading first level fill       
		foreach my $playerA_dau (@dayghters_playerA) {        
				my @daughters_playerB_dau=$playerA_dau->daughters;
		
	if (@daughters_playerB_dau) {
                                } else {
				my $name=$playerA_dau->name;
				print "$name\n";
				}

            foreach my $playerB_dau (@daughters_playerB_dau) {	

			    my @daughters_playerC_dau=$playerB_dau->daughters;
			       
            if (@daughters_playerC_dau){	} else {
				my $name=$playerB_dau->name;
				print "$name\n";
				}				
			   
              foreach my $playerC_dau (@daughters_playerC_dau) {	
					    my $name=$playerC_dau->name; 
                        print "$name\n";					
				}						
            }
		}
    }
