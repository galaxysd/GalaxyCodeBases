#!/usr/bin/perl
###########################################################################
##  ViralFusionSeq
##  Software for discovering and annotating viral integration event and
##    fusion transcript
##  
##  Version 1.0 -- November 15, 2012
##  
##  Copyright (C) 2012 by Jing-Woei Li & Raymond Wan, All rights reserved.
##  Contact:  marcoli@cuhk.edu.hk, rwan@cuhk.edu.hk
##  Organization:  Hong Kong Bioinformatics Centre, School of Life Sciences, The
##                 Chinese University of Hong Kong, Shatin, NT,
##                 Hong Kong SAR
##  
##  This file is part of ViralFusionSeq.
##  
##  ViralFusionSeq is free software; you can redistribute it and/or 
##  modify it under the terms of the GNU General Public License 
##  as published by the Free Software Foundation; either version 
##  3 of the License, or (at your option) any later version.
##  
##  ViralFusionSeq is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public 
##  License along with ViralFusionSeq; if not, see 
##  <http://www.gnu.org/licenses/>.
###########################################################################


##  $LastChangedDate: 2012-11-27 17:32:45 +0800 (Tue, 27 Nov 2012) $
##  $LastChangedRevision: 1092 $

use List::Util qw( min max );
use strict;
use FileHandle;
$| = 1; #disable buffer

my $usage = qq{Usage: $0 <infected chromsome> <output file of SC method> <output file of RP method>
	           Make sure you have indexed your infected chromosome
};
die($usage) if (@ARGV == 0);
print "Executing $0\n";
print "Infected chromosome $ARGV[0]\n";
print "Output file of SC method: $ARGV[1]\n";
print "Output file of RP method: $ARGV[2]\n";

#####################
## Path to programs #
#####################
my $BLAST_bin = "blastall";
my $formatdb = "formatdb";
my $thread = 8;
my $evalue = 5;

#########################
## Capture user's input #
#########################

my $Infected_chr = $ARGV[0];
my $SC_output_file = $ARGV[1];
my $RP_output_file = $ARGV[2];

my @SC_a = @{file_2_array($SC_output_file)};
my @RP_a = @{file_2_array($RP_output_file)};

my @files_for_serch = ();
#format BLASTdb
my $cmd = "$formatdb -p F -i $Infected_chr";
system ($cmd);
	
#convert CSm output to fasta for location search
open my $CSm_wFH, '>', "$SC_output_file.LOCs.fa";
for (@SC_a){
	my @s = split (/\t/, $_);
	my ($v_read, $h_read, $full_read) = ($s[12], $s[13], $s[14]);
	my @readnameCORE_split = split (/_/, $s[0]);
	my $readnameCORE = $readnameCORE_split[0];
	my $v_read_wh_str = ">".$readnameCORE."_v"."\n".$v_read."\n";
	my $h_read_wh_str = ">".$readnameCORE."_h"."\n".$h_read."\n";
	my $full_read_wh_str = ">".$readnameCORE."_f"."\n".$full_read."\n";
	print $CSm_wFH $v_read_wh_str;
	print $CSm_wFH $h_read_wh_str;
	print $CSm_wFH $full_read_wh_str;
}
close($CSm_wFH);
push (@files_for_serch, "$SC_output_file.LOCs.fa");

#convert RPm output to readID, original position of reads on infected chromosome is already stored in header of simulated reads (no need to search)
open my $RPm_wFH, '>', "$RP_output_file.LOCs.fa";
open my $RPm_readIDonly_wFH , '>', "$RP_output_file.readID";
for (@RP_a){
	my @s = split (/\t/, $_);
	my $v_readname_raw = $s[3];
	my $v_seq = $s[-21];
	my $v_len = length($v_seq);
	my @v_readname_split = split (/\//, $v_readname_raw);
	my $h_readname_raw = $s[-17];
	my $h_seq = $s[-1];
	my $h_len = length($h_seq);
	#my @h_readname_split = split (/\//, $h_readname_raw);
	print $RPm_wFH ">".$v_readname_raw."_$v_len"."\n".$v_seq."\n".">".$h_readname_raw."_$v_len"."\n".$h_seq."\n";
	print $RPm_readIDonly_wFH $v_readname_split[0]."\n";
}
close ($RPm_wFH);
close($RPm_readIDonly_wFH);
push (@files_for_serch, "$RP_output_file.LOCs.fa");

## Do the search and parse for Evlauation script
for my $curr_file_4_search (@files_for_serch){
	my $cmd = "$BLAST_bin -p blastn -d $Infected_chr -m 9 -i $curr_file_4_search -a $thread -e $evalue >$curr_file_4_search.parsed";
	print "$cmd\n";
	system ($cmd);
	my @result_a = @{file_2_array("$curr_file_4_search.parsed")};
	if ($curr_file_4_search =~ /CSm./){
		print "Working $curr_file_4_search.parsed\n";
		open my $CSparsed_FULLoutFH, '>',"$curr_file_4_search.parsed.full.out";
		open my $CSparsed_EVALoutFH, '>',"$curr_file_4_search.parsed.eval.out";
		open my $CSparsed_WGSIMCKoutFH, '>',"$curr_file_4_search.parsed.wgsim_check.out";
		my %m = ();
		my @readname = ();
		my %readname_full_to_evalues = ();
		my %readcore_to_fullreadname = ();
		for (@result_a){
			my @s = split (/\t/, $_);
			my $readIDraw = $s[0];
			my @sRAWread = split (/_/, $readIDraw);
			my $readname_core = $sRAWread[0];
			push (@{$m{$readIDraw}{$s[-2]}}, $_); #readname evalue full_line
			push (@{$readname_full_to_evalues{$readIDraw}}, $s[-2]);
			push (@readname, $readname_core); #order of reads
			push (@{$readcore_to_fullreadname{$readname_core}}, $readIDraw);
		}
		#remove redundant record in readname array, in order
		my %redundnat = ();
		my @rename_non_redundnat = ();
		for (@readname){
			unless (defined $redundnat{$_}){
				push (@rename_non_redundnat, $_);
				$redundnat{$_} = "";
			}
		}
		for my $readname_core(@rename_non_redundnat){ #Parse in read pair
			my @readname_full_two = @{$readcore_to_fullreadname{$readname_core}};
			my ($viral_segment, $human_segment, $full_read);
			for (@readname_full_two){
				if (/_v/){
					$viral_segment = $_
				}
				if (/_h/){
					$human_segment = $_;
				}
				if (/_f/){
					$full_read = $_;
				}
			}
			my $At_least_one_best_match = 0;
			my @VreadID_w_lowest_evalue = ();
			my @HreadID_w_lowest_evalue = ();
			my @FreadID_w_lowest_evalue = ();
			if (($readname_full_to_evalues{$viral_segment})&&($readname_full_to_evalues{$human_segment})&&($readname_full_to_evalues{$full_read})){
				##viral segment
				my $Vevalue_aref = $readname_full_to_evalues{$viral_segment};
			    my @Vevalues = @{$Vevalue_aref};
				my $Vmin_evalue = min(@Vevalues);
				@VreadID_w_lowest_evalue = @{$m{$viral_segment}{$Vmin_evalue}};
				my $VreadID_w_lowest_evalue_count = @VreadID_w_lowest_evalue;
				if ($VreadID_w_lowest_evalue_count>1){
					warn "$viral_segment has more than 1 best match\n";
				}
				##human segment
				my $Hevalue_aref = $readname_full_to_evalues{$human_segment};
			    my @Hevalues = @{$Hevalue_aref};
				my $Hmin_evalue = min(@Hevalues);
				@HreadID_w_lowest_evalue = @{$m{$human_segment}{$Hmin_evalue}};
				my $HreadID_w_lowest_evalue_count = @HreadID_w_lowest_evalue;
				if ($HreadID_w_lowest_evalue_count>1){
					warn "$human_segment has more than 1 best match\n";
				}
				##Full read
				my $Fevalue_aref = $readname_full_to_evalues{$full_read};
			    my @Fevalues = @{$Fevalue_aref};
				my $Fmin_evalue = min(@Fevalues);
				@FreadID_w_lowest_evalue = @{$m{$full_read}{$Fmin_evalue}};
				my $FreadID_w_lowest_evalue_count = @FreadID_w_lowest_evalue;
				if ($FreadID_w_lowest_evalue_count>1){
					#warn "$full_read has more than 1 best match\n";
				}
				##Leverage Full read position to infer viral-human junction
				for my $v (@VreadID_w_lowest_evalue){
					for my $h (@HreadID_w_lowest_evalue){
						for my $f (@FreadID_w_lowest_evalue){
							print $CSparsed_FULLoutFH "$v\t$h\t$f\n";
							##Parse for junction
							my @vs = split (/\t/, $v);
							my @hs = split (/\t/, $h);
							my @fs = split (/\t/, $f);
							my ($v_st, $v_ed) = ($vs[8], $vs[9]);
							my ($h_st, $h_ed) = ($hs[8], $hs[9]);
							my ($f_st, $f_ed) = ($fs[8], $fs[9]);
							#print "$v_st\t$v_ed\t$h_st\t$h_ed\t$f_st\t$f_ed\n";
							if (($v_st-$h_ed) == 1){ ## Human--Viral
								my $str = $h_ed++;
								print $CSparsed_EVALoutFH "$readname_core\t$str\n";
								$At_least_one_best_match = 1;
							}
							if (($h_st-$v_ed) == 1){ ## Viral--Human
								print $CSparsed_EVALoutFH "$readname_core\t$h_st\n";
								$At_least_one_best_match = 1;
							}
						}
					}
				}
			}
			if ($At_least_one_best_match == 0){
				#warn "no match\n";
				for my $v (@VreadID_w_lowest_evalue){
					for my $h (@HreadID_w_lowest_evalue){
						for my $f (@FreadID_w_lowest_evalue){
							#print "$v\t$h\t$f\n";
							my @vs = split (/\t/, $v);
							my @hs = split (/\t/, $h);
							my @fs = split (/\t/, $f);
							my ($v_st, $v_ed) = ($vs[8], $vs[9]);
							my ($h_st, $h_ed) = ($hs[8], $hs[9]);
							my ($f_st, $f_ed) = ($fs[8], $fs[9]);
							print $CSparsed_WGSIMCKoutFH "$readname_core\t$v_st\t$v_ed\t$h_st\t$h_ed\t$f_st\t$f_ed\n";
						}
					}
				}
			}
		}
		close ($CSparsed_FULLoutFH);
		close ($CSparsed_EVALoutFH);
	}
	else { #RPm
		##select read and location with best location
		print "Working $curr_file_4_search.parsed\n";
		open my $RPparsed_FULLoutFH, '>',"$curr_file_4_search.parsed.full.out";
		open my $RPparsed_EVALoutFH, '>',"$curr_file_4_search.parsed.eval.out";
		my %m = ();
		my @readname = ();
		my %readname_full_to_evalues = ();
		my %readcore_to_fullreadname = ();
		for (@result_a){
			my @s = split (/\t/, $_);
			my $readIDraw = $s[0];
			my @sRAWread = split (/\//, $readIDraw);
			my @sStrandANDlength = split (/_/, $sRAWread[1]);
			my ($readname_core, $strand, $len) = ($sRAWread[0], $sStrandANDlength[0], $sStrandANDlength[1]);
			push (@{$m{$readIDraw}{$s[-2]}}, $_); #readname evalue full_line
			push (@{$readname_full_to_evalues{$readIDraw}}, $s[-2]);
			push (@readname, $readname_core); #order of reads
			push (@{$readcore_to_fullreadname{$readname_core}}, $readIDraw);
		}
		#remove redundant record in readname array, in order
		my %redundnat = ();
		my @rename_non_redundnat = ();
		for (@readname){
			unless (defined $redundnat{$_}){
				push (@rename_non_redundnat, $_);
				$redundnat{$_} = "";
			}
		}
		for my $readname_core(@rename_non_redundnat){ #Parse in read pair
			my @readname_full_two = @{$readcore_to_fullreadname{$readname_core}};
			my ($forward_readname, $reverse_readname);
			for (@readname_full_two){
				if (/\/1/){
					$forward_readname = $_
				}
				if (/\/2/){
					$reverse_readname = $_;
				}
			}
			#Forward read
			my @FreadID_w_lowest_evalue = ();
			my @RreadID_w_lowest_evalue = ();
			if (defined $readname_full_to_evalues{$forward_readname}){
				my $Fevalue_aref = $readname_full_to_evalues{$forward_readname};
			    my @Fevalues = @{$Fevalue_aref};
				my $Fmin_evalue = min(@Fevalues);
				@FreadID_w_lowest_evalue = @{$m{$forward_readname}{$Fmin_evalue}};
				my $FreadID_w_lowest_evalue_count = @FreadID_w_lowest_evalue;
				if ($FreadID_w_lowest_evalue_count>1){
					warn "$forward_readname has more than 1 best match\n";
				}
			}
			#Reverse read
			if (defined $readname_full_to_evalues{$reverse_readname}){
				my $Revalue_aref = $readname_full_to_evalues{$reverse_readname};
			    my @Revalues = @{$Revalue_aref};
				my $Rmin_evalue = min(@Revalues);
				@RreadID_w_lowest_evalue = @{$m{$reverse_readname}{$Rmin_evalue}};
				my $RreadID_w_lowest_evalue_count = @RreadID_w_lowest_evalue;
				if ($RreadID_w_lowest_evalue_count>1){
					warn "$reverse_readname has more than 1 best match\n";
				}
			}
			##compute all combinations
			if ((defined $readname_full_to_evalues{$forward_readname}) && (defined $readname_full_to_evalues{$reverse_readname})){
				for my $F_result (@FreadID_w_lowest_evalue){
					for my $R_result (@RreadID_w_lowest_evalue){
						print $RPparsed_FULLoutFH "$F_result\t$R_result\n";
						my @Fs = split (/\t/, $F_result);
						my @Rs = split (/\t/, $R_result);
						my $readname_original_in = $Fs[0];
						my @readname_original_inS = split (/\//, $readname_original_in);
						my $readname_original = $readname_original_inS[0];
						my ($Fst, $Fed) = ($Fs[8], $Fs[9]);
						my ($Rst, $Red) = ($Rs[8], $Rs[9]);
						if ($Rst > $Fed){ #Forward --- Reverse
							print $RPparsed_EVALoutFH "$readname_original\t$Fed\t$Rst\n";
						} else { #Reverse -- Forward
							print $RPparsed_EVALoutFH "$readname_original\t$Red\t$Fst\n";
						}
					}
				}
			}
	    }
	    close ($RPparsed_FULLoutFH);
	    close ($RPparsed_EVALoutFH);
	}
}


sub file_2_array {
	my $curr_func = (caller(0))[3];
	my $count = @_;
	my $file_in = $_[0];
	my @out = ();
	open my $fh, '<', $file_in or die "Can't open $file_in for reading : $!";
	while (<$fh>){
		chomp;
		unless (/^\s*$/){ #push unless blank line
			unless (/^#/){ # remove comment line
				push (@out, $_);
			}
		}
	}
	close ($fh);
	return \@out;
}