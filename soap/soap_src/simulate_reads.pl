#! /usr/bin/perl -w

use strict;
use Getopt::Long;
srand (time ^ $$);

############################## to be edit ####################################
my %options;
GetOptions (\%options, "para=s", "output=s", "snp=s", "indel=s", "fasta=s", "change=s", "type=s", "error=s", "depth=s", "qual=s");
$options{para} and $options{fasta} and $options{type} 
and $options{depth} or &usage();

sub usage {
    print <<HELPOUT;

usage: $0 [options] <output file>
       
options: 
	 -p|--para       fasta sequence parameter [seq_name](:[seq_length])
	 -s|--snp        snp parameters [R:r,0.2] or [R:n,2000] or [D:dbfile]
	 -i|--indel      indel parameters [indel_num:indel_maxlength]
	 -f|--fasta      fasta input file
	 -t|--type       single or pair
	 -e|--error      sequence error: 1 or 0
	 -d|--depth      coverage of reads
	 -q|--qual       add quality to fastq: 1 or 0

HELPOUT
exit 0;
}

###################### you can modify these parameters ###################
## copy
my $known_snp_percent = 0.5;
my $three_snp_rate = 0.5;
my $extend_insert_prob = 0.5;
my $extend_delete_prob = 0.5;
my $insert_percentage = 0.5;
#my $multi_snp_rate = 0.2;

## read
my $read_len = 32;
my $pair_dist = 500;
my $pair_dist_sd = 4;
my $sequence_error_rate = 0.001;
my $quality_error_rate = 0.01;

################## main ####################################
my ($seq, $seq_file, $seq_para);
$seq_file = $options{fasta};
$seq_para = $options{para};

print "begining to read sequence!\n";
my ($start_pos, $seq_length, $seq_name) = &read_seq (\$seq, $seq_file, $seq_para);
print "reading sequence is over!\n";  # debug

print "begining to generate copy!\n";  #debug
my $snp = $options{snp};
my $known_snp = "knownSnp.tab";
my $indel = $options{indel};  
my $copy;
my %indel;
&generate_copy (\$copy, \%indel, $seq, $seq_name, $seq_length, $start_pos, $snp, $known_snp, $indel);
print "generating copy is over!\n";   # debug

print "begining to generate reads!\n";  # debug
my $type = $options{type};  # single or pair
my $sequence_error = $options{error};
my $coverage = $options{depth};
my $outfile = "reads.fq";
my $real_align = "right_align.tab";
my $qual = $options{qual};
&generate_reads ($seq_name, $start_pos, $seq_length, \$copy, \%indel,
		 $type, $sequence_error, $qual, $coverage, $outfile, $real_align);
print "generating reads is over!\n";   # debug

##############################################################################

### 
sub read_seq {
    my ($seq_pt, $seq_file, $seq_para) = @_;
# read seq_name and seq_length

    my ($seq_name, $seq_length, $start_pos);
    if ($seq_para =~ /:/) {
	($seq_name, $seq_length) = split /:/, $seq_para;
	if ($seq_length =~ /^(\d+)([mMkK])$/) {
	    if ($2 eq 'm' or $2 eq 'M') {
		$seq_length = $1 * 1000000;
	    } else {
		$seq_length = $1 * 1000;
	    }
	}
    } elsif ($seq_para =~ /^\S+$/) {
	$seq_name = $seq_para;
	$seq_length = -1;
    } else {
	die "error in sequence parameters: must [seq_name](:[seqlength])\n";
    }

# fetch sequence    
    open SEQFILE, $seq_file or
	die "error in open $seq_file: $!\n";
    my $to_read = 0;
    my $seq_whole_len = 0;
    $$seq_pt = '';
    while (<SEQFILE>) {
	chomp;
	if (/^>/) {
	    $to_read = 0;
	    if (/^>$seq_name/) {
		$to_read = 1;
	    }
	} elsif ($to_read == 1) {
	    $$seq_pt .= $_;
	    $seq_whole_len += length ($_);
	}
    }

# random fetch a certain length sub-sequence    
    if ($seq_length != -1) {
	my $random_max = $seq_whole_len - $seq_length + 1;
	if ($random_max < 1) {
	    die "the sequence length is bigger than whole length!";
	}

	$start_pos = &rand_num ($random_max);
	$$seq_pt = substr $$seq_pt, $start_pos, $seq_length;
    } else {
	$start_pos = 0;
	$seq_length = $seq_whole_len;
    }

    return (($start_pos + 1), $seq_length, $seq_name);
}

###
sub generate_copy {
    my ($copy_pt, $indel_hash_pt, $seq, $seq_name, $seq_length, $start_pos, $snp, $known_snp, $indel) = @_;
    my (@change_pos, $flag);

######################
# deal with snp

    my ($snp_num, %snp_known, %snp_unknown);
    
    if ($snp) {
	my @snp = split /:/, $snp;

	if ($snp[0] eq 'R') {
	    my @r_snp = split /,/, $snp[1];
	    if ($r_snp[0] eq 'n') {
		$snp_num = $r_snp[1];
	    } else {
		$snp_num = int ($seq_length * $r_snp[1]);
	    }

	    my $known_snp_num = int ($snp_num * $known_snp_percent);
	    foreach (1..$known_snp_num) {
		my $randnum = &rand_num ($seq_length);
		$flag = 0;
		foreach (@change_pos) {
		    if ($randnum - $_ > -32 and $randnum - $_ < 32) {
			$flag = 1;
			last;
		    }
		}
		if ($flag) {
		    redo;
		} else {
		    $snp_known{$randnum} = -1;
		    push @change_pos, $randnum;
		}

	    }
	    foreach (($known_snp_num+1)..$snp_num) {
		my $randnum = &rand_num ($seq_length);
		$flag = 0;
		foreach (@change_pos) {
		    if ($randnum - $_ > -32 and $randnum - $_ < 32) {
			$flag = 1;
			last;
		    }
		}
		if ($flag) {
		    redo;
		} else {
		    $snp_unknown{$randnum} = -1;
		    push @change_pos, $randnum;
		}
	    }

	    &set_snp (\$seq, $copy_pt, $seq_length, \%snp_known, \%snp_unknown);
	    
	} elsif ($snp[0] eq 'D') {
	    open DBSNP, $snp[1] or
		die "error open $snp[1]: $!\n";
	    my $relative_pos;

	    while (<DBSNP>) {
		chomp;
		my @temp = split /\t/, $_;
		if (@temp == 22 and $temp[21] eq "reference" and $temp[6] eq $seq_name) {
		    if ($temp[11] =~ /^\d+$/) {

			$relative_pos = $temp[11] - $start_pos;
			next if ($temp[16] eq '0');
			if ($temp[13] ne '0' and $temp[13] !~ /^\s*$/) {
			    next if (&rand_num_float(1) > $temp[13]);
			} else {
			    next;
			}
			next if ($relative_pos < 0 or $relative_pos > $seq_length);
			next if (exists $snp_known{$relative_pos});
			$snp_known{$relative_pos} = -1;
			push @change_pos, $relative_pos;
		    } else {
			next;
		    }
		}
	    }

## must add some random SNP;
	    my $unknown_snp_num = ((scalar (keys %snp_known)) / $known_snp_percent);
	    $unknown_snp_num = (scalar (keys %snp_known)) - $unknown_snp_num;
	    foreach (1..$unknown_snp_num) {
		my $randnum = &rand_num ($seq_length);
		next if (exists $snp_known{$randnum});
		$snp_unknown{$randnum} = -1;
		push @change_pos, $randnum;
	    }

	    &set_snp (\$seq, $copy_pt, $seq_length, \%snp_known, \%snp_unknown);
	}
	
	$seq = $$copy_pt;
    }  # end SNP

    print "SNP is over!\n";  #debug
	
# deal with indel
    
    if ($indel) {
	my ($indel_num, $indel_maxlength) = split /:/, $indel;
	my @indel_pos;
	my $insert_num = int ($indel_num * $insert_percentage);
	my $delete_num = $indel_num - $insert_num;
	$$copy_pt = '';

	foreach (1..$insert_num) {
	    my $rand_num = &rand_num ($seq_length + 1);
	    $flag = 0;
	    foreach (@change_pos) {
		if ($rand_num - $_ > -32 + $indel_maxlength  and $rand_num - $_ < 32 - $indel_maxlength) {
		    $flag = 1;
		    last;
		}
	    }
	    if ($flag) {
		redo;
	    } else {
		$$indel_hash_pt{$rand_num} = '+';
		push @change_pos, $rand_num;
	    }
	}

	foreach (1..$delete_num) {
	    my $rand_num = &rand_num ($seq_length - $indel_maxlength);
	    $flag = 0;
	    foreach (@change_pos) {
		if ($rand_num - $_ > -32 + $indel_maxlength  and $rand_num - $_ < 32 - $indel_maxlength) {
		    $flag = 1;
		    last;
		}
	    }
	    if ($flag) {
		redo;
	    } else {
		$$indel_hash_pt{$rand_num} = '-';
		push @change_pos, $rand_num;
	    }
	}
	
	@indel_pos = sort {$a<=>$b} (keys %$indel_hash_pt);

	my ($random_length_prob, $random_length, $begin, $me, $end, $check_case);
	my $last = 0;
	my @length_array = (2..$indel_maxlength);
	@length_array = reverse @length_array;

	foreach my $pos (@indel_pos) {
	    $$copy_pt .= substr $seq, $last, ($pos - $last);

	    if ($$indel_hash_pt{$pos} eq '+') {
		$random_length_prob = &rand_num_float (1);
		$random_length = 1;

		foreach (@length_array) {
		    if ($random_length_prob <= ($extend_insert_prob ** ($_ - 1))) {
			$random_length = $_;
			last;
		    }
		}
		
		$me = &rand_seq ($random_length);
		$$copy_pt .= $me;
		$$indel_hash_pt{$pos} .= ((length $me ) . ":$me");
		$last = $pos;

	    } elsif ($$indel_hash_pt{$pos} eq '-') {
		$random_length_prob = &rand_num_float (1);
		$random_length = 1;
		foreach (@length_array) {
		    if ($random_length_prob <= ($extend_delete_prob ** ($_ - 1))) {
			$random_length = $_;
			last;
		    }
		}
		
		$me = substr $seq, $pos, $random_length;
		$me = uc ($me);
		$$indel_hash_pt{$pos} .= ($random_length . ":$me");
		$last = $pos + $random_length;
	    }
	}
	
	$$copy_pt .= substr $seq, $last, ($seq_length - $last);
    } # end indel

######### print change_file
    open OUT, ">$known_snp" or
	die "error write $known_snp: $!\n";
    
    my $index = 0;
    @change_pos = sort {$a<=>$b} @change_pos;
    foreach (@change_pos) {
		if (exists $snp_known{$_}) {
			$index ++;
			print OUT "$index\t$seq_name\t$_\t$snp_known{$_}\t0\t0\n";
		}
	}
    close OUT;
}
	
###
sub generate_reads {
    my ($seq_name, $start_pos, $seq_len, $copy_pt, $indel_hash_pt,
	$type, $sequence_error, $qual, $coverage, $outfile, $real_align) = @_;

    my $copy_length = length ($$copy_pt);
    my @indel_pos;

    @indel_pos = sort {$a<=>$b} (keys %$indel_hash_pt);

    if ($type eq 'single') {
	open READS, ">$outfile" or
	    die "error write $outfile: $!\n";
	my ($read, $read_pos, $read_qual, $read_change);
	foreach (1 .. ($coverage * $copy_length)) {
	    $read_pos = &rand_num ($copy_length - $read_len + 1);
	    $read = substr $$copy_pt, $read_pos, $read_len;
	    if ($sequence_error) {
		($read, $read_change) = &seq_error ($read);
	    }
	    
	    my $tmp_length;
	    foreach my $tmp_pos (@indel_pos) {
		if ($read_pos > $tmp_pos) {
		    if ($$indel_hash_pt{$tmp_pos} =~ /^([+-]\d+):/) {
			$tmp_length = $1;
		    }
		    my $tmp_read_pos = $read_pos - $tmp_length;
		    if ($tmp_read_pos < $tmp_pos) {
			$read_pos = $tmp_pos;
			last;
		    } else {
			$read_pos = $tmp_read_pos;
		    }
		} else {
		    last;
		}
	    }
	    
	    print READS "\@$seq_name:$start_pos:$seq_len\__ref:+:$read_pos\__$read_change\n$read\n";
	    if ($qual) {
		$read_qual = &read_qual ();
		print READS "+$seq_name:$start_pos:$seq_len\__ref:+:$read_pos\__$read_change\n$read_qual\n";
#		print READS "$read_qual\n";
	    }
	}

    } elsif ($type eq 'pair') {
	
	my $index = 0;

	open ALIGN, ">$real_align" or
		die "error write $real_align: $!\n";

	open READSA, ">$outfile.a" or
	    die "error write $outfile.a:$!\n";
	open READSB, ">$outfile.b" or
	    die "error write $outfile.b:$!\n";
	my ($reada, $reada_pos, $reada_qual, $changea, $readb, $readb_pos, $readb_qual, $changeb);
	my $loop_num = int ($coverage * $copy_length/$read_len / 2);

	foreach (1 .. $loop_num) {
	    $reada_pos = &rand_num ($copy_length - $pair_dist - 2 * $read_len);
	    $reada = substr $$copy_pt, $reada_pos, $read_len;
	    $readb_pos = int (&norm_rand ($pair_dist, $pair_dist_sd) + $reada_pos);
	    if ($readb_pos + $read_len <= $copy_length) {
		$readb = substr $$copy_pt, $readb_pos, $read_len;
	    } else {
		$readb_pos = $copy_length - $read_len;
		$readb = substr $$copy_pt, $readb_pos, $read_len;
	    }
	    
	    my $pos_diff = 0;
	    my $tmp_length;
	    my ($reada_pos_tmp, $readb_pos_tmp) = ($reada_pos, $readb_pos);

# 	    foreach my $tmp_pos (@indel_pos) {
# 		if ($reada_pos_tmp - $pos_diff > $tmp_pos) {
# 		    if ($$indel_hash_pt{$tmp_pos} =~ /^([+-]\d+):/) {
# 			$tmp_length = $1;
# 		    } else {
# 			print STDERR "error in $tmp_pos of indel\n";
# 		    }
# 		    $pos_diff += $tmp_length;
# 		    if ($reada_pos_tmp - $pos_diff < $tmp_pos) {
# 			$reada_pos = $tmp_pos;
# 		    } else {
# 			$reada_pos = $reada_pos_tmp - $pos_diff;
# 		    }
# 		} elsif ($readb_pos_tmp - $pos_diff > $tmp_pos) {
# 		    if ($$indel_hash_pt{$tmp_pos} =~ /^([+-]\d+):/) {
# 			$tmp_length = $1;
# 		    } else {
# 			print STDERR "error in $tmp_pos of indel\n";
# 		    }
# 		    $pos_diff += $tmp_length;
# 		    if ($readb_pos_tmp - $pos_diff < $tmp_pos) {
# 			$readb_pos = $tmp_pos;
# 		    } else {
# 			$readb_pos = $readb_pos_tmp - $pos_diff;
# 		    }
# 		}
# 	    }

	    foreach my $tmp_pos (@indel_pos) {
		if ($reada_pos > $tmp_pos) {
		    if ($$indel_hash_pt{$tmp_pos} =~ /^([+-]\d+):/) {
			$tmp_length = $1;
		    }
		    my $tmp_read_pos = $reada_pos - $tmp_length;
		    if ($tmp_read_pos < $tmp_pos) {
			$reada_pos = $tmp_pos;
			last;
		    } else {
			$reada_pos = $tmp_read_pos;
		    }
		} else {
		    last;
		}
	    }

	    foreach my $tmp_pos (@indel_pos) {
		if ($readb_pos > $tmp_pos) {
		    if ($$indel_hash_pt{$tmp_pos} =~ /^([+-]\d+):/) {
			$tmp_length = $1;
		    }
		    my $tmp_read_pos = $readb_pos - $tmp_length;
		    if ($tmp_read_pos < $tmp_pos) {
			$readb_pos = $tmp_pos;
			last;
		    } else {
			$readb_pos = $tmp_read_pos;
		    }
		} else {
		    last;
		}
	    }

		
	    if ($sequence_error) {
		($reada, $changea) = &seq_error ($reada);
		($readb, $changeb) = &seq_error ($readb);
	    }
	    

		$readb = reverse $readb;
		$readb =~ tr/atgcnATGCN/tacgnTACGN/;
		print READSA "\@$index\n$reada\n";
		print READSB "\@$index\n$readb\n";
		my $loc_a = $start_pos+$reada_pos;
		my $loc_b = $start_pos+$readb_pos;
		print ALIGN "$index\ta\t+\t$seq_name\t$loc_a\t$changea\n";
		print ALIGN "$index\tb\t-\t$seq_name\t$loc_b\t$changeb\n";
		print ALIGN "\n";

		if ($qual) {
 		    $reada_qual = &read_qual();
 		    $readb_qual = &read_qual();
 		    print READSA "+$index\n$reada_qual\n";
 		    print READSB "+$index\n$readb_qual\n";
 		}
		$index++;
	}
	}
}

###
sub rand_num {
    my ($max) = @_;
    return (int (rand ($max)));
}

###
sub rand_num_float {
    my ($max) = @_;
    return (rand ($max));
}

###
sub set_snp {
    my ($seq_pt, $copy_pt, $seq_length, $snp_known_pt, $snp_unknown_pt) = @_;
    my @snp_pos = ((keys %$snp_known_pt), (keys %$snp_unknown_pt));
    @snp_pos = sort {$a<=>$b} @snp_pos;
    my @base = ('A', 'T', 'G', 'C');
    my %base = ('A' => 0,
		'T' => 1,
		'G' => 2,
		'C' => 3,
		);
    my $last = 0;
    $$copy_pt = '';
    foreach my $pos (@snp_pos) {
	$$copy_pt .= substr $$seq_pt, $last, ($pos - $last);
	my $me = substr $$seq_pt, $pos, 1;
	my $me_new = uc $me;
	my $rand_num;
	if (exists $base{$me_new}) {
	    $rand_num = ($base{$me_new} + &rand_num (3) + 1) % 4;
	    $me_new = $base[$rand_num];
	} else {
	    $rand_num = &rand_num (4);
	    $me_new = $base[$rand_num];
	}

	if (exists $$snp_known_pt{$pos}) {
	    $$snp_known_pt{$pos} = "$me/$me_new";
	} else {
	    $$snp_unknown_pt{$pos} = "$me/$me_new";
	}
	
	$$copy_pt .= $me_new;
	$last = $pos + 1;
    }
    $$copy_pt .= substr $$seq_pt, $last, ($seq_length - $last);
}

###
sub rand_seq {
    my ($len) = @_;
    my %base = (0 => 'A',
		1 => 'T',
		2 => 'G',
		3 => 'C',
		);
    my $string = '';
    foreach (1..$len) {
	$string .= $base{&rand_num (4)};
    }
    return $string;
}

###
sub norm_rand {
    my ($mean, $sd) = @_;
    my $sum = 0;
    foreach (0..11) {
	$sum += &rand_num_float (1);
    }
    return ($mean + $sd * ($sum - 6));
}

###
sub seq_error {
    my ($read) = @_;
    my $prob = $read_len * $sequence_error_rate;
    my $change = "no_error";
    if (&rand_num_float (1) < $prob) {
	my $error_pos = &rand_num ($read_len);
	my $me = substr $read, $error_pos, 1;
	my $begin = substr $read, 0, $error_pos;
	my $end = substr $read, ($error_pos + 1), ($read_len - $error_pos - 1);
	my @base = ('A', 'T', 'G', 'C');
	my %base = ('A' => 0,
		    'T' => 1,
		    'G' => 2,
		    'C' => 3,
		    );
	my $me_new = uc ($me);
	my $rand_num;
	if (exists $base{$me_new}) {
	    $rand_num = ($base{$me_new} + &rand_num (3) + 1) % 4;
	    $me_new = $base[$rand_num];
	} else {
	    $rand_num = &rand_num (4);
	    $me_new = $base[$rand_num];
	}

	$read = $begin . $me_new . $end;
	$change = "$error_pos:$me/$me_new";
    }
    return ($read, $change);
}

###
sub read_qual {
    my $temp;
    foreach (1..$read_len) {
	if (&rand_num_float (1) < $quality_error_rate) {
	    $temp .= '0';
	} else {
	    $temp .= '9';
	}
    }
    return $temp;
}

