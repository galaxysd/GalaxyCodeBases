#!/bin/env perl
#use lib "/ifs1/ST_ASMB/USER/huxuesong/public/lib";
use lib '/export/data0/gentoo/tmp';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
our $opts='o:b';
our($opt_o, $opt_b);

#our $desc='';
our $help=<<EOH;
\t-o output prefix (./stat).{alistat,insert}
\t-b No pause for batch runs
EOH
our $ARG_DESC='soap_files{,.gz,.bz2}';

ShowHelp();
$opt_o='./stat' if ! $opt_o;
die "[x]No input files found !\n" unless @ARGV;
#die "[!]Max 252 files supported.\n" if @ARGV>252;

print STDERR "From [@ARGV] to [$opt_o]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <STDIN>;}

my $start_time = [gettimeofday];
#BEGIN
sub combineC($) {
	my $href=$_[0];
	if ($href and %$href) {
		my (@str,$m);
		$m = (sort {$a<=>$b} keys %$href)[-1];
		for (1..$m) {
			#$$href{$_} += 0;
			push @str,join(':',$_,$$href{$_}||0);
		}
		return \join(',',@str);
	} else {return \'.';}
}

sub combineJ($) {
	my $href=$_[0];
	if ($href and %$href) {
		my @str;
		for (sort {$a<=>$b} keys %$href) {
			push @str,join(':',$_,$$href{$_});
		}
		return \join(',',@str);
	} else {return \'.';}
}

sub getRealpos($$$$) {
	my ($len,$strand,$realpos,$trim)=@_;
	if ($strand eq '-') {	# Negative
		$realpos += $len;	# should be $len-1. So, starting 0. (+ & -, never meets.)
		if ($trim =~ /(\d+)S$/) {
			$realpos += $1;	# $1 only reset after next /()/
		}
	} elsif ($strand eq '+') {	# Positive
		if ($trim =~ /^(\d+)S/) {
			$realpos -= $1;
		}
	} else {
		$realpos=-1;
	}
	return $realpos;
}
=pod
case 'R': opt->FR = 0;

alnSeq[0] = mseqs->seqList+i;
alnSeq[1] = mseqs->seqList+i+1;
pe_aux.len = alnSeq[0]->len;

fprintf(stderr,"%c\n", "+-"[strain]);

int strain1 = (p->info >> 24)&1;
int strain2 = (q->info >> 24)&1;
if(o->FR){	# FR
	if(!strain1 && q->pos-p->pos + o->len >= o->min_ins && q->pos-p->pos + o->len <= o->max_ins) return TRUE;
	else if(strain1 && p->pos-q->pos+o->len >= o->min_ins && p->pos-q->pos+o->len <= o->max_ins) return TRUE;
	else{ return FALSE;};
}
else if(!o->FR){	# RF, (-R)
	if(strain1 && q->pos - p->pos >= o->min_ins && q->pos - p->pos + o->len <= o->max_ins) return TRUE;
	else if(!strain1 && p->pos-q->pos >= o->min_ins && p->pos-q->pos+o->len <= o->max_ins) return TRUE;
	else{ return FALSE;}
}

Well, InsertSize should always be defined from the 2 terminals of a sequencing fragument.
Thus, for FR, the same as what from the gel; for RF, should be less than the gel result by the length of what actually been sequenced.
=cut

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.bz2$/) {
	    open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}
sub checkfiletype($) {
    my $fhc=$_[0];
    my $type;
    my $head=<$fhc>;
    if ($head =~ /^\[bwa_sai2sam_pe_core\]/) {
        $type='samlog';
    } elsif ($head =~ /^\@SQ/) {
        $type='sam';
    } elsif ($head =~ /[ATCG]/) {
        $type='soap';
    } elsif ($head =~ /^$/) {
=pod
00000000  0a 42 65 67 69 6e 20 50  72 6f 67 72 61 6d 20 53  |.Begin Program S|
00000010  4f 41 50 61 6c 69 67 6e  65 72 2f 73 6f 61 70 32  |OAPaligner/soap2|
00000020  0a 46 72 69 20 4a 75 6c  20 31 35 20 30 39 3a 32  |.Fri Jul 15 09:2|
00000030  39 3a 34 32 20 32 30 31  31 0a 52 65 66 65 72 65  |9:42 2011.Refere|
00000040  6e 63 65 3a 20 2e 2f 72  65 66 2f 68 75 6d 61 6e  |nce: ./ref/human|
00000050  2e 66 61 2e 69 6e 64 65  78 0a 51 75 65 72 79 20  |.fa.index.Query |
00000060  46 69 6c 65 20 61 3a 20  2e 2f 30 6e 6f 6d 61 73  |File a: ./0nomas|
=cut
        $type='soaplog';
    }
    return $type;
}

sub sumsoapdata($$) {
    my ($hit,$len,$chr,$types,$trim,$mistr)=@$_[0];
    my $dathref=$_[1];
    my ($BPOut,$ReadsOut,$MisSum,$TrimedBP,$TrimedReads,%Hit9r,%Hit9bp,%misMatch,%Indel)=(0,0,0,0,0);
    #my (%chrBPOut,%chrReadsOut,%chrMisSum,%chrTrimedBP,%chrTrimedReads,%chrHit9r,%chrHit9bp,%chrmisMatch,%chrIndel);
    my (@trims,$missed);
		$BPOut += $len;#	$chrBPOut{$chr} += $len;
		++$ReadsOut;#	++$chrReadsOut{$chr};
		$missed=$mistr=~tr/ATCGatcg//;
		++$misMatch{$missed};
		#++$chrmisMatch{$chr}{$missed};
		$MisSum += $missed;
		#$chrMisSum{$chr} += $missed;

		$hit=4 if $hit>4;	# max to count 3, then >=4. Ancient Chinese wisdom, and a bit more ...
		++$Hit9r{$hit};
		#++$chrHit9r{$chr}{$hit};
		$Hit9bp{$hit} += $len;
		#$chrHit9bp{$chr}{$hit} += $len;
=pod
		if ($types < 100) {
			#++$misMatch{$types};
			#++$chrmisMatch{$chr}{$types};
		} elsif ($types < 200) {	# '3S33M9D39M', '32M1D14M29S' exists ?
			++$Indel{$types-100};
			++$chrIndel{$chr}{$types-100};
		} else {
			++$Indel{200-$types};
			++$chrIndel{$chr}{200-$types};
		}
=cut
		if ($types > 200) {
			++$Indel{200-$types};
			#++$chrIndel{$chr}{200-$types};
		} elsif ($types > 100) {
			++$Indel{$types-100};
			#++$chrIndel{$chr}{$types-100};
		}

		@trims = $trim =~ /(\d+)S/;
		if (@trims) {
			++$TrimedReads;
			#++$chrTrimedReads{$chr};
			for ( @trims ) {
				$TrimedBP += $_;
				#$chrTrimedBP{$chr} += $_;
			}
		}
    $dathref->{'BPOut'} += $BPOut;
    $dathref->{'ReadsOut'} += $ReadsOut;
    $dathref->{'MisSum'} += $MisSum;
    $dathref->{'TrimedReads'} += $TrimedReads;
    $dathref->{'TrimedBP'} += $TrimedBP;
    $dathref->{'misMatch'}{$_} += $misMatch{$_} for keys %misMatch;
    $dathref->{'Indel'}{$_} += $Indel{$_} for keys %Indel;
    $dathref->{'Hit9r'}{$_} += $Hit9r{$_} for keys %Hit9r;
    $dathref->{'Hit9bp'}{$_} += $Hit9bp{$_} for keys %Hit9bp;
}
sub statsoap($) {
    my $fh=$_[0];
    my (%datsum);
    $datsum{'PEuniqPairs'}=0;
    my ($BadLines,$PESE,$pairs,$lastpos,$line1,$line2,$pp,$pn,$calins,%insD)=(0,'PE',0);
	while ($line1=<$fh>) {
		my ($id1, $n1, $len1, $f1, $chr1, $x1, $types1, $m1, $mistr1)
		 = (split "\t", $line1)[0,3,5,6,7,8,9,-2,-1];
		unless ($types1) {    # soap2 output always more than 10 columes.
		    ++$BadLines;
		    last;
		}
		$line2=<$fh>;
		last unless $line2;
		my ($id2, $n2, $len2, $f2, $chr2, $x2, $types2, $m2, $mistr2)
		 = (split "\t", $line2)[0,3,5,6,7,8,9,-2,-1];
		unless ($types2) {    # soap2 output always more than 10 columes.
		    ++$BadLines;
		    last;
		}
			#($soapid,$hit,$len,$strand,$chr,$pos,$trim) = (split(/\t/))[0,3,5,6,7,8,-2];
			#        ($hit,$len,        $chr,$types,$trim,$mistr) = @lines[3,5,7,9,-2,-1];
		$id1 =~ s/\/[12]$//;
		$id2 =~ s/\/[12]$//;
		&sumsoapdata([$n1, $len1, $chr1, $types1, $m1, $mistr1],\%datsum);
		&sumsoapdata([$n2, $len2, $chr2, $types2, $m2, $mistr2],\%datsum);
		if (($PESE eq 'SE') or ($id1 ne $id2)){	# single
			$PESE='SE';
			next;
		}
		next if $n1+$n2>2 or $chr1 ne $chr2;
=pod
		if ($f1 eq '+') {
			($pp,$pn)=($x1,$x2);
		} else {
			($pp,$pn)=($x2,$x1);
		}
		if ($ins > 1500) {	# FR => +.pos < -.pos; RF => -.pos < +.pos
			next if $pp < $pn;
		} else { next if $pp > $pn; }
=cut
		++$pairs;
		$line1=&getRealpos($len1, $f1, $x1, $m1);	# $len,$strand,$realpos,$trim
		$line2=&getRealpos($len2, $f2, $x2, $m2);	# Well, $line{1,2} is recycled.
		$calins=abs($line1-$line2);	# -.starting=0
		++$insD{$calins};
	}
	if ($PESE eq 'PE') {
    	$datsum{'PEuniqPairs'}=$pairs;
	}
	return [$PESE,[0,0,0],\%datsum,\%insD];
}
sub statsoaplog($) {
    my $fh=$_[0];
    my ($Pairs,$Paired,$Singled,$Reads,$Alignment)=(0,0,0,0,0);
    while(<$fh>) {
	    $Pairs = (split)[-2] if /^Total Pairs:/;
	    $Paired = (split)[1] if /^Paired:/;
	    $Singled = (split)[1] if /^Singled:/;
	    $Reads = (split)[-1] if /^Total Reads/;
	    $Alignment = (split)[1] if /^Alignment:/;
    }
    my @RET;
    if ($Reads) {
        @RET=['SE',[$Reads,0,$Alignment]];
    } else {
        @RET=['PE',[$Pairs*2,$Paired*2,$Singled]];
    }
    return \@RET;
=pod
Total Pairs: 34776407 PE
Paired:      17719335 (50.95%) PE
Singled:     30467634 (43.81%) SE

IN:	34776407 reads x 2 fq file = 69552814 reads from _1 & _2.fq
-o:	17719335 pairs = 17719335 x 2 =35438670 lines in .soap
-2:	30467634 lines in .single

17719335/34776407 = 0.50952172833726037310294878939046
30467634/69552814 = 0.43805034257851882168275750856033


Total Reads: 25
Alignment:   22 (88.00%)
=cut
}
sub statsam($) {
    my $fh=$_[0];
    return 'a';
}
my %dostat=(
    'sam' => \&statsam,
    'soap' => \&statsoap,
#    'samlog' => sub {},
    'soaplog' => \&statsoaplog,
);

sub mergehash($$) {
    my ($inhref,$intohref)=@_;
}

my $files=0;
while($_=shift @ARGV) {
    ++$files;
    my $infile;
    $infile=openfile($_);
    my $type=checkfiletype($infile);
    close $infile;
    print STDERR "$files\t[$type] $_ ...";
    $infile=openfile($_);
    if (exists $dostat{$type}) {
        my $ret=&{$dostat{$type}}($infile);
        print "$ret\n";
    } else {
        print STDERR "\b\b\bskipped."
    }
    close $infile;
    print STDERR "\n";
}
warn "[!]Files Read: $files\n";


__END__

my ($PESE,$rlins,$G,$opath,$Ref,$fqextpath,$fqname,$u) = @ARGV;
#print '[',join('] [',@ARGV),"]\n";
#print '[',join('] [',$PESE,$ins,$G,$Ref,$fqextpath,@fqnames),"]\n";
my ($readlen,$ins)=split ',',$rlins;
my ($ext,$path)=split ',',$fqextpath;
my @fqnames=split ',',$fqname;

}
my ($n,$fqcmd,$tmpcmd,$avg,$std,$Lsd,$Rsd,$max_y,$max_x)=(0);
$ext .= '.gz' unless -s "$path$fqnames[0]$ext";	# Well, sometime we gzip those .fqs later on and ...
if ($PESE eq 'PE') {	# PE
###
	#open TP,'<',"$opath$fqnames[0].tp";
	open( TP,"-|","gzip -dc $opath$fqnames[0].tp.gz") or die "Error: $!\n";
	my ($pairs,$lastpos,$line1,$line2,$pp,$pn,$calins,%insD)=(0);
	while ($line1=<TP>) {
		last if eof TP;
		$lastpos=tell TP;
		my ($id1, $n1, $len1, $f1, $chr1, $x1, $m1) = (split "\t", $line1)[0,3,5,6,7,8,-2];
		$line2=<TP>;
		my ($id2, $n2, $len2, $f2, $chr2, $x2, $m2) = (split "\t", $line2)[0,3,5,6,7,8,-2];
			#($soapid,$hit,$len,$strand,$chr,$pos,$trim) = (split(/\t/))[0,3,5,6,7,8,-2];
		$id1 =~ s/\/[12]$//;
		$id2 =~ s/\/[12]$//;
		if ($id1 ne $id2){	# single
			seek (TP, $lastpos, 0);
			next;
		}
		next if $n1+$n2>2 or $chr1 ne $chr2;
		if ($f1 eq '+') {
			($pp,$pn)=($x1,$x2);
		} else {
			($pp,$pn)=($x2,$x1);
		}
		if ($ins > 1500) {	# FR => +.pos < -.pos; RF => -.pos < +.pos
			next if $pp < $pn;
		} else { next if $pp > $pn; }
		++$pairs;
		$line1=&getRealpos($len1, $f1, $x1, $m1);	# $len,$strand,$realpos,$trim
		$line2=&getRealpos($len2, $f2, $x2, $m2);	# Well, $line{1,2} is recycled.
		$calins=abs($line1-$line2);	# -.starting=0
		++$insD{$calins};
	}
	close TP;
	my ($sum,$sum2,$v)=(0,0);
	open O,'>',"$opath$fqnames[0].insD";
	for my $k (sort {$a <=> $b} keys %insD) {
		$v=$insD{$k};
		print O "$k\t$v\n";
		$sum += $k * $v;
		$n += $v;
		$sum2 += $k*$k * $v;
	}
	$avg = $sum/$n;
	$std = sqrt($sum2/$n-$avg*$avg);
	print O "# $avg ± $std\n";
	($max_y,$max_x)=(-1,0);
	for my $k ($avg-$std .. $avg+$std) {	# 68.27%
		next unless $v=$insD{$k};
		if ($max_y < $v) {
			$max_x = $k;
			$max_y = $v;
		}
	}
	my $cutoff = $max_y / 1000;
	$cutoff = 3 if $cutoff<3;
	my ($diff, $Lc, $Rc);
	for my $k (keys %insD) {
		$v=$insD{$k};
		next if $v < $cutoff;
		$diff = $k - $max_x;
		if ($diff < 0) {
			$Lsd += $v * $diff * $diff;
			$Lc += $v;
		}
		elsif ($diff > 0) {
			$Rsd += $v * $diff * $diff;
			$Rc += $v;
		}
	}
	$Lsd = sqrt($Lsd/$Lc);
	$Rsd = sqrt($Rsd/$Rc);
	print O "# +$Lsd -$Rsd\n";
	close O;
	($min,$max)=(int($avg-$Lsd*2.576-$soap2patch),int($avg+$Rsd*2.576+$soap2patch+.5));	# 99%
}
my ($Pairs,$Paired,$Singled,$Reads,$Alignment);

open NFO,'>',"$opath$fqnames[0].nfo" or die "[x]Error opening $opath$fqnames[0].nfo: $!\n";
if ($max==0) {
	print NFO "#fmtS\tTotal_Reads\tAlignment\n";
	print NFO "Summary\t",join("\t",$Reads,$Alignment),"\n";
	@ARGV=("$opath$fqnames[0].se");
} else {
	my $p=sprintf "%.2f",100*$max_y/$n;
	$avg=int($avg*10+.5)/10;
	$Lsd=int($Lsd*100+.5)/100;
	$Rsd=int($Rsd*100+.5)/100;
	print NFO "#fmtS\tTotal_Pairs\tPaired\tSingled\tMode(p%),Lsd,Rsd,InsAvg,STD\tInsMin,InsMax\n";
	print NFO "Summary\t",join("\t",$Pairs,$Paired,$Singled,"$max_x($p %),$Lsd,$Rsd,$avg,$std","$min,$max"),"\n";
	#@ARGV=("$opath$fqnames[0].soap", "$opath$fqnames[0].single");
	open(OUTS,"-|","gzip -dc $opath$fqnames[0].soap.gz $opath$fqnames[0].single.gz") or die "Error: $!\n";
}
my ($BadLines,$BPOut,$ReadsOut,$MisSum,$TrimedBP,$TrimedReads,%Hit9r,%Hit9bp,%misMatch,%Indel)=(0,0,0,0,0,0);
my (%chrBPOut,%chrReadsOut,%chrMisSum,%chrTrimedBP,%chrTrimedReads,%chrHit9r,%chrHit9bp,%chrmisMatch,%chrIndel);
# for 46999 ChrIDs, one hash took 12m RES for {}=$ and 125m VIRT for {}{10}=$, thus fine to work with scaffolds. ( 1 gb VIRT max ? )
my (@lines,$hit,$len,$chr,$types,$trim,$mistr,$missed);
while (<OUTS>) {
	@lines = split /\t/;
	if (@lines > 10) {	# soap2 output always more than 10 columes.
		($hit,$len,$chr,$types,$trim,$mistr) = @lines[3,5,7,9,-2,-1];
		$BPOut += $len;	$chrBPOut{$chr} += $len;
		++$ReadsOut;	++$chrReadsOut{$chr};
		$missed=$mistr=~tr/ATCGatcg//;
		++$misMatch{$missed};
		++$chrmisMatch{$chr}{$missed};
		$MisSum += $missed;
		$chrMisSum{$chr} += $missed;

		$hit=4 if $hit>4;	# max to count 3, then >=4. Ancient Chinese wisdom, and a bit more ...
		++$Hit9r{$hit};
		++$chrHit9r{$chr}{$hit};
		$Hit9bp{$hit} += $len;
		$chrHit9bp{$chr}{$hit} += $len;
=pod
		if ($types < 100) {
			#++$misMatch{$types};
			#++$chrmisMatch{$chr}{$types};
		} elsif ($types < 200) {	# '3S33M9D39M', '32M1D14M29S' exists ?
			++$Indel{$types-100};
			++$chrIndel{$chr}{$types-100};
		} else {
			++$Indel{200-$types};
			++$chrIndel{$chr}{200-$types};
		}
=cut
		if ($types > 200) {
			++$Indel{200-$types};
			++$chrIndel{$chr}{200-$types};
		} elsif ($types > 100) {
			++$Indel{$types-100};
			++$chrIndel{$chr}{$types-100};
		}

		@lines = $trim =~ /(\d+)S/;
		if (@lines) {
			++$TrimedReads;
			++$chrTrimedReads{$chr};
			for ( @lines ) {
				$TrimedBP += $_;
				$chrTrimedBP{$chr} += $_;
			}
		}
	} else {	# rare event ...
		++$BadLines;
		next;
	}
}
close OUTS;

print NFO "\n#fmtC\tReadsOut\tBPOut\tMisSum\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\tBadLines\n";
print NFO join("\t",'ALL',$ReadsOut,$BPOut,$MisSum,$TrimedReads,$TrimedBP,
	${&combineC(\%misMatch)},${&combineJ(\%Hit9r)},${&combineJ(\%Hit9bp)},${&combineJ(\%Indel)},$BadLines),"\n\n";
print NFO "#fmtP\tReadsOut\tBPOut\tMisSum\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\n";
print NFO join("\t",";$_",$chrReadsOut{$_},$chrBPOut{$_},$chrMisSum{$_}||0,$chrTrimedReads{$_}||0,$chrTrimedBP{$_}||0,
	${&combineC(\%{$chrmisMatch{$_}})},${&combineJ(\%{$chrHit9r{$_}})},${&combineJ(\%{$chrHit9bp{$_}})},${&combineJ(\%{$chrIndel{$_}})}),"\n" for sort keys %chrReadsOut;
close NFO;

END:
#system("bzip2 -9 $opath$fqnames[0].unmap");
