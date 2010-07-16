#!/bin/env perl
use strict;
use warnings;

my $bin='/nas/RD_09C/resequencing/soft/bin/soap/soap2.20';
my $arg0='-p 6 -t -s 40 -l 32';
my $soap2patch=5;

unless (@ARGV>6){
	print "perl $0 PESE ReadLen,insSize gap outpath/ ref.index fqext,path/ fqname[,fqname2] [u]\n";
	exit;
}

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

my ($PESE,$rlins,$G,$opath,$Ref,$fqextpath,$fqname,$u) = @ARGV;
#print '[',join('] [',@ARGV),"]\n";
#print '[',join('] [',$PESE,$ins,$G,$Ref,$fqextpath,@fqnames),"]\n";
my ($readlen,$ins)=split ',',$rlins;
my ($ext,$path)=split ',',$fqextpath;
my @fqnames=split ',',$fqname;
my $mismatch=$readlen>70?3:1;
$mismatch = 5 if $readlen >= 90;
my ($min,$max)=(140,1200);
if ($ins > 1500) {
	$arg0 .= ' -R';
	($min,$max)=(1000,$ins*2-1000);
}
my ($n,$fqcmd,$tmpcmd,$avg,$std,$Lsd,$Rsd,$max_y,$max_x)=(0);
if ($PESE eq 'PE') {	# PE
	$fqcmd="-a $path$fqnames[0]$ext -b $path$fqnames[-1]$ext -o $opath$fqnames[0].soap -2 $opath$fqnames[0].single";
	$tmpcmd="-a $path$fqnames[0]$ext -b $path$fqnames[-1]$ext -o $opath$fqnames[0].tp -2 $opath$fqnames[0].ts";
	system("$bin $tmpcmd -D $Ref $arg0 -m $min -x $max -v $mismatch -g $G 2>$opath$fqnames[0].tlog")==0 or die "[x]system soap failed: $?";
###
	open TP,'<',"$opath$fqnames[0].tp";
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
	#system("bzip2 -9 $opath$fqnames[0].tp $opath$fqnames[0].ts");
	unlink qq($opath$fqnames[0].tp $opath$fqnames[0].ts);
###
	$arg0 .= " -m $min -x $max";
} else {	# SE
	$fqcmd="-a $path$fqnames[0]$ext  -o $opath$fqnames[0].se";
}
$fqcmd .= " -u $opath$fqnames[0].unmap";# if $G or $u && $u eq 'u';
# always put -abDo2 in the first for the poor case without break ...
my $sh="$bin $fqcmd -D $Ref $arg0 -v $mismatch -g $G 2>$opath$fqnames[0].log";
print "[$sh]\n";

my ($Pairs,$Paired,$Singled,$Reads,$Alignment);
if (-s "$opath$fqnames[0].nfo") {
	system("mv -f $opath$fqnames[0].archive $opath$fqnames[0].archive.0"); #if (-e "$opath$fqnames[0].archive");
	system("mv -f $opath$fqnames[0].nfo $opath$fqnames[0].nfo.0") unless -s "$opath$fqnames[0].nfo.0";
	# Check both .nfo and .log, whether goto END;
} else {
	RUN:
	open OUT,'>',"$opath$fqnames[0].archive" or warn "[!]Error opening $opath$fqnames[0].archive: $!\n";
	print OUT "#!/bin/sh\n$sh\n";
	close OUT;
	system($sh)==0 or die "[x]system [$sh] failed: $?";
}
open LOG,'<',"$opath$fqnames[0].log" or (warn "[x]Error opening $opath$fqnames[0].log: $!\n" and goto RUN);
while (<LOG>) {
	$Pairs = (split)[-2] if /^Total Pairs:/;
	$Paired = (split)[1] if /^Paired:/;
	$Singled = (split)[1] if /^Singled:/;
	$Reads = (split)[-1] if /^Total Reads/;
	$Alignment = (split)[1] if /^Alignment:/;
}
close LOG;
unless ($Pairs or $Reads) {
	system("mv -f $opath$fqnames[0].log $opath$fqnames[0].log.0");
	goto RUN;
}
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
	@ARGV=("$opath$fqnames[0].soap", "$opath$fqnames[0].single");
}
my ($BadLines,$BPOut,$ReadsOut,$MisSum,$TrimedBP,$TrimedReads,%Hit9r,%Hit9bp,%misMatch,%Indel)=(0,0,0,0,0,0);
my (%chrBPOut,%chrReadsOut,%chrMisSum,%chrTrimedBP,%chrTrimedReads,%chrHit9r,%chrHit9bp,%chrmisMatch,%chrIndel);
# for 46999 ChrIDs, one hash took 12m RES for {}=$ and 125m VIRT for {}{10}=$, thus fine to work with scaffolds. ( 1 gb VIRT max ? )
my (@lines,$hit,$len,$chr,$types,$trim,$mistr,$missed);
while (<>) {
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

		$hit=10 if $hit>10;	# max to count 9, then >=10
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
print NFO "\n#fmtC\tReadsOut\tBPOut\tMisSum\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\tBadLines\n";
print NFO join("\t",'ALL',$ReadsOut,$BPOut,$MisSum,$TrimedReads,$TrimedBP,
	${&combineC(\%misMatch)},${&combineJ(\%Hit9r)},${&combineJ(\%Hit9bp)},${&combineJ(\%Indel)},$BadLines),"\n\n";
print NFO "#fmtP\tReadsOut\tBPOut\tMisSum\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\n";
print NFO join("\t",";$_",$chrReadsOut{$_},$chrBPOut{$_},$chrMisSum{$_}||0,$chrTrimedReads{$_}||0,$chrTrimedBP{$_}||0,
	${&combineC(\%{$chrmisMatch{$_}})},${&combineJ(\%{$chrHit9r{$_}})},${&combineJ(\%{$chrHit9bp{$_}})},${&combineJ(\%{$chrIndel{$_}})}),"\n" for sort keys %chrReadsOut;
close NFO;

END:
system("bzip2 -9 $opath$fqnames[0].unmap");
