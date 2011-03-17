#!/bin/env perl
use strict;
use warnings;

unless (@ARGV>0){
	print "perl $0 base\n";
	exit;
}
my ($soapbase)=@ARGV;

### Erf < ###
sub gam1($) {
	my $x1 = $_[0];
	my $t = 0;
	my @a = (0.0000677106,-0.0003442342,0.0015397681,-0.0024467480,
					 0.0109736958,-0.0002109075,0.0742379071,0.0815782188,
					 0.4118402518,0.4227843370,1.0);
	unless($x1 > 0) {
		print "x must greater than 0!\n";
		exit 0;
	}
	if($x1 < 1) {
		$t = 1/($x1*($x1+1));
		$x1 += 2;
	} elsif($x1 <= 2) {
		$t = 1/$x1;
		$x1 += 1;
	} elsif($x1 <= 3) {
		$t = 1;
	} else {
		$t = 1;
		while($x1 > 3) {
			$x1 -= 1;
			$t *= $x1;
		}
	}
	my $s1 = $a[0];
	my $u = $x1-2;
	foreach (@a) {
		$s1 = $s1*$u+$_;
	}
	$s1 *= $t;
}

sub gam2($$) {
	my ($a_gam,$x_gam) = ($_[0],$_[1]);
	my ($n_gam,$p_gam,$q_gam,$d_gam,$s_gam,$s1_gam,$p0_gam,$q0_gam,$p1_gam,$q1_gam,$qq_gam);
	if(($a_gam <= 0)||($x_gam < 0)) {
		if($a_gam <= 0) {
			print "err**a<=0!\n";
		}
		if($x_gam < 0) {
			print "err**x<0!\n";
		}
		exit 0;
	}
	if($x_gam == 0) {
		return 0;
	}
	if($x_gam > 1e35) {
		return 1;
	}
	$q_gam = $a_gam*log($x_gam);
	$qq_gam = exp($q_gam);
	if($x_gam < 1+$a_gam) {
		$p_gam = $a_gam;
		$s_gam = $d_gam = 1/$a_gam;
		for($n_gam = 1; $n_gam <= 100; $n_gam++) {
			$p_gam++;
			$d_gam = $d_gam*$x_gam/$p_gam;
			$s_gam += $d_gam;
			 if(abs($d_gam) < abs($s_gam)*1e-7) {
				$s_gam=$s_gam*exp(-$x_gam)*$qq_gam/&gam1($a_gam);
				return $s_gam;
			}
		}
	} else {
		$s_gam = 1/$x_gam;
		$p0_gam = 0;
		$q0_gam = $p1_gam=1;
		$q1_gam = $x_gam;
		for($n_gam = 1; $n_gam <= 100; $n_gam++) {
			$p0_gam = $p1_gam+($n_gam-$a_gam)*$p0_gam;
			$q0_gam = $q1_gam+($n_gam-$a_gam)*$q0_gam;
			$p_gam = $x_gam*$p0_gam+$n_gam*$p1_gam;
			$q_gam = $x_gam*$q0_gam+$n_gam*$q1_gam;
			if(abs($q_gam) != 0) {
				$s1_gam = $p_gam/$q_gam;
				$p1_gam = $p_gam;
				$q1_gam = $q_gam;
				if(abs(($s1_gam-$s_gam)/$s1_gam) < 1e-7) {
					$s_gam = $s1_gam*exp(-$x_gam)*$qq_gam/&gam1($a_gam);
					return (1-$s_gam);
				}
				$s_gam = $s1_gam;
			}
			$p1_gam = $p_gam;
			$q1_gam = $q_gam;
		}
	}
	print "a too large !\n";
	$s_gam = 1-$s_gam*exp(-$x_gam)*$qq_gam/&gam1($a_gam);
}

sub erf($) {
	my $x_erf = $_[0];
	my $y_erf;
	if($x_erf >= 0) {
		$y_erf=&gam2(0.5,$x_erf**2);
	} else {
		$y_erf=-&gam2(0.5,$x_erf);
	}
}
### Erf > ###


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

my ($Pairs,$Paired,$Singled,$Reads,$Alignment);
open LOG,'<',"$soapbase.log" or die "[x]Error opening $soapbase.log: $!\n";
while (<LOG>) {
	$Pairs = (split)[-2] if /^Total Pairs:/;
	$Paired = (split)[1] if /^Paired:/;
	$Singled = (split)[1] if /^Singled:/;
	$Reads = (split)[-1] if /^Total Reads/;
	$Alignment = (split)[1] if /^Alignment:/;
}
close LOG;

open NFO,'>',"$soapbase.nfo" or die "[x]Error opening $soapbase.nfo: $!\n";

	print NFO "#fmtS\tTotal_Pairs\tPaired\tSingled\n";
	print NFO "Summary\t",join("\t",$Pairs,$Paired,$Singled),"\n";
	#open(OUTS,"-|","gzip -dc $soapbase.soap.gz $soapbase.single.gz") or die "Error: $!\n";
	open(TP,"-|","gzip -dc $soapbase.soap.gz") or die "Error: $!\n";
	open(TS,"-|","gzip -dc $soapbase.single.gz") or die "Error: $!\n";

my ($BadLines,$BPOut,$ReadsOut,$MisSum,$TrimedBP,$TrimedReads,%Hit9r,%Hit9bp,%misMatch,%Indel)=(0,0,0,0,0,0);
my (%chrBPOut,%chrReadsOut,%chrMisSum,%chrTrimedBP,%chrTrimedReads,%chrHit9r,%chrHit9bp,%chrmisMatch,%chrIndel);
# for 46999 ChrIDs, one hash took 12m RES for {}=$ and 125m VIRT for {}{10}=$, thus fine to work with scaffolds. ( 1 gb VIRT max ? )
#my (@lines,$hit,$len,$chr,$types,$trim,$mistr,$missed);

sub getStatNfo($) {
	my ($hit,$len,$chr,$types,$trim,$mistr) = @{$_[0]}[3,5,7,9,-2,-1];
	my $missed=$mistr=~tr/ATCGatcg//;
	$BPOut += $len;	$chrBPOut{$chr} += $len;
	++$ReadsOut;	++$chrReadsOut{$chr};
	++$misMatch{$missed};
	++$chrmisMatch{$chr}{$missed};
	$MisSum += $missed;
	$chrMisSum{$chr} += $missed;

	$hit=4 if $hit>4;	# max to count 3, then >=4. Ancient Chinese wisdom, and a bit more ...
	++$Hit9r{$hit};
	++$chrHit9r{$chr}{$hit};
	$Hit9bp{$hit} += $len;
	$chrHit9bp{$chr}{$hit} += $len;
	if ($types > 200) {
		++$Indel{200-$types};
		++$chrIndel{$chr}{200-$types};
	} elsif ($types > 100) {
		++$Indel{$types-100};
		++$chrIndel{$chr}{$types-100};
	}

	my @Trimed = $trim =~ /(\d+)S/;
	if (@Trimed) {
		++$TrimedReads;
		++$chrTrimedReads{$chr};
		for ( @Trimed ) {
			$TrimedBP += $_;
			$chrTrimedBP{$chr} += $_;
		}
	}

}

my ($pairs,$lastpos,$line1,$line2,$pp,$pn,$calins,%insD,@lines)=(0);
while ($line1=<TP>) {
	last if eof TP;
	$lastpos=tell TP;
	@lines = split /\t/,$line1;
	if (@lines < 11) {
		++$BadLines;
		next;
	}
	&getStatNfo(\@lines);
	my ($id1, $n1, $len1, $f1, $chr1, $x1, $m1) = @lines[0,3,5,6,7,8,-2];
	$line2=<TP>;
	@lines = split /\t/,$line2;
	if (@lines < 11) {
		++$BadLines;
		next;
	}
	my ($id2, $n2, $len2, $f2, $chr2, $x2, $m2) = @lines[0,3,5,6,7,8,-2];
	#($soapid,$hit,$len,$strand,$chr,$pos,$trim) = (split(/\t/))[0,3,5,6,7,8,-2];
	$id1 =~ s/\/[12]$//;
	$id2 =~ s/\/[12]$//;
	if ($id1 ne $id2){	# single
		seek (TP, $lastpos, 0);
		next;
	}
	&getStatNfo(\@lines);
	next if $n1+$n2>2 or $chr1 ne $chr2;
	if ($f1 eq '+') {
		($pp,$pn)=($x1,$x2);
	} else {
		($pp,$pn)=($x2,$x1);
	}
=pod
	if ($ins > 1500) {	# FR => +.pos < -.pos; RF => -.pos < +.pos
		next if $pp < $pn;
	} else { next if $pp > $pn; }
=cut
	next if $pp > $pn;
	++$pairs;
	$line1=&getRealpos($len1, $f1, $x1, $m1);	# $len,$strand,$realpos,$trim
	$line2=&getRealpos($len2, $f2, $x2, $m2);	# Well, $line{1,2} is recycled.
	$calins=abs($line1-$line2);	# -.starting=0
	++$insD{$calins};
}
close TP;
while (<TS>) {
	@lines = split /\t/;
	if (@lines > 10) {	# soap2 output always more than 10 columes.
		&getStatNfo(\@lines);
	} else {	# rare event ...
		++$BadLines;
		next;
	}
}
close TS;

print NFO "\n#fmtC\tReadsOut\tBPOut\tMisSum\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\tBadLines\n";
print NFO join("\t",'ALL',$ReadsOut,$BPOut,$MisSum,$TrimedReads,$TrimedBP,
	${&combineC(\%misMatch)},${&combineJ(\%Hit9r)},${&combineJ(\%Hit9bp)},${&combineJ(\%Indel)},$BadLines),"\n\n";
print NFO "#fmtP\tReadsOut\tBPOut\tMisSum\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\n";
print NFO join("\t",";$_",$chrReadsOut{$_},$chrBPOut{$_},$chrMisSum{$_}||0,$chrTrimedReads{$_}||0,$chrTrimedBP{$_}||0,
	${&combineC(\%{$chrmisMatch{$_}})},${&combineJ(\%{$chrHit9r{$_}})},${&combineJ(\%{$chrHit9bp{$_}})},${&combineJ(\%{$chrIndel{$_}})}),"\n" for sort keys %chrReadsOut;
close NFO;

my ($n,$fqcmd,$tmpcmd,$avg,$std,$Lsd,$Rsd,$max_y,$max_x)=(0);
	my ($sum,$sum2,$v)=(0,0);
	open O,'>',"$soapbase.insD";
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

END:
#system("bzip2 -9 $opath$fqnames[0].unmap");
