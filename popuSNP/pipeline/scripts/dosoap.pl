#!/bin/env perl
use strict;
use warnings;

my $bin='/nas/RD_09C/resequencing/soft/bin/soap/soap2.20';
my $arg0='-p 6 -t -s 40 -l 32';

unless (@ARGV>6){
	print "perl $0 PESE ReadLen,insMIN,insMAX gap outpath/ ref.index fqext,path/ fqname[,fqname2] [u]\n";
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

my ($PESE,$ins,$G,$opath,$Ref,$fqextpath,$fqname,$u) = @ARGV;
#print '[',join('] [',@ARGV),"]\n";
#print '[',join('] [',$PESE,$ins,$G,$Ref,$fqextpath,@fqnames),"]\n";
my ($readlen,$min,$max)=split ',',$ins;
my ($ext,$path)=split ',',$fqextpath;
my @fqnames=split ',',$fqname;
my $mismatch=$readlen>70?3:1;
$mismatch = 5 if $readlen >= 90;
$arg0 .= ' -R' if $min > 1500;
my $fqcmd;
unless ($max==0) {	# PE
	$fqcmd="-a $path$fqnames[0]$ext -b $path$fqnames[-1]$ext -o $opath$fqnames[0].soap -2 $opath$fqnames[0].single";
} else {	# SE
	$fqcmd="-a $path$fqnames[0]$ext  -o $opath$fqnames[0].se";
}
$fqcmd .= " -u $opath$fqnames[0].unmap" if $G or $u && $u eq 'u';
# always put -abDo2 in the first for the poor case without break ...
my $sh="$bin $fqcmd -D $Ref $arg0 -m $min -x $max -v $mismatch -g $G 2>$opath$fqnames[0].log";
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
	print NFO "#fmtS\tTotal_Pairs\tPaired\tSingled\n";
	print NFO "Summary\t",join("\t",$Pairs,$Paired,$Singled),"\n";
	@ARGV=("$opath$fqnames[0].soap", "$opath$fqnames[0].single");
}
my ($BadLines,$BPOut,$ReadsOut,$TrimedBP,$TrimedReads,%Hit9r,%Hit9bp,%misMatch,%Indel)=(0,0,0,0,0);
my (%chrBPOut,%chrReadsOut,%chrTrimedBP,%chrTrimedReads,%chrHit9r,%chrHit9bp,%chrmisMatch,%chrIndel);
# for 46999 ChrIDs, one hash took 12m RES for {}=$ and 125m VIRT for {}{10}=$, thus fine to work with scaffolds. ( 1 gb VIRT max ? )
my (@lines,$hit,$len,$chr,$types,$trim);
while (<>) {
	@lines = split /\t/;
	if (@lines > 10) {	# soap2 output always more than 10 columes.
		($hit,$len,$chr,$types,$trim) = @lines[3,5,7,9,-2];
		$BPOut += $len;	$chrBPOut{$chr} += $len;
		++$ReadsOut;	++$chrReadsOut{$chr};

		$hit=10 if $hit>10;	# max to count 9, then >=10
		++$Hit9r{$hit};
		++$chrHit9r{$chr}{$hit};
		$Hit9bp{$hit} += $len;
		$chrHit9bp{$chr}{$hit} += $len;

		if ($types < 100) {
			++$misMatch{$types};
			++$chrmisMatch{$chr}{$types};
		} elsif ($types < 200) {	# '3S33M9D39M', '32M1D14M29S' exists
			++$Indel{$types-100};
			++$chrIndel{$chr}{$types-100};
		} else {
			++$Indel{200-$types};
			++$chrIndel{$chr}{200-$types};
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
print NFO "\n#fmtC\tReadsOut\tBPOut\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\tBadLines\n";
print NFO join("\t",'ALL',$ReadsOut,$BPOut,$TrimedReads,$TrimedBP,
	${&combineC(\%misMatch)},${&combineJ(\%Hit9r)},${&combineJ(\%Hit9bp)},${&combineJ(\%Indel)},$BadLines),"\n\n";
print NFO "#fmtP\tReadsOut\tBPOut\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\n";
print NFO join("\t",";$_",$chrReadsOut{$_},$chrBPOut{$_},$chrTrimedReads{$_},$chrTrimedBP{$_},
	${&combineC(\%{$chrmisMatch{$_}})},${&combineJ(\%{$chrHit9r{$_}})},${&combineJ(\%{$chrHit9bp{$_}})},${&combineJ(\%{$chrIndel{$_}})}),"\n" for sort keys %chrReadsOut;
close NFO;

END:
