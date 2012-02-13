#!/bin/env perl
use strict;
use warnings;

unless (@ARGV>1){
	print "perl $0 soap_stderr soap [single]\n";
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

my $soaplog = shift @ARGV;

my ($Pairs,$Paired,$Singled,$Reads,$Alignment);

open LOG,'<',"$soaplog" or die "[x]Error opening $soaplog: $!\n";
while (<LOG>) {
	$Pairs = (split)[-2] if /^Total Pairs:/;
	$Paired = (split)[1] if /^Paired:/;
	$Singled = (split)[1] if /^Singled:/;
	$Reads = (split)[-1] if /^Total Reads/;
	$Alignment = (split)[1] if /^Alignment:/;
}
close LOG;

my ($BadLines,$BPOut,$ReadsOut,$MisSum,$TrimedBP,$TrimedReads,%Hit9r,%Hit9bp,%misMatch,%Indel)=(0,0,0,0,0,0);
my (%chrBPOut,%chrReadsOut,%chrMisSum,%chrTrimedBP,%chrTrimedReads,%chrHit9r,%chrHit9bp,%chrmisMatch,%chrIndel);
# for 46999 ChrIDs, one hash took 12m RES for {}=$ and 125m VIRT for {}{10}=$, thus fine to work with scaffolds. ( 1 gb VIRT max ? )
my (@lines,$hit,$len,$chr,$types,$trim,$mistr,$missed);

my @files=@ARGV;
while (my $fq = shift) {
	if ($fq =~ /\.gz$/) {
		open( OUTS,'-|',"gzip -dc $fq") or die "[x]Error on $fq: $!\n";
	} elsif ($fq =~ /\.bz2$/) {
		open( OUTS,'-|',"bzip2 -dc $fq") or die "[x]Error on $fq: $!\n";
	} else {open( OUTS,'<',$fq) or die "[x]Error on $fq: $!\n";}

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
}
print STDERR "\n#fmtC\tReadsOut\tBPOut\tMisSum\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\tBadLines\n";
print STDERR join("\t",'ALL',$ReadsOut,$BPOut,$MisSum,$TrimedReads,$TrimedBP,
	${&combineC(\%misMatch)},${&combineJ(\%Hit9r)},${&combineJ(\%Hit9bp)},${&combineJ(\%Indel)},$BadLines),"\n\n";
print STDERR "#fmtP\tReadsOut\tBPOut\tMisSum\tTrimedReads\tTrimedBP\tmisMatchReads\tReads\@Hit\tBP\@Hit\tIndelReads\n";
print STDERR join("\t",";$_",$chrReadsOut{$_},$chrBPOut{$_},$chrMisSum{$_}||0,$chrTrimedReads{$_}||0,$chrTrimedBP{$_}||0,
	${&combineC(\%{$chrmisMatch{$_}})},${&combineJ(\%{$chrHit9r{$_}})},${&combineJ(\%{$chrHit9bp{$_}})},${&combineJ(\%{$chrIndel{$_}})}),"\n" for sort keys %chrReadsOut;
