#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
use Galaxy::ShowHelp;
use Time::HiRes qw ( gettimeofday tv_interval );
use Data::Dump qw(ddx);

$main::VERSION=0.0.1;

our $opts='i:o:n:g:v:r:c:l:p:bqd';
our($opt_i, $opt_n, $opt_g, $opt_r, $opt_c, $opt_l, $opt_o, $opt_v, $opt_b, $opt_q, $opt_d, $opt_h, $opt_p);

our $help=<<EOH;
\t-i Blast filted list (./markerblastf.lst)
\t-l linkage map list (./linkagemap.lst)
\t-p marker dat (./markerpospa64.dat)
\t-n N zone file (Pa64.Nzone)
\t-o Output file
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='./markerblastf.lst' if ! $opt_i;
$opt_p='./markerpospa64.dat' if ! $opt_p;
$opt_l='./linkagemap.lst' if ! $opt_l;
$opt_n='Pa64.Nzone' if ! $opt_n;
$opt_v=0 if ! $opt_v;

print STDERR "From [$opt_i] to [$opt_o] refer to [$opt_n][$opt_l][$opt_p]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my %Scaffords;	# \%( $chr->[$sPosRel,$scM,$sWeight] )
my %LinkageMap;

open L,'<',$opt_l or die "Error: $!\n";
while (<L>) {
	chomp;
	my ($ChrID,$File)=split /\t/;
	open LM,'<',$File or die "Error: $!\n";
	while (<LM>) {
		next if /^#/;
		my ($Chr,$Pos,$cM)=split /\t/;
		die "[x][$ChrID] not match [$File] !\n" if $Chr ne $ChrID;
		$LinkageMap{$Chr}{$Pos}=$cM;
	}
}
close L;

my (%cMPos,%cMlst);
open L,'<',$opt_p or die "Error: $!\n";
while (<L>) {
	next if /^#/;
	chomp;
	my ($Marker,$Chr,$cM,$Pos,$strand,$weight)=split /\t/;
	$cMPos{$Chr}{$cM}=$Pos;
}
close L;
for my $chr (keys %cMPos) {
	$cMlst{$chr} = [sort {$a <=> $b} keys %{$cMPos{$chr}}];
}
my $Err=50/135;	# +- 100 * 0.5/135
sub cMtoPos($$) {
	my ($Chr,$cM,$Strand)=@_;
	($Chr,$Strand)=split /\t/,$Chr;
	my ($thecMlower,$thecMUpper)=(0,0);
	if ($cM > $cMlst{$Chr}->[0]) {
		for (@{$cMlst{$Chr}}) {
			$thecMUpper=$_;
			last if $_ > $cM;
			next if $_ == $cM;
			$thecMlower=$thecMUpper;
		}
	} else {
		$thecMUpper=$cMlst{$Chr}->[1];
	}
	my ($thePoslower,$thePosUpper)=($cMPos{$Chr}{$thecMlower},$cMPos{$Chr}{$thecMUpper});
	die "[x]Pos Equal @ $Chr,$cM,$thecMlower,$thecMUpper" if $thePosUpper==$thePoslower;
	my $BPtocM=($thePosUpper-$thePoslower)/($thecMUpper-$thecMlower);
	die "[x]Pos Order Error @ $Chr,$thecMlower,$thecMUpper" if $BPtocM <= 0;	# So, deal with unorder mannually ...
	my $Pos=($cM-$thecMlower)*$BPtocM+$thePoslower;
	my $ErrR=$Err*$BPtocM;
	my @R=($Pos,$ErrR,$thePoslower,$thePosUpper);
	$_=int(0.5+10*$_)/10 for @R;
	return [@R,$thecMlower,$thecMUpper];
}

sub getRel($$$$$) {
	my ($Qid,$Qs,$Qe,$Ss,$Se,$BTOP)=@_;
	my ($mChr,$mPos,$mSiderLen)=@{&splitMarkerid($Qid)};
	my $cM=$LinkageMap{$mChr}{$mPos};
	unless (defined $cM) {
		warn "[!]Cannot find Marker $Qid @ $mChr,$mPos !\n";
		$cM=0;
	}
	my $LeftBPalongQ=$mSiderLen-$Qs+1;
	my $WalkingOn=$LeftBPalongQ;
	my @btop0=split /(\D+)/,$BTOP;
	# $ perl -le '$a="45YT9-Ac-c-c-11TC4";@b=split /(\D+)/,$a;print "[",join("|",@b),"]"'
	# [45|YT|9|-Ac-c-c-|11|TC|4]
	my @btop=();
	for (@btop0) {
		if (/\d+/) {
			push @btop,$_;
		} else {
			my @bin=split /([\w-]{2})/;
			# $ perl -le '$a="-Ac-c-c-";@b=split /([\w-]{2})/,$a,0;print "[",join("|",@b),"]"'
			# [|-A||c-||c-||c-]
			for (@bin) {
				next unless $_;
				push @btop,$_;
			}
		}
	}
	my $LeftBPalongS=0;
	for (@btop) {
		if ($WalkingOn <= 0) {
			print STDERR "$Qid [$_] " if $opt_v and $WalkingOn == 0;
			last;
		}
		print STDERR "-$Qid [$_]-" if $opt_v>1;
		if (/\d/) {
			$LeftBPalongS += $_;
			$WalkingOn -= $_;
		} else {
			my @op=split //;
			if (/-/) {
				--$WalkingOn if $op[1] eq '-' and $op[0] ne '-';
				++$LeftBPalongS if $op[0] eq '-' and $op[1] ne '-';
			} else {
				--$WalkingOn;
				++$LeftBPalongS;
			}
		}
		print STDERR "-$WalkingOn $LeftBPalongS\n" if $opt_v>1;
		my $NOP;
	}
	my $strand;
	if ($Ss < $Se) {
		$strand=1;
		$mChr .= "\t+";
	} else {
		$strand=-1;
		$mChr .= "\t-";
	}
	my $Spos=$Ss + $strand*$LeftBPalongS;
	warn "$mChr,$Spos,$cM\n" if $opt_v;
	return [$mChr,$Spos,$cM];
}
sub splitMarkerid($) {
	my $MarkerID=$_[0];
	my ($mChr,$mPos,$mSiderLen)=split /[_m]/,$MarkerID,3;
	return [$mChr,$mPos,$mSiderLen];
}
sub pushScaf($) {
	my ($Qid,$Sid,$Pidentity,$AlnLen,$identical,$mismatches,$Gap,$Qs,$Qe,$Ss,$Se,$E,$bitScore,$BTOP,$Hit)=@{$_[0]};
	my ($sPosRel,$scM,$sWeight,$Count)=(0,0,0,0);
	my ($mChr,$mPos,$mSiderLen)=@{&splitMarkerid($Qid)};
	if ($Ss < $Se) {
		$mChr .= "\t+";
	} else {
		$mChr .= "\t-";
	}
	my $weight=-log($E)/$Hit;
	my %Dat;
	if (exists $Scaffords{$Sid}) {
		%Dat=%{$Scaffords{$Sid}};
		if (exists $Scaffords{$Sid}{$mChr}) {
			($sPosRel,$scM,$sWeight,$Count)=@{$Dat{$mChr}};
		}
	} else {
		%Dat=($mChr=>[$sPosRel,$scM,$sWeight,$Count]);
	}
	my ($chr,$pos,$cM)=@{&getRel($Qid,$Qs,$Qe,$Ss,$Se,$BTOP)};
	$sPosRel += $pos*$weight;
	$scM += $cM*$weight;
	$sWeight += $weight;
	++$Count;
	$Dat{$mChr}=[$sPosRel,$scM,$sWeight,$Count];
	$Scaffords{$Sid}=\%Dat;
}

open L,'<',$opt_i or die "Error: $!\n";
while (<L>) {
	chomp;
	my ($ChrID,$File)=split /\t/;
	open IN,'<',$File or die "Error: $!\n";
	while (<IN>) {
		next if /^#/;
	# Fields: query id, subject id, % identity, alignment length, identical, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, BTOP
	#-outfmt '6 qseqid sseqid pident length nident mismatch gapopen qstart qend sstart send evalue bitscore btop'
	#Chr01_457584m45 Scaffold011460  95.45   88      87      4       0       1       88      1020    933     1e-34    147    16CA2WA25MA6WA35        1
		chomp;
		my @Dat=split /\t/;
		my ($Qid,$Sid,$Pidentity,$AlnLen,$identical,$mismatches,$Gap,$Qs,$Qe,$Ss,$Se,$E,$bitScore,$BTOP,$Hit)=@Dat;
		next if $Sid =~ /^chr/i;	# Just skip ...
		&pushScaf(\@Dat);
	}
	close IN;
}
ddx \%Scaffords if $opt_v>1;
open O,'>',$opt_o or die "Error: $!\n";
print O "#ScaffordID\tScaffordAnchorOffect\tChrAnchored\tStrand\tcM\tPosT\tPosErr\tWeight,Count,thecMlower,thecMUpper,thePoslower,thePosUpper\n";
for my $ScaffID (keys %Scaffords) {
	my @dat=();	# [key,value(weight)]
	for (keys %{$Scaffords{$ScaffID}}) {
		push @dat,[$_,${$Scaffords{$ScaffID}{$_}}[2]];
	}
	@dat = sort {$b->[1] <=> $a->[1]} @dat;
	my ($sPosRel,$scM,$sWeight,$Count)=@{$Scaffords{$ScaffID}{$dat[0][0]}};
	$sPosRel=int(0.5+10*$sPosRel/$sWeight)/10;
	$scM=int(0.5+1000*$scM/$sWeight)/1000;
	$Count=-$Count;	# just flag
	my ($pos,$err,$thePoslower,$thePosUpper,$thecMlower,$thecMUpper)=@{&cMtoPos($dat[0][0],$scM)};
	$Scaffords{$ScaffID}=[$dat[0][0],$sPosRel,$scM,$sWeight,$Count];
	print O join("\t",$ScaffID,$sPosRel,$dat[0][0],$scM,$pos,$err,"$sWeight,$Count,$thecMlower,$thecMUpper,$thePoslower,$thePosUpper"),"\n";
}
close O;
ddx \%Scaffords if $opt_v>1;

__END__
grep -h caff ./out/f3545Chr*.pa64 > f3545ChrScaff.pa64
./relocate.pl -bi f3545ChrScaff.pa64 -o f3545ChrScaff.pos
