#!/bin/env perl
use strict;
use warnings;
use Data::Dump qw(dump ddx);

unless (@ARGV){
	print "perl $0 <soaps.lst> <soaps.stat>\n";	# soaps.nfo can store the file size of soap. Too late, useless.
	exit;
}

my ($fqlst,$statout) = @ARGV;
my (%DATrbrf,%nfo);

open LST,'<',$fqlst or die "[x]Error opening $fqlst: $!\n";
while (<LST>) {
	chomp;
	my ($PESE,$sample,$lib,$FL,$ReadLen,$nfofpath)=split /\t/;
	$nfo{$sample}{$lib}{$FL}=[$nfofpath,$PESE,$ReadLen];
	$DATrbrf{Lane}{$sample}{$lib}{$FL}={ALL => [0,0,0,0,{},{},{},{},0], Summary => [0,0,0], };
	# "ReadsOut, BPOut, TrimedReads, TrimedBP, misMatchReads, Reads@Hit, BP@Hit, IndelReads, BadLines"
	# "Total_Reads(=Total_Pairs*2), Paired, Singled/Alignment"
	$DATrbrf{Lib}{$sample}{$lib}={ALL => [0,0,0,0,{},{},{},{},0], Summary => [0,0,0], } unless exists $DATrbrf{Lib}{$sample}{$lib};	# may be faster ?
	$DATrbrf{Sample}{$sample}={ALL => [0,0,0,0,{},{},{},{},0], Summary => [0,0,0], } unless exists $DATrbrf{Sample}{$sample};
}

sub combineC($) {
	my $href=$_[0];
	if ($href and %$href) {
		my (@str,$m);
		$m = (sort {$a<=>$b} keys %$href)[-1];
		for (1..$m) {
			push @str,$$href{$_}||0;
		}
		return \join('|',@str);
	} else {return \'.';}
}
sub combineLineA($) {
	my $ref=$_[0];
	return join ',',@{$$ref{Summary}},@{$$ref{ALL}}[0..3],${&combineC($$ref{ALL}[4])},${&combineC($$ref{ALL}[5])},${&combineC($$ref{ALL}[6])},${&combineC($$ref{ALL}[7])},$$ref{ALL}[8];
}
sub combineLine($) {
	my $ref=$_[0];
	return join ',',@$ref[0..3],${&combineC($$ref[4])},${&combineC($$ref[5])},${&combineC($$ref[6])},${&combineC($$ref[7])};
}
sub sumcsv ($$) {
	my ($href,$stref)=@_;
	return if $$stref eq '.';
	my %new = map {split /:/} split(/,/,$$stref);
	$$href{$_} += $new{$_} for (keys %new);
}
sub sumup ($$$) {
	my ($PESE,$sum,$item)=@_;
	my ($a,$m)=($$sum{Summary},$$sum{ALL});
	my ($b,$n)=($$item{Summary},$$item{ALL});
	if ($PESE eq 'PE') {
		$$a[0] += $$b[0]+$$b[0];
		$$a[1] += $$b[1];
		$$a[2] += $$b[2];
	} else {
		$$a[0] += $$b[0];
		$$a[2] += $$b[1];
	}
	$$m[-1] += $$n[-1];
	for my $chr (keys %$item) {
		next if $chr eq 'Summary';
		#$$sum{$_}=[0,0,0,0,{},{},{},{}] unless $$sum{$_};
		($m,$n)=($$sum{$chr},$$item{$chr});
		$$sum{$chr}=$m=[0,0,0,0,{},{},{},{}] unless $m;
		$$m[$_] += $$n[$_] for (0..3);
		&sumcsv($$m[$_],\$$n[$_]) for (4..7);
#warn "$chr\n";
	}
}

my %NFO;
for my $sample (sort keys %nfo) {
	for my $lib (keys %{$nfo{$sample}}) {
		for my $FL (keys %{$nfo{$sample}{$lib}}) {
			my ($nfofpath,$PESE,$ReadLen)=@{$nfo{$sample}{$lib}{$FL}};
			#print "[$_]\n";
			open NFO,'<',"$nfofpath" or (warn "[!]Error opening $nfofpath: $!\n" and next);
warn "$nfofpath\n";
			while (<NFO>) {
				next if /^(#|$)/;
				chomp;
				my ($key,@values)=split /\t/;
				$key='ALL' if $key eq '_ALL_';	# for old format
				$NFO{$key}=\@values;
			}
			close NFO;
			&sumup($PESE,$DATrbrf{Sample}{$sample},\%NFO);
			&sumup($PESE,$DATrbrf{Lib}{$sample}{$lib},\%NFO);
			&sumup($PESE,$DATrbrf{Lane}{$sample}{$lib}{$FL},\%NFO);
		}
	}
}
#ddx \%DATrbrf;
my %Chr;
open O,'>',$statout or die "[x]Error opening $statout: $!\n";
print O "#Summary\tSubItemOrder: Total_Reads,Aligned_Pairs,Aligned_Single,ReadsOut,BPOut,TrimedReads,TrimedBP,misMatchReads|ASC,Reads\@Hit|ASC,BP\@Hit|ASC,IndelReads|ASC,BadLines\n";
my ($Rsample,$Rlib,$RFL);
for my $sample (sort keys %nfo) {
	$Rsample=&combineLineA($DATrbrf{Sample}{$sample});#join ',',@{$DATrbrf{Sample}{$sample}{Summary}},@{$DATrbrf{Sample}{$sample}{ALL}}[0..3],${&combineC($DATrbrf{Sample}{$sample}{ALL}[4])},${&combineC($DATrbrf{Sample}{$sample}{ALL}[5])},${&combineC($DATrbrf{Sample}{$sample}{ALL}[6])},${&combineC($DATrbrf{Sample}{$sample}{ALL}[7])},$DATrbrf{Sample}{$sample}{ALL}[8];
	for (keys %{$DATrbrf{Sample}{$sample}}) {
		next if /(Summary|ALL)/;
		++$Chr{$_};
	}
	for my $lib (sort keys %{$nfo{$sample}}) {
		$Rlib=&combineLineA($DATrbrf{Lib}{$sample}{$lib});
		for my $FL (sort keys %{$nfo{$sample}{$lib}}) {
			$RFL=&combineLineA($DATrbrf{Lane}{$sample}{$lib}{$FL});
			print O join("\t",'ALL',$sample,$Rsample,$lib,$Rlib,$FL,$RFL),"\n";
		}
	}
}

print O "\n#ByChr\tSubItemOrder: ReadsOut,BPOut,TrimedReads,TrimedBP,misMatchReads|ASC,Reads\@Hit|ASC,BP\@Hit|ASC,IndelReads|ASC\n";
for my $Chr (sort keys %Chr) {
	for my $sample (sort keys %nfo) {
		$Rsample=&combineLine($DATrbrf{Sample}{$sample}{$Chr});
		for my $lib (sort keys %{$nfo{$sample}}) {
			$Rlib=&combineLine($DATrbrf{Lib}{$sample}{$lib}{$Chr});
			for my $FL (sort keys %{$nfo{$sample}{$lib}}) {
				$RFL=&combineLine($DATrbrf{Lane}{$sample}{$lib}{$FL}{$Chr});
				my $ChrID=$Chr;
				$ChrID =~ s/^;//;
				print O join("\t",$ChrID,$sample,$Rsample,$lib,$Rlib,$FL,$RFL),"\n";
			}
		}
	}
}
close O;
