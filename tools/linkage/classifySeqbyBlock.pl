#!/usr/bin/perl -w
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
use Time::HiRes qw ( gettimeofday tv_interval );
use GalaxyXS::ChromByte 1.02;# ':all';
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;

our $opts='i:o:v:r:l:c:p:bd';
our($opt_i, $opt_o, $opt_v, $opt_b, $opt_d, $opt_l, $opt_c, $opt_p);

our $help=<<EOH;
\t-i soaps.lst (2soap/soaps.lst)
\t-l block file lst (./block.lst)
\t-c Chromosome NFO file (chr.nfo) in format: /^ChrID\\s+ChrLen\\s?.*\$/
\t-o Output path (./cfq) for [Sample/F_L_Lib.dat.gz], will mkdir if not exist
\t-p Parallel id (0)
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='./2soap/soaps.lst' if ! $opt_i;
die "[x]-i $opt_i not exists !\n" unless -f $opt_i;
$opt_l='./block.lst' if ! $opt_l;
$opt_o='./ffq' if ! $opt_o;
$opt_c='chr.nfo' if ! $opt_c;

no warnings;
$opt_v=int $opt_v;
$opt_p=int $opt_p;
use warnings;

$opt_o=~s#/+$##;

print STDERR "From [$opt_i],$opt_p [$opt_l] to [$opt_o] refer to [$opt_c]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (%ChrMem,%ChrLen);
open ChrNFO,'<',$opt_c or die "[x]Error opening $opt_c: $!\n";
while (<ChrNFO>) {
	next if /^#/;
	my ($chrid,$len)=split /\t/;
	$ChrLen{$chrid}=$len;
}
close ChrNFO;

system('mkdir','-p',$opt_o);

open LST,'<',$opt_l or die "[x]Error opening $opt_l: $!\n";
my (%BlockD,%BlockS,%BlockDO);
while (<LST>) {
	chomp;
	my ($chrid,$block)=split /\t/;
	open B,'<',$block or die "[x]Error opening $block: $!\n";
	print STDERR "[!] Block $chrid: ";
	chomp($_=<B>);
	s/^#//;
	my @Sample=split /\s+/;
	system('mkdir','-p',$opt_o.'/'.$_) for @Sample;
	$BlockS{$chrid}=\@Sample;
	while (<B>) {
		chomp;
		my ($pos,@Dat)=split /\s+/;
		next if $pos !~ /\d+/;
		$BlockD{$chrid}{$pos}=join('',@Dat);
	}
	print STDERR ".\b";
	my @PosOrder=sort {$a<=>$b} keys %{$BlockD{$chrid}};
	$BlockDO{$chrid}=\@PosOrder;
	$ChrMem{$chrid}=initchr($ChrLen{$chrid});
	warn "$ChrLen{$chrid}.\n";
}
close LST;
my $fileNFO;

my %varToGT=(
	0   => 1,
	1   => 2,
	'-' => 0,
	0.5 => 3,
);
sub readGTtoMEM($) {
	my ($sample)=@_;
	for my $chrid (keys %BlockS) {
		my $sp=0;	# Starts from 0
		for (@{$BlockS{$chrid}}) {
			last if $_ eq $sample;
			++$sp;
		}
		setbases($ChrMem{$chrid},0,$ChrLen{$chrid},0);	# Zero it
		my $lastPos=1;
		for my $pos (@{$BlockDO{$chrid}}) {
			my $CurGT=substr($BlockD{$chrid}{$pos},$sp,1);
			my $v=0;
			$v=$varToGT{$CurGT} if exists $varToGT{$CurGT};
			setbases($ChrMem{$chrid},$lastPos,$pos,$v);
			$lastPos=$pos;
		}
	}
}
sub getGT($$) {
	my ($chrid,$pos)=@_;
	return 0 unless exists $ChrMem{$chrid};
	return 0 if $pos > $ChrLen{$chrid} or $pos <= 0;
	return getbase($ChrMem{$chrid},$pos);
}

sub getFH($) {
	my ($file)=@_;
	my $FH;
	my ($PESE,$sample,$lib,$FL)=@$fileNFO;
	unless (-s $file) {
		open $FH,'-|',"gzip -dc $file.gz" or die "[x]Error opening $PESE [$file.gz] with gzip: $!\n";
		print STDERR '[!]',join(', ',$sample,$lib,$FL,"$file.gz"),"\n";
	} else {
		open $FH,'<',$file or die "[x]Error opening $PESE [$file]: $!\n";
		print STDERR '[!]',join(', ',$sample,$lib,$FL,$file),"\n";
	}
	return $FH;
}

sub readLine($) {
	my ($FDH)=@_;
	my ($FH,$Dat,$CurtLine,$NextLine)=@{$FDH};
	my ($line,$ID, $Chr, $Pos, $len);
	$CurtLine=$NextLine;
	if ($line=<$FH>) {
		($ID, $Chr, $Pos, $len)=(split /\t/,$line)[0,7,8,5];
		$NextLine=[$ID, $Chr, $Pos, $len];
	} else {
		$NextLine=[-2,'',0,0];
	}
	$FDH->[2] = $CurtLine;
	$FDH->[3] = $NextLine;
}
sub nextPair($) {
	my ($FDH)=@_;
	my ($FH,$Dat,$CurtLine,$NextLine)=@{$FDH};
	if ( $NextLine->[0] == -2 or $CurtLine->[0] < $NextLine->[0] ) {
		unless ($NextLine->[0] == -2) {
			$Dat=[$CurtLine->[0],$CurtLine->[1],$CurtLine->[2],$CurtLine->[2]+$CurtLine->[3]-1,'',0,0];
			&readLine($FDH);
		} else {
			$Dat=[-2,'',0,0,'',0,0];
		}
	} elsif ( $CurtLine->[0] == $NextLine->[0] ) {
		$Dat=[$CurtLine->[0],$CurtLine->[1],$CurtLine->[2],$CurtLine->[2]+$CurtLine->[3]-1,$NextLine->[1],$NextLine->[2],$NextLine->[2]+$NextLine->[3]-1];
		&readLine($FDH);
		&readLine($FDH);
	} else { die "[x]Not a raw SOAP2 file !"; }
	$FDH->[1] = $Dat;
}

sub getSeq($$) {
	my ($FDH,$lastID)=@_;
	my @nFDH;
	FILES: for my $fdh (@$FDH) {
		while ($fdh->[1]->[0] <= $lastID) {
			if ($fdh->[1]->[0] == -2) {
				#$fdh=undef;
				next FILES;
			}
			&nextPair($fdh);	# Sync
		}
		push @nFDH,$fdh;
	}
	$_[0]=\@nFDH;
	@nFDH = sort { $a->[1]->[0] <=> $b->[1]->[0] } @nFDH;	# ASC
	if (@nFDH) {
		return $nFDH[0]->[1];
	} else {
		return [-2,'',0,0,'',0,0];
	}
}

open LST,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
my @LSTdat=<LST>;
close LST;
if ($opt_p>0) {
	--$opt_p;
	@LSTdat=($LSTdat[$opt_p]);
}
for (@LSTdat) {
	chomp;
	my ($PESE,$sample,$lib,$FL,$ReadLen,$nfofpath)=split /\t/;
	readGTtoMEM($sample);
	$fileNFO=[$PESE,$sample,$lib,$FL];
	$nfofpath =~ s/\.nfo$//;
	my @Files;
	if ($PESE eq 'PE') {
		@Files=("$nfofpath.soap","$nfofpath.single");
	} else {
		@Files=("$nfofpath.soap");
		next;
	}
	my $shfile="$nfofpath.archive";
	if (-s $shfile) {
		open SH,'<',$shfile or die "[x]Error opening $shfile: $!\n";
	} elsif (-s "$shfile.0") {
		open SH,'<',"$shfile.0" or die "[x]Error opening $shfile.0: $!\n";
	} else {
		warn "[!]$shfile not found, skipped [$PESE,$sample,$lib,$FL].\n";
		next;
	}
	open SH,'<',"$nfofpath.archive" or die "[x]Error opening $nfofpath.archive: $!\n";
	my %FQ;
	while (<SH>) {
		next unless / -[ab] /;
		while (/ -([ab]) (\S+)\b/g) {
			$FQ{$1}=$2 if $2;
		}
	}
	open OUT,'|-',"gzip -9c - >$opt_o/$sample/${FL}_$lib.dat.gz" or die "[x]Error opening [$opt_o/$sample/${FL}_$lib.dat.gz] with gzip: $!\n";
	print OUT "#$FQ{'a'},$FQ{'b'}\n#$PESE,$sample,$lib,$FL,$ReadLen,$nfofpath\n";
	my $tdat=[-1];
	my $tbuf=[-1,'',0,0];
	my @FHD;
	push @FHD,[&getFH($_),$tdat,$tbuf,$tbuf] for @Files;
	# Handle, [id,chr1,pos1,chr2,pos2],bufferCur[id,chr,pos,ab],bufferNext[id,chr,pos,ab]
	my $lastID=-1;
	while ($lastID > -2) {
		my ($id,$chr1,$pos1,$end1,$chr2,$pos2,$end2)=@{&getSeq(\@FHD,$lastID)};
		my $GTa=&getGT($chr1,$pos1);
		my $GTb=&getGT($chr2,$pos2);
		my $GTae=&getGT($chr1,$end1);
		my $GTbe=&getGT($chr2,$end2);
		print OUT join("\t",$id,$GTa|$GTb|$GTae|$GTbe,"$GTa,$GTae;$GTb,$GTbe","$chr1,$pos1-$end1; $chr2,$pos2-$end2"),"\n" if $id >=0;
		$lastID=$id;
	}
	close OUT;
	warn "[!]$opt_o/$sample/${FL}_$lib.dat done.\n",'-' x 75,"\n";
}
#close LST;

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
$ ./classifySeqbyBlock.pl -bi ./to9311g10a/2soap/soaps.lst

less -S ffq/BI10-1/FC704U5AAXX_L5_ORYqzpRALDIAAPEI-4.dat
gzip -dc ./to9311g10a/2soap/BI10-1/ORYqzpRALDIAAPEI-4/100506_I328_FC704U5AAXX_L5_ORYqzpRALDIAAPEI-4_1.soap.gz | awk '{print $1"\t"$8"\t"$9"\t"$5}' |less -S
gzip -dc ./to9311g10a/2soap/BI10-1/ORYqzpRALDIAAPEI-4/100506_I328_FC704U5AAXX_L5_ORYqzpRALDIAAPEI-4_1.single.gz | awk '{print $1"\t"$8"\t"$9"\t"$5}' |less -S

PR  NI  VIRT  RES  SHR S %CPU %MEM    TIME+  COMMAND
16   0  681m 605m 1760 S    1  3.8   3:17.54 classifySeqbyBl

#!/bin/sh
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#$ -cwd -r y -l vf=950m
#$ -o /dev/null -e /dev/null
#$ -S /bin/bash -t 1-195
perl ./classifySeqbyBlock.pl -bi ./to9311g10a/2soap/soaps.lst -p $SGE_TASK_ID -o ./ffq2  2> ./log/classifySeqbyBlock.$SGE_TASK_ID.log
