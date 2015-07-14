#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );

use Data::Dump qw(ddx);

die "Usage: $0 <Ref> <Virus> <Outprefix>\n" if @ARGV <3;

my ($Reff,$Virf,$outp)=@ARGV;

$Reff='/share/users/huxs/git/toGit/perl/readsim/Chr1.fa';
$Virf='/share/users/huxs/work/bsvir/HBV.AJ507799.2.fa';

my $SampleCnt = 100;
my $Depth = 50;
my $PEinsertLen=200;
my $SeqReadLen=90;
my $VirFragMax = 500;
my $VirFragMin = 20;
my $RefBorder = $PEinsertLen + 1000;
my $RefNratioMax = 0.02;

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.xz$/) {
		open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

sub getRefChr1st($) {
	my $GENOME = $_[0];
	while (<$GENOME>) {
		s/^>//;
		/^(\S+)/ or next;
		my $seqname = $1;
		print STDERR ">$seqname ...";
		$/=">";
		my $genome=<$GENOME>;
		chomp $genome;
		$genome=~s/\s//g;
		$/="\n";
		my $thelength = length $genome;
		print STDERR "\b\b\b", $thelength, ".\n";
		return $genome;
	}
}

my $Refh = openfile($Reff);
my $Refstr = getRefChr1st($Refh);
my $RefLen = length $Refstr;
close $Refh;
my $Virfh = openfile($Virf);
my $Virstr = getRefChr1st($Virfh);
my $VirLen = length $Virstr;
close $Virfh;
$Virstr .= $Virstr;	# circle


my (@Refticks,@Virticks);
sub getticks($$$) {
	my ($RefBorder,$Refstr,$RefLen) = @_;
	my @theticks;
	while (@theticks < 100) {
		my $pos0 = int(rand($RefLen-(2*$RefBorder)))+$RefBorder;
		my $str0 = substr $Refstr,($pos0-$PEinsertLen),2*$PEinsertLen;
		my $seq = $str0;
		my $N = $seq=~tr/Nn//;
		next if $N > 2*$PEinsertLen*$RefNratioMax;
		push @theticks,$pos0
	}
	@theticks = sort {$a<=>$b} @theticks;
	return \@theticks;
}

@Refticks = @{getticks($RefBorder,$Refstr,$RefLen)};
@Virticks = @{getticks($VirFragMax,$Virstr,$VirLen)};
#ddx \@Refticks,\@Virticks;

sub revcom($) {
    my $str = $_[0];
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    my $rev = reverse $str;
    $rev    =~ tr/[](){}<>/][)(}{></;
    return $rev;
}
open O,'>',"$outp.Ref.fa";
open R1,'>',"$outp.1.fq";
open R2,'>',"$outp.2.fq";
for my $pRef (@Refticks) {
	my $seqR1 = substr $Refstr,($pRef-$PEinsertLen),$PEinsertLen;
	my $seqR2 = substr $Refstr,$pRef,$PEinsertLen;
	my $pVir = shift @Virticks;
	my $LenV = int(rand($VirFragMax-$VirFragMin))+$VirFragMin;
	my $startV = $pVir-int(0.5*$LenV);
	my $seqV = substr $Virstr,$startV,$LenV;
	my $isReverse = int(rand(2));
	my $strand = '+';
	if ($isReverse) {
		$seqV = revcom($seqV);
		$strand = '-';
	}
	my $newSeq = join('',uc $seqR1,lc $seqV,uc $seqR2);
	my $tID = join('_','Ref',$pRef-$PEinsertLen,$pRef,$pRef+$PEinsertLen,'Vir',$strand,$startV,$startV+$LenV);
	print O '>',$tID,"\n$newSeq\n\n";
	my $maxP = length($newSeq) - $PEinsertLen;
	my $step = int(2*$SeqReadLen/$Depth) or die "2* $SeqReadLen /$Depth => 0\n";
	my $p = 0;
	while ($p <= $maxP) {
		my $PE = substr $newSeq,$p,$PEinsertLen;
		my $R1 = substr $PE,0,$SeqReadLen;
		my $R2 = substr $PE,$PEinsertLen-$SeqReadLen,$SeqReadLen;
		#my $revR1 = revcom($R1);
		my $revR2 = revcom($R2);
		my $Qual = 'e' x $SeqReadLen;
		print R1 "\@sf${p}_${tID}/1\n$R1\n+\n$Qual\n";
		print R2 "\@sf${p}_${tID}/2\n$revR2\n+\n$Qual\n";
		#print R1 "\@sr${p}_${tID}/1\n$revR1\n+\n$Qual\n";
		#print R2 "\@sr${p}_${tID}/2\n$R2\n+\n$Qual\n";
		$p += $step;
	}
}
close O;
close R1; close R2;

__END__
./virusinserts.pl /share/users/huxs/git/toGit/perl/readsim/Chr1.fa /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa test
