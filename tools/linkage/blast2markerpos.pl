#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
#use GalaxyXS::ChromByte;
use Data::Dump qw(ddx);

unless (@ARGV > 0) {
    print "perl $0 <marker linkage> <marker blast_m6c rd> <Out_prefix>\n";
    exit 0;
}

my ($markerf,$blastf,$out)=@ARGV;
my $opt_v=0;
my (%MarkercM,%Scaffords);
#open C,'<',$chrnfof or die "Error:[$chrnfof] $!\n";
#while (<C>) {
#	next if /^#/;
#	my ($chrid,$len,$efflen)=split /\t/;
#	$ChrLen{$chrid}=$len;
#}
#close C;
open C,'<',$markerf or die "Error:[$markerf] $!\n";
while (<C>) {
	next if /^#/;
	my ($chrid,$pos,$cm)=split /\t/;
	$MarkercM{$chrid}{$pos}=$cm;
}
close C;

sub getcM($$) {
	my ($chrid,$pos)=@_;
	if (exists $MarkercM{$chrid}{$pos}) {
		return $MarkercM{$chrid}{$pos};
	} else {
		return "-";	# too slow ...
		my @Pos=sort {$a <=> $b} keys %{$MarkercM{$chrid}};
		my $i;
		for ($i=0;$i<=$#Pos;$i++) {
			last if $Pos[$i] > $pos;
		}
		my ($left,$right)=($Pos[$i-1],$Pos[$i]);
		if ($pos-$left > $right-$pos) {
			return $MarkercM{$chrid}{$right};
		} else {
			return $MarkercM{$chrid}{$left};
		}
	}
}
sub splitMarkerid($) {
	my $MarkerID=$_[0];
	my ($mChr,$mPos,$mSiderLen)=split /[_m]/,$MarkerID,3;
	my $mcM=&getcM($mChr,$mPos);
	return [$mChr,$mPos,$mSiderLen,$mcM];
}
sub GetRelPos(@) {
	my ($Qid,$Sid,$Pidentity,$AlnLen,$identical,$mismatches,$Gap,$Qs,$Qe,$Ss,$Se,$E,$bitScore,$BTOP,$Hit)=@_;
	my ($mChr,$mPos,$mSiderLen,$mcM)=@{&splitMarkerid($Qid)};
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
	} else {
		$strand=-1;
	}
	my $Spos=$Ss + $strand*$LeftBPalongS;
	warn "$mChr,$Spos,$strand\n" if $opt_v;
	return [$Spos,$strand,$mcM];
}

my ($CountI,$CountO)=(0,0);
open B,'<',$blastf or die "Error:[$blastf] $!\n";
open O,'>',$out or die "Error:[$out] $!\n";
print O "#Markerid\tMarkercM\tSid\tpos\tstrand\tPidentity\tE\tBTOP\n";
while (<B>) {
	next if /^#/;
# Fields: query id, subject id, % identity, alignment length, identical, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, BTOP
#-outfmt '6 qseqid sseqid pident length nident mismatch gapopen qstart qend sstart send evalue bitscore btop'
#Chr01_457584m45 Scaffold011460  95.45   88      87      4       0       1       88      1020    933     1e-34    147    16CA2WA25MA6WA35        1
	chomp;
	my @Dat=split /\t/;
	my ($Qid,$Sid,$Pidentity,$AlnLen,$identical,$mismatches,$Gap,$Qs,$Qe,$Ss,$Se,$E,$bitScore,$BTOP,$Hit)=@Dat;
	#next if $Sid =~ /^chr/i;	# Just skip ...
	++$CountI;
	next if $Dat[14] > 1;
	my ($pos,$strand,$mcM)=@{&GetRelPos(@Dat)};
	#print "$pos,$strand\t@Dat\n";
	print O join("\t",$Qid,$mcM,$Sid,$pos,$strand,$Pidentity,$E,$BTOP),"\n";
	++$CountO;
}
close B;
print O "# $CountI,$CountO\n";
print "[!]$blastf $CountI->$CountO\t$out\n";
close O;

__END__
/ifs1/POPULATION/Rice/gapfilling/f3545ChrScaff.pos.n3

/ifs1/POPULATION/Rice/gapfilling/denovo20110113/Rice_PA64_63Km.scafSeq.chr.nfo
/ifs1/POPULATION/Rice/gapfilling/denovo20110113/out20110113/ex45Chr01.mTns6.rd
/ifs1/POPULATION/Rice/allfq/linkage/markerl/cmChr01.18

cat ../chrorder|while read a;do ./readblast.pl 15 ./out20110113/ex45$a.mTns6.rd ./out20110113/ex45$a.mTns6;done
cat ../chrorder|while read a;do ./readblast.pl 15 ./out20110113/ex45$a.mToc6.rd ./out20110113/ex45$a.mToc6;done

./blast2markerpos.pl ./markerl/cmChr01.18 ./out20110113/ex45Chr01.mTns6.rd Rice_PA64_63Km.scafSeq.chr.nfo t.out

cat ../chrorder|while read a;do ./blast2markerpos.pl ./markerl/cm$a.18 ./out20110113/ex45$a.mTns6.rd ./markerpos/m2s$a.pos;done
cat ../chrorder|while read a;do ./blast2markerpos.pl ./markerl/cm$a.18 ./out20110113/ex45$a.mToc6.rd ./markerpos/m2c$a.pos;done
