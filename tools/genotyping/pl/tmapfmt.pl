#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
#use GalaxyXS::ChromByte;
use Data::Dump qw(ddx);

unless (@ARGV > 0) {
    print "perl $0 <In_file> <Out_prefix> <ChrID>\n";
    exit 0;
}
=pod
my %varToGT=(
	0   => 'A',
	1   => 'B',
	'-' => '-',
	2   => 'H',
);
=cut
my ($Count,$population)=(0);
my (%DAT,%deDup);
open I,'<',$ARGV[0] or die $!;
my $t=<I>;
$t=~s/^#//;
while (<I>) {
    chomp;
    my ($Pos,$D)=split /\s+/,$_,2;
	$DAT{$Pos}=$D;
	push @{$deDup{$D}},$Pos;
	#++$Count;
}
close I;
no warnings;
$population=split(/\s+/,(values %DAT)[0]);
use warnings;

#ddx %deDup;
#warn $population;
#die;

for my $D (keys %deDup) {
	my ($min,$max)=sort {$a<=>$b} @{$deDup{$D}};
	$deDup{$D}=$max?[$min,$max]:[$min];
}
my @OrderD=sort {$deDup{$a}->[0] <=> $deDup{$b}->[0]} keys %deDup;
#ddx %deDup;
my ($lastPos,$CurrPos)=0;
for my $D (@OrderD) {
	$CurrPos = ($#{$deDup{$D}}==1)?($deDup{$D}->[1]):($deDup{$D}->[0]);
	warn "[!]$lastPos covers [",join(',',@{$deDup{$D}}),"].\n" if $CurrPos < $lastPos;
	$lastPos=$CurrPos;
}
$Count=keys %deDup;

open O,'>',$ARGV[1].'.intercross' or die $!;
open OO,'>',$ARGV[1].'.order' or die $!;
print O "data type f2 intercross\n$population $Count\n";
#for my $pos (sort {$a<=>$b} keys %DAT) {
for my $D (@OrderD) {
	my $pos=join('-',@{$deDup{$D}});
	my $str=$D;
	$str =~ tr/012/ABH/;
	print O $ARGV[2],"m$pos\t",join("\t",split /\s+/,$str),"\n";
	print OO '*',$ARGV[2],"m$pos\n";
}
close OO;
close O;

__END__
./tmapfmt.pl ./dat2/Chr01_num Chr01.ic Chr01
cat chrorder | while read a; do ./tmapfmt.pl ./dat2/${a}_num de$a $a;done &

./tmap-1.0/pedmerge ril1206.Chr01.intercross tt
cat chrorder | while read a; do ./tmap-1.0/pedmerge ./marker/de$a.intercross ./marker/de$a.pedigree;done &
(Name beginning with Chr12m21531432- too long, truncating)

$ cat doim2.sh
#!/bin/sh
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH,BLASTDB
#$ -cwd -r y -l vf=3g
#$ -S /bin/bash -t 1-12
#$ -o /dev/null -e /dev/null
SEEDFILE=./chrorder
a=$(sed -n -e "$SGE_TASK_ID p" $SEEDFILE)
eval ./tmap-1.0/tmap -i ./marker/de$a.pedigree ./marker/de$a.order > ./marker/de$a.imo 2> ./marker/de$a.imo.log
