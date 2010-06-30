#!/usr/bin/perl -w
#use strict;
#use Time::HiRes qw ( gettimeofday tv_interval );
use Getopt::Long;
my ($depth_file,$fa_file,$out_file,$win_size,$help);
GetOptions(
	"depth:s"=>\$depth_file,
	"fa:s"=>\$fa_file,
	"output:s"=>\$out_file,
	"win:i"=>\$win_size,
	"help"=>\$help
);
if(!$depth_file || !$fa_file || !$out_file || !$win_size || $help) {
	$|=1;
	print "\t -depth: the depth file\n";
	print "\t -fa: the .fa file\n";
	print "\t -output: the output file\n";
	print "\t -win: the window size\n";
	print "\t -help: if you need the help\n";
	exit 1;
}

#getDepth
print "Begin to count average of depth\n";
my %depth=();
my $cnt=0;
my $sum=0;
my $chr;
open(IN,'<',"$depth_file")or die "[x]Can not open $depth_file :$!\n";
while(<IN>) {
	chomp;
	if($_ !~ /^>/) {
		while($_=~/(\d+)/g) {
			my $num=$1;
			if($num==65535){$num=0;}
			$sum+=$num;
			$cnt++;
			if($cnt==$win_size) {
				my $aver=$sum/$cnt;
				push @{$chr},$aver;
				$sum=0;
				$cnt=0;
			}
		}
	} else {
		$chr=$_;
		$depth{$chr}=$chr;
		print "$chr\n";
	}
}
close IN;
close OUT;
print "Finish the depth counting\n";

#countGC
$chr=();
$cnt=0;
$sum=0;
print "Begin to count the gc\n";
open(IN,"<$fa_file")or die "can not open $fa_file:$!";
while(<IN>) {
	chomp;
	if($_!~/^>/) {
		my $num=s/[gc]/Y/ig;
		$sum+=$num;
		$cnt+=length($_);
		if($cnt==$win_size) {
			my $aver=int($sum/$cnt*100);
			push @{$chr},$aver;
			$sum=0;
			$cnt=0;
		} elsif($cnt>$win_size) {
			my $dis=$cnt-$win_size;
			my @chr=split(//,$_);
			my $len=@chr;
			my $count=0;
			for(my $i=$len-$dis;$i<$len;$i++) {
				if($chr[$i]=~/[gc]/i) {
					$count++;
				}
			}
			my $aver=int(($sum-$count)/$win_size*100);
			push @{$chr},$aver;
			$sum=$count;
			$cnt=$dis;
		} else {next;}
	} else {
		$chr="gc$_";
		print "$_\n";
	}
}
close IN;

foreach(keys %depth) {
	my $chrname=$_;
	foreach(@{"gc$chrname"}) {
		my $dep=shift @{$chrname};
		push @{$_},$dep;
	}
}

my $str=0;
my $end=90;
my $na=0;
my $line=$end-$str;
my $max_depth=0;
my $cnt_depth=0;
open(OUT,">$out_file.dat")or die "$!\n";
for(my $i=$str;$i<=$end;$i++) {
	print OUT "$i\t";
}
print OUT "\n";
while($na<=$line) {
	$na=0;
	for(my $j=$str;$j<=$end;$j++) {
		my $dep=shift @{$j};
		if(!defined $dep){$dep="NA";$na++;}
		#else{$max_depth+=$dep;$cnt_depth++;}
		print OUT "$dep\t";
	}
	print OUT "\n";
}
close OUT;

#R plot
open(OUT,">$out_file.R")or die;
print OUT 'a=read.delim("',"$out_file.dat",'");',"\n";
print OUT 'pdf("',"$out_file.pdf",'",width=8,height=6,paper="special");',"\n";
print OUT 'b=boxplot(a,plot=F);
m=ceiling(max(b$stats,na.rm=T));
bxp(b,axes=F,outline=F,xlab=\'GC%\',ylab=\'Depth\',ylim=c(0,m));
axis(1,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=c(0,10,20,30,40,50,60,70,80,90,100));
axis(2,las=2);
dev.off();
q();
n
';
close OUT;
system("Rscript $out_file.R");
