#!/usr/bin/perl -w
use strict;
#use Time::HiRes qw ( gettimeofday tv_interval );
use Getopt::Long;
my ($depth_path,$fa_path,$out_file,$win_size,$help,$chrlst);
GetOptions(
	"chrlst:s"=>\$chrlst,
	"depth:s"=>\$depth_path,
	"fa:s"=>\$fa_path,
	"output:s"=>\$out_file,
	"win:i"=>\$win_size,
	"help"=>\$help
);
if(!$depth_path || !$fa_path || !$out_file || !$win_size || $help) {
	$|=1;
	print "\t -chrlst: the ChrID list\n";
	print "\t -depth: the depth path\n";
	print "\t -fa: the .fa path\n";
	print "\t -output: the output file\n";
	print "\t -win: the window size\n";
	print "\t -help: if you need the help\n";
	exit 1;
}

my ($windepcnt,%DEPGCcnt,%wincvgcnt)=(0);
open L,'<',$chrlst or die "[x]Can not open $chrlst :$!\n";
while (<L>) {
	chomp;
	my $chrid=$_;
	my ($title,$seqname,$genome);
	open FA,'<',$fa_path.'/'.$chrid.'.fa' or die "can not open $fa_path/$chrid.fa:$!";
	while (<FA>) {
		s/^>//;
		$title = $_;
		$seqname = $1 if($title =~ /^(\S+)/);
		print STDERR "loading >$seqname ...\t";
		$/=">";
		$genome=<FA>;
		chomp $genome;
		$genome=~s/\s//g;
		$/="\n";
	}
	close FA;
	my $len=length $genome;
	print "$len -> ";
	open DP,'<',$depth_path.'/'.$chrid.'.coverage' or die "can not open $depth_path/$chrid.coverage:$!";
	$title=<DP>;
	$title = $1 if($title =~ /^>(\S+)/);
	die "\n[!]ChrID not match for [$seqname] and [$title] !\n" if $title ne $seqname;
	my ($sum,$pos,$tpos)=(0,0,0);
	while (<DP>) {
		#chomp;
		#my $w=0;
		while($_=~/(\d+)/g) {
			my $num=$1;
			if ($num != 65535) {
				$sum+=$num;
				#++$w;
			}
			++$pos;
			if ( ! ($pos % $win_size) ) {
				if ($sum) {
					#my $averD=$sum/$w;
					my $seq=substr $genome,$pos-1,$win_size;
					my @STR=split //,$seq;
					my ($gccnt,$ncnt)=(0,0);
					for (@STR) {
						++$gccnt if $_ =~ /[gc]/i;
						++$ncnt if $_ =~ /n/i;
					}
					#warn "[!]N count not match @$pos:$seq [$w!=",$win_size-$ncnt,"] !\n" if $w != $win_size-$ncnt;
					my $w=$win_size-$ncnt;
					next if $w == 0;
					my $averD=$sum/$w;
					my $gcR=int(2*$gccnt/$w)/2;
					++$wincvgcnt{$gcR};
					++$windepcnt;
					++$DEPGCcnt{$averD}{$gcR};
					$tpos=$pos;
				}
				$sum=0;	#$w=0;
			}
		}	# the last half window will be skipped ...
	}
	close DP;
	print "$tpos\n";
}
close L;

open O,'>',$out_file.'.dat' or die "[x]Can not open $out_file.dat :$!\n";
open OD,'>',$out_file.'.db' or die "[x]Can not open $out_file.db :$!\n";
my @GC;
print OD '#';
for my $gc (sort {$a<=>$b} keys %wincvgcnt) {
	print OD "$gc:",$wincvgcnt{$gc},"\t";
	if ($wincvgcnt{$gc} >= 10) {
		push @GC,$gc;
	}
}
#@GC=sort {$a<=>$b} @GC;
print O join("\t",@GC),"\n";
print OD "\n#Total:$windepcnt\n\n",join("\t",@GC),"\n";
for my $dep (sort { $a<=>$b } keys %DEPGCcnt) {
	my (@out,$tout)=();
	for my $gc (@GC) {
		if (exists $DEPGCcnt{$dep}{$gc}) {
			$tout = $windepcnt*$DEPGCcnt{$dep}{$gc}/$wincvgcnt{$gc};
			print OD $DEPGCcnt{$dep}{$gc},"\t";
		} else {
			$tout = 'NA';
			print OD "NA\t";
		}
		push @out,$tout;
	}
	print O join("\t",@out),"\n";
	print OD "\n";
}
close O;
close OD;

#R plot
open(OUT,">$out_file.R")or die;
print OUT 'a=read.delim("',"$out_file.dat",'");',"\n";
print OUT 'pdf("',"$out_file.pdf",'",width=8,height=6,paper="special");',"\n";
print OUT 'b=boxplot(a,plot=F);
m=ceiling(max(b$stats,na.rm=T));
bxp(b,axes=F,outline=F,xlab=\'GC%\',ylab=\'Std Depth\',ylim=c(0,m));
axis(1,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=c(0,10,20,30,40,50,60,70,80,90,100));
axis(2,las=2);
dev.off();
q();
n
';
close OUT;
system("Rscript $out_file.R");

__END__
./gcdepchrlst.pl -chrlst chr.lst -depth ./depth/ -fa ./ref/chromFa/ -win 100 -output da100

