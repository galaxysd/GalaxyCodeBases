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
	print STDERR "$len -> ";
	open DP,'<',$depth_path.'/'.$chrid.'.coverage' or die "can not open $depth_path/$chrid.coverage:$!";
	$title=<DP>;
	$title = $1 if($title =~ /^>(\S+)/);
	die "\n[!]ChrID not match for [$seqname] and [$title] !\n" if $title ne $seqname;
	my ($sum,$pos,$tpos,$cnt)=(0,0,0,0);
	while (<DP>) {
		#chomp;
		while($_=~/(\d+)/g) {
			my $num=$1;
			if ($num != 0 and $num != 65535) {
				$sum+=$num;
				++$cnt;
			}
			++$pos;
			if ( ! ($pos % $win_size) ) {
				my $seq=substr $genome,$pos-1,$win_size;
				my @STR=split //,$seq;
				my ($gccnt,$ncnt)=(0,0);
				for (@STR) { 
					++$gccnt if $_ =~ /[gc]/i;
					++$ncnt if $_ =~ /n/i;
				}
				my $w=$win_size-$ncnt;
				next unless $w;
				$windepcnt += $sum;
				#my $gcR=int(200*$gccnt/$w)/2;
				my $gcR=int($gccnt/$w);	# later plotting needs 0..100
				++$wincvgcnt{$gcR};
				if ($w>=3) {	# at least 3 points needed.
					#my $averD=$sum/$w;
					#warn '[!]N count not match @',"$pos:$seq [$w<$cnt] !\n" if $cnt > $w;
					#next if $w < 3;
					my $averD=$sum/$w;
					++$DEPGCcnt{$averD}{$gcR};
					$tpos=$pos;
#print "$cnt $w,";
				}
				$sum=0;	$cnt=0;
			}
		}	# the last half window will be skipped ...
	}
	close DP;
	print STDERR "$tpos\t$windepcnt\n";
}
close L;

open O,'>',$out_file.'.dat' or die "[x]Can not open $out_file.dat :$!\n";
open OD,'>',$out_file.'.db' or die "[x]Can not open $out_file.db :$!\n";
open OXY,'>',$out_file.'.xy' or die "[x]Can not open $out_file.xy :$!\n";
my ($gcsum,@GC)=(0);
print OD '#';
for my $gc (sort {$a<=>$b} keys %wincvgcnt) {
	print OD "$gc:",$wincvgcnt{$gc},"\t";
	#if ($wincvgcnt{$gc} >= 30) {
		push @GC,$gc;
		$gcsum += $wincvgcnt{$gc};
	#}
}
print OXY "GC\tR\n";
for my $gc (@GC) {
	print OXY "$gc\t",$wincvgcnt{$gc}/$gcsum,"\n";
}
close OXY;
#@GC=sort {$a<=>$b} @GC;
@GC=0..100;
print O join("\t",@GC),"\n";
print OD "\n#Total:$windepcnt\n\n",join("\t",@GC),"\n";
for my $dep (sort { $a<=>$b } keys %DEPGCcnt) {
	my (@out,$tout)=();
	for my $gc (@GC) {
		if (exists $DEPGCcnt{$dep}{$gc}) {
			$tout = $DEPGCcnt{$dep}{$gc}*$gcsum/($wincvgcnt{$gc}*$windepcnt);
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
print OUT 'a=read.delim("',"$out_file.db",'");',"\n";
print OUT 'l=read.delim("',"$out_file.xy",'");',"\n";
print OUT 'pdf("',"$out_file.pdf",'",width=8,height=6,paper="special");',"\n";
print OUT 'b=boxplot(a,plot=F);
m=ceiling(max(b$stats,na.rm=T));
n=floor(max(b$stats[3,],na.rm=T)/max(l$R))
bxp(b,axes=F,outline=F,xlab=\'GC%\',ylab=\'Std Depth\',ylim=c(0,m));
lines(l$GC,l$R*n,type=\'l\',col=\'blue\');
axis(1,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=c(0,10,20,30,40,50,60,70,80,90,100));
axis(2,las=2);
axis(4,at=0.01*n*c(0,2,4,6,8,10,15,20,25,30,50,75,100),labels=c(0,2,4,6,8,10,15,20,25,30,50,75,100));
dev.off();
q();
n
';
close OUT;
system("Rscript $out_file.R");

__END__
./gcdepchrlst.pl -chrlst chr.lst -depth ./depth/ -fa ./ref/chromFa/ -win 100 -output da100

