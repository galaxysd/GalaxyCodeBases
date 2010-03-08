#!/bin/env perl
use strict;
#输入snp文件_wild（已经trim好）
#输入snp文件_cultivar（已经trim好）
#输入非unique文件
use Getopt::Long;
#use Data::Dumper;

###############
my %opts;

GetOptions(\%opts,"snp_w=s","snp_c=s","snpdb=s","o=s","n=s","bin=s","h");

#&help()if(defined $opts{h});
if(!defined($opts{snp_w}) || !defined($opts{snp_c}) ||!defined($opts{n}) ||!defined($opts{snpdb})||!defined($opts{bin})||defined($opts{h}) ){
	my $ver="1.0";
	print <<"	Usage End.";
	Description:
		calculate pi、theta、D、diff 、Fst in different groups
		Version: $ver 0.1

	Usage:
		-snp_w/snp_c   snp file        Must be given, multi snp result
		-snpdb	selected snps		groups snps (including copy number infor)
		-o   output       Must be given output file names
		-n   window size
		-bin   sliding bin       for example: 10%(0.1) for moving
		-h    Help document

	out put format :
		location pi-1 theta-1 jajima's_D-1 pi-2 theta-2 jajima's_D-2 difference Fst length;

	Usage End.

	exit;
}

#open FILT,$opts{ratio};
#if (!defined $opts{cutoff}){$opts{cutoff}=2;}
#my %pval;
#while(<FILT>){
#	chomp;
#	my @inf=split(/\t/);
#	if ($inf[12]<=$opts{cutoff}){
#		$pval{$inf[1]}=$inf[12];
#	}
#}
#print "a\n";
open SC,"$opts{snp_c}";
open SW,"$opts{snp_w}";open DB,"$opts{snp_w}";
my $dbn=<DB>;
my @st=split/\s+/,$dbn;

my $tmp=<SC>;
my @t=split/\s+/,$tmp;
my $nc=$t[7]+$t[8];
#$nc = 38;
$tmp=<SW>;
@t=split/\s+/,$tmp;
my $nw=$t[7]+$t[8];
#$nw = 16;
close SC;
close SW;
open SNP,"$opts{snpdb}";
open OUT,">$opts{o}";
my $window=$opts{n};
my $bin=$opts{bin};
#print "$bin\n";
my $n;my @slid;
open SW,"$opts{snp_w}";
open SC,"$opts{snp_c}";
my $locname=0;
#print "$locname\n";
for (my $i=0; ;$i+=$window*$bin) {
	my $length=0;my $pi_c=0;my $theta_c=0;my $pi_w=0;my $theta_w=0;my $dif=0;
	#length 需要去掉任何没有mapping到（depth《20或》250）、copy number大于1.5的长度。
	#pi 值计算符合snp要求的位点，即在群体里符合depth 20-250 quality》=1.5  copy number小于等于1.5 theta也如此
	#先读入总体snp文件，对每一个碱基判断是否需要计算长度、是否需要计算snp，不管是否需要计算，都读入两个snp文件，
	#如果不符合，则不计算，如果符合，则进行计算
	#由于需要做sliding，所以需要记住前面的结果，如果一开始等于while里面的计数，则将该数字前面45000的结果都加起来
	#这个数组只有45000+5000个，如果往里面新加，在while里面，从5000开始替换，一直替换到while结束，每次while一开始计数都是0；

	my $s=0;
	#my $lengthe=0;my $pi_ce=0;my $theta_ce=0;my $pi_we=0;my $theta_we=0;my $dife=0;
	for (my $k=$window*$bin-1;$k<$window ;$k++) {
		if (!exists $slid[$k]{leth}) {next;	}


#				$pi_we+=$slid[$k]{piw};
#				$pi_ce+=$slid[$k]{pic};
#				$dife+=$slid[$k]{dif};
#				$theta_we+=($slid[$k]{tw}>0);
#				$theta_ce+=($slid[$k]{tc}>0);
#循环替换，向前移

				$slid[$s]{piw}=$slid[$k]{piw};
				$slid[$s]{pic}=$slid[$k]{pic};
				$slid[$s]{dif}=$slid[$k]{dif};
				$slid[$s]{tw}=$slid[$k]{tw};
				$slid[$s]{tc}=$slid[$k]{tc};







		#$lengthe+=$slid[$k]{leth};
		$slid[$s]{leth}=$slid[$k]{leth};$s++;		}

	my $m=0;
	while (my $block=<SNP>) {
		chomp $block;
		my @inf=split/\s+/,$block;
		$slid[$s]{piw}=0;
		$slid[$s]{pic}=0;
		$slid[$s]{dif}=0;
		$slid[$s]{tw}=0;
		$slid[$s]{tc}=0;
		$slid[$s]{leth}=0;



		#print "@inf\n";exit;
	# 20,250,1.5 need to be set !
		if ($inf[3]<20||$inf[3]>95||$inf[10]>1.16) {$slid[$s]{leth}=0;<SC>;<SW>;}#不符合
		if ($inf[3]>=20&&$inf[3]<=95&&$inf[10]<=1.16) {

				my $wilds=<SW>;
				my $culs=<SC>;
				chomp $wilds;chomp $culs;
				my @line_w=split/\s+/,$wilds;my @line_c=split/\s+/,$culs;
				if ($line_w[3]>0&&$line_c[3]>0) {$length++;$slid[$s]{leth}=1;
#				$slid[$s]{piw}="";
#				$slid[$s]{pic}="";
#				$slid[$s]{dif}="";
#				$slid[$s]{tw}="";
#				$slid[$s]{tc}="";
#
				}#可以统计的位点


				$locname=$line_w[1];
			#if ($inf[9]>=15&&($line_w[3]>0&&$line_c[3]>0)&&(exists $pval{$inf[1]})) {
			if ($inf[9]>=15&&($line_w[3]>0&&$line_c[3]>0)) {
				if ($line_w[1]!=$inf[1]||$line_c[1]!=$inf[1]) {print OUT "$inf[1]\t$line_w[1]\t$line_c[1]\terror\n";exit;				}


				my $piw=$line_w[8]*$line_w[7]*2/(($line_w[8]+$line_w[7])*($line_w[8]+$line_w[7]-1));
				$pi_w+=$piw;
				my %fw;my %fc;
				$fw{$line_w[6]}=$line_w[8]/($line_w[8]+$line_w[7]); $fw{$line_w[5]}=1-$fw{$line_w[6]};

				my $pic=$line_c[8]*$line_c[7]*2/(($line_c[8]+$line_c[7])*($line_c[8]+$line_c[7]-1));
					$pi_c+=$pic;
				$fc{$line_c[6]}=$line_c[8]/($line_c[8]+$line_c[7]); $fc{$line_c[5]}=1-$fc{$line_c[6]};
				my $diff=1-($fw{$inf[5]}*$fc{$inf[5]}+(1-$fw{$inf[5]})*(1-$fc{$inf[5]}));
				$dif+=$diff;
				if ($line_w[8]>0) {$theta_w++;			}
				if ($line_c[8]>0) {$theta_c++;				}
				$slid[$s]{piw}=$piw;
				$slid[$s]{pic}=$pic;
				$slid[$s]{dif}=$diff;
				$slid[$s]{tw}=($line_w[8]>0);
				$slid[$s]{tc}=($line_c[8]>0);



			}#有效snp位点



		}
		$s++;
		$n++;
		if ($n>$i+$window) {last;		}
		$m=1;

	}
	#print "$m\n";
	if ($m!=1) {exit;	}#判断是否有文件读入，如果没有读入，跳出
	my $lengthe=0;my $pi_ce=0;my $theta_ce=0;my $pi_we=0;my $theta_we=0;my $dife=0;

#	$length+=$lengthe;$pi_c+=$pi_ce;$theta_c+=$theta_ce;$pi_w+=$pi_we;$theta_w+=$theta_we;$dif+=$dife;
	for (my $in=0;$in<$window ;$in++) {
				$pi_we+=$slid[$in]{piw};
				$pi_ce+=$slid[$in]{pic};
				$dife+=$slid[$in]{dif};
				$theta_we+=($slid[$in]{tw}>0);
				$theta_ce+=($slid[$in]{tc}>0);
				$lengthe+=$slid[$in]{leth};


	}
	if ($lengthe<=($window*0.3)||$pi_ce==0) {print OUT "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";next;	}#如果能够mapping的区域小于window的4/10则不输出结果

	my $diff = $dife / $lengthe;
	my $Dc=tajima($nc,$pi_ce,$theta_ce);
	my $Dw=tajima($nw,$pi_we,$theta_we);
	my $fst=fst($pi_ce,$nc,$pi_we,$nw,$dife);
	my $piavc=$pi_ce/$lengthe;my $piavw=$pi_we/$lengthe;
	print  OUT "$i\t$piavc\t$theta_ce\t$Dc\t$piavw\t$theta_we\t$Dw\t$fst\t$diff\t$lengthe\n";

}

warn "[!]Done.\n";


sub tajima{
	my $num=$_[0];
	my $pi=$_[1];
	my $t=$_[2];

	my $a1c;
	for (my $i=1;$i<=$num-1;$i++) {
		$a1c+=1/$i;
	}

	my $a2c;

	for (my $i=1;$i<=$num-1;$i++) {
		$a2c+=1/($i*$i);
	}

	my $b1c=($num+1)/(3*($num-1));

	my $b2c=2*($num*$num+$num+3)/(9*$num*($num-1));

	my $c1c=$b1c-(1/$a1c);

	my $c2c=$b2c-(($num+2)/($a1c*$num))+$a2c/($a1c*$a1c);

	my $e1c=$c1c/$a1c;
	my $e2c=$c2c/($a1c*$a1c+$a2c);
	my $Dc = 0;
	if ($e1c*$t+$e2c*$t*($t-1) != 0){
		$Dc=($pi-($t/$a1c))/sqrt($e1c*$t+$e2c*$t*($t-1));
	}
	return $Dc;
}

sub fst{
	#should input the parwise comparision of each group and parwise comparsion between two groups
	#should input the number of each group
	#format fst(group1 pi,group1 number,group2 pi,group2 number,difference between two groups)
	my $p1=$_[0];my $n1=$_[1];my $p2=$_[2];my $n2=$_[3];my $di=$_[4];
	my $hw=($n1*($n1-1)*$p1+$n2*($n2-1)*$p2)/($n1*($n1-1)+$n2*($n2-1));
	my $hb=$di;
	my $fst=1-($hw)/($hb);
	return $fst;
}

__END__
cat chrorder | while read a;do echo "#$ -N \"${a}_poly\"" >./outpoly/$a.sh;echo "#$ -cwd -r y -l vf=1M,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e ./outpoly/$a.err" >> ./outpoly/$a.sh;echo ./polymorphism.pl -snp_w ./wild/$a.add_cn -snp_c ./cultivate/$a.add_cn -snpdb ./Add/$a.add_cn -n 5000 -bin 0.1 -o ./outpoly/$a.polymorphism >> ./outpoly/$a.sh; done
