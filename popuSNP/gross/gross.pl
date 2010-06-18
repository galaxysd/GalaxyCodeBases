#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
#输入snp文件_wild（已经trim好）
#输入snp文件_cultivar（已经trim好）
#输入非unique文件
#use Getopt::Long;
#use Data::Dumper;

$main::VERSION=0.0.1;

our $opts='i:c:w:n:l:o:f:vb';
our ($opt_i, $opt_o, $opt_c, $opt_v, $opt_b, $opt_l, $opt_f, $opt_w, $opt_n);

our $help=<<EOH;
\t-i Population SNP file (./pop.add_cn.filter)
\t-w Wild/Original Group PSNP file (./wild.add_cn)
\t-c Cultural/Evolutionized Group PSNP file (./cul.add_cn)
\t   For multiple ChrID in single file, sort by position ASC
\t-n window size (20000) bp
\t-l step length (2000) bp
\t-f filter out smaller regins, 30% (0.3) of window size
\t-o Output File (./result.polymorphism)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./result.polymorphism' if ! defined $opt_o;
$opt_i='./pop.add_cn.filter' if ! $opt_i;
$opt_c='./cul.add_cn' if ! $opt_c;
$opt_w='./wild.add_cn' if ! $opt_w;
$opt_n=20000 if ! $opt_n;
$opt_l=2000 if ! $opt_l;
$opt_f=0.3 if ! $opt_f;
my $minCountP=int(.5+$opt_n*$opt_f);

print STDERR "From [$opt_i] to [$opt_o] for [$opt_w]->[$opt_c], with [$opt_n][$opt_l][$opt_f -> $minCountP]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}
###############
#起点	Pi(A)	Theta_Pi(A)	Tajima's D(A)	Pi(B)	Theta_Pi(B)	Tajima's D(B)	difference	Fst	length	chr_id
my ($n,@slid,$ChrID);
open SC,'<',$opt_c or die "[x]Error opening $opt_c: $!\n";
open SW,'<',$opt_w or die "[x]Error opening $opt_w: $!\n";
open SNP,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";

my ($nc,$nw,@t);
my $locname=0;
@t=split /\t/,scalar <SC>;
$nc=$t[7]+$t[8];
seek SC,0,0;
@t=split /\t/,scalar <SW>;
$nw=$t[7]+$t[8];
seek SW,0,0;
warn "[!]Haploid NFO: Original,$nw -> Evolutionized,$nc\n";

sub CalGROSS($$) {
	my ($len,$datp)=@_;
}

my (%SNP,%LastPos);
while (<SNP>) {
	my ($chr,$pos,undef,$depth,$cnr,$base1,$base2,$nbase1,$nbase2,$Q,$cn)=split /\t/;
## format_of_copyNumRST5p_v2.pl.output
#chr	loc	refBase	depth	cpnum	base1	base2	numOfBase1	numOfBase2	quality	copynum	5bp_edge	no_edge	Uniq_reads	ration_check	RST_check
#chromosome01	173	C	95	0.594737	C	T	98	4	18	1.17	27	68	84	0.00	0.364752
#seg01	62	G	26	0.500000	A	G	5	7	44	1.46	13	13	18	0.75	0.448375
	$SNP{$chr}{$pos}{P}=[$depth,$base1,$base2,$nbase1,$nbase2,$Q];
	$LastPos{$chr} = $pos unless (exists $LastPos{$chr} && $LastPos{$chr} > $pos);
}
while (<SC>) {
	my ($chr,$pos,undef,$depth,$cnr,$base1,$base2,$nbase1,$nbase2,$Q)=split /\t/;
	if ($cnr > 0) {
		$SNP{$chr}{$pos}{C}=[$depth,$base1,$base2,$nbase1,$nbase2,$Q];
		#$LastPos{$chr} = $pos unless (exists $LastPos{$chr} && $LastPos{$chr} > $pos);
	}
}
while (<SW>) {
	my ($chr,$pos,undef,$depth,$cnr,$base1,$base2,$nbase1,$nbase2,$Q)=split /\t/;
	if ($cnr > 0) {
		$SNP{$chr}{$pos}{W}=[$depth,$base1,$base2,$nbase1,$nbase2,$Q];
		#$LastPos{$chr} = $pos unless (exists $LastPos{$chr} && $LastPos{$chr} > $pos);
	}
}

my %CountP;
for my $chr (keys %SNP) {
	for (my $i=1;$i<=$LastPos{$chr}-$opt_l+1;$i++) {
		$CountP{$chr}=0;
		for (my $pos=$i;$pos<$i+$opt_n;$pos++) {
			++$CountP{$chr} if exists $SNP{$chr}{$pos}{P};# && exists $SNP{$chr}{$pos}{C} && exists $SNP{$chr}{$pos}{W};
		}
		print $CountP{$chr},"\n";
	}
}



__END__
open OUT,'>',$opt_o or die "[x]Error opening $opt_o: $!\n";
my ($i,$LastPosSNP,$LastPosSC,$LastPosSW)=($opt_l-$opt_n,0,0,0);
my ($wstart,$wend,%DAT,%CountP,$res);
CYCLE: while (1) {	# last on file end
	#($wstart,$wend)=($i+1,$i+$opt_n);

	for my $chr (keys %DAT) {
		for my $pos (1,$opt_l) {
			delete $DAT{$chr}{$pos};
		}
		for my $pos ($opt_l+1,$opt_n) {
			$DAT{$chr}{$pos-$opt_l}=$DAT{$chr}{$pos};
			delete $DAT{$chr}{$pos};
		}
	}

	SNP: seek SNP,$LastPosSNP,0;
	while (<SNP>) {
		my ($chr,$pos,undef,$depth,$cnr,$base1,$base2,$nbase1,$nbase2,$Q,$cn)=split /\t/;
## format_of_copyNumRST5p_v2.pl.output
#chr	loc	refBase	depth	cpnum	base1	base2	numOfBase1	numOfBase2	quality	copynum	5bp_edge	no_edge	Uniq_reads	ration_check	RST_check
#chromosome01	173	C	95	0.594737	C	T	98	4	18	1.17	27	68	84	0.00	0.364752
#seg01	62	G	26	0.500000	A	G	5	7	44	1.46	13	13	18	0.75	0.448375
		next if $pos-$i < $opt_n-$opt_l;
		last if $pos-$i > $opt_n;
		push @{$DAT{$chr}{$pos-$i}},[$depth,$base1,$base2,$nbase1,$nbase2,$Q];
		$LastPosSNP=tell SNP;
#warn "SNP $pos\t$LastPosSNP\n";
	}
	seek SW,$LastPosSW,0;
	while (<SW>) {
		my ($chr,$pos,undef,$depth,undef,$base1,$base2,$nbase1,$nbase2,$Q)=split /\t/;
		next if $pos-$i < $opt_n-$opt_l;
		last if $pos-$i > $opt_n;
		push @{$DAT{$chr}{$pos-$i}},[$depth,$base1,$base2,$nbase1,$nbase2,$Q];
		$LastPosSW=tell SW;
#warn "SW $pos\t$LastPosSW\n";
	}
	;
	seek SC,$LastPosSC,0;
	while (<SC>) {
		my ($chr,$pos,undef,$depth,undef,$base1,$base2,$nbase1,$nbase2,$Q)=split /\t/;
		next if $pos-$i < $opt_n-$opt_l;
		last if $pos-$i > $opt_n;
		push @{$DAT{$chr}{$pos-$i}},[$depth,$base1,$base2,$nbase1,$nbase2,$Q];
		$LastPosSC=tell SC;
#warn "SC $pos\t$LastPosSC\n";
	}
	%CountP=();
	for my $chr (keys %DAT) {
		$CountP{$chr}=0;
		for my $winpos (keys %{$DAT{$chr}}) {
			++$CountP{$chr} if $DAT{$chr}{$winpos} and scalar @{$DAT{$chr}{$winpos}} == 3;
		}
		#next if $CountP{$chr} < $minCountP;
warn "$i $chr $CountP{$chr} [",join ',',(sort {$a <=> $b} keys %{$DAT{$chr}})[0,-1],"]\n";
		$res=&CalGROSS($CountP{$chr},$DAT{$chr});
		print OUT "";
	}





	$i+=$opt_l;
}
close OUT;

__END__

for (my $i=0; ;$i+=$opt_l) {
	my $length=0;my $pi_c=0;my $theta_c=0;my $pi_w=0;my $theta_w=0;my $dif=0;
	#length 需要去掉任何没有mapping到（depth《20或》250）、copy number大于1.5的长度。
	#pi 值计算符合snp要求的位点，即在群体里符合depth 20-250 quality》=1.5  copy number小于等于1.5 theta也如此
	#先读入总体snp文件，对每一个碱基判断是否需要计算长度、是否需要计算snp，不管是否需要计算，都读入两个snp文件，
	#如果不符合，则不计算，如果符合，则进行计算
	#由于需要做sliding，所以需要记住前面的结果，如果一开始等于while里面的计数，则将该数字前面45000的结果都加起来
	#这个数组只有45000+5000个，如果往里面新加，在while里面，从5000开始替换，一直替换到while结束，每次while一开始计数都是0；

	my $s=0;	# s serves as a position indicatior
	#my $lengthe=0;my $pi_ce=0;my $theta_ce=0;my $pi_we=0;my $theta_we=0;my $dife=0;
	for (my $k=$opt_l-1;$k<$opt_n ;$k++) {
		if (!exists $slid[$k]{leth}) {next;	}
#循环替换，向前移
		$slid[$s]{piw}=$slid[$k]{piw};
		$slid[$s]{pic}=$slid[$k]{pic};
		$slid[$s]{dif}=$slid[$k]{dif};
		$slid[$s]{tw}=$slid[$k]{tw};
		$slid[$s]{tc}=$slid[$k]{tc};
		#$lengthe+=$slid[$k]{leth};
		$slid[$s]{leth}=$slid[$k]{leth};$s++;
	}

	my $m=0;	# m is the counter of the window.
	my $flag=1;
	my ($wilds,$culs);
	READSNP: while (my $block=<SNP>) {
		chomp $block;
		my @inf=split/\s+/,$block;
		$ChrID=$inf[0];
		$n=$inf[1]-1;
		$s=$n-$i;
		$slid[$s]{piw}=0;
		$slid[$s]{pic}=0;
		$slid[$s]{dif}=0;
		$slid[$s]{tw}=0;
		$slid[$s]{tc}=0;
		$slid[$s]{leth}=0;

		#print "@inf\n";exit;
	# 20,250,1.5 need to be set !
		#if ($inf[3]<20||$inf[3]>95) {$slid[$s]{leth}=0;<SC>;<SW>;}#不符合
		#if ($inf[3]>=20&&$inf[3]<=95) {
			$wilds=<SW> if $flag;
			$culs=<SC> if $flag;
			READSUB: chomp $wilds;chomp $culs;
			my @line_w=split/\s+/,$wilds;my @line_c=split/\s+/,$culs;
			my $maxP=(sort {$a <=> $b} ($line_w[1],$line_c[1]))[1];
			if ($maxP > $line_w[1]) {
				defined($wilds=<SW>) or last;
#				warn "[w]$s\t$inf[1]\t$line_w[1]\t$line_c[1]\n";
				goto READSUB;
			} elsif ($maxP > $line_c[1]) {
				defined($culs=<SC>) or last;
#				warn "[c]$s\t$inf[1]\t$line_w[1]\t$line_c[1]\n";
				goto READSUB;
			}
			if ($maxP > $inf[1]) {
				$flag=0;
				$slid[$s]{leth}=1;
				#++$s;
#				warn "[p]$s\t$inf[1]\t$line_w[1]\t$line_c[1]\n";
				goto READSNP;
			}
			$flag=1;
			if ($line_w[3]>0&&$line_c[3]>0) {
				$length++;$slid[$s]{leth}=1;
#				$slid[$s]{piw}="";
#				$slid[$s]{pic}="";
#				$slid[$s]{dif}="";
#				$slid[$s]{tw}="";
#				$slid[$s]{tc}="";
				warn "[p]$inf[1]\t$s\t$m\t$n\n";
			}#可以统计的位点

			$locname=$line_w[1];
			if ($inf[9]>=15&&($line_w[3]>0&&$line_c[3]>0)) {
			#if ($inf[9]>=15&&($line_w[3]>0&&$line_c[3]>0)) {
				if ($line_w[1]!=$inf[1]||$line_c[1]!=$inf[1]) {print OUT "$inf[1]\t$line_w[1]\t$line_c[1]\terror\n";exit;}

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
		#}
#		$s++;
#		$n++;
		if ($n>$i+$opt_n) {last;		}
		$m=1;

	}
	#print "$m\n";
	if ($m!=1) {exit;	}#判断是否有文件读入，如果没有读入，跳出
	my $lengthe=0;my $pi_ce=0;my $theta_ce=0;my $pi_we=0;my $theta_we=0;my $dife=0;

#	$length+=$lengthe;$pi_c+=$pi_ce;$theta_c+=$theta_ce;$pi_w+=$pi_we;$theta_w+=$theta_we;$dif+=$dife;
	for (my $in=0;$in<$opt_n ;$in++) {
		$pi_we+=$slid[$in]{piw};
		$pi_ce+=$slid[$in]{pic};
		$dife+=$slid[$in]{dif};
		$theta_we+=($slid[$in]{tw}>0);
		$theta_ce+=($slid[$in]{tc}>0);
		$lengthe+=$slid[$in]{leth};
	}
	if ($lengthe<=($opt_l)||$pi_ce==0) {print OUT "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";next;	}#如果能够mapping的区域小于window的4/10则不输出结果
	# 0.1 (<0.3) still makes sense, right ?

	my $diff = $dife / $lengthe;
	my $Dc=tajima($nc,$pi_ce,$theta_ce);
	my $Dw=tajima($nw,$pi_we,$theta_we);
	my $fst=fst($pi_ce,$nc,$pi_we,$nw,$dife);
	my $piavc=$pi_ce/$lengthe;my $piavw=$pi_we/$lengthe;
	print OUT "$i\t$piavc\t$theta_ce\t$Dc\t$piavw\t$theta_we\t$Dw\t$fst\t$diff\t$lengthe\t$ChrID\n";
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
cat chrorder | while read a;do echo "#$ -N \"${a}_poly\"" >./shell/${a}_PO.sh;echo "#$ -cwd -r y -l vf=1M,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e ./outpoly/_$a.err" >> ./shell/${a}_PO.sh;echo ./polymorphism.pl -snp_w ./wild/$a.add_cn -snp_c ./cultivate/$a.add_cn -snpdb ./Add/$a.add_cn -ratio ./population/$a.population.snp -n 5000 -bin 0.1 -cutoff 3 -chr $a -cn 1.16 -f 0.1 -o ./outpoly/$a.polymorphism >> ./shell/${a}_PO.sh; done

cat chrorder | while read a;do echo "#$ -N \"${a}_polyf\"" >./shell/${a}_POf.sh;echo "#$ -cwd -r y -l vf=1M,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e ./outpoly/_$a.err" >> ./shell/${a}_POf.sh;echo ./polymorphism.pl -snp_w ./wild/$a.add_cn -snp_c ./cultivate/$a.add_cn -snpdb ./Add/$a.add_cn -ratio ./population/$a.population.snp.f -n 5000 -bin 0.1 -cutoff 3 -chr $a -cn 1.16 -f 0.1 -o ./outpoly/$a.polymorphism >> ./shell/${a}_POf.sh; done

rm -f ./outpoly/mix.mpoly && find ./outpoly/*.polymorphism | xargs cat >> ./outpoly/mix.mpoly

perl -lane 'print unless /^NA/' outpoly/*.polymorphism |les
perl -lane '@a=split /\t/;print unless (/^NA/ or $#a < 5)' outpoly/*.polymorphism |les

perl -lane 'print unless /^NA/' outpoly/*.polymorphism > ./outpoly/mix.mpoly



cat chrorder | while read a;do echo "#$ -N \"${a}_poly\"" >./shell/${a}_wPO.sh;echo "#$ -cwd -r y -l vf=10M,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e /dev/null" >> ./shell/${a}_wPO.sh;echo ./polymorphism.pl -snp_w ./wildadd/$a.add_cn -snp_c ./culadd/$a.add_cn -snpdb ./alladd/$a.add_cn -ratio ./finalsnp/$a.population.snp -n 20000 -bin 0.01 -cutoff 3 -chr $a -cn 1.5 -f 0.1 -o ./outpoly_20k_10_15/$a.polymorphism >> ./shell/${a}_wPO.sh; done
cat chrorder | while read a;do echo "#$ -N \"${a}_poly\"" >./shell/${a}_qPO.sh;echo "#$ -cwd -r y -l vf=10M,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e /dev/null" >> ./shell/${a}_qPO.sh;echo ./polymorphism.pl -snp_w ./wildadd/$a.add_cn -snp_c ./culadd/$a.add_cn -snpdb ./alladd/$a.add_cn -ratio ./finalsnp/$a.population.snp -n 20000 -bin 0.01 -cutoff 3 -chr $a -cn 1.5 -f 0.2 -o ./outpoly_20k_20_15/$a.polymorphism >> ./shell/${a}_qPO.sh; done
