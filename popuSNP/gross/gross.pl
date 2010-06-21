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
\t-i Chr.nfo (Chr.nfo)
\t-w Wild/Original Group PSNP file (./wild.snp)
\t-c Cultural/Evolutionized Group PSNP file (./cul.snp)
\t-n window size (20000) bp
\t-l step length (2000) bp
\t-f filter out smaller regins, 30% (0.3) of window size
\t-o Output File (./result.polymorphism)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./result.polymorphism' if ! defined $opt_o;
$opt_i='Chr.nfo' if ! $opt_i;
$opt_c='./cul.snp' if ! $opt_c;
$opt_w='./wild.snp' if ! $opt_w;
$opt_n=20000 if ! $opt_n;
$opt_l=2000 if ! $opt_l;
$opt_f=0.3 if ! $opt_f;
my $minCountP=int(.5+$opt_n*$opt_f);

print STDERR "From [$opt_i] to [$opt_o] for [$opt_w]->[$opt_c], with [$opt_n][$opt_l][$opt_f -> $minCountP]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}
###############
#起点	Pi(A)	Theta_Pi(A)	Tajima's D(A)	Pi(B)	Theta_Pi(B)	Tajima's D(B)	difference	Fst	length	chr_id
my (@ChrIDs,%ChrLen);
open CHRLEN,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
while (<CHRLEN>) {
	chomp;
	s/\r//g;
	my ($chr,$len)=split /\s+/;
	push @ChrIDs,$chr;
	$ChrLen{$chr}=$len;
}
close CHRLEN;

open SC,'<',$opt_c or die "[x]Error opening $opt_c: $!\n";
open SW,'<',$opt_w or die "[x]Error opening $opt_w: $!\n";

my ($nc,$nw,@t);
my $locname=0;
@t=split /\t/,scalar <SC>;
$nc=$t[7]+$t[8];
seek SC,0,0;
@t=split /\t/,scalar <SW>;
$nw=$t[7]+$t[8];
seek SW,0,0;
warn "[!]Haploid NFO: Original,$nw -> Evolutionized,$nc\n";

sub Calpre($$) {
	my ($base1,$base2,$nbase1,$nbase2)=@_;
	my $pi=$nbase2*$nbase1*2/(($nbase2+$nbase1)*($nbase2+$nbase1-1));
	my %f;
	$f{$base2}=$nbase2/($nbase2+$nbase1);
	$f{$base1}=1-$f{$base2};
	return [$pi,\%f];
}

sub CheckPos($$) {
	my ($chr,$pos)=@_;
	;
	return 1;	# True
}
my %SNP;
## format_of_copyNumRST5p_v2.pl.output
#chr	loc	refBase	depth	cpnum	base1	base2	numOfBase1	numOfBase2	quality	copynum	5bp_edge	no_edge	Uniq_reads	ration_check	RST_check
#chromosome01	173	C	95	0.594737	C	T	98	4	18	1.17	27	68	84	0.00	0.364752
#seg01	62	G	26	0.500000	A	G	5	7	44	1.46	13	13	18	0.75	0.448375

while (<SC>) {
	my ($chr,$pos,$ref,$depth,$cnr,$base1,$base2,$nbase1,$nbase2,$Q)=split /\t/;
	next unless &CheckPos($chr,$pos);
	${$SNP{$chr}{$pos}}[0]=$ref;
	${$SNP{$chr}{$pos}}[2]=&Calpre($base1,$base2,$nbase1,$nbase2);
}
while (<SW>) {
	my ($chr,$pos,$ref,$depth,$cnr,$base1,$base2,$nbase1,$nbase2,$Q)=split /\t/;
	${$SNP{$chr}{$pos}}[0]=$ref;
	${$SNP{$chr}{$pos}}[1]=&Calpre($base1,$base2,$nbase1,$nbase2);
}
for my $chr (sort keys %SNP) {
	for my $pos (sort {$a <=> $b} keys %{$SNP{$chr}}) {
		my ($ref,$SW,$SC)=@{$SNP{$chr}{$pos}};
		my (%fw,%fc);
		if ($SW) {
			%fw=%{$$SW[1]};
		} else { %fw=($ref => 1,); }
		if ($SC) {
			%fc=%{$$SC[1]};
		} else { %fc=($ref => 1,); }
		my (%t,$ta);
		++$t{$_} for keys %fw;
		++$t{$_} for keys %fc;
		($ta)=sort {$t{$b} <=> $t{$a}} keys %t;
		$fw{$ta}=0 unless defined $fw{$ta};
		$fc{$ta}=0 unless defined $fc{$ta};
		my $diff=1-($fw{$ta}*$fc{$ta}+(1-$fw{$ta})*(1-$fc{$ta}));
		${$SNP{$chr}{$pos}}[3]=$diff;
	}
}
warn "[!]loading done !\n";

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

open OUT,'>',$opt_o or die "[x]Error opening $opt_o: $!\n";
my ($posA, $pi_we, $pi_ce, $dife, $theta_we, $theta_ce, $lengthe)=(0);
for my $chr (sort keys %SNP) {
	for(my $pos0 = 1; $pos0 <= $ChrLen{$chr}-$opt_n; $pos0+=$opt_l) {	# $pos0 < $ChrLen{$chr}-$opt_n+1; +1 is useless if $opt_l > 1
		($pi_we, $pi_ce, $dife, $theta_we, $theta_ce, $lengthe)=(0,0,0,0,0,0);
		for (my $pos = $pos0; $pos < $pos0+$opt_n; $pos++) {
			#last if $pos > $ChrLen{$chr};
			++$lengthe if &CheckPos($chr,$pos);
			next unless defined $SNP{$chr}{$pos};
			my ($ref,$SW,$SC,$diff)=@{$SNP{$chr}{$pos}};
			if ($SW) {
				++$theta_we;
				$pi_we += $$SW[0];
			}
			if ($SC) {
				++$theta_ce;
				$pi_ce += $$SC[0];
			}
			$dife += $diff;
		}
		if ($lengthe<=($opt_n*$opt_f)||$pi_ce==0) {
			print OUT "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";next;
		}#如果能够mapping的区域小于window的4/10则不输出结果
		my $diff = $dife / $lengthe;
		my $Dc=tajima($nc,$pi_ce,$theta_ce);
		my $Dw=tajima($nw,$pi_we,$theta_we);
		my $fst=fst($pi_ce,$nc,$pi_we,$nw,$dife);
		my $piavc=$pi_ce/$lengthe;my $piavw=$pi_we/$lengthe;
		#print OUT $posA+$pos0,"\t$pos0\t$piavc\t$theta_ce\t$Dc\t$piavw\t$theta_we\t$Dw\t$fst\t$diff\t$lengthe\t$chr";
		print OUT "$pos0\t$piavc\t$theta_ce\t$Dc\t$piavw\t$theta_we\t$Dw\t$fst\t$diff\t$lengthe\t$chr\n";
	}
	#$posA += int(.5+$ChrLen{$chr}/1000)*1000;
}
close OUT;
warn "[!]All done !\n";


__END__
cat chrorder | while read a;do echo "#$ -N \"${a}_poly\"" >./shell/${a}_PO.sh;echo "#$ -cwd -r y -l vf=1M,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e ./outpoly/_$a.err" >> ./shell/${a}_PO.sh;echo ./polymorphism.pl -snp_w ./wild/$a.add_cn -snp_c ./cultivate/$a.add_cn -snpdb ./Add/$a.add_cn -ratio ./population/$a.population.snp -n 5000 -bin 0.1 -cutoff 3 -chr $a -cn 1.16 -f 0.1 -o ./outpoly/$a.polymorphism >> ./shell/${a}_PO.sh; done

cat chrorder | while read a;do echo "#$ -N \"${a}_polyf\"" >./shell/${a}_POf.sh;echo "#$ -cwd -r y -l vf=1M,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e ./outpoly/_$a.err" >> ./shell/${a}_POf.sh;echo ./polymorphism.pl -snp_w ./wild/$a.add_cn -snp_c ./cultivate/$a.add_cn -snpdb ./Add/$a.add_cn -ratio ./population/$a.population.snp.f -n 5000 -bin 0.1 -cutoff 3 -chr $a -cn 1.16 -f 0.1 -o ./outpoly/$a.polymorphism >> ./shell/${a}_POf.sh; done

rm -f ./outpoly/mix.mpoly && find ./outpoly/*.polymorphism | xargs cat >> ./outpoly/mix.mpoly

perl -lane 'print unless /^NA/' outpoly/*.polymorphism |les
perl -lane '@a=split /\t/;print unless (/^NA/ or $#a < 5)' outpoly/*.polymorphism |les

perl -lane 'print unless /^NA/' outpoly/*.polymorphism > ./outpoly/mix.mpoly



cat chrorder | while read a;do echo "#$ -N \"${a}_poly\"" >./shell/${a}_wPO.sh;echo "#$ -cwd -r y -l vf=10M,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e /dev/null" >> ./shell/${a}_wPO.sh;echo ./polymorphism.pl -snp_w ./wildadd/$a.add_cn -snp_c ./culadd/$a.add_cn -snpdb ./alladd/$a.add_cn -ratio ./finalsnp/$a.population.snp -n 20000 -bin 0.01 -cutoff 3 -chr $a -cn 1.5 -f 0.1 -o ./outpoly_20k_10_15/$a.polymorphism >> ./shell/${a}_wPO.sh; done
cat chrorder | while read a;do echo "#$ -N \"${a}_poly\"" >./shell/${a}_qPO.sh;echo "#$ -cwd -r y -l vf=10M,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e /dev/null" >> ./shell/${a}_qPO.sh;echo ./polymorphism.pl -snp_w ./wildadd/$a.add_cn -snp_c ./culadd/$a.add_cn -snpdb ./alladd/$a.add_cn -ratio ./finalsnp/$a.population.snp -n 20000 -bin 0.01 -cutoff 3 -chr $a -cn 1.5 -f 0.2 -o ./outpoly_20k_20_15/$a.polymorphism >> ./shell/${a}_qPO.sh; done
