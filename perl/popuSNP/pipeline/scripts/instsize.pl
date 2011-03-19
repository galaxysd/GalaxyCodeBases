#!/bin/env perl
use strict;
use warnings;
#use File::Basename;

unless (@ARGV){
	print "perl $0 <fq list file> <max read len> <ref> <output prefix>\n";
	exit;
}

my ($fqlst,$readlen,$ref,$out)=@ARGV;
my $bin='/panfs/GAG/huxuesong/scripts/soap2.20';
my $arg0='-p 8 -x 1000 -t -s 40 -l 32';# -m 10 -g 5 -r 2';
my $OKvalue=0.45;	# Peak is OK is use as it is main peak(dropped < 0.45)
my $WARNvalue=0.4;	# soap 10-1000 is OK as single < 0.4
my $minlen;
open LEN,'<',$readlen or die "[x]Error opening $readlen: $!\n";
$readlen = <LEN>;
$minlen = <LEN>;
chomp $readlen;
chomp $minlen;
close LEN;
my $mismatch=$readlen>70?3:1;
$minlen -= $mismatch;
open LST,'<',$fqlst or die "[x]Error opening $fqlst: $!\n";
my $count=0;
system('mkdir','-p',$out);
my $sh;
while (<LST>) {
	chomp;
	my ($type,$fq1,$fq2)=split /\t/;
	if ($type eq 'PE') {
		++$count;
		#my $outsh=$out.'/insoap_'.$count.'.sh';
		my $outbase=$fq1;
		$outbase =~ s/(_1)?\.fq$//i;
		$sh="$bin -a $fq1 -b $fq2 -D $ref -o $outbase.soap -2 /dev/null $arg0 -m $minlen -v $mismatch 2>$outbase.log";
	TEST:
		if (-s "$outbase.log") {
			system("mv -f ${outbase}_insoap.sh ${outbase}_insoap.oldsh") if (-e "${outbase}_insoap.sh");
		} else {
			open OUT,'>',"${outbase}_insoap.sh" or warn "[!]Error opening ${outbase}_insoap.sh: $!\n";
			print OUT "#!/bin/sh\n$sh\n";
			close OUT;
			system($sh)==0 or die "[x]system [$sh] failed: $?";
		}
		my ($Pairs,$Paired,$Singled);
		open LOG,'<',"$outbase.log" or die "[x]Error opening $outbase.log: $!\n";
		while (<LOG>) {
			$Pairs = (split)[-2] if /^Total Pairs:/;
			$Paired = (split)[1] if /^Paired:/;
			$Singled = (split)[1] if /^Singled:/;
		}
		close LOG;
		unless ($Pairs) {
			system("mv -f $outbase.log $outbase.log.0");
			goto TEST;
		}
		my $SEratio=$Singled/(2*$Pairs);
		my $WARN=($SEratio<$WARNvalue)?0:1;
		open SOAP,'<',"$outbase.soap" or die "[x]Error opening $outbase.soap: $!\n";
my $total_map_reads = 0;
my $total_single = 0;
my $total_repeat_pair = 0;
my $total_uniq_pair = 0;
my $total_uniq_low_pair = 0;
my $total_uniq_normal_pair = 0;
my ($insert,%insert,%insertN);
		while (<SOAP>) {
			$total_map_reads+=2;
			my ($id1, $n1, $len1, $f1, $chr1, $x1, $m1) = (split "\t")[0,3,5,6,7,8,-2];
			if (eof SOAP) {
				--$total_map_reads;
				last;
			}
			$_=<SOAP>;
			my ($id2, $n2, $len2, $f2, $chr2, $x2, $m2) = (split "\t")[0,3,5,6,7,8,-2];
			#$id1 =~ s/\/[12]$//;
			#$id2 =~ s/\/[12]$//;
			if ($id1 ne $id2){	#single
				seek (SOAP, -length($_),1);
				--$total_map_reads;
				++$total_single;
			} elsif ($chr1 ne $chr2) {	#single
				++$total_single;
			} elsif ($n1!=1 or $n2!=1) {	#repeat
				++$total_repeat_pair;
			} elsif ($f1 eq '+' && $f2 eq '-') {
				$insert = $x2 - $x1 + $len2;
				$insert{$insert}++;
			} elsif ($f2 eq '+' && $f1 eq '-') {
				$insert = $x1 - $x2 + $len1;
				$insert{$insert}++;
			}
		}
		close SOAP;
		%insertN=%insert;
		my ($value,$peak,$sum0,$delta)=@{&cut(\%insert,\%insertN)};
		my $dropratio=$delta/$sum0;
		my $OK=($dropratio<$OKvalue)?0:1;


		my ($k, $v);
		open O1,'>',"$outbase.o1" or die "[x]Error opening $outbase.o1: $!\n";
		print O1 "$_\t$insert{$_}\n" for ( sort {$a <=> $b} keys %insert);
		close O1;
		open O2,'>',"$outbase.o2" or die "[x]Error opening $outbase.o2: $!\n";
		print O2 "$_\t$insertN{$_}\n" for ( sort {$a <=> $b} keys %insertN);
		close O2;
		open O1,'>',"$outbase.tt" or die "[x]Error opening $outbase.tt: $!\n";
		print O1 "$dropratio = $delta / $sum0\nTotal Pairs:\t$Pairs\nPaired:\t$Paired\nSingled:\t$Singled\nWarn: $WARN\t$SEratio\nOK: $OK\t$dropratio\n";
		close O1;

	}
}

sub cut($$) {
	my ($insertref,$insertNref)=@_;
	my $minvalue=160;	# min Lib insert size is 200 bp
	my @keys=sort {$a <=> $b} keys %$insertref;
	my ($peak,$value,$sum0,$delta)=(0,0);
	for my $k (@keys) {
		my $t=$$insertref{$k};
		$sum0 += $t;
		next if $k < $minvalue;
		if ($t > $value) {
			$value=$t;$peak=$k;
		}
	}
	my $curval=$value;
	for my $k (@keys) {
		next if $k < $peak;
		my $t=$$insertNref{$k};
		if ($t <= $curval) {
			$curval=$t;
		} else {
			$delta += $t-$curval;
			$$insertNref{$k}=$curval;
		}
	}
	$curval=$value;
	for my $k (reverse @keys) {
		next if $k > $peak;
		my $t=$$insertNref{$k};
		if ($t <= $curval) {
			$curval=$t;
		} else {
			$delta += $t-$curval;
			$$insertNref{$k}=$curval;
		}
	}
	return [$value,$peak,$sum0,$delta];
}












sub calerf {
	my ($arg, $result, $jint) = @_;
	local $[ = 1;
#------------------------------------------------------------------
#
# This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
#   for a real argument  x.  It contains three FUNCTION type
#   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
#   and one SUBROUTINE type subprogram, CALERF.  The calling
#   statements for the primary entries are:
#
#                   Y=ERF(X)     (or   Y=DERF(X)),
#
#                   Y=ERFC(X)    (or   Y=DERFC(X)),
#   and
#                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
#
#   The routine  CALERF  is intended for internal packet use only,
#   all computations within the packet being concentrated in this
#   routine.  The function subprograms invoke  CALERF  with the
#   statement
#
#          CALL CALERF(ARG,RESULT,JINT)
#
#   where the parameter usage is as follows
#
#      Function                     Parameters for CALERF
#       call              ARG                  Result          JINT
#
#     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
#     ERFC(ARG)     ABS(ARG) < XBIG           ERFC(ARG)         1
#     ERFCX(ARG)    XNEG < ARG < XMAX         ERFCX(ARG)        2
#
#   The main computation evaluates near-minimax approximations
#   from "Rational Chebyshev approximations for the error function"
#   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
#   transportable program uses rational functions that theoretically
#   approximate  erf(x)  and  erfc(x)  to at least 18 significant
#   decimal digits.  The accuracy achieved depends on the arithmetic
#   system, the compiler, the intrinsic functions, and proper
#   selection of the machine-dependent constants.
#
#*******************************************************************
#*******************************************************************
#
# Explanation of machine-dependent constants
#
#   XMIN   = the smallest positive floating-point number.
#   XINF   = the largest positive finite floating-point number.
#   XNEG   = the largest negative argument acceptable to ERFCX;
#            the negative of the solution to the equation
#            2*exp(x*x) = XINF.
#   XSMALL = argument below which erf(x) may be represented by
#            2*x/sqrt(pi)  and above which  x*x  will not underflow.
#            A conservative value is the largest machine number X
#            such that   1.0 + X = 1.0   to machine precision.
#   XBIG   = largest argument acceptable to ERFC;  solution to
#            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
#            W(x) = exp(-x*x)/[x*sqrt(pi)].
#   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
#            machine precision.  A conservative value is
#            1/[2*sqrt(XSMALL)]
#   XMAX   = largest acceptable argument to ERFCX; the minimum
#            of XINF and 1/[sqrt(pi)*XMIN].
#
#   Approximate values for some important machines are:
#
#                          XMIN       XINF        XNEG     XSMALL
#
#  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
#  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
#  IEEE (IBM/XT,
#    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
#  IEEE (IBM/XT,
#    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
#  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
#  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
#  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
#  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
#
#
#                          XBIG       XHUGE       XMAX
#
#  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
#  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
#  IEEE (IBM/XT,
#    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
#  IEEE (IBM/XT,
#    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
#  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
#  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
#  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
#  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
#
#*******************************************************************
#*******************************************************************
#
# Error returns
#
#  The program returns  ERFC = 0      for  ARG >= XBIG;
#
#                       ERFCX = XINF  for  ARG < XNEG;
#      and
#                       ERFCX = 0     for  ARG >= XMAX.
#
#
# Intrinsic functions required are:
#
#     ABS, AINT, EXP
#
#
#  Author: W. J. Cody
#          Mathematics and Computer Science Division
#          Argonne National Laboratory
#          Argonne, IL 60439
#
#  Latest modification: March 19, 1990
#
# Translation to Perl by Peter J. Acklam, December 3, 1999
#
#------------------------------------------------------------------
	my ($i);
	my ($x, $del, $xden, $xnum, $y, $ysq);
#------------------------------------------------------------------
#  Mathematical constants
#------------------------------------------------------------------
	my ($four, $one, $half, $two, $zero) = (4, 1, 0.5, 2, 0);
	my $sqrpi = 5.6418958354775628695e-1;
	my $thresh = 0.46875;
	my $sixten = 16;
#------------------------------------------------------------------
#  Machine-dependent constants
#------------------------------------------------------------------
	my ($xinf, $xneg, $xsmall) = (1.79e308, -26.628, 1.11e-16);
	my ($xbig, $xhuge, $xmax) = (26.543, 6.71e7, 2.53e307);
#------------------------------------------------------------------
#  Coefficients for approximation to  erf  in first interval
#------------------------------------------------------------------
	my @a = (3.16112374387056560e00, 1.13864154151050156e02,
			 3.77485237685302021e02, 3.20937758913846947e03,
			 1.85777706184603153e-1);
	my @b = (2.36012909523441209e01, 2.44024637934444173e02,
			 1.28261652607737228e03, 2.84423683343917062e03);
#------------------------------------------------------------------
#  Coefficients for approximation to  erfc  in second interval
#------------------------------------------------------------------
	my @c = (5.64188496988670089e-1, 8.88314979438837594e0,
			 6.61191906371416295e01, 2.98635138197400131e02,
			 8.81952221241769090e02, 1.71204761263407058e03,
			 2.05107837782607147e03, 1.23033935479799725e03,
			 2.15311535474403846e-8);
	my @d = (1.57449261107098347e01, 1.17693950891312499e02,
			 5.37181101862009858e02, 1.62138957456669019e03,
			 3.29079923573345963e03, 4.36261909014324716e03,
			 3.43936767414372164e03, 1.23033935480374942e03);
#------------------------------------------------------------------
#  Coefficients for approximation to  erfc  in third interval
#------------------------------------------------------------------
	my @p = (3.05326634961232344e-1, 3.60344899949804439e-1,
			 1.25781726111229246e-1, 1.60837851487422766e-2,
			 6.58749161529837803e-4, 1.63153871373020978e-2);
	my @q = (2.56852019228982242e00, 1.87295284992346047e00,
			 5.27905102951428412e-1, 6.05183413124413191e-2,
			 2.33520497626869185e-3);
#------------------------------------------------------------------
	$x = $arg;
	$y = abs($x);
	if ($y <= $thresh) {
#------------------------------------------------------------------
#  Evaluate  erf  for  |X| <= 0.46875
#------------------------------------------------------------------
		$ysq = $zero;
		if ($y > $xsmall) { $ysq = $y * $y }
		$xnum = $a[5]*$ysq;
		$xden = $ysq;
		for (my $i = 1 ; $i <= 3 ; ++$i) {
			$xnum = ($xnum + $a[$i]) * $ysq;
			$xden = ($xden + $b[$i]) * $ysq;
		}
		$$result = $x * ($xnum + $a[4]) / ($xden + $b[4]);
		if ($jint != 0) { $$result = $one - $$result }
		if ($jint == 2) { $$result = exp($ysq) * $$result }
		goto x800;
#------------------------------------------------------------------
#  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
#------------------------------------------------------------------
	} elsif ($y <= $four) {
		$xnum = $c[9]*$y;
		$xden = $y;
		for (my $i = 1 ; $i <= 7 ; ++$i) {
			$xnum = ($xnum + $c[$i]) * $y;
			$xden = ($xden + $d[$i]) * $y;
		}
		$$result = ($xnum + $c[8]) / ($xden + $d[8]);
		if ($jint != 2) {
			$ysq = int($y*$sixten)/$sixten;
			$del = ($y-$ysq)*($y+$ysq);
			$$result = exp(-$ysq*$ysq) * exp(-$del) * $$result;
		}
#------------------------------------------------------------------
#  Evaluate  erfc  for |X| > 4.0
#------------------------------------------------------------------
	} else {
		$$result = $zero;
		if ($y >= $xbig) {
			if (($jint != 2) || ($y >= $xmax)) { goto x300 }
			if ($y >= $xhuge) {
				$$result = $sqrpi / $y;
				goto x300;
			}
		}
		$ysq = $one / ($y * $y);
		$xnum = $p[6]*$ysq;
		$xden = $ysq;
		for (my $i = 1 ; $i <= 4 ; ++$i) {
			$xnum = ($xnum + $p[$i]) * $ysq;
			$xden = ($xden + $q[$i]) * $ysq;
		}
		$$result = $ysq *($xnum + $p[5]) / ($xden + $q[5]);
		$$result = ($sqrpi -  $$result) / $y;
		if ($jint != 2) {
			$ysq = int($y*$sixten)/$sixten;
			$del = ($y-$ysq)*($y+$ysq);
			$$result = exp(-$ysq*$ysq) * exp(-$del) * $$result;
		}
	}
#------------------------------------------------------------------
#  Fix up for negative argument, erf, etc.
#------------------------------------------------------------------
  x300:
	if ($jint == 0) {
		$$result = ($half - $$result) + $half;
		if ($x < $zero) { $$result = -$$result }
	} elsif ($jint == 1) {
		if ($x < $zero) { $$result = $two - $$result }
	} else {
		if ($x < $zero) {
			if ($x < $xneg) {
				$$result = $xinf;
			} else {
				$ysq = int($x*$sixten)/$sixten;
				$del = ($x-$ysq)*($x+$ysq);
				$y = exp($ysq*$ysq) * exp($del);
				$$result = ($y+$y) - $$result;
			}
		}
	}
  x800:
	return 1;
#---------- Last card of CALERF ----------
}

sub erf {
	my $x = @_ ? $_[0] : $_;
#--------------------------------------------------------------------
#
# This subprogram computes approximate values for erf(x).
#   (see comments heading CALERF).
#
#   Author/date: W. J. Cody, January 8, 1985
#
# Translation to Perl by Peter J. Acklam, December 3, 1999
#
#--------------------------------------------------------------------
	my $result;
	my $jint = 0;
	calerf($x, \$result, $jint);
	my $erf = $result;
	return $erf;
#---------- Last card of ERF ----------
}
