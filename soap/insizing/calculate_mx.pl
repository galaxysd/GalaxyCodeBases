#!/usr/bin/perl -w
use strict;

unless(@ARGV>0)
{
	print "perl $0 <soap file name> <insert size result> <frequence statistics> <drop rate(0~1)>\n";
	exit 0;
}

my ($sum,$count,$sum_2,@insertdata,%f);
my $soapfile = $ARGV[0];
my $insert = $ARGV[1];
my $frequence = $ARGV[2];
my $drop = $ARGV[3];
#open SOAP, "<", "$soapfile";
open( SOAP,"-|","gzip -dc $soapfile");
open INSERT, ">", "$insert";
open FREQUENCE, ">", "$frequence";

my $prograss;
my $filetotal = `ls -l $soapfile`;
my @tmp = split/\s+/,$filetotal;
$filetotal = $tmp[4];

print "\nProcessing......\n\n";
while(my $line1 = <SOAP>)
{
	my @arry1 = split/\s+/,$line1;
	my $line2 = <SOAP>;
	my @arry2 = split/\s+/,$line2;
	unless($arry1[0] eq $arry2[0] && $arry1[7] eq $arry2[7])
	{
		$line2 = <SOAP>;
		@arry2 = split/\s+/,$line2;
	}
	if($arry1[6] eq '-')
	{
		my $result1 = $arry1[8]+$arry1[5]-$arry2[8];
		print INSERT $result1,"\n";
		$sum += $result1;
		$sum_2 += ($result1**2);
		$f{$result1}++;
	}
	else
	{
		my $result2 = $arry2[8]+$arry2[5]-$arry1[8];
		print INSERT $result2,"\n";
		$sum += $result2;
		$sum_2 += ($result2**2);
		$f{$result2}++;
	}
	$count++;
	if($count%640 == 0)
	{
		&prograss_bar((tell SOAP),$filetotal);
	}
}
close INSERT;

sub prograss_bar
{
	my ($i,$total) = ($_[0],$_[1]);
	local $| = 1;
	print "\r [".("=" x int(($i/$total)*25)).">".(" " x (24 - int(($i/$total)*25)))."]";
	printf("%2.1f %%",$i/$total*100);
	local $| = 0;
}

my @maxfr = (0,0);
foreach (sort keys %f)
{
	if($f{$_} > $maxfr[1])
	{
		$maxfr[0] = $_;
		$maxfr[1] = $f{$_};
	}
	print FREQUENCE $_,"\t";
	print FREQUENCE $f{$_},"\n";
}

print "\n\n";
print "Process Done!\n\n";
my $average = $sum/$count;
print "The average value is : ", $average,".\n";
print "The number of sample is : ", $count,".\n";
my $sd = sqrt(($sum_2-$count*($average**2))/($count-1));
print "The SD value is : ", $sd,".\n\n";

print "The Max Frequence Value is : ", $maxfr[0], " ", $maxfr[1], "\n";
#print "Now Please Input the presentage you want to drop (0~1): ";
#chomp(my $drop = <STDIN>);
my $k = 0;
$sum = $sum_2 = $count = 0;
open INSERT, "<", "$insert";
while(<INSERT>)
{
	chomp;
	$insertdata[$k] = $_;
	$k++;
}
foreach (sort keys %f)
{
	if($f{$_} <= $maxfr[1]*$drop)
	{
		delete $f{$_};
	}
}
foreach (@insertdata)
{
	if(exists $f{$_})
	{
		$sum += $_;
		$sum_2 += ($_**2);
		$count++;
	}
}

sub gam1{
	my $x1 = $_[0];
	my $t = 0;
	my @a = (0.0000677106,-0.0003442342,0.0015397681,-0.0024467480,
					 0.0109736958,-0.0002109075,0.0742379071,0.0815782188,
					 0.4118402518,0.4227843370,1.0);
	unless($x1 > 0)
	{
		print "x must greater than 0!\n";
		exit 0;
	}
	if($x1 < 1)
	{
		$t = 1/($x1*($x1+1));
		$x1 += 2;
	}
	elsif($x1 <= 2)
	{
		$t = 1/$x1;
		$x1 += 1;
	}
	elsif($x1 <= 3)
	{
		$t = 1;
	}
	else
	{
		$t = 1;
		while($x1 > 3)
		{
			$x1 -= 1;
			$t *= $x1;
		}
	}
	my $s1 = $a[0];
	my $u = $x1-2;
	foreach (@a)
	{
		$s1 = $s1*$u+$_;
	}
	$s1 *= $t;
}

sub gam2
{
	my ($a_gam,$x_gam) = ($_[0],$_[1]);
  my ($n_gam,$p_gam,$q_gam,$d_gam,$s_gam,$s1_gam,$p0_gam,$q0_gam,$p1_gam,$q1_gam,$qq_gam);
  if(($a_gam <= 0)||($x_gam < 0))
  {
  	if($a_gam <= 0)
  	{
  		print "err**a<=0!\n";
  	}
    if($x_gam < 0)
    {
    	print "err**x<0!\n";
    }
    exit 0;
	}
  if($x_gam == 0)
  {
  	return 0;
  }
  if($x_gam > 1e35)
  {
  	return 1;
  }
  $q_gam = $a_gam*log($x_gam);
  $qq_gam = exp($q_gam);
  if($x_gam < 1+$a_gam)
  {
  	$p_gam = $a_gam;
  	$s_gam = $d_gam = 1/$a_gam;
  	for($n_gam = 1; $n_gam <= 100; $n_gam++)
    {
    	$p_gam++;
    	$d_gam = $d_gam*$x_gam/$p_gam;
    	$s_gam += $d_gam;
	    if(abs($d_gam) < abs($s_gam)*1e-7)
      {
      	$s_gam=$s_gam*exp(-$x_gam)*$qq_gam/&gam1($a_gam);
        return $s_gam;
      }
    }
  }
  else
  {
  	$s_gam = 1/$x_gam;
  	$p0_gam = 0;
  	$q0_gam = $p1_gam=1;
  	$q1_gam = $x_gam;
  	for($n_gam = 1; $n_gam <= 100; $n_gam++)
  	{
  		$p0_gam = $p1_gam+($n_gam-$a_gam)*$p0_gam;
  		$q0_gam = $q1_gam+($n_gam-$a_gam)*$q0_gam;
  		$p_gam = $x_gam*$p0_gam+$n_gam*$p1_gam;
  		$q_gam = $x_gam*$q0_gam+$n_gam*$q1_gam;
  		if(abs($q_gam) != 0)
  		{
  			$s1_gam = $p_gam/$q_gam;
  			$p1_gam = $p_gam;
  			$q1_gam = $q_gam;
  			if(abs(($s1_gam-$s_gam)/$s1_gam) < 1e-7)
  			{
  				$s_gam = $s1_gam*exp(-$x_gam)*$qq_gam/&gam1($a_gam);
  				return (1-$s_gam);
  			}
  			$s_gam = $s1_gam;
  		}
  		$p1_gam = $p_gam;
  		$q1_gam = $q_gam;
  	}
  }
  print "a too large !\n";
  $s_gam = 1-$s_gam*exp(-$x_gam)*$qq_gam/&gam1($a_gam);
}

sub erf
{
	my $x_erf = $_[0];
	my $y_erf;
	if($x_erf >= 0)
	{
		$y_erf=&gam2(0.5,$x_erf**2);
	}
	else
	{
		$y_erf=-&gam2(0.5,$x_erf);
	}
}

my $average_drop = $sum/$count;
print "The average value after modification is : ", $average_drop,".\n";
print "The number of sample after modification is : ", $count,".\n";
my $pa = 3.1415926536;
my $variance = ($sum_2-$count*($average_drop**2))/($count-1);
print "The SD value after modification is : ",sqrt($variance), "\n";
my $sd_drop = sqrt($variance*&erf(sqrt((-1)*log($drop)))/(&erf(sqrt((-1)*log($drop)))-2*$drop*sqrt((-1)*log($drop))/sqrt($pa)));
print "The estimated SD value based on the modified sample is : ", $sd_drop,".\n";
close SOAP;
close INSERT;