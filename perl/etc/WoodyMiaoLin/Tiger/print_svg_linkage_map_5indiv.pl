#!/usr/bin/perl
use strict;
use warnings;

die "Usage: $0 <input_bcf> <scaffold> <output_svg>\n" if @ARGV<3;

my @samples = qw/ JHH001 BHX011 BHX019 GZXJ36 GZXJ37 /;
my $ns = @samples;
my $in = shift;
my $sc = shift;
my $out = shift;

open I, "-|", "bcftools view -I $in";
open O, ">", "$out";

# Read genotype from bcf file, and push into @gt;
my @gt;
while (<I>) {
	next unless /^$sc\t/;
	chomp;
	my @line = split /\t/;
	splice @line, -17, 8;
	splice @line, -7, 4;
	splice @line, -1, 1;
	next if $line[5] < 20;
	$line[7] =~ /;FQ=([0-9\-.]+)/;
	next unless $1 > 0;
	my $a = 0;
	my @a;
	foreach (9 .. ($ns+8)) {
		my @sample = split /:/, $line[$_];
		if ($sample[4] >= 20 and $sample[2] >= 1) {
			++$a;
			push @a, $sample[0];
		}
	}
	next if $a < $ns;
	next if (!grep {$_ ne $a[0]} @a);
	my @b;
	$b[0] = $line[0];
	$b[1] = $line[1];
	my $b0 = 0;
	my $b1 = 0;
	foreach (0 .. ($ns-1)) {
		my ($a1, $a2) = split /\//, $a[$_];
		if (($a1 eq "0") and ($a2 eq "0")) {
			++$b0;
			push @b, $line[3];
		} elsif (($a1 eq "1") and ($a2 eq "1")) {
			++$b1;
			push @b, $line[4];
		} else {
			push @b, "$line[3]$line[4]";
		}
	}
	if ($b0 > $b1) {
		push @b, $line[3];
	} elsif ($b0 < $b1) {
		push @b, $line[4];
	} else {
		push @b, "equal";
	}
	push @gt, \@b;
}
close I;

# SVG attribute definitions;
my $w = 30*$ns + 130;
my $h = 10*@gt + 20;
print O "<?xml version=\"1.0\"?>\n";
print O "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$w\" height=\"$h\">\n";
print O "<g transform=\"translate(10,10)\">\n\n";

# Print title;
my $x0 = 112;
foreach (@samples) {	
	print O "<text x=\"$x0\" y=\"0\" font-family=\"Courier\" font-size=\"7\" fill=\"black\">$_</text>\n";
	$x0 += 30;
}

# Print genotype;
my $yy;
foreach my $a (@gt) {
	++$yy;
	my $y = 10*$yy;
	print O "<text x=\"0\" y=\"$y\" font-family=\"Courier\" font-size=\"8\" fill=\"black\">${$a}[0]</text>\n";
	print O "<text x=\"70\" y=\"$y\" font-family=\"Courier\" font-size=\"8\" fill=\"black\">${$a}[1]</text>\n";
	foreach (2 .. ($ns+1)) {
		my $xt = 30*$_ + 60;
		my $yt = 10*$yy;
		my $xr = 30*$_ + 50;
		my $yr = 10*$yy - 8;
		my $h = length ${$a}[$_];
		my $color;
		if ($h == 2) {
			$color = "green";
		} elsif ($h == 1) {
			if (${$a}[$ns+2] eq "equal") {
				$color = "pink";
			} elsif (${$a}[$ns+2] eq ${$a}[$_]) {
				$color = "yellow";
			} else {
				$color = "blue";
			}
		} else {
			$color = "white";
		}
	print O "<rect x=\"$xr\" y=\"$yr\" width=\"30\" height=\"10\" fill=\"$color\" stroke-width=\"0.5\" stroke=\"black\"/>\n";
	print O "<text x=\"$xt\" y=\"$yt\" font-family=\"Courier\" font-size=\"8\" fill=\"black\">${$a}[$_]</text>\n";    
	}
}
print O "\n</g>\n</svg>\n";
close O;
