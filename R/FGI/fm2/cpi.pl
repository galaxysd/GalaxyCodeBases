$_=<DATA>;
chomp;
my @id = split /\t/,$_;
#print "@id\n";
my $cpi;

sub s2a($) {
	my @x = split '',@_[0];
	my %hash;
	$hash{$_}++ for (@x);
	@ks = sort keys %hash;
	$r = join('',@ks);
	return $r,length($r);
}
my $cpi = 1;
while(<DATA>) {
	chomp;
	s/\///g;
	my ($gF,$gM,$gC) = split /\t/,$_;
	my ($F,$Fl) = s2a($gF);
	my ($M,$Ml) = s2a($gM);
	my ($C,$Cl) = s2a($gC);
	my $pi=10**(-5);
	my $type='NA';
	if ($F eq $M and $M eq $C and $Fl==1) {
		$pi = 2;
		$type = "PPP 1/p";
	} elsif ($F ne $M and $Fl+$Ml==2 and $Cl==2 and $C =~ /$F/ and $C=~/$M/) {
		$pi = 2;
		$type = "P PQ Q 1/q";
	} elsif ($C eq $M and $Cl+$Ml==2 and $Fl==2 and $F=~/$M/) {
		$pi = 1;
		$type = "P P PQ 1/2p";
	} elsif ($Ml==1 and $Cl==2 and $Fl==2 and $C=~/$M/) {
		$pi = 1;
		$type = "P PQ PQ/QR 1/2q";
	} elsif ($Ml==2) {
		if ($Cl==1 and $M=~/$C/) {
			if ($Fl==1 and $F=~/$C/) {
				$pi = 2;
				$type = "PQ QQ QQ 1/q";
			} elsif ($Fl==2 and $F=~/$C/) {
				$pi = 1;
				$type = "PQ QQ QR 1/2q";
			}
		} elsif ($Cl==2) {
			if ($M eq $C) {
				if ($F eq $C or $Fl==1) {
					$pi = 1;
					$type = "PQ PQ PP/QQ/PQ 1/(p+q)";
				} elsif ($Fl==2) {
					$pi = 1/2;
					$type = "PQ PQ PR 1/2(p+q)";
				}
			} else {
				if ($Fl==1 and $F=~/$C/) {
					$pi = 2;
					$type = "PQ QR RR 1/r";
				} elsif ($Fl==2 and $F=~/$C/) {
					$pi = 1;
					$type = "PQ QR RS 1/2r";
				} elsif ($F eq $C) {
					$pi = 1;
					$type = "PQ PR PR 1/2r";
				}
			}
		}
	}
	print "($M,$C,$F)\t$pi=$type\n";
	$cpi *= $pi;
}
print $cpi,"\n";

__DATA__
G2F	G2M	G2S
G/G	G/G	G/G
G/G	A/G	G/G
C/C	G/C	G/C
T/C	C/C	C/C
A/G	A/G	A/G
A/A	A/A	A/A
C/T	C/C	C/C
A/G	G/G	G/G
G/A	G/G	G/G
G/A	G/A	A/A
T/T	T/T	T/T
G/G	A/G	A/G
T/T	T/T	T/T
T/C	T/C	T/C
A/A	T/T	T/A
C/C	T/C	T/C
C/T	C/T	C/T
G/A	A/A	A/A
A/T	T/T	T/T
G/G	G/G	G/G
G/G	A/A	A/G
G/C	C/C	G/C
A/G	A/G	A/G
A/A	A/A	A/A
T/T	T/T	T/T
C/C	C/C	C/C
A/A	A/A	A/A
G/A	G/G	G/G
A/A	T/T	T/A
T/C	T/C	T/T
G/G	G/G	G/G
T/C	T/T	T/C
G/G	G/G	G/G
G/T	G/T	G/G
T/T	T/T	T/T
T/T	C/C	C/T
G/G	C/G	G/G
C/C	C/C	C/C
C/C	C/C	C/C
A/G	A/G	A/A
T/T	T/T	T/T
A/A	A/A	A/A
A/A	A/C	A/C
C/A	C/C	C/C
T/T	C/T	T/T
T/T	G/G	T/T
G/T	G/G	G/T
G/G	G/G	G/G
G/A	G/G	G/G
G/G	G/G	G/G
T/A	T/A	T/T
G/T	G/G	G/G
T/T	T/T	T/T
A/G	A/A	A/A
G/A	G/A	A/A
T/C	T/C	C/C
G/A	G/A	A/A
G/G	G/G	G/G
G/T	T/T	T/T
A/A	A/A	A/A
G/G	T/G	G/G
A/G	A/A	A/A
C/G	G/G	C/G
G/G	G/G	G/G
G/G	G/G	G/G
A/G	G/G	G/G
G/G	G/G	G/G
G/A	A/A	G/A
C/C	G/C	C/C
G/A	A/A	A/A
T/C	T/C	T/T
C/T	C/C	C/C
C/T	C/C	C/C
C/A	A/A	C/A
A/A	A/A	A/A
A/A	A/A	A/A
C/T	C/T	C/C
C/C	C/C	C/C
C/C	A/C	A/C
T/C	C/C	T/C
C/C	C/C	C/C
C/C	C/C	C/C
C/C	C/C	C/C
C/C	C/T	C/T
C/C	T/T	T/C
A/G	G/G	A/G
C/C	A/A	A/A
G/G	A/A	G/A
A/G	A/A	A/G
A/A	C/C	C/A
A/T	T/T	T/T
A/G	A/A	A/G
T/C	T/C	C/C
A/T	A/T	A/T
G/G	C/C	C/G
A/G	A/A	A/G
T/C	T/C	C/C
C/C	T/C	C/C
T/T	T/C	T/T
C/C	C/C	C/C
T/T	T/C	T/T
A/G	A/G	G/G
A/C	A/C	A/A
T/C	T/C	T/T
A/A	G/A	A/A
G/A	G/A	G/A
G/A	G/A	G/A
C/C	C/C	C/C
T/T	T/T	T/T
T/T	T/T	T/T
C/C	T/T	T/C
T/C	T/C	C/C
G/G	G/G	G/G
C/T	T/T	C/T
G/G	C/G	G/G
T/C	C/C	C/C
T/T	T/C	T/C
A/G	A/G	A/A
T/A	A/A	T/A
G/A	G/A	G/A
G/G	A/A	G/A
G/G	G/A	G/A
T/A	T/T	T/A
C/C	C/T	C/C
C/C	C/C	C/C
A/A	A/A	A/A
T/T	T/T	T/T
C/T	T/T	T/T
T/C	T/C	T/T
C/C	C/T	C/T
T/C	T/T	T/T
A/A	A/A	A/A
G/G	A/A	A/G
A/A	A/G	A/G
C/T	T/T	C/T
C/C	T/C	C/C
T/A	A/A	T/A
C/G	C/G	C/G
T/T	T/T	T/T
C/C	C/T	C/C
A/A	T/A	A/A
