package Galaxy::IO::FASTAQ;
use strict;
use warnings;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(readfq);
our @EXPORT_OK   =qw();
our $VERSION   = v1.0.0;

sub readfq($$) {
	my ($fh, $aux) = @_;
	@$aux = [undef, 0] if (!defined(@$aux));
	return if ($aux->[1]);
	if (!defined($aux->[0])) {
		while (<$fh>) {
			chomp;
			if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
				$aux->[0] = $_;
				last;
			}
		}
		if (!defined($aux->[0])) {
			$aux->[1] = 1;
			return;
		}
	}
	my $name = /^.(\S+)/? $1 : '';
	my $seq = '';
	my $c;
	$aux->[0] = undef;
	while (<$fh>) {
		chomp;
		$c = substr($_, 0, 1);
		last if ($c eq '>' || $c eq '@' || $c eq '+');
		$seq .= $_;
	}
	$aux->[0] = $_;
	$aux->[1] = 1 if (!defined($aux->[0]));
	return ($name, $seq) if ($c ne '+');
	my $qual = '';
	while (<$fh>) {
		chomp;
		$qual .= $_;
		if (length($qual) >= length($seq)) {
			$aux->[0] = undef;
			return ($name, $seq, $qual);
		}
	}
	$aux->[1] = 1;
	return ($name, $seq);
}

1;

__END__

https://github.com/lh3/readfq/blob/master/readfq.pl

my @aux = undef;
my ($name, $seq, $qual);
my ($n, $slen, $qlen) = (0, 0, 0);
while (($name, $seq, $qual) = readfq(\*STDIN, \@aux)) {
	++$n;
	$slen += length($seq);
	$qlen += length($qual) if ($qual);
}
print join("\t", $n, $slen, $qlen), "\n";
