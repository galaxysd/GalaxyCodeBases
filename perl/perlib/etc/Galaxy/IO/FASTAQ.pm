package Galaxy::IO::FASTAQ;
use strict;
use warnings;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(readfq);
our @EXPORT_OK   =qw(readfq getQvaluesFQ);
our $VERSION   = v1.1.0;

sub readfq($$) {
	my ($fh, $aux) = @_;
	$aux = [undef, 0] unless (@$aux);
	return 0 if ($aux->[1]);
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
			return 0;
		}
	}
	#my $name = /^.(\S+)/? $1 : '';
	my ($name,$comment);
	if (/^.(\S+)( (.+)$)?/) {
		$name = $1;
	} else { $name = ''; }
	$comment = $3? $3 : '';
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
	return [$name, $comment, $seq] if ($c ne '+');
	my $qual = '';
	while (<$fh>) {
		chomp;
		$qual .= $_;
		if (length($qual) >= length($seq)) {
			$aux->[0] = undef;
			return [$name, $comment, $seq, $qual];
		}
	}
	$aux->[1] = 1;
	return [$name, $comment, $seq];
}

sub getQvaluesFQ(@) {
	my @Qstr=split //,$_[0];
	my $qbase = $_[1]?$_[1]:33;
	my @Qvalue=();
	push @Qvalue,ord($_)-$qbase for @Qstr;
	return \@Qvalue;
}

1;

__END__

https://github.com/lh3/readfq/blob/master/readfq.pl

use Galaxy::IO::FASTAQ qw(readfq getQvaluesFQ);

my @aux = undef;
my ($name, $comment, $seq, $qual, $ret, @QV);
my ($n, $slen, $qlen) = (0, 0, 0);
while ( $ret = &readfq(\*STDIN, \@aux) ) {
	($name, $comment, $seq, $qual) = @$ret;
	++$n;
	$slen += length($seq);
	if ($qual) {
		$qlen += length($qual);
		@QV = @{getQvaluesFQ($qual)};
	}
	print join("\t", $name, $comment, $seq, $qual, join(',',@QV)), "|\n";
}
print join("\t", $n, $slen, $qlen), "\n";
