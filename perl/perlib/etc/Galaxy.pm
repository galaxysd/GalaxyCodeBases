package Galaxy;
use strict;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(commify);
our @EXPORT_OK   =qw();
our $VERSION   = v1.0.0;

sub commify {	# http://www.perlmonks.org/?node_id=653
	my $input = shift;
	$input = reverse $input;
	$input =~ s<(\d\d\d)(?=\d)(?!\d*\.)><$1,>g;
	$input = reverse $input;
	return $input;
}
