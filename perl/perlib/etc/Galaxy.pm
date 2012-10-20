package Galaxy;
use strict;
use 5.010;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(commify colormap);
our @EXPORT_OK   =qw(colormap);
our $VERSION   = v1.0.0;

sub commify {	# http://www.perlmonks.org/?node_id=653
	my $input = shift;
	$input = reverse $input;
	$input =~ s<(\d\d\d)(?=\d)(?!\d*\.)><$1,>g;
	$input = reverse $input;
	return $input;
}

sub colormap {	# http://cresspahl.blogspot.com/2012/03/expanded-control-of-octaves-colormap.html
	my ($value,$MR,$MG,$MB) = @_;
	$MR //= [[0,0],[0.02,0.3],[0.3,1],[1,1]];
	$MG //= [[0,0],[0.3,0],[0.7,1],[1,1]];
	$MB //= [[0,0],[0.7,0],[1,1]];
    # http://www.effectiveperlprogramming.com/blog/704
	my $r = int( 0.5 + 255*_getColorMap($value,$MR) );
	my $g = int( 0.5 + 255*_getColorMap($value,$MG) );
	my $b = int( 0.5 + 255*_getColorMap($value,$MB) );
	my @color = map { sprintf "%02x",$_; } ($r,$g,$b);
	return '#'.join('',@color);
}
sub _getColorMap {
	my ($v,$matrix) = @_;
	my ($oldvalue,$oldout);
	for (@$matrix) {
		my ($value,$out) = @$_;
		if ($value<=$v) {
			$oldvalue = $value;
			$oldout = $out;
			return $out if $v == $value;
		} else {
			my $k = ($out - $oldout)/($value - $oldvalue);
			my $b = $out - $k * $value;
			return $k*$v + $b;
		}
	}
}


1;
