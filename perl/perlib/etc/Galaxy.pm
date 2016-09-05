package Galaxy;
use strict;
use 5.010;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(commify colormap alphanum);
our @EXPORT_OK   =qw(colormap);
our $VERSION   = v1.0.0;

sub commify {	# http://www.perlmonks.org/?node_id=653 给数字加逗号
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

# http://www.davekoelle.com/alphanum.html
# usage:
#my @sorted = sort { alphanum($a,$b) } @strings;
sub alphanum {
  # split strings into chunks
  my @a = chunkify($_[0]);
  my @b = chunkify($_[1]);
  # while we have chunks to compare.
  while (@a && @b) {
    my $a_chunk = shift @a;
    my $b_chunk = shift @b;
    my $test =
        (($a_chunk =~ /\d/) && ($b_chunk =~ /\d/)) ? # if both are numeric
            $a_chunk <=> $b_chunk : # compare as numbers
            $a_chunk cmp $b_chunk ; # else compare as strings
    # return comparison if not equal.
    return $test if $test != 0;
  }
  # return longer string.
  return @a <=> @b;
}
# split on numeric/non-numeric transitions
sub chunkify {
  my @chunks = split m{ # split on
    (?= # zero width
      (?<=\D)\d | # digit preceded by a non-digit OR
      (?<=\d)\D # non-digit preceded by a digit
    )
  }x, $_[0];
  return @chunks;
}

1;

__END__

# Tips:

## http://stackoverflow.com/questions/2957879/perl-map-need-to-map-an-array-into-a-hash-as-arrayelement-array-index
%hash = map { $arr[$_] => $_ } 0..$#arr;

## http://stackoverflow.com/questions/2700302/how-do-i-get-a-slice-from-an-array-reference/2702033#2702033
my @slice =   @   array   [1,3,2];
my @slice =   @ { $aref } [1,3,2];
my $slice_ref = [ @$aref[1,3,2] ];
