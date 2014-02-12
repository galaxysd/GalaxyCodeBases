package Galaxy::ChromString;
#package main;
use strict;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(ChrStrInit ChrStrGetBase ChrStrSetBase ChrStrANDBase ChrStrORBase);
our @EXPORT_OK   =qw();
our $VERSION   = v1.0.0;

=head1 Purpose
Provide some common tools.
=cut

sub ChrStrInit($$) {
	my ($aChr,$aLen)=@_;
	my (%ChrLen,%ChrStr);
	@ChrLen{@$aChr} = @$aLen;
	for (@$aChr) {
		$ChrStr{$_} = "\0" x $ChrLen{$_};
	}
	return \%ChrStr;
}

sub ChrStrSetBase($$$$) {
	my ($aChrStr,$Chr,$Pos,$Value) = @_;
	my $v = pack('C',$Value);
	substr $$aChrStr{$Chr}, $Pos-1, 1, $v;
	return $v;
}

sub ChrStrANDBase($$$$) {
	my ($aChrStr,$Chr,$Pos,$Value) = @_;
	my $str = substr $$aChrStr{$Chr}, $Pos-1, 1;
	my $v = unpack('C',$str);
	$v &= $Value;
	$v = pack('C',$v);
	substr $$aChrStr{$Chr}, $Pos-1, 1, $v;
	return $v;
}
sub ChrStrORBase($$$$) {
	my ($aChrStr,$Chr,$Pos,$Value) = @_;
	my $str = substr $$aChrStr{$Chr}, $Pos-1, 1;
	my $v = unpack('C',$str);
	$v |= $Value;
	$v = pack('C',$v);
	substr $$aChrStr{$Chr}, $Pos-1, 1, $v;
	return $v;
}

sub ChrStrGetBase($$$) {
	my ($aChrStr,$Chr,$Pos) = @_;
	my $str = substr $$aChrStr{$Chr}, $Pos-1, 1;
	return unpack 'C',$str;
}

1;
