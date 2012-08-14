package Galaxy::SeqTools;
#package main;
use strict;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(revcom);
our @EXPORT_OK   =qw(com translate);
our $VERSION   = v1.0.0;

=head1 Purpose
Provide some common tools.
=cut

#http://doc.bioperl.org/releases/bioperl-1.4/Bio/Tools/SeqPattern.html#CODE6
sub revcom($) {
    my $str = $_[0];
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    my $rev = reverse $str;
    $rev    =~ tr/[](){}<>/][)(}{></;
    return $rev;
}

sub com($) {
    my $str = $_[0];
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    # tr/acgtACGT/tgcaTGCA/
#print "$str\n";
    return $str;
}

sub translate($$) {
	my ($seq,$codonH) = @_;
	$seq =~ s/\s//g;
	my $len = length $seq;
	# get a list of codon start locations
	my @codon_starts =
	  map { 3 * $_ } ( 0 .. ( int($len / 3) - 1 ) );
	my $peptide;
	for (@codon_starts) {
		my $code = substr $seq,$_,3;
		if (exists $codonH->{$code}) {
			$peptide .= $codonH->{$code};
		} else {
			$peptide .= 'x';
			#print "[$code $_]\n";
		}
	}
	return $peptide;
}
