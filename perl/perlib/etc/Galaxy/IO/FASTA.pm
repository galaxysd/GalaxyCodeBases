package Galaxy::IO::FASTA;
use strict;
use warnings;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(FastaReadNext);
our @EXPORT_OK   =qw(FastaReadNext);
our $VERSION   = v1.0.0;

sub FastaReadNext($) {
	my $fh = $_[0];
	my ($seqname,$genome);
	while (<$fh>) {
		s/^>//;
		/^(\S+)/ or next;
		$seqname = $1;
		$/=">";
		$genome=<$fh>;
		chomp $genome;
		$genome=~s/\s//g;
		$/="\n";
		return [$seqname,$genome];
	}
	$genome='';
	return 0;
}

1;
