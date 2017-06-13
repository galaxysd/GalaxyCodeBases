package Galaxy::IO::FASTA;
use strict;
use warnings;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(FastaReadNext);
our @EXPORT_OK   =qw(FastaReadNext FastaReadNextA);
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

sub FastaReadNextA($) {
	my $fh = $_[0];
	my ($seqname,$genome,$seqdesc);
	while (<$fh>) {
		s/^>//;
		/^(\S+) ?(.*)?/ or next;
		$seqname = $1;
		$seqdesc = $2?$2:'';
		$/=">";
		$genome=<$fh>;
		chomp $genome;
		$genome=~s/\s//g;
		$/="\n";
		return [$seqname,$genome,$seqdesc];
	}
	$genome='';
	return 0;
}

1;
