package Galaxy::IO::FASTA;
use strict;
use warnings;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(readwholefa);
our @EXPORT_OK   =qw(readwholefa);
our $VERSION   = v1.0.0;

sub readwholefa($) {
	my $fh = $_[0];
	my %Genome;
	while (<$fh>) {
		s/^>//;
		/^(\S+)/ or next;
		my $seqname = $1;
		print STDERR " >$seqname ...";
		$/=">";
		my $genome=<$fh>;
		chomp $genome;
		$genome=~s/\s//g;
		$/="\n";
		$Genome{$seqname}=$genome;
		print STDERR "\b\b\b",length $Genome{$seqname},".\n";
		$genome='';
	}
	return \%Genome;
}

1;
