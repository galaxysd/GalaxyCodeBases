use strict;
use warnings;
use File::Basename;
#package main;

# Static Var.
our $DEBUG;
our ($minLen,%Genome,%ChrLen);

# Functions
sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.xz$/) {
		open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

sub getsamChrLen($) {
	my $in = $_[0];
	if ( ($DEBUG - int($DEBUG)) < 0.85 or ($in =~ /HBV\.AJ507799\.2\.fa$/) ) {
		my $GENOME = openfile($in);
		while (<$GENOME>) {
			s/^>//;
			/^(\S+)/ or next;
			my $seqname = $1;
			print STDERR " >$seqname ...";
			$/=">";
			my $genome=<$GENOME>;
			chomp $genome;
			$genome=~s/\s//g;
			$/="\n";
			$Genome{$seqname}=$genome;
			my $thelength = length $Genome{$seqname};
			print STDERR "\b\b\b", $thelength, ".\n";
			$genome='';
			$ChrLen{$seqname} = $thelength;
		}
		close $GENOME;
	}
}

sub revcom($) {
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	my $rev = reverse $str;
	$rev =~ tr/[](){}<>/][)(}{></;
	return $rev;
}

1;
