package Galaxy::IO;
#package main;
use strict;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(openfile openpipe);
our @EXPORT_OK   =qw(openfile openpipe);
our $VERSION   = v1.0.0;

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

sub openpipe($$) {
	my ($cmd,$filename)=@_;
	my $infile;	# undef
	if ( length($filename) ) {
		open( $infile,"-|","$cmd $filename") or die "Error opening [$cmd],[$filename]: $!\n";
	}
	return $infile;
}
