package Galaxy::SeqTools;
#package main;
use strict;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(com revcom);
our @EXPORT_OK   =qw();
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
