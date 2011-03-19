package Galaxy::SGE::DoneMark;

#use 5.008008;
use strict;
use warnings;
use File::Basename;
use Exporter 'import';

our @EXPORT_OK = qw(MakeJobSH);
our @EXPORT = qw(transname);

our $VERSION = '0.01';

sub transname($) {
	my $str=$_[0];
	return '.' if $str !~ /\w/;
	my ($filename, $directories, $suffix) = fileparse($str, qr/\.[^.]*/);
	$filename=join('_',split(//,$filename));
	$filename =~ s/_\._/\./g;
	$suffix=join('-',split(//,$suffix));
	$suffix =~ s/\.-/\./;
	my $retstr=$directories.'.d-'.$filename.$suffix;
	return $retstr;
}

# Preloaded methods go here.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Galaxy::SGE::DoneMark - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Galaxy::SGE::DoneMark;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Galaxy::SGE::DoneMark, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

A. U. Thor, E<lt>huxuesong@localE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by A. U. Thor

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.


=cut
