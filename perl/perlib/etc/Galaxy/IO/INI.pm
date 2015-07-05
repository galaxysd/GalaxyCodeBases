package Galaxy::IO::INI;
# I need a hash with order.
use 5.004;
use strict;
use warnings;

our $VERSION = '0.02';
$__PACKAGE__::errstr = '';

# Create an empty object
sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = {']' => []};	# ini section can never be ']'
	tie %{$self},'INIHash';
	return bless $self, $class;
}
# Create an object from a file
sub read {
	my $class = shift;

	# Check the file
	my $file = shift or return $class->_error( 'You did not specify a file name' );
	return $class->_error( "File '$file' does not exist" )              unless -e $file;
	return $class->_error( "'$file' is a directory, not a file" )       unless -f _;
	return $class->_error( "Insufficient permissions to read '$file'" ) unless -r _;

	# Slurp in the file
	local $/ = undef;
	open CFG, $file or return $class->_error( "Failed to open file '$file': $!" );
	my $contents = <CFG>;
	close CFG;

	$class->read_string( $contents );
}

# Create an object from a string
sub read_string {
	my $self = shift;

	# Parse the file
	my $ns = '_';
	my $counter = 0;
	foreach ( split /(?:\015{1,2}\012|\015|\012)/, shift ) {
		$counter++;

		# Skip comments and empty lines
		next if /^\s*(?:\#|\;|$)/;

		# Remove inline comments
		s/\s\;\s.+$//g;
#print "$_\n";
		# Handle section headers
		if ( /^\s*\[\s*(.+?)\s*\]\s*$/ ) {
			# Create the sub-hash if it doesn't exist.
			# Without this sections without keys will not
			# appear at all in the completed struct.
			$self->{$ns = $1} ||= {'=' => []};	# ini key can never be '='
			push @{$$self{']'}},$ns unless exists $$self{$ns};
			next;
		}

		# Handle properties
		if ( /^\s*([^=]+?)\s*=\s*(.*?)\s*$/ ) {
			push @{$$self{$ns}{'='}},$1 unless exists $$self{$ns}{$1};
			$self->{$ns}->{$1} = $2;
			next;
		}

		return $self->_error( "Syntax error at line $counter: '$_'" );
	}
	return $self;
}

# Save an object to a file
sub write {
	my $self = shift;
	my $file = shift or return $self->_error(
		'No file name provided'
		);

	# Write it to the file
	open( CFG, '>' . $file ) or return $self->_error(
		"Failed to open file '$file' for writing: $!"
		);
	print CFG $self->write_string;
	close CFG;
}

# Save an object to a string
sub write_string {
	my $self = shift;
	my $contents = '';
	foreach my $section (@{$$self{']'}}) {
		my $block = $self->{$section};
		$contents .= "\n" if length $contents;
		$contents .= "[$section]\n" unless $section eq '_';
		my %Properties = map { $_ => 1 } keys %{$block};
		delete $Properties{'='};
		foreach my $property ( @{$$block{'='}} ) {
			$contents .= "$property=$block->{$property}\n";
			delete $Properties{$property};
		}
		foreach my $property ( sort keys %Properties ) {
			$contents .= "$property=$block->{$property}\n";
		}
	}
	$contents;
}

# Error handling
sub errstr { $__PACKAGE__::errstr }
sub _error { $__PACKAGE__::errstr = $_[1]; undef }

1;


=pod

=head1 NAME

Galaxy::IO::INI - Read/Write .ini style files with as little code as possible

=head1 SYNOPSIS

    # In your configuration file
    rootproperty=blah

    [section]
    one=twp
    three= four
    Foo =Bar
    empty=

    # In your program
    use Galaxy::IO::INI;

    # Create a config
    my $Config = Galaxy::IO::INI->new();

    # Open the config
    $Config->read( 'file.conf' );

    # Reading properties
    my $rootproperty = $Config->{_}->{rootproperty};
    my $one = $Config->{section}->{one};
    my $Foo = $Config->{section}->{Foo};

    # Changing data
    $Config->{newsection} = { this => 'that' }; # Add a section
    $Config->{section}->{Foo} = 'Not Bar!';     # Change a value
    delete $Config->{_};                        # Delete a value or section

    # Save a config
    $Config->write( 'file.conf' );

=head1 DESCRIPTION

C<Galaxy::IO::INI> is a perl class to read and write .ini style configuration
files with as little code as possible, reducing load time and memory
overhead. Most of the time it is accepted that Perl applications use a lot
of memory and modules. The C<::Tiny> family of modules is specifically
intended to provide an ultralight alternative to the standard modules.

This module is primarily for reading human written files, and anything we
write shouldn't need to have documentation/comments. If you need something
with more power move up to L<Config::Simple>, L<Config::General> or one of
the many other C<Config::> modules. To rephrase, L<Config::Tiny> does B<not>
preserve your comments, whitespace, or the order of your config file.

Well, L<Config::Tiny> WILL preserve the order of sections and keys in your
config file. That is why I write this.

The order stores as @{$Config->{']'}} and @{$Config->{section}->{'='}}, so
do NOT parse them as hash ref.

=head1 CONFIGURATION FILE SYNTAX

Files are the same format as for windows .ini files. For example:

	[section]
	var1=value1
	var2=value2

If a property is outside of a section at the beginning of a file, it will
be assigned to the C<"root section">, available at C<$Config-E<gt>{_}>.

Lines starting with C<'#'> or C<';'> are considered comments and ignored,
as are blank lines.

When writing back to the config file, all comments, custom whitespace,
and the ordering of your config file elements is discarded. If you need
to keep the human elements of a config when writing back, upgrade to
something better, this module is not for you.

=head1 METHODS

=head2 new

The constructor C<new> creates and returns an empty C<Galaxy::IO::INI> object.

=head2 read $filename

The C<read> constructor reads a config file, and returns a new
C<Galaxy::IO::INI> object containing the properties in the file.

Returns the object on success, or C<undef> on error.

When C<read> fails, C<Galaxy::IO::INI> sets an error message internally
you can recover via C<<Galaxy::IO::INI->errstr>>. Although in B<some>
cases a failed C<read> will also set the operating system error
variable C<$!>, not all errors do and you should not rely on using
the C<$!> variable.

=head2 read_string $string;

The C<read_string> method takes as argument the contents of a config file
as a string and returns the C<Galaxy::IO::INI> object for it.

=head2 write $filename

The C<write> method generates the file content for the properties, and
writes it to disk to the filename specified.

Returns true on success or C<undef> on error.

=head2 write_string

Generates the file content for the object and returns it as a string.

=head2 errstr

When an error occurs, you can retrieve the error message either from the
C<$Galaxy::IO::INI::errstr> variable, or using the C<errstr()> method.

=head2 property_string

This method is called to produce the string used to represent the property in a
section.  It is passed the section name and property name.

=head2 set

This is a convenience is called to set a value found in the parsed config string.  It is
passed the section name, property name, and value.

=head1 SUPPORT

Bugs should be reported the author.

=head1 AUTHOR

Adam Kennedy E<lt>galaxy001@gmail.comE<gt>

=head1 ACKNOWLEGEMENTS

This package is based on Adam Kennedy's Config::Tiny v2.12

Thanks to Sherzod Ruzmetov E<lt>sherzodr@cpan.orgE<gt> for
L<Config::Simple>, which inspired this module by being not quite
"simple" enough for me :)

=head1 SEE ALSO

L<Config::Tiny>, L<Config::Simple>, L<Config::General>, L<ali.as>

http://search.cpan.org/~rsavage/Config-Tiny/lib/Config/Tiny.pm

=head1 COPYRIGHT

Copyright by Simba Galaxy.

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

=cut

package INIHash;
use Carp;
require Tie::Hash;

@INIHash::ISA = qw(Tie::StdHash);

sub STORE {
	if ($_[1] eq ']') {
		carp "[!]INI section can never be ']'";
		return;
	}
	#$_[0]->{$_[1]} = $_[2];
	push @{$_[0]->{']'}},$_[1] unless exists $_[0]->{$_[1]};
	for (keys %{$_[2]}) {
		next if $_ eq '=';
		push @{$_[0]->{$_[1]}->{'='}},$_ unless exists $_[0]->{$_[1]}->{$_};
		$_[0]->{$_[1]}->{$_}=$_[2]->{$_};
	}
	$_[0]->{$_[1]}->{'='};	# Why ?
}

1;

__END__
