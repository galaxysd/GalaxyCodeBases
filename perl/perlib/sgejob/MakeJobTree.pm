package Galaxy::SGE::MakeJobTree;

#use 5.008008;
use strict;
use warnings;
use Exporter 'import';

our @EXPORT_OK = qw(MakeJobTree);
our $VERSION = '0.01';

sub MakeJobTree() { __PACKAGE__ }

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = {
		name => 'batch',
		vf => '1G',
		env => 'PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH',
		sh => '/bin/bash',
		cwd => '1',
		rerun => '1',
		nodotoe => '1',
		cmd => 'sleep 1',
		waiter => '/share/raid010/resequencing/user/huxs/release/tools/waiter.pl',
		marker => '/share/raid010/resequencing/user/huxs/release/tools/marker.pl',
		@_,	# 覆盖以前的属性
	};
	return bless $self, $class;
}

sub clone {
	my $model = shift;
	my $self = $model->new(%$model, @_);
	return $self;	# 前面被 ->new 赐福过了
}

#sub DESTROY {
#	my $self = shift;
#	my $fh = $self->{mailhandle};
#	my $id = $self->{name};
	#print $fh "\n$id is signing off at " . localtime( ) . "\n";
	#close $fh;   # 关闭mailer的管道
#}

for my $field (qw(cmd name vf req sh env cwd rerun nodotoe waiter marker waitopt markopt)) {
	no strict 'refs';	# 这样指向类型团的符号引用就可以用了
	*$field = sub {
		my $self = shift;
		$self->{$field} = shift if @_;
		return $self->{$field};
	};
}











sub as_string {
	my $self = shift;
	my $str;
	for (sort keys %$self) { $str .= "$_=[$$self{$_}], " }
	$str =~ s/, $//;
	return $str;
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Galaxy::SGE::MakeJobSH - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Galaxy::SGE::MakeJobSH;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Galaxy::SGE::MakeJobSH, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.

$cmd can end with STD, ERR, LOG, NONE
default is NONE

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
