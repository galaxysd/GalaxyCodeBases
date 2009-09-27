package Galaxy::SGE::MakeJobSH;

#use 5.008008;
use strict;
use warnings;
use Exporter 'import';

our @EXPORT_OK = qw(MakeJobSH);
our $VERSION = '0.06';

#INIT {
#	our $Shell=$ENV{SHELL};	# fix if you are in a prefix env., which breaks `pwd`
#	#print "[$Shell] got.\n";
#}
#our $Shell;	# Well, I prefer $ sh /path_to_file/do.sh. Thus shebang not needed.

sub MakeJobSH() { __PACKAGE__ }

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

sub format {
	my $self = shift;
	my %translate=(
		name => ['aq','-N '],
		vf => ['a','-l vf='],
		req => ['a','-l '],
		sh => ['a','-S '],
		env => ['a','-v '],
		cwd => ['b','-cwd'],
		rerun => ['b','-r y'],
		waitopt => ['m','waiter'],	# 'waiter' is the key of %self
		markopt => ['m','marker'],
		nodotoe => ['b','-o /dev/null -e /dev/null']
	);

	my $fixhead=<<'EOH';
ARC=lx26-amd64
QSUB=$SGE_ROOT/bin/$ARC/qsub
SLEEP=20

cmd="$QSUB $0 $0"

# started by SGE or manually
if [ "$JOB_ID" = "" ]; then
	echo "submitting $0 ..."
	$cmd
	while [ "x$?" != "x0" ]; do
	echo "$0: qsub failed - retrying .." >&2
		sleep $SLEEP
		$cmd
	done
else
	if [ "$1" = "" ]; then
		MAIN=${JOB_NAME}_${JOB_ID}${TASK_ID}
	else
		MAIN=$1
	fi
	if [ "$RESTARTED" = "0" ]; then
		echo \# Begin  @ `date` >${MAIN}.err
		cat /dev/null > ${MAIN}.out
	else
		echo \#-------------------- >>${MAIN}.err
		echo \#-------------------- >>${MAIN}.out
		echo \#Restart @ `date` [$RESTARTED]>>${MAIN}.err
	fi
	echo \#$ENVIRONMENT $JOB_NAME of $QUEUE @ Host:$HOSTNAME as Job:[$JOB_ID],Task:[$TASK_ID] >>${MAIN}.err
##### job starts #####
EOH
# well, ends with "\n"

my $fixtail=<<'EOT';	# well, starts with "\n"

##### job  ends  #####
	ENDVALUE=$?
	cat <<EOFSTAT >> ${MAIN}.err
#  End  @ `date`
# Used $SECONDS sec.
#Job ended as $ENDVALUE
#\$PWD is $PWD

#\$SGE_O_WORKDIR is [$SGE_O_WORKDIR]
#PATH is [$PATH]
#LD_LIBRARY_PATH is [$LD_LIBRARY_PATH]
#PERL5LIB is [$PERL5LIB]
#PYTHONPATH is [$PYTHONPATH]

#`qstat -j $JOB_ID | grep usage`
#I am [$0]
#I was [$@]
#All done !
EOFSTAT
fi
EOT

	my ($file,$cmd,$str,$key,%DMcmd)=("#!/bin/sh\n");	# Global symbol "$file" requires explicit package name
	#$file="#!/bin/env bash\n";
	#$file="#!$Shell\n";
	my %Functions=(
		a => sub {$file .= '#$ '.$str.$_[0]."\n";},
		aq => sub {$_[0] =~ s/\"//g;$file .= '#$ '.$str.'"'.$_[0]."\"\n";},
		#as => sub {$_[0] =~ s/\'/\\\"/g;$file .= '#$ '.$str.$_[0]."\n";},	# -v cannot contail ['] in SGE 6.0u8
		b => sub {$file .= '#$ '.$str."\n" if $_[0] !~ /^[0fn]$/i;},
		m => sub {
				return unless $_[0];
				chomp $_[0];
				if ($str eq 'marker') {$DMcmd{$str} .= "\n### $str starts ###\n";}
				$DMcmd{$str} .= $self->{$str};
				$DMcmd{$str} .= ' '.$_[0].' 1>&2 2>>${MAIN}.err';	# ERR4ALL
				if ($str eq 'waiter') {$DMcmd{$str} .= "\n### $str ends ###\n";}
				},
	);
	no strict 'refs';	# better to be outside of the cycle
	for $key (keys %$self) {
		if (defined $translate{$key}) {	# will exists faster than defined ?
			($cmd,$str)=@{$translate{$key}};
		} else { next }
		if ($Functions{$cmd}) {$Functions{$cmd}->($self->{$key})}
		 else { warn "Operate:$cmd not defined." }	# never happens
	}
	my (@newcmd,$thecmd);	# empty each time, ready to be push in
	for (split /\n/,$$self{cmd}) {
		next if /^\s*$/;
		if (s/ LOG$//) {
			$thecmd=$_ . ' >>${MAIN}.out 2>>${MAIN}.err';
			push @newcmd,$thecmd;
		} elsif (s/ ERR$//) {
			$thecmd=$_ . ' 2>>${MAIN}.err';
			push @newcmd,$thecmd;
		} elsif (s/ ERR4ALL$//) {
			$thecmd=$_ . ' 1>&2 2>>${MAIN}.err';
			push @newcmd,$thecmd;
		} elsif (s/ STD4ALL$//) {
			$thecmd=$_ . ' 2>&1 >>${MAIN}.out';
			push @newcmd,$thecmd;
		} elsif (s/ STD$//) {
			$thecmd=$_ . ' >>${MAIN}.out';
			push @newcmd,$thecmd;
		} elsif (s/ NONE$//) {
			push @newcmd,$_;
		} else { push @newcmd,$_; }
	}
	$thecmd = join "\n",@newcmd;
	$thecmd = $DMcmd{waiter}.$thecmd if $DMcmd{waiter};
	$thecmd .= $DMcmd{marker} if $DMcmd{marker};
	return $file.$fixhead.$thecmd.$fixtail;
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
