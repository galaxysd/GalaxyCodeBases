#package Galaxy::ShowHelp;
package main;
use strict;
use Getopt::Std;
#use Exporter 'import';
#our @EXPORT = qw(ShowHelp);
our $VERSION = '0.06';

$Getopt::Std::STANDARD_HELP_VERSION=1;

sub main::ShowHelp() {
#	my ($opt,$help)=@_;
	if (@main::ARGV == 0) {
		&main::VERSION_MESSAGE();
		&main::HELP_MESSAGE();
		die "\n";
		#return 0;
	}
	getopts($main::opts);
}

sub main::HELP_MESSAGE() {
    #my ($scr) = ($0 =~ m,([^/\\]+)$,);
	$main::help =~ s|\[|[\033[0;0m|g;
	$main::help =~ s|\]|\033[32;1m]|g;
	$main::help =~ s|\(|(\033[0;1m|g;
	$main::help =~ s|\)|\033[32;1m)|g;
	$main::help =~ s|:(\s*\n?\s*)(\S)|:$1\033[0;1m$2|g;

	$main::help =~ s|\\\[\033\[0;0m|[|g;
	$main::help =~ s|\\\033\[32;1m\]|]|g;
	$main::help =~ s|\\\(\033\[0;1m|(|g;
	$main::help =~ s|\\\033\[32;1m\)|)|g;
	$main::help =~ s|\\:(\s*\n?\s*)\033\[0;1m|:$1|g;

	$main::help =~ s|\n|\033[32;1m\n|g;
	$main::ARG_DESC='[PROGRAM_ARG1 ...]' unless $main::ARG_DESC;
	print STDERR <<EOH;
\nUsage: \033[0;1m$0\033[0;0m [-OPTIONS [-MORE_OPTIONS]] [--] $main::ARG_DESC

The following single-character options are accepted:
\033[32;1m$main::help\033[0;0mOptions may be merged together.  -- stops processing of options.
Space is not required between options and their arguments.
EOH
}
sub main::VERSION_MESSAGE() {
	my $perlv = $];
	$perlv = sprintf "%vd", $^V if $] >= 5.006;
	my $ver = sprintf "%vd", $main::VERSION;
	my ($scr) = ($0 =~ m,([^/\\]+)$,);
	if ($main::desc) {
		print STDERR <<EOH;
\033[32;1m$main::desc\033[0;0m ($scr) version \033[0;1m$ver\033[0;0m,
 running under Perl version $perlv.
EOH
	} else {
		print STDERR <<EOH;
\033[32;1m$scr\033[0;0m version \033[0;1m$ver\033[0;0m, running under Perl version $perlv.
EOH
	}
}

1;
