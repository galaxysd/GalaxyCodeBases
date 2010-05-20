#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use File::Basename;
use Galaxy::ShowHelp;
use FindBin qw($Bin);

$main::VERSION=0.0.1;
my $SCRIPTS='$Bin/../scripts';

=pod
Note:
If a PE Lib come with only one .adapter.list file, all will use that one. If more than 2, only (sort)[0,-1] is used as _1 and _2.

./0rawfq , the path to fq file cannot contain more than 1 "_{$LibName}", since it is searched directly on fullpath.

$ perl -MTree::Suffix -e 'my $tree = Tree::Suffix->new(qw(stringxxxssx stringyx1xxssx axxxssxstring));my @lcs = $tree->lcs;print "[$_]\n" for @lcs;'
[string]

$ perl -MTree::Suffix -e 'my $tree = Tree::Suffix->new(qw(zzzzzzxxaaaaaastringaxxssx tyaaaaaastringzzzzzzyaxxssx zzzzzz1aaaaaaaaxxssxstring));my @lcs = $tree->lcs;print "[$_]\n" for @lcs;'
[zzzzzz]
[aaaaaa]
[axxssx]
[string]

=cut
our $opts='i:o:s:bvqd';
our($opt_i, $opt_s, $opt_o, $opt_v, $opt_b, $opt_q, $opt_d);

our $desc='1.filter fq, 2.stats';
our $help=<<EOH;
\t-i FASTAQ file path (./0rawfq) with *.adapter.list and *.fq or *.fq.gz
\t-s Sample list (sample.lst) in format: /^Sample\\tLib\\tInsertSize\$/
\t-o Project output path (./1fqfilted), will mkdir if not exist
\t-q run qsub automatically
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='./0rawfq' if ! $opt_i;
$opt_s='sample.lst' if ! $opt_s;
$opt_o='./1fqfilted' if ! $opt_o;

# `readlink -f` will be blank if target not exists.
system('mkdir','-p',$opt_o);

die "[x]-i $opt_i not exists !\n" unless -d $opt_i;

print STDERR "From [$opt_i] to [$opt_o] refer to [$opt_s]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
