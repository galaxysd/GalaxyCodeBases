#!/bin/env perl
#use lib "/ifs1/ST_ASMB/USER/huxuesong/public/lib";
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
our $opts='o:';
our($opt_o);

#our $desc='';
our $help=<<EOH;
\t-o output prefix (./matrix).{mcount,mratio}
EOH
our $ARG_DESC='{sam,soap}pe_files';

ShowHelp();
$opt_o='./mchrpos.lst' if ! $opt_o;

print STDERR "From [@ARGV] to [$opt_o]\n";
my $start_time = [gettimeofday];
#BEGIN
my %Dat;
while(<>) {
    chomp;
    my ($chr,$pos)=split /\s+/;
    ++$Dat{$chr}{$pos};
}
open O,'>',$opt_o or die "Error: $!\n";
my %Stat;
for my $chr (sort keys %Dat) {
    for my $pos (sort {$a<=>$b} keys %{$Dat{$chr}}) {
        print O "$chr\t$pos\t$Dat{$chr}{$pos}\n";
        ++$Stat{$Dat{$chr}{$pos}};
    }
}
print "Stat:\n";
print "$_\t$Stat{$_}\n" for sort {$a<=>$b} keys %Stat;
#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

