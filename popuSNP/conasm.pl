#!/bin/env perl
use strict;
use warnings;
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;

our $opts='i:o:c:bvne';
our($opt_i, $opt_o, $opt_c, $opt_v, $opt_b, $opt_n, $opt_e);

our $desc='Genome creater for population analyse';
our $help=<<EOH;
\t-c Consensus file with refbase
\t-i filtered Indel file
\t-o Output file prefix for the sample, directories must exist
\t-n mask noncoverage bases to N, otherwise little cased
\t-e skip hete indels, otherwise both homo and hete(little cased) used
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

warn "[x]Must specify -c xx.consensus.txt !\n" unless $opt_c;
warn "[x]Must specify -i xx.indel.txt.filter !\n" unless $opt_i;
warn "[x]Must specify -o ./xx !\n" unless $opt_o;
exit 1 unless $opt_i and $opt_o and $opt_c;

print STDERR "From [$opt_c][$opt_i] to [$opt_o]\n ";
print STDERR "noncoverage mask to n.   " if $opt_n;
print STDERR "hete indels will be skiped." if $opt_e;
warn "\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

open CNS,'<',$opt_c or die "Error opening $opt_c: $!\n";
my (%Genome,%Depth);
my ($chr,$pos,$best,$depth,$ref,$indel,$bases,$homhet,$strand,$t);
while (<CNS>) {
    chomp;
    ($chr,$pos,undef,$best,undef,undef,$depth,$ref) = split /\t/;
    unless ($ref) {   # if file not completed.
        print STDERR '.';
        next;
    }
    if ($depth==0) {
        unless ($opt_n) {$best=lc $ref;}
         else {$best='n';}
    } else {$Depth{$chr}->{$pos}=chr $depth;} # string is smaller than double
# WARING: Something might be wrong if $depth > 255 !
# However, people should not have such much money.
#$best = $ref; # DEBUG ONLY !!!
    $Genome{$chr}->[$pos]=$best;
    #$Depth{$chr}->[$pos]=$depth;
}
close CNS;
warn "[!]CNS done.\n";

open INDEL,'<',$opt_i or die "Error opening $opt_i: $!\n";
while (<INDEL>) {
    ($chr,$pos,$indel,$bases,$strand,$homhet,undef,$depth,$t) = split /\t/;
    unless ($t) {   # if file not completed.
        print STDERR '\'';
        next;
    }
    next unless exists $Genome{$chr};
    if ($homhet eq 'hete') {
        unless ($opt_e) {$bases = lc $bases;}
         else {next;}
    }
    if ($indel =~ /^([ID])(\d+)$/) {
        if ($1 eq 'I') {
            $Genome{$chr}->[$pos] .= $bases;  # if length($Genome{$chr}->[$t]) == 1
            $Depth{$chr}{$pos} .= chr($depth) x $2;
        } else {
            my $i=$2-1;
            #++$pos;    # deletion including it.
            if ($opt_v) {
                print "[D]$indel:${1}_${2}@","$pos $strand\t";
                #print $Genome{$chr}->[$_] for ($pos..$pos+$i);
                print(($t=$Genome{$chr}->[$_])?$t:'.') for ($pos-5..$pos-1);
                print '<';
                print(($t=$Genome{$chr}->[$_])?$t:'.',' ') for ($pos..$pos+$i);
                print '>';
                print(($t=$Genome{$chr}->[$_])?$t:'.') for ($pos+$i+1..$pos+$i+5);
                print ' - ',$bases,"\n";
            }
            for $t ($pos..$pos+$i) {
                $Genome{$chr}->[$t]='';
                delete $Depth{$chr}{$pos};# if exists $Depth{$chr}{$pos};
            }
        }
    } else {print STDERR '|';next;}
}
close INDEL;
warn "[!]INDEL done.\n";

my ($i,$dep);
for $chr (keys %Genome) {
    $Genome{$chr}->[0]='';  # no longer undef
    $t=$opt_o.'.'.$chr;
    open OUT,'>',$t.'.fa' or die "Error opening ${t}.fa: $!\n";
    print OUT ">$chr\n";
    open DEP,'>',$t.'.dep' or die "Error opening ${t}.dep: $!\n";
    print DEP ">$chr\n";
    $t=0;
    $i=-1;
    for (@{$Genome{$chr}}) {
        ++$i;
        next if $_ eq '';
        ++$t;
        print OUT $_;
        $dep=$Depth{$chr}{$i} or $dep="\0";
        $dep = join(' ',map ord,split //,$dep);
        #print DEP '[',$dep,']';
        if ($t > 80) {
            $t=0;
            print OUT "\n";
            print DEP "$dep\n";
        } else {print DEP $dep,' ';}
    }
    close OUT;
    close DEP;
}
warn "[!]Output done.\n";

__END__
vf=11.9g for 1.3G Feb 20 07:23 consensus/QRS9.Gm03.txt, 47,781,076
