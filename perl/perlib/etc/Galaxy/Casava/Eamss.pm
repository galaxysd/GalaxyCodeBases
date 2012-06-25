
=pod
CASAVA_v1.8.2/src/perl/lib/Casava/Common/Eamss.pm @ Apr 15  2011

Copyright (c) 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

=cut

#package Casava::Common::Eamss;
package Galaxy::Casava::Eamss;

use strict;
use warnings;
use Data::Dumper;
use Carp;
use List::Util qw[min max];

BEGIN {
    use Exporter();
    use vars qw (@ISA @EXPORT @EXPORT_OK);
    @ISA    = qw(Exporter);
    @EXPORT = qw();
    @EXPORT_OK = qw(&maskQvalsByEamss);
}

my $lowScore = 1;
my $mediumScore = 0;
my $highScore = -2;
my $mediumThreshold = ord('O');
my $highThreshold = ord('^');
my $minScore = 1;
my @motifList = ();

sub eamss($);
sub findStr($$$$);

sub maskQvalsByEamss($$) #(std::string &qValues, const std::string &baseCalls) const
{
    my ($qValues, $baseCalls) = @_;
    my ($score, $position) = eamss($qValues);
    return if ($score < $minScore);
    # look for troublemaker motifs
    my $extendedPosition = $position;
    foreach my $motif (@motifList)
    {
        my $troublemakerPosition = findStr($baseCalls, $motif, max(0, $position - 15), min(scalar(@$baseCalls), $position));
        $extendedPosition = min($extendedPosition, $troublemakerPosition) if (-1 < $troublemakerPosition);
    }
    $position = min($position, $extendedPosition);
    # look for extended run of polyG, with at least 90% G bases
    my $numG = 0;
    my $numBases = 0;
    my $maskPolyG = 0;
    my $maskStart = $position;
    for (my $curPos = $position; 0 <= $curPos; $curPos--)
    {
        $numG++ if (ord('G') == ord($baseCalls->[$curPos]));
        $numBases++;
        next if (10 > $numBases);
        my $gFrac = $numG/$numBases;
        if ($gFrac >= 0.9 && ord('G') == ord($baseCalls->[$curPos])) # only start masking at G
        {
            $maskPolyG = 1;
            $maskStart = $curPos;
        }
        elsif ($gFrac < 0.9)
        {
            last;
        }
    }
    $position = $maskStart if ($maskPolyG);
    for (my $idx = $position; scalar(@$qValues) > $idx; $idx++)
    {
        $qValues->[$idx] = 'B';
    }
}

sub eamss($)
{
    my ($qValues) = @_;
    my $curScore = 0;
    # initialize the bestscore to something lower than the first value of curScore
    my $bestScore = min($highScore, $mediumScore, $lowScore) - 1;
    my $bestPosition = -1;
    for (my $idx = $#{$qValues}; 0 <= $idx; $idx--)
    {
        if (ord($qValues->[$idx]) >= $highThreshold)
        {
            $curScore += $highScore;
        }
        elsif (ord($qValues->[$idx]) >= $mediumThreshold)
        {
            $curScore += $mediumScore;
        }
        else
        {
            $curScore += $lowScore;
        }
        if ($curScore >= $bestScore)
        {
            $bestScore = $curScore;
            $bestPosition = $idx;
        }
    }
    return ($bestScore, $bestPosition);
}

sub findStr($$$$)
{
    my ($targetString, $queryString, $start, $stop) = @_;
    my $substring = substr($targetString, $start, $stop - $start + 1);
    return index($substring, $queryString);
}

1;
__END__
