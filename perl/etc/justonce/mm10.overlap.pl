#!/usr/bin/env perl
use strict;
use warnings;
#use Data::Dump qw(ddx);

my %Dat;
while(<>) {
  chomp;
  my ($chr,undef,$type,$start,$end,undef,$strand,undef,undef,$geneid) = split;
  next if $type ne 'transcript';
  $geneid =~ s/[\";\,]//g;
  #$geneid =~ s/(\.\d+$)//g;
  $Dat{"$chr\t$strand"}{$start}{$end} = $geneid;
  #push @{$Dat{"$chr\t$strand"}{$start}},[$end,$geneid];
}

my ($chrStrand,$lastChrStrand);
my ($lastStart,$lastEnd,$laseGene) = (0,0,'');
for $chrStrand (sort keys %Dat) {
  my $startDat = $Dat{$chrStrand};
  my $start;
  for $start (sort {$a <=> $b} keys %{$startDat}) {
    my $endGeneDat = $startDat->{$start};
    #my $end;
    #for $end (sort {$b <=> $a} keys %{$endGeneDat}) {
    #  my $gene = $endGeneDat->{$end};
    #  print join("\t",$chrStrand,$start,$end,$gene),"\n";
    #}
    my @ends = sort {$b <=> $a} keys %{$endGeneDat};
    my $thisGene = $endGeneDat->{$ends[0]};
    if ($lastEnd >= $start and $thisGene ne $laseGene and $lastChrStrand eq $chrStrand) {
      print join("\t",$chrStrand,$start,$ends[0],$thisGene,'<',$lastStart,$lastEnd,$laseGene,'=',$lastEnd-$start,$lastEnd-$ends[0]),"\n";
      #ddx $Dat{$chrStrand}{$start};
    }
    $lastEnd = $ends[0];
    $laseGene = $endGeneDat->{$lastEnd};
    $lastStart = $start;
    $lastChrStrand = $chrStrand;
  }
}
