#!/bin/env perl
use strict;
use warnings;
use Bio::AlignIO;
use Data::Dump qw(ddx);

die "Usage: $0 <input msf> <output>\n" if @ARGV < 2;
my ($inf,$outf)=@ARGV;

open OUT,'>',$outf or die;

my $in = Bio::AlignIO->new(-file   => $inf,
                           -format => "msf" );

my $alnObj = $in->next_aln(); # get entire alignment data

my (%Seq,$cns,$combined,@CNS,$t,@ids);
foreach my $seqObj ($alnObj->each_seq) {
	print join(',',$seqObj->display_id), "\n";
	push @ids,$seqObj->display_id;
	my $seq = $seqObj->seq;
	$Seq{$seqObj->display_id}=[-1,0,$seq];
	my @bases = split //,$seq;
	$t=0;
	for my $b (@bases) {
		if ($b ne '.') {
			if (exists $CNS[$t]) {
				++$CNS[$t]->{$b};
			} else {
				$CNS[$t] = {$b => 1};
			}
			if ($Seq{$seqObj->display_id}->[0] == -1) {
				$Seq{$seqObj->display_id}->[0] = $t + 1;
			}
		}
		++$t;
	}
	$seq =~ s/\.+$//g;
	$Seq{$seqObj->display_id}->[1] = length($seq);
}

ddx \@CNS;
#ddx \%Seq;
print "@ids\n";

#/share/users/huxs/git/toGit/perl/perlib/etc/Galaxy/Data.pm
our %REV_IUB = (A       => 'A',
                T       => 'T',
                C       => 'C',
                G       => 'G',
                AC      => 'M',
                AG      => 'R',
                AT      => 'W',
                CG      => 'S',
                CT      => 'Y',
                'GT'    => 'K',
                ACG     => 'V',
                ACT     => 'H',
                AGT     => 'D',
                CGT     => 'B',
                ACGT=> 'N',
                N       => 'N'
                );

my %ATGtoBIN = ( A => 1,C => 2,G => 4,T => 8 );

$cns=$combined='';
my @CNSbases;
for my $base (@CNS) {
	my %thebase = %{$base};
	$t = (sort { $thebase{$b} <=> $thebase{$a} } keys %thebase)[0];
	push @CNSbases,$t;
	$cns .= $t;
	$t = join '',(sort { $a cmp $b } keys %thebase);
	$combined .= $REV_IUB{$t} or die;
}

print "$cns\n$combined\n";

my $Main = (split /\*/,$ids[0])[0];
print OUT <<HEAD;
<?xml version="1.0" encoding="UTF-8"?>
<reference schema="350">
<name>${Main}_Alleles</name>
<library>${Main}</library>
<comments>No Comment</comments>
<version>2012_Tang</version>
<start>0</start>
<max_deletion>5</max_deletion>
<consensus>$cns</consensus>
<combined>$combined</combined>

<allele-list>
HEAD

for my $id (@ids) {
	my ($start,$stop,$seq) = @{$Seq{$id}};
	my @bases = split //,$seq;
	$t=0;
	my @variants=();
	for my $b (@bases) {
		if ($b ne '.' and $b ne $CNSbases[$t]) {
			push @variants,join(' ',$t+1,$ATGtoBIN{$b});
			print "$id:$t $b ne $CNSbases[$t]\n";
		}
		++$t;
	}
	$t = join(' ',@variants);
	print OUT <<ITEM;
	<allele>
		<name>$id</name>
		<type>0</type>
		<start>$start</start>
		<stop>$stop</stop>
		<variants>$t</variants>
	</allele>
ITEM
}

print OUT '</allele-list>
</reference>
';
close OUT;
