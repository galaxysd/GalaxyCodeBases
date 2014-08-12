#!/usr/bin/perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my @patterns = qw[
atag
attc
cagt
ccgg
tagc
aact
gaga
gcac
tggc
gcga
acgt
aaa
ata
gga
cga
tgc
taa
tac
gac
ggc
cca
cac
];
sub revcom($) {
    my $str = $_[0];
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    my $rev = reverse $str;
    $rev    =~ tr/[](){}<>/][)(}{></;
    return $rev;
}
my $maxLen = 500;
my $minLen = 50;
my $flanking = 200;
print "[@patterns]\n";
my (@bothpatterns,%t);
push @bothpatterns,revcom($_) for @patterns;
++$t{$_} for (@patterns,@bothpatterns);
@bothpatterns = sort { length($a) <=> length($b) || $a cmp $b } keys %t;
print "[@bothpatterns]\n";

my $in = shift;
#my $out = shift;
my $out;
$in =~ /_(\w+)\.mfa\.gz$/;
$out =$1.'.lst';
print "From [$in] to [$out]\n";

open IN, "-|", "zcat $in";
open O,'>',$out;

sub getpattern($$) {
	my ($seq,$pattLen) = @_;
	#my $pattLen = length $patt;
	my $mincpy = int($minLen/$pattLen);
	my $maxcpy = int(1+$maxLen/$pattLen);
	my (%count, $size);
	$size=$pattLen;
#$mincpy=2;
#$flanking = 10;
	$count{uc(substr($seq, $_, $size))}++ for 0..(length($seq) - $size);
#ddx \%count;
	for (sort {$count{$b} <=> $count{$a} || $a cmp $b} keys %count) {
	    last if $count{$_} < $mincpy;
		next if /^N*$/i;
		my $pattern = $_;
		my (@poses,@groups,$lastpos,$groupid);
		while ($seq =~ /$_/ig) {
			push @poses,pos($seq)-$size;
		}
		$lastpos = shift @poses;
		$groupid = 0;
		for (@poses) {
			if ($_-$lastpos == $size) {
				push @{$groups[$groupid]},$_;	# lacking the 1st one.
			} else {
				++$groupid;
			}
			$lastpos = $_;
		}
		for (@groups) {
			next unless defined $_;
			my @eachposes = @$_;
			if (@eachposes > $mincpy-1) {
				my $start = $eachposes[0]-$size;
				my $left = substr($seq, $start-$flanking, $flanking);
				my $right = substr($seq, $eachposes[-1]+$size,$flanking);
				my $ssrlen = (1+@eachposes) * $size;
				my $ssrseq = substr($seq, $start, $ssrlen);
				my $all = "$left ($ssrseq) $right";
				print O join("\t",$pattern, $ssrlen, 1+$start,$all),"\n";
			}
		}
	    #print "$_ :";
	    #print ' '.(pos($seq)-$size) while $seq =~ /$_/g;
	    #print " ($count{$_} matches)\n";
	}
=pod
	while($seq=~m/($patt){$mincpy,$maxcpy}/ig){
	#while($seq=~m/([ATCG]+)(\1){$mincpy,$maxcpy}/ig){
		#push(@index,$-[2]);
		#pos($seq)=$-[0];
		my $len = length $&;
		print O join("\t",$patt,$len,$-[0],$&,substr($seq, $-[0]-$flanking, $flanking),substr($seq, $+[0],$flanking) ),"\n";# if $len >= $minLen and $len <= $maxLen;
	}
=cut
}

{
	$/ = ">";
	my $junk = <IN>;
	while ( my $record = <IN> ) {
		chomp $record; # Remove the ">" from the end of $record, and realize that the ">" is already gone from the begining of the record
		my ($defLine, @seqLines) = split /\n/, $record;
		my $sequence = join('',@seqLines); # Concatenates all elements of the @seqLines array into a single string.
		$sequence =~ s/\s//g;
		print "$defLine\n"; # Print your definition; remember the ">" has already been removed. Remember to print a newline.
		my $seqlen = length($sequence);
		print "Seq Length: $seqlen\n"; # Print the sequence length and a newline
		print O ">$defLine $seqlen\n";
		#getpattern($sequence,$_) for (@bothpatterns);
		getpattern($sequence,$_) for (3,4);
		print O "\n";

	}
	$/ = "\n";
}
print "\n";

close IN;
close O;


__END__
scp galaxy@192.168.0.83:/Users/Galaxy/t/ssr.pl .

find *.mfa.gz|while read a;do perl ssr.pl "$a";done
awk '{print $2}' *.lst|sort -n|uniq -c

wc -l *.lst|grep -v ' 2 '|grep lst|awk '{print $2}'|xargs tar -czvf ssr.tgz
