#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dump;

open IA,'<','crep_all_tsv_new.txt.up2.txt' or die $!;
open IB,'<','crep_all_tsv_new.txt.up2.rev.txt' or die $!;
open IR,'<','resLOC.rpratio' or die $!;

open OA,'>','ricexpro.up2.clu.txt' or die $!;
open OB,'>','ricexpro.up2.rev.clu.txt' or die $!;
open OS,'>','ricexpro.up2.stat.txt' or die $!;

my (%dA,%dB,%dFlag,%cntdflag);
while (<IA>) {
	my ($id,$cnt,$value) = split /\t/;
	$id =~ s/^LOC_Os/t/;
	$dA{$id} = log($value) / log(2);
	++$dFlag{$id};
}
close IA;
while (<IB>) {
	my ($id,$cnt,$value) = split /\t/;
	$id =~ s/^LOC_Os/t/;
	$dB{$id} = log($value) / log(2);
	$dFlag{$id} += 100;
}
close IB;

for my $k (keys %dFlag) {
	push @{ $cntdflag{$dFlag{$k}} },$k;
}
for my $k (keys %cntdflag) {
	print "$k -> ",scalar @{ $cntdflag{$k} },': ';
	print "\n";	# join(',',@{ $cntdflag{$k} }),
}
# 1 -> 5332:
# 101 -> 531:
# 100 -> 8173:
my ($max,$min,$sum,$n,$ss)=(-999,999,0,0,0);
for my $id (@{ $cntdflag{'101'} }) {
	#print OS "# $id: $dA{$id}, $dB{$id}\n";
	$sum += $dA{$id}+$dB{$id};
	$ss += $dA{$id}*$dA{$id} + $dB{$id}*$dB{$id};
	$n += 2;
	$max = $dA{$id} if $max < $dA{$id};
	$min = $dA{$id} if $min > $dA{$id};
	$max = $dB{$id} if $max < $dB{$id};
	$min = $dB{$id} if $min > $dB{$id};
}
my $mean = $sum/$n;
print "Range:[$min,$max] Mean: $mean SD:",sqrt($ss/$n - $mean*$mean)," Cnt:$n\n";
=pod
1 -> 5332:
101 -> 531:
100 -> 8173:
Range:[-5.28337053983209,5.82964609302319] Mean: -0.575308843011672 SD:1.62278420169048 Cnt:1062
=cut
my $outHead = <IR>;
print OS "$outHead";
print OA "$outHead";
print OB "$outHead";
while (<IR>) {
	my ($id,$acc,@dat) = split /\t/;
	my $LOCid = $id;
	$id =~ s/^LOC_Os/t/;
	my $tmp = join('', map {my $a; if (defined $_) {$a = int($_); $a="+$a" if $_>=0; $a="-$a" if $_<0 and $a==0; $a; } } ($dA{$id},$dB{$id}));
	my $newLOCid = $LOCid . $tmp;
	my $newid = "${id}_${acc}$tmp";
	if ($dFlag{$id} == '101') {
		print OS join("\t",$newid,   $dA{$id} .';'. $dB{$id}   ,@dat);
	} elsif ($dFlag{$id} == '1') {
		print OA join("\t",$newid,$newLOCid,@dat);
	} elsif ($dFlag{$id} == '100') {
		print OB join("\t",$newid,$newLOCid,@dat);
	} else { die; }
}
close IR;
close OS;
close OA;
close OB;

__END__
$ wc -l crep_all_tsv_new.txt.up2.*txt ricexpro.up2.* resLOC.rpratio
   8704 crep_all_tsv_new.txt.up2.rev.txt
   5863 crep_all_tsv_new.txt.up2.txt
   4756 ricexpro.up2.clu.txt
   7344 ricexpro.up2.rev.clu.txt
    401 ricexpro.up2.stat.txt
  12499 resLOC.rpratio
  39567 total

./cluster -f ricexpro.up2.clu.txt -l -cg a -ng -g 2 -e 2
./cluster -f ricexpro.up2.rev.clu.txt -l -cg a -ng -g 2 -e 2
