#!/bin/sh

samtools view -F4 $1 |perl -lane 'BEGIN {my ($cnt,$sum);} ++$cnt; my $a;$a+=$1 while ($F[5]=~m/(\d+)M/g); $sum+=$a; END {print "F4\t$cnt -> $sum, ",int(.5+10*$sum/$cnt)/10;}' > $1.rF4

samtools view -f4 $1 |perl -lane 'BEGIN {my ($cnt,$sum);} ++$cnt; my $a;$a+=$1 while ($F[5]=~m/(\d+)M/g); $sum+=$a; END {print "f4\t$cnt -> $sum, ",int(.5+10*$sum/$cnt)/10;}' > $1.f4

#samtools view -F4 Blood-MDA.sort.bam|perl -lane 'BEGIN {my ($cnt,$sum);} ++$cnt; my $a;$a+=$1 while ($F[5]=~m/(\d+)M/g); $sum+=$a; END {print "$cnt -> $sum";}' > Blood-MDA.sort.F4

#echo samtools view -F4 $1 \|perl -lane \''BEGIN {my ($cnt,$sum)=(1,2);} ++$cnt; my $a;$a+=$1 while ($F[5]=~m/(\d+)M/g); $sum+=$a; END {print "F4\t$cnt -> $sum";}'\' \> $1.rF4

