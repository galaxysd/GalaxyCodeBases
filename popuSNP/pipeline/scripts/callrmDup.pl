#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <chrorder> <soap_list> <out prefix> <dp>\n";
	exit;
}
# $out.soap, $out.single, $out.log
my ($chrO,$ilst,$out,$dp)=@ARGV;
$dp = $dp?'-dp':'';
#warn $DM;
my @o=split /\//,$out;
my $prefix=pop @o;
my $path=join '/',@o;
my $bin='/nas/RD_09C/resequencing/soft/pipeline/popSNP/scripts/rmdSM';
my $sh="$bin -chr $chrO $dp -out_dir $path -prefix $prefix -in_list $ilst >$out.log 2>$out.err";
# if you need monoploid calling mode, -m must be just in front of -I
open OUT,'>',"${out}_rmd.sh.archive" or warn "[!]Error opening ${out}_rmd.sh.archive: $!\n";
print OUT "#!/bin/sh\n$sh\n";
close OUT;
system($sh) or system("echo done ! > $out.tag");

__END__
 /nas/RD_09C/resequencing/soft/bin/rmdSM -chr chrorder -dp -out_dir ./ttt/ -prefix ORYqzpRFSDIAAPEI-5 -in_list ORYqzpRFSDIAAPEI-5.ll >s.log 2>s.err &
