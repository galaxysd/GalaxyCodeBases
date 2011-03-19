#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <ref> <soap list> <out prefix>\n";
	exit;
}
# $out.soap, $out.single, $out.log
my ($ref,$soaplst,$out)=@ARGV;
my $bin='/nas/RD_09C/resequencing/soft/pipeline/popSNP/bin/soap.coverage';

my $sh="$bin -cvg -refsingle $ref -il $soaplst -o $out.coverage -depthsingle $out.depth >$out.log 2>$out.err";
open OUT,'>',"${out}_cvg.sh.archive" or warn "[!]Error opening ${out}_cvg.sh.archive: $!\n";
print OUT "#!/bin/sh\n$sh\n";
close OUT;
system($sh) or system("echo done ! > $out.tag");

__END__
/ifs1/GAG/population/huxuesong/watermelon/v4/2soap/depth/draw/subBin/soap.coverage -cvg -refsingle /share/raid010/resequencing/resequencing/tmp/pub/Genome/watemelon/watermelon_v4_2509/wmp.merge.fa -il /ifs1/GAG/population/huxuesong/watermelon/v4/2soap/depth/GS-5/soap.l -o /ifs1/GAG/population/huxuesong/watermelon/v4/2soap/depth/GS-5/coverage/total_coverage.info -depthsingle /ifs1/GAG/population/huxuesong/watermelon/v4/2soap/depth/GS-5/depth/total_depthsingle -addn /ifs1/GAG/population/huxuesong/watermelon/v4/2soap/depth/GS-5/depth/GS-5.ngap
