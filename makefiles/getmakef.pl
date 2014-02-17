#!/bin/env perl
use strict;
#use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

die "Usage: $0 <fq.gz path> <out path>\n" if @ARGV < 2;
my ($inp,$outp)=@ARGV;

my $Ref = 'Dasnov3';
my $BwaSai = "bwa aln -l 17 -q 10 $Ref";
my $BwaSam = "bwa sampe -a 800 $Ref";
my $BwaSort = 'perl pickU_sort.pl';
my $BwaRmDup = 'samtools rmdup';
my $readsCounter = '/share/users/huxs/git/toGit/c_cpp/faststater/readsCounter';

opendir my($dh), $inp or die "Couldn't open dir '$inp': $!";
my @files = readdir $dh;
closedir $dh;
# http://perlmeme.org/faqs/file_io/directory_listing.html

@files = grep(/\.f(ast|)q\b/i, @files);

#ddx \@files;

my (%Pairs,@AllTargets);
for ( @files ) {
	m~^([^/]+?)([\W_]*?)([12])?(\.f(ast|)q(\.gz)?)$~i or die $_;
	my ($m,$rp,$r,$ext) = ($1,$2,$3,$4);
	$m =~ s/\W//g and die $m;
	print "$1, $2, $3, $4, $5, $6, $7\t$m,$rp,$r,$ext, $_\n";
die unless ($r == 1) or ($r==2);	# no SE now
	push @{$Pairs{$m}},[$_,$r];
}
@AllTargets = sort keys %Pairs;

ddx \%Pairs;
ddx \@AllTargets;

open M,'>','Makefile' or die $!;
print M 'all: ',join(' ',@AllTargets,'_AdditionalTG_'),"\n";
my @AdditionalTG;

mkdir 'sai';
mkdir 'sam';
mkdir 'bam';

for my $Target ( @AllTargets ) {
	my @Files = @{$Pairs{$Target}};
die if @Files != 2;	# no SE now
	my (@fqFiles,@newFiles);
	for (@Files) {
		my ($fqname,$read12) = @$_;
		$newFiles[$read12-1] = "${Target}_$read12";
		$fqFiles[$read12-1] = $fqname;
		print M "
sai/${Target}_$read12.sai: $inp$fqname
\t$BwaSai $inp$fqname >sai/${Target}_$read12.sai 2>sai/${Target}_$read12.logsai
\tdate >> sai/${Target}_$read12.logsai

sai/${Target}_$read12.fqstat: $inp$fqname
\t$readsCounter -o sai/${Target}_$read12.fqstat $inp$fqname
";
		push @AdditionalTG, "sai/${Target}_$read12.fqstat";
	}
	print M "
sam/${Target}.sam.gz: ",join(' ',map { "sai/$_.sai" } @newFiles),"
\t$BwaSam ",join(' ',(map { "sai/$_.sai" } @newFiles),(map { "$inp$_" } @fqFiles) )," 2>sam/${Target}.logsam |gzip -9c >sam/${Target}.sam.gz
\tdate >> sam/${Target}.logsam
";
	print M "
sam/${Target}.sort.bam: sam/${Target}.sam.gz
\t$BwaSort sam/${Target}.sam.gz sam/${Target}.sort.bam
";
	print M "
bam/${Target}.rmdup.bam: sam/${Target}.sort.bam
\t$BwaRmDup sam/${Target}.sort.bam bam/${Target}.rmdup.bam
";
	print M "${Target}: bam/${Target}.rmdup.bam\n";
}

	print M join(' ','_AdditionalTG_:',@AdditionalTG),"\n";
close M;

__END__
perl getmakefln.pl fq/ .
