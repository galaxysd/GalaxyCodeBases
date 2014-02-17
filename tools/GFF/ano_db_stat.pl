#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Vcf;
#use GTF;
#use DBI;
use Galaxy::IO;
use Galaxy::ChromString;
use Data::Dump qw(ddx);

die "Usage: $0 <db_file> <bwa ann file> <out>\n" if @ARGV<3;
my $dbfile = shift;
my $annfile = shift;
my $outfs = shift;

warn "From [$dbfile] to [$outfs] with [$annfile]\n";

open D,'<',$dbfile or die;
open A,'<',$annfile or die;
open O,'>',$outfs or die;
print O "# From [$dbfile] to [$outfs] with [$annfile]\n";

warn "[!]Begin $annfile.\n";
my (%ChrLen,$ChrDat,@IDs,@Lens);
<A>;
while (<A>) {
	my (undef,$id) = split / /;	# gi|477502308|ref|NW_004457745.1|
	if ($id =~ /ref\|([^|]+)\|/) {
		$id = $1;
	}
	$_ = <A>;
	my (undef,$len) = split / /;
	$ChrLen{$id} = $len;
	push @IDs,$id;
	push @Lens,$len;
#print "$id, $len\n";
}
#$ChrDat = ChrStrInit(\@IDs,\@Lens);
close A;

my %bitflag = (
		'CDS' => 1,
		'mRNA' => 2
	);

warn "[!]Begin $dbfile.\n";
my $secname = ']';
my (%GeneDat, $strand, $seqname, $flag);
while (<D>) {
	next if /^(#|((\s)*$))/;
	my ($primary, $start, $end, $frame);
	if (/^\[([^]]*)\] ([+-]) ([^ ]+) (.+)$/) {
#ddx $Gene2Chr{$secname};
		$secname = $1;
		$strand = $2;
		$seqname = $3;	# NW_004465862.1
		$flag = $4;
		#print "[$secname, $strand, $seqname, $flag]\n";
		die "[$_]" if length $secname == 0;
	} else {
		chomp;
		($primary, $start, $end, $frame) = split /\t/,$_;
		if ( $primary eq 'CDS' or $primary eq 'mRNA' ) {
			push @{$GeneDat{$seqname}},[$secname,$primary,$strand,$start, $end, $frame];
		}
	}
	#ddx $GeneDat{$seqname};
}
close D;

warn "[!]Begin Stat.\n";
my (%Results);
my @ChrList = sort keys %GeneDat;
my $ChrCount = scalar @ChrList;
my $currChr=0;
for my $chr (@ChrList) {
	++$currChr;
	my %Stat;
	for my $v ( @{$GeneDat{$chr}} ) {
		for my $pos ( $$v[3] .. $$v[4] ) {
			++$Stat{$pos}{$$v[1]};
		}
	}
	for my $pos (keys %Stat) {
		if (exists $Stat{$pos}{'CDS'}) {
			++$Results{'CDS'};
			++$Results{'mCDS'} if $Stat{$pos}{'CDS'} > 1;
		}
		if (exists $Stat{$pos}{'mRNA'}) {
			++$Results{'mRNA'};
			++$Results{'mmRNA'} if $Stat{$pos}{'mRNA'} > 1;
		}
	}
	print STDERR '[!] ',int(10000*$currChr/$ChrCount)/100," %   \r";
}
print STDERR "\n";
ddx \%Results;
for my $k (sort keys %Results) {
	print O join("\t",$k,$Results{$k}),"\n";
}

close O;

__END__
perl ano_db_stat.pl Dasnov3.db Dasnov3.ann Dasnov3.db.stat
