#!/usr/bin/perl -w
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;

our $opts='i:o:n:g:v:r:c:f:h:bqd';
our($opt_i, $opt_n, $opt_g, $opt_r, $opt_c, $opt_f, $opt_o, $opt_v, $opt_b, $opt_q, $opt_d, $opt_h);

our $help=<<EOH;
\t-i soaps.lst
\t-c Chromosome NFO file (chr.nfo) in format: /^ChrID\\s+ChrLen\\s?.*\$/
\t-o Output path (./ffq) for [ChrID/Sample.ffq], will mkdir if not exist
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

die "[x]-i $opt_i not exists !\n" unless -f $opt_i;
$opt_c='chr.nfo' if ! $opt_c;
$opt_o='./ffq' if ! $opt_o;

no warnings;
$opt_v=int $opt_v;
use warnings;

$opt_o=~s#/+$##;

print STDERR "From [$opt_i] to [$opt_o] refer to [$opt_c]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
system('mkdir','-p',$opt_o);
my %ChrOutFH;
open ChrNFO,'<',$opt_c or die "[x]Error opening $opt_c: $!\n";
while (<ChrNFO>) {
	next if /^#/;
	my ($chrid)=split /\t/;
	open $ChrOutFH{$chrid},'|-',"gzip -9c - >$opt_o/$chrid.ffq.gz" or die "[x]Error opening $opt_o/$chrid.ffq with gzip: $!\n";
}
close ChrNFO;
open $ChrOutFH{'__UnKnown__'},'|-',"gzip -9c - >$opt_o/__UnKnown__.ffq.gz" or die "[x]Error opening $opt_o/__UnKnown__.ffq with gzip: $!\n";

open LST,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
while (<LST>) {
	chomp;
	my ($PESE,$sample,$lib,$FL,$ReadLen,$nfofpath)=split /\t/;
	$nfofpath =~ s/\.nfo$//;
	my @Files;
	if ($PESE eq 'PE') {
		@Files=("$nfofpath.soap","$nfofpath.single");
	} else {
		@Files=("$nfofpath.soap");
	}
	unless (-s "$nfofpath.unmap") {
		open FQ,'-|',"gzip -dc $nfofpath.unmap.gz" or warn "[x]Error opening $PESE [$nfofpath.unmap.gz] with gzip: $!\n" and next;
	} else { open FQ,'<',"$nfofpath.unmap" or warn "[x]Error opening $PESE [$nfofpath.unmap]: $!\n" and next; }
	my ($name,$seq,$Qstr);
	while (<FQ>) {
		chomp;
		s/^[@>]//;
		$name=$_;
		chomp($seq=<FQ>);
		<FQ>;
		chomp($Qstr=<FQ>);
		print ${$ChrOutFH{'__UnKnown__'}} join("\t",$FL,$name,'.','0','0',$seq,$Qstr,$sample,$lib),"\n";
	}
	close FQ;
	for my $file (@Files) {
		unless (-s $file) {
			open SP,'-|',"gzip -dc $file.gz" or warn "[x]Error opening $PESE [$file.gz] with gzip: $!\n" and next;
			print STDERR '[!]',join(', ',$sample,$lib,$FL,"$file.gz");
		} else {
			open SP,'<',$file or warn "[x]Error opening $PESE [$file]: $!\n" and next;
			print STDERR '[!]',join(', ',$sample,$lib,$FL,$file);
		}
		while (<SP>) {
			my ($soapid,$seq,$Qstr,$hit,$strand,$chr,$pos) = (split(/\t/))[0,1,2,3,6,7,8];
			next unless defined $pos;
			$chr='__UnKnown__' unless exists $ChrOutFH{$chr};
			if ($hit == 1) {
				print ${$ChrOutFH{$chr}} join("\t",$FL,$soapid,$chr,$pos,$hit,$seq,$Qstr,$sample,$lib),"\n";
			} else {
				print ${$ChrOutFH{'__UnKnown__'}} join("\t",$FL,$soapid,$chr,$pos,$hit,$seq,$Qstr,$sample,$lib),"\n";
			}
		}
		warn ", done.\n";
	}
}
close LST;

close $ChrOutFH{$_} for keys %ChrOutFH;

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
