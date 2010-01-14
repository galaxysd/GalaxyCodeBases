#!/usr/bin/perl -w
#use threads 1.73;
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
#use Galaxy::Data;
#use Fcntl qw(:DEFAULT :flock);

$main::VERSION=0.1.1;

our $opts='i:o:bv';
our($opt_i, $opt_o, $opt_v,$opt_b);

our $help=<<EOH;
\t-i Annotation files list (./anno.lst) [ID\\tPath_to_file\\n]
\t-o Output mixed file (mixed_annot.txt)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH
# \t-d dbSNP SQLite data file (_tdbSNP.sqlite)

ShowHelp();

$opt_i='./anno.lst' if ! defined $opt_i;
$opt_o='mixed_annot.txt' if ! $opt_o;

print STDERR "From [$opt_i] to [$opt_o]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

my (@Annots,%Keys);
open( SAMP,'<',$opt_i) or die "Error: $!\n";
open( OUT,'>',$opt_o) or die "Error: $!\n";
print OUT 'Gene';
while (<SAMP>) {
	chomp;
	my ($id,$file)=split /\s+/;
	my %dat;
	getHash($file,\%dat);
	++$Keys{$_} for keys %dat;
	push @Annots,[$id,\%dat];
	print OUT "\t$id\t${id}_nfo";
}
close SAMP;
print OUT "\n";

sub getHash {
	my( $file ,$hash) = @_;
	open (F,"$file") or die $!;
	while (<F>) {
		chomp;
		my($key,$a,$b)=split /\s+/,$_,3;
		$key=~s/\.\d+$//;
		$a=~s/\t/ /g;
		$b=~s/\t/ /g;
		my $value = "$a\t$b";
	#print "$key\t$value\n";
		$$hash{$key} = $value;
	}
	close F;
}


for my $gene (sort keys %Keys) {
	print OUT $gene;
	for (@Annots) {
		my ($id,$dathash)=@$_;
		my $v=(defined $$dathash{$gene})?$$dathash{$gene}:".\t.";
		print OUT "\t$v";
	}
	print OUT "\n";
}

close OUT;

my $stop_time = [gettimeofday];
#$|=1;
print STDERR "\n Time Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s)\n";

print STDERR "\033[32;1mgrep -i 'transcription factor' mixed_annot.txt > out\033[0;0m\n";

__END__
