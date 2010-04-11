#!/usr/bin/perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.1.2;

our $opts='i:o:bv';
our($opt_i, $opt_o, $opt_v,$opt_b);

our $help=<<EOH;
\t-i Add_ref files list (./indsnp.lst)
\t-o Output PWM file (snp.pwm)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH
# \t-d dbSNP SQLite data file (_tdbSNP.sqlite)

ShowHelp();

$opt_i='./indsnp.lst' if ! defined $opt_i;
$opt_o='snp.pwm' if ! $opt_o;

print STDERR "From [$opt_i] to [$opt_o]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

#http://doc.bioperl.org/bioperl-live/Bio/Tools/IUPAC.html#BEGIN1
my %IUB = ( A => [qw(A)],
	     C => [qw(C)],
	     G => [qw(G)],
	     T => [qw(T)],
	     U => [qw(U)],
	     M => [qw(A C)],
	     R => [qw(A G)],
	     W => [qw(A T)],
	     S => [qw(C G)],
	     Y => [qw(C T)],
	     K => [qw(G T)],
	     V => [qw(A C G)],
	     H => [qw(A C T)],
	     D => [qw(A G T)],
	     B => [qw(C G T)],
	     X => [qw(G A T C)],
	     N => [qw(G A T C)]
	     );

my $start_time = [gettimeofday];

my @files;
open L,'<',$opt_i or die "Error opening $opt_i: $!\n";
while (<L>) {
	chomp;
	if (-s $_) {push @files,$_;}
	 else {warn "[!]$_ not available.\n";}
}
close L;
warn '[!]Total: ',scalar @files," files.\n";

my $snpcount;
open O,'>',$opt_o or die "Error opening $opt_o: $!\n";
for my $file (@files) {
	open I,'<',$file or die "Error opening $file: $!\n";
	while (<I>) {
		chomp;
		my ($chr,$pos,$bases) = split /\t/;
		++$snpcount;
		my @base = split / /,$bases;
		print O "$chr\t$pos\t$base[0]\t";
		#if ($opt_v) {
		#	$bases =~ s/ //g;
		#	print ">$chr\t$pos\t$base[0]\t$bases\n";
		#}
		print ">$chr\t$pos\t" if $opt_v;
		my ($d,$i,%counter,%stat);
		$stat{1}=$stat{2}=$stat{3}=$stat{4}=0;
		for (@base) {
			my $arr=$IUB{$_} or next;	# So, $i != @base;
			$d = @$arr;
			print "$_$d" if $opt_v;
			++$stat{$d};
			$d = 1/$d;
			for (@$arr) {
				$counter{$_} += $d;
				#$i += $d;
			}
			++$i;
		}
		print O "$i,$stat{1},$stat{2},$stat{3},$stat{4}\t";
		@base=();
		for (sort {$counter{$b} <=> $counter{$a}} keys %counter) {
			push @base,sprintf('%s:%f',"$_",$counter{$_}/$i);
		}
		print O join(',',@base),"\n";
		print "\n$i,$stat{1},$stat{2},$stat{3},$stat{4}\t",join(',',@base),"\n\n" if $opt_v;
	}
	close I;
}
close O;

warn "[!]Total: $snpcount SNPs.\n[!]Done !\n";

my $stop_time = [gettimeofday];
#$|=1;
print STDERR "\n Time Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s)\n";
