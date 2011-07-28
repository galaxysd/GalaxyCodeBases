#!/usr/bin/perl -w
use strict;
use warnings;
#use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.3.1;

our $opts='i:o:g:w:bv';
our ($opt_i, $opt_o, $opt_g, $opt_w, $opt_v, $opt_b);

our $DEPTH_SINGLE_LENGTH_PER_ROW=100;

our $help=<<EOH;
\t-i Input depth file (total_depthsingle)
\t-g Genome file (human.fa)
\t-w Window size / $DEPTH_SINGLE_LENGTH_PER_ROW (10) [10 for win_size=1000, must use integer]
\t  depthsingle file must contain $DEPTH_SINGLE_LENGTH_PER_ROW values pre line except for the last.
\t-o Output SQLite data file (data.tsv)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='total_depthsingle' if ! defined $opt_i;
$opt_o='data.tsv' if ! $opt_o;
$opt_w=10 if ! $opt_w;
$opt_g='human.fa' if ! $opt_g;

my $win_s=$opt_w * $DEPTH_SINGLE_LENGTH_PER_ROW;

print STDERR "From [$opt_i][$opt_g] to [$opt_o] with [$opt_w]x100\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}
$|=1;
my $start_time = [gettimeofday];
#BEGIN
sub calGConek($) {
	my $seq_ref=$_[0];
	my ($i,$gccount,$str,@STR,@GC)=(0);
	my $len=length $$seq_ref;
print "$len -> ";
	while ($i <= $len-$win_s) {
		use integer;
		$gccount=0;
		$str = substr $$seq_ref,$i,$win_s;
#print '.';
		@STR=split //,$str;
#print 'o';
		for (@STR) {
			++$gccount if $_ =~ /[gc]/i;
		}
		no integer;
		push @GC,($gccount/$win_s);
		$i += $win_s;
#print " $gccount\t"; sleep 1;
	}
		$gccount=0;
		$str = substr $$seq_ref,$i;
		@STR=split //,$str;
		for (@STR) {
			++$gccount if $_ =~ /[gc]/i;
		}
		push @GC,($gccount/$win_s);
print 1+$#GC,"\n";
	return [@GC];
}

sub cal_Depth_1k($) {
	my $arref=$_[0];
	my ($dep,$len);
	$len=1+$#$arref;
	return -1 if $len == 0;	# should never happens
	for (@$arref) {
		$_ = 0 if $_ == 65535;
		$dep += $_;
	}
	return $dep/$len;
}

my (%GC,%DEPTH,%RESULT,@line,@win,$seqname,$dep,$count,$countchr);
$countchr=$count=0;
open DEPTH,'<',$opt_i or die "Error: $!\n";
while (<DEPTH>) {
	#s/^>//;
	my $title = $_;
	if($title =~ /^>(\S+)/) {
		if ($#win > -1) {	# tailing ...
			$dep=cal_Depth_1k(\@win);
			push @{$DEPTH{$seqname}},$dep;
			print STDERR "$countchr -> ",1+$#{$DEPTH{$seqname}},"\n";
		}
		@win=();$countchr=$count=0;
		$seqname = $1;
		$seqname =~ s/^chr
		(?>
			((?<=^chr)o)?
			((?<=^chro)m)?
			((?<=^chrom)o)?
			((?<=^chromo)s)?
			((?<=^chromos)o)?
			((?<=^chromoso)m)?
			((?<=^chromosom)e)?
		)//xi;
		print STDERR "loading Depth > $seqname ...\t";
		next;
	}
	chomp;
	@line=split / /;
	++$countchr;
	if ($count < $opt_w) {
		push @win,@line;
#print $#line,' ',$line[-1],"\n"; sleep 1;
		++$count;
		next;
	} else {
		$dep=cal_Depth_1k(\@win);
		push @{$DEPTH{$seqname}},$dep;
		@win=@line;$count=1;	# never miss the last one
	}

=pod
	$/=">";
	my $genome=<DEPTH>;
	chomp $genome;
	$genome=~s/\n//g;
	$/="\n";
	$DEPTH{$seqname}=cal_Depth_1k(\$genome);
	$genome='';
=cut
#for (@{$DEPTH{$seqname}}) {
#	print ">$_\t";
#	sleep 1;
#}
}
if ($#win > -1) {	# last line ...
	$dep=cal_Depth_1k(\@win);
	push @{$DEPTH{$seqname}},$dep;
	print STDERR "$countchr -> ",1+$#{$DEPTH{$seqname}},"\n";
}
@win=();$count=0;
close DEPTH;

#print $#{$DEPTH{1}},"\t";
#print "$_ " for @{$DEPTH{1}};
#print "\n";die;

open GENOME,'<',$opt_g or die "Error: $!\n";
while (<GENOME>) {
	s/^>//;
	my $title = $_;
	my $seqname = $1 if($title =~ /^(\S+)/);
	$seqname =~ s/^chr
		(?>
			((?<=^chr)o)?
			((?<=^chro)m)?
			((?<=^chrom)o)?
			((?<=^chromo)s)?
			((?<=^chromos)o)?
			((?<=^chromoso)m)?
			((?<=^chromosom)e)?
		)//xi;
print STDERR "loading Seq. > $seqname ...\t";
	$/=">";
	my $genome=<GENOME>;
	chomp $genome;
	$genome=~s/\s//g;
	$/="\n";
#print "$seqname done.",substr($genome,0,10)," \n";
	$GC{$seqname} = &calGConek(\$genome);
#print $GC{$seqname},"-\n";
	$genome='';
#print ${$GC{$seqname}}[0],"\t",${$GC{$seqname}}[1],"\t",${$GC{$seqname}}[-1],"\n";
#for (@{$GC{$seqname}}) {
#	print ">$_\t";
#	sleep 1;
#}

}
close GENOME;


for my $chr (sort keys %DEPTH) {
	for (@{$DEPTH{$chr}}) {
		my $gcvalue=shift @{$GC{$chr}};
#print "$gcvalue\t";
		$gcvalue=int($gcvalue*100);
#print "$gcvalue\n";
		push @{$RESULT{$gcvalue}},$_;
	}
}

open OUT,'>',$opt_o or die "Error: $!\n";
for my $gcv (sort {$a <=> $b} keys %RESULT) {
	print OUT $gcv,"\t";
	for (@{$RESULT{$gcv}}) {
		print OUT "$_\t";
	}
	print OUT "\n";
}
close OUT;
#END
my $stop_time = [gettimeofday];

$|=1;
print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";

