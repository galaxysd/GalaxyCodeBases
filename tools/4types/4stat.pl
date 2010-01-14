#!/bin/env perl
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use lib 'E:\BGI\toGit\perlib\etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
use GalaxyXS::ChromByte 1.02;
use DBI;

$main::VERSION=0.0.1;

our $opts='i:o:c:d:l:t:bv';
our ($opt_i, $opt_o, $opt_c, $opt_v, $opt_b, $opt_d, $opt_l, $opt_t);

our $help=<<EOH;
\t-i GFF files list (gff.lst) [SampleName\\tPath_to_file\\n](8 lines max !)
\t-c Chromosome length list (chr.len) [ChrID\\tLen\\n]
\t-l Minimal overlap length (10)
\t-o Output Stat (stat.txt)
\t-d Details dump to (details.lst)
\t-t tmpfile for swap, will be removed (./._tmp_)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='gff.lst' if ! defined $opt_i;
$opt_o='stat.txt' if ! $opt_o;
$opt_d='details.lst' if ! $opt_d;
$opt_c='chr.len' if ! $opt_c;
$opt_t='./._tmp_' if ! $opt_t;
$opt_l=10 if ! $opt_l;
$opt_l=int($opt_l);
$opt_l=1 if $opt_l<1;

print STDERR "From [$opt_i] with [$opt_c][$opt_l][$opt_t] to [$opt_o] [$opt_d]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

open( T,'>',$opt_t) or die "[x]Error creating swap file: $!\n";
close T;
unlink $opt_t;

sub unescape {
  my $v = shift;
  $v =~ tr/+/ /;
  $v =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $v;
}
sub read_gene_gff($$$$) {
	my ($file,$bit,$genehash,$sth)=@_;
	my @dat;
	while (<$file>) {
		next if /^\s+$/ or /^;+/;
		chomp;
		s/\r//g;	# There DO be some file with \r. DAMN the MS and APPLE !
		my ($seqname, $source, $primary, $start, $end,
		$score, $strand, $frame, $groups) = split /\t/;	# well, reading file no need to Optimize
		my @groups = split(/\s*;\s*/, $groups);
		my (%groups,$name);
		for my $group (@groups) {
			my ($tag,$value) = split /=/,$group;
			next unless defined $value;
			$tag = unescape($tag);
			my @values = map {unescape($_)} split /,/,$value;
			$groups{$tag}=\@values;	# patch for those alter-splices
		}
		my @name_order=qw/ID Target/;
		for (@name_order) {
			if ($groups{$_}) {$name=$groups{$_};last;}
		}
		for (@$name) {
			push @{$$genehash{$seqname}{$bit}},[$_,$start,$end];
			$sth->execute( $bit,$seqname,$start,$end,$_ );	# samplebit,chrid,start,end,name
			print "$_\t$seqname,$start,$end,$primary,$strand,$frame\n" if $opt_v;
		}
	}
	1;
}

my $start_time = [gettimeofday];

my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0
);
my $dbh = DBI->connect('dbi:SQLite:dbname='.$opt_t,'','',\%attr) or die $DBI::errstr;

my ($i,$bit,%ChrLen,%Genes,@Samples,%SampBit,%BitSamp,%Combines,%Log2,%Summary,%Groups,%Bits,%Tables)=(0);
my $sql=q/
CREATE TABLE IF NOT EXISTS gff
(  samplebit INTEGER,
   chrid TEXT,
   start INTEGER,
   end INTEGER,
   name TEXT );
/;
for (split /;/,$sql) {
	next if /^\s*$/;
	$dbh->do($_) or die $dbh->errstr;
}
$dbh->commit;
my $sth = $dbh->prepare( "INSERT INTO gff ( samplebit,chrid,start,end,name ) VALUES ( ?,?,?,?,? )" );
my $gth = $dbh->prepare( "SELECT name FROM gff WHERE samplebit=:1 AND chrid=:2 AND start <= :3 AND end >= :4" );

open( C,'<',$opt_c) or die "[x]Error: $!\n";
while (<C>) {
	chomp;
	my ($chr,$len)=split /\s+/;
	$ChrLen{$chr}=$len;
}
close C;

open( D,'>',$opt_d) or die "[x]Error creating details file: $!\n";
print D ";Group Dump file\n[SampleID]\n";
open( L,'<',$opt_i) or die "[x]Error: $!\n";
while (<L>) {
	use integer;
	chomp;
	my ($id,$file)=split /\s+/;
	push @Samples,$id;
	$bit=1;
	$bit *= 2 for (1..$i);
	warn "$i  $id\t$file\n";
	print D chr(65+$i)," = $id\n";
	$Log2{$bit}=$i;
	$SampBit{$id}=$bit;
	$BitSamp{$bit}=$id;
	my $infile;
	open( $infile,"<",$file) or die "[x]Error: $!\n";
	read_gene_gff($infile,$bit,\%Genes,$sth);
	close $infile;
	++$i;
	if ($i>8) {
		warn "[!]Only first 8 samples used.\n";
		last;
	}
	#++$Keys{$_} for keys %dat;
	#push @Annots,[$id,\%dat];
	#print OUT "\t$id\t${id}_nfo";
}
close L;

$dbh->do('CREATE INDEX IF NOT EXISTS scse ON gff(samplebit,chrid,start,end);') or die $dbh->errstr;
$dbh->commit;
my $read_time = [gettimeofday];

$bit *= 2;
--$bit;
for my $v (1..$bit) {
	use integer;
	for my $b (keys %Log2) {
		my $n=$Log2{$b};
		my $t=$v & $b;
		next if $t == 0;
		next if $v >> $n+1 == 0;
		my $w=($v << $i-$n) & $bit;
		next if $w != 0;
		push @{$Combines{$b}},$v;
		$Summary{$v}=0;
	}
	my @tt=split //,(sprintf "%b",$v);	# "%.${i}b"
	my ($ii,$tt,@Tab)=0;
	for (reverse @tt) {
		if ($_ == 1) {
			++$tt;
			push @Tab,chr(65+$ii);
		}
		++$ii;
	}
	$Bits{$v}=$tt;
	$Tables{$v}=join '',@Tab;
	#warn "$v\t@tt\t$tt\n";
}
#@AllCombines = sort { $Bits{$b} <=> $Bits{$a} } @AllCombines;
#print D "\n[Data]\nChr";
#for (@AllCombines) {
#	print D "\t",$Tables{$_};
#}
#print D "\n";

for my $chr (sort keys %Genes) {
	use integer;
	#%Groups=();
	print STDERR ">$chr   $ChrLen{$chr}\t";
	my $handle=&initchr($ChrLen{$chr});
	my $Gffs=$Genes{$chr};
	for my $b (keys %Log2) {
		for (@{$$Gffs{$b}}) {
			&orbase($handle,$_,$b) for ($$_[1]..$$_[2]);
			#my $bb=getbase($handle,$$_[2]-1);warn "$bb";
		}
	}
	print STDERR "Filled.\t";
	for my $bit (keys %Log2) {
		for (@{$$Gffs{$bit}}) {
			my $id=$$_[0];
			my (%Contiune,%Flag,%Dat);
			for my $pos ($$_[1]..$$_[2]+1) {	# the extra 1 makes ++$Flag{$_} work at the end.
				my $x=getbase($handle,$pos);
				#warn "   $pos\n";
				for (@{$Combines{$bit}}) {
					my $t=$x & $_;
					#warn "$x,$_" if $x != $_ and $x != 15;
					if ($x == $_) {
						++$Contiune{$_};
						#warn "   $pos,$x\t$Contiune{$_}\n" if $x != 15;
					} else {
						if ($Contiune{$_} and $Contiune{$_} >= $opt_l) {
							++$Flag{$_};
							push @{$Dat{$_}},[$pos-$Contiune{$_},$pos-1];
							#warn "$pos\t$_\t$Contiune{$_}\t$Flag{$_}\t@{${$Dat{$_}}[-1]}\n" if $opt_v;
							warn "\n[!]More hits for $id @ $chr:$pos as $_ for $Contiune{$_}\n" if $#{$Dat{$_}} > 0;
						}
						$Contiune{$_}=0;
					}
				}
			}
			# search done.
			if (%Flag) {
				my @R = sort { $Bits{$b} <=> $Bits{$a} } keys %Flag;
				my $Rid=$R[0];
				my $Range=$Dat{$Rid}->[0];
				warn "> $Rid\t@$Range\t$id\t",$$Range[1]-$$Range[0]+1,"\n" if $opt_v;
				my @otherNames=($id);
				for my $bitf (keys %Log2) {
					my $x = $Rid & $bitf;
					next if $x == 0;
					my $name;
					$gth->execute($bitf,$chr,@$Range);
					my $qres = $gth->fetchall_arrayref;
					if ($#$qres == 0) {
						$name=$$qres[0][0];
					} elsif ($#$qres > 0) {
						warn "\n$#$qres more hit(s) for $bitf,$chr,@$Range !\n";
						$name=$$qres[0][0];	# map join if EP
					} else {
						warn "\n[x]No info. for $R[0] -> $bitf,$chr,@$Range\t$x !\n";
						$name='.';
					}
					push @otherNames,$name;
				}
				++$Summary{$R[0]};
				push @{$Groups{$R[0]}{$chr}},\@otherNames;
			}
		}
	}
	warn "Parsed.\n";
	freechr($handle);
=pod
	print D "$chr\t";
	for (@AllCombines) {
		unless ($Groups{$_}) {
			print D "$chr\t";
		}
		print D ;
	}
=cut
}
for (sort { $Bits{$b} <=> $Bits{$a} } keys %Groups) {
	print D "\n[$Tables{$_}]\n";
	my $hash=$Groups{$_};
	for my $chr (sort keys %{$hash}) {
		print D $chr,"\t",join(',',@$_),"\n" for @{$$hash{$chr}};
	}
}

close D;
open( O,'>',$opt_o) or die "[x]Error creating Summary file: $!\n";
for (sort { $Bits{$b} <=> $Bits{$a} } keys %Summary) {
	print O "$Tables{$_}\t$Summary{$_}\n";
}
close O;
my $work_time = [gettimeofday];

$dbh->commit;
$dbh->disconnect;
#unlink $opt_t;
my $stop_time = [gettimeofday];

$|=1;
print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time ),
	" second(s).\n   INSERT used:\t",tv_interval( $start_time, $read_time ),
	" second(s).\n   Parse used:\t",tv_interval( $read_time, $work_time ),
	" second(s).\n   Output used:\t",tv_interval( $work_time, $stop_time )," second(s).\n";
