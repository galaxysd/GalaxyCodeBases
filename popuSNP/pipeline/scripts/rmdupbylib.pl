#!/usr/bin/perl -w
#use threads;
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

######
=pod
Changelog:
0.3.10->9	Li Jun 3 found few days before that the types of mismatch contains only count in the seed.
	However, GuoXs confirmed from LiYR that soapSNP will count mismatch from sequence itself. Thus no correction needed.
	And, soapSNP is able to read soap2 out directly, so in future versions, transformat to soap1 will be skipped.
0.3.9	Guoxs said, tirm will be only 1 at 5', others are 3'
0.3.8	change input format
0.3.7	add cache option
0.3.6	filter out by col.10 in PE for soapSNP mis-read with soap2 -g
0.3.5	$sth4 SELECTs strand to ensure 2 lines out
0.3.4	only -b -o will trim
=cut
######

$main::VERSION=0.3.9;

our $opts='i:o:c:bvmf';
our($opt_i, $opt_o, $opt_c, $opt_v, $opt_b, $opt_m, $opt_f);

our $desc='SoapSort library PCR PE Duplication Remover & Merger (Atom Edition)';
our $help=<<EOH;
\t-i SOAP result list, in format: /^[PS]E\\tid\\t/path/to/soap.file\\n\$/
\t-c Chromosome name, the single one to parse
\t-o Merged output file prefix for a sample, directories must exist
\t   Output file(with path) will be [{-o}.{-c}]
\t-m cache soap file so that no more file reading at output
\t-f Filter out reads with gaps on [soap2 -g]
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

warn "[x]Must specify -i xx.soaplist !\n" unless $opt_i;
warn "[x]Must specify -c 'Chromosome Name' !\n" unless $opt_c;
warn "[x]Must specify -o ./xx !\n" unless $opt_o;
exit 1 unless $opt_i and $opt_o and $opt_c;

my $outfile=${opt_o}.'.'.$opt_c;
print STDERR "From [$opt_i] to [$outfile]\n";
print STDERR "-f on\t" if $opt_f;
print STDERR "-m on\t" if $opt_m;
print STDERR "\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
unless (-s $opt_i) {die "[x]Soaplist [$opt_i] is nothing !\n";}
my (%FILES,%FID);
open LST,'<',$opt_i or die "Error opening $opt_i: $!\n";
while (<LST>) {
	chomp;
	my ($pe,$fid,$file)=split /\t/;
	push @{$FILES{$pe}},$file;
	$FID{$file}=$fid;
}
close LST;
#if ($opt_v) {}
#system('mkdir','-p',$opt_o);
open OUT,'>',$outfile or die "Cannot create $outfile: $!\n";	# better to check before

### Begin SQL
my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0 );
my $dbh = DBI->connect('dbi:SQLite:dbname=:memory:','','',\%attr) or die $DBI::errstr;

my $sql=q/
CREATE TABLE IF NOT EXISTS se
(  soapid TEXT,
   pos INTEGER,
   fileid INTEGER,
   offset INTEGER  );
CREATE TABLE IF NOT EXISTS soap
(  soapid TEXT,
   len INTEGER,
   hit INTEGER,
   sumQasc INTEGER,
   pos INTEGER,
   strand TEXT,
   realpos INTEGER,
   fileid INTEGER,
   offset INTEGER,
   isTrim INTEGER,
   chosen INTEGER  );
/;
# The rowid value can be accessed using one of the special names "ROWID", "OID", or "_ROWID_".
for (split /;/,$sql) {
	next if /^\s*$/;
	$dbh->do($_) or die $dbh->errstr;
}
$dbh->commit;

sub doindex() {
	my $sql=q/
CREATE INDEX IF NOT EXISTS rps ON soap(realpos);
CREATE INDEX IF NOT EXISTS sf ON soap(soapid,fileid);
/;
	for (split /;/,$sql) {
		next if /^\s*$/;
		$dbh->do($_) or warn $dbh->errstr;
	}
	$dbh->commit;
}

my $sthpe = $dbh->prepare( 'INSERT INTO soap ( soapid,pos,fileid,offset,len,hit,sumQasc,strand,realpos,isTrim ) VALUES ( ?,?,?,?,?,?,?,?,?,? )' );
my $sthse = $dbh->prepare( 'INSERT INTO se ( soapid,pos,fileid,offset ) VALUES ( ?,?,?,? )' );
my $sth0 = $dbh->prepare( 'SELECT DISTINCT soapid,fileid FROM soap;' );
my $sth1 = $dbh->prepare( 'SELECT DISTINCT realpos,strand FROM soap WHERE soapid=? AND fileid=? ORDER BY strand ASC;' );	# '+', which comes first, < '-'
my $sth2 = $dbh->prepare( 'SELECT DISTINCT soapid,fileid FROM soap WHERE realpos=? AND strand=?;' );
my $sth3 = $dbh->prepare( 'UPDATE soap SET chosen=1 WHERE soapid=? AND fileid=?;' );
my $sth4 = $dbh->prepare( 'SELECT DISTINCT len,hit,sumQasc,strand FROM soap WHERE soapid=? AND fileid=? ORDER BY strand ASC;' );	# strand must be here as the other maybe same.
my $sth5 = $dbh->prepare( 'SELECT DISTINCT soapid,fileid,offset,\'PE\',isTrim,pos FROM soap WHERE chosen=1 UNION
SELECT DISTINCT soapid,fileid,offset,\'SE\',0,pos FROM se
 ORDER BY pos ASC;' );	# order to make a sorted out
#EXPLAIN QUERY PLAN:
#0|0|TABLE soap
#0|0|TABLE se
# cannot use index on single table

#my $sth5 = $dbh->prepare( 'SELECT DISTINCT soapid,fileid,offset,pos FROM soap WHERE chosen=1 ORDER BY pos ASC;' );	# order to make a sorted out
#my $sth6 = $dbh->prepare( 'SELECT DISTINCT soapid,fileid,offset,pos FROM se ORDER BY pos ASC;' );
# no need to CREATE INDEX IF NOT EXISTS c ON soap(chosen) since not so many dups.

$|=1;
#if ($opt_m) {
my ($o_f,$o_fs,%PSE)=(0,0);
#}
my ($fileid,%FH,$Offset,$the_time,$isTrim)=(0);
my ($allDroped,$maxDroped,$itemsInP,$itemsOutP0,$itemsOutP,$itemsOutS,$itemsInS,$PEtrimed,$SEtrimed)=(0,0,0,0,0,0,0,0,0);	# $itemsIn = $allDroped + $itemsOut
my ($bpsIP,$bpsOP,$bpsIS,$bpsOS)=(0,0,0,0);
my $count;
for my $file (@{$FILES{PE}}) {
	$fileid=$FID{$file};
	$the_time = [gettimeofday];
	printf STDERR ">Parsing\033[32;1m PE %3u @ %11.6f\033[0;0m sec.\n[%s] ",$fileid,tv_interval($start_time, $the_time),$file;
	open $FH{$fileid},'<',$file or die "Error opening PE [$file]: $!\n";
	$Offset=tell $FH{$fileid};
	#print STDERR 'p',$fileid;
	my ($sumQasc,$soapid,$Qstr,$hit,$len,$strand,$chr,$pos,$realpos,@Q,$trim,$trimed,$types);
	$count=0;
	while ($_=readline $FH{$fileid}) {
# soapSNP 2.20 cannot run with -g, so filter out by col. 10
		($soapid,$Qstr,$hit,$len,$strand,$chr,$pos,$types,$trim) = (split(/\t/))[0,2,3,5,6,7,8,9,-2];
		if ($chr ne $opt_c) {
			$Offset=tell $FH{$fileid};
			next;
		}
		if ($opt_f and length $types > 2) {
			print "[v]Skipped:[${fileid}_$soapid]\@Chr:$chr,$pos,$strand\t${types}\n" if $opt_v;
			$Offset=tell $FH{$fileid};
			++$o_f;
			next;
		}
		++$count;	++$itemsInP;	$bpsIP += $len;
		$trimed=$sumQasc=0;	# must reset to 0 cyclely !
		$realpos=$pos;
		if ($strand eq '-') {
			$realpos += $len;	# should be $len-1. So, starting 0. (+ & -, never meets.)
			if ($trim =~ /(\d+)S$/) {
				$realpos += $trimed = 1;	# $1;	# $1 only reset after next /()/
			}
		} elsif ($strand eq '+') {
			if ($trim =~ /^(\d+)S/) {
				$realpos -= $trimed = 1;	# $1;	# Gxs said tirm will be only 1 at 5', others are 3'
			}
		} else {	# elsif ($strand ne '+')
			$realpos=0;
			print STDERR "[!][${fileid}_$soapid] @ Chr:$chr $pos with wrong strand [$strand]\n";
		}
		@Q = unpack('C*', $Qstr);
		$sumQasc += $_ for @Q;
		#@Q = split //,$Qstr;
		#for (@Q) {	$sumQasc += ord;	}
		if ($trimed > 0) {
			++$PEtrimed;
			$isTrim=1;
			print "[v]Trimed:[${fileid}_$soapid]\@Chr:$chr,$pos,$strand\t${trim}->\tΔ:$trimed,l:$len,$realpos:",$pos-$realpos,"\tQ:$sumQasc\t$PEtrimed\n" if $opt_v;
		} else { $isTrim=0; }
		if ($opt_m) {
			$PSE{$fileid}{$Offset}=$_;
		}
		$sthpe->execute($soapid,$pos,$fileid,$Offset,$len,$hit,$sumQasc,$strand,$realpos,$isTrim);
		$Offset=tell $FH{$fileid};
	}
	print STDERR $count,"\n";
	#close $FH{$fileid};
	$dbh->commit;
}
for my $file (@{$FILES{SE}}) {
	$fileid=$FID{$file};
	$the_time = [gettimeofday];
	printf STDERR ">Parsing\033[32;1m SE %3u @ %11.6f\033[0;0m sec.\n[%s] ",$fileid,tv_interval($start_time, $the_time),$file;
	open $FH{$fileid},'<',$file or die "Error opening SE [$file]: $!\n";
	$Offset=tell $FH{$fileid};
	#print STDERR 's',$fileid;
	my ($soapid,$len,$chr,$pos);
	$count=0;
	while ($_=readline $FH{$fileid}) {
		($soapid,$len,$chr,$pos) = (split(/\t/))[0,5,7,8];
		if ($chr ne $opt_c) {
			$Offset=tell $FH{$fileid};
			next;
		}
#		if ($trim =~ /\dS/) { $isTrim=1;	++$SEtrimed; }
#		 else { $isTrim=0; }
		++$count;	++$itemsInS;	$bpsIS += $len;
		if ($opt_m) {
			$PSE{$fileid}{$Offset}=$_;
		}
		$sthse->execute($soapid,$pos,$fileid,$Offset);
		$Offset=tell $FH{$fileid};
	}
	print STDERR $count,"\n";
	$dbh->commit;
}

$the_time = [gettimeofday];
printf STDERR "-Finish Loading\033[32;1m @ %11.6f\033[0;0m sec.\n",tv_interval($start_time, $the_time);
print STDERR 'L';
##FINISH Loading soaps

print STDERR "i\b";
&doindex;	print STDERR 'I';

print STDERR "p\b";
$sth0->execute;	print STDERR 'P';

print STDERR "d\b";
my ($soapid,%mergeids,$aa,$bb,$rres,$mergeid,$qres);	# "my" variable $fileid masks earlier declaration
while ($rres = $sth0->fetchrow_arrayref) {
	($soapid,$fileid)=@$rres;	# soapid,fileid
	$mergeid=$fileid.'_'.$soapid;
	++$mergeids{$mergeid};	# init mark
	next if $mergeids{$mergeid} > 1;	# check mark
	$sth1->execute($soapid,$fileid);	# realpos,strand	'+'<'-'
	$qres = $sth1->fetchall_arrayref;
	if ($#$qres == 1) {	# Normal PE
		($aa,$bb)=@$qres;	# just mid. value of [realpos,strand]; aa->'+', bb->'-'
		$sth2->execute($$aa[0],$$aa[1]);
		$aa=$sth2->fetchall_arrayref;	# soapid,fileid
		$sth2->execute($$bb[0],$$bb[1]);
		$bb=$sth2->fetchall_arrayref;
	} elsif ($#$qres > 1) {	# more than PE, only if file with -t was merged, next
		print "[!]Reads:[",1+$#$qres,"], not a PE soap for [$mergeid] ! Ignored.(Maybe soap file with -t was merged with cat or cp)\n";
		next;
	} elsif ($#$qres == 0) {	# SE, impossable, next
		if ($opt_f) {
			print "[!]Reads:[1], is a SE soap for [$mergeid] at strand:$$qres[0][1] ! Ignored. (This is OK since -f used.)\n";
			++$o_fs;
		} else {
			print "[!]Reads:[1], is a SE soap for [$mergeid] at strand:$$qres[0][1] ! Ignored.(Please use given program to generate input list file, ONLY those soap -b -o are PE.)\n";	# Well, we can INSERT INTO se. But, is robust so important when I am tired ?
		}
		next;
	} else { die "[x]System Error. (Memory or DBI)"; }	# -1, no such soapid

	my (%c,@dup,$dmergeid);	#=();
	for (@$aa,@$bb) {
		$dmergeid=$$_[1].'_'.$$_[0];
		next if $dmergeid eq $mergeid;
		++$c{$dmergeid};
	}
	for (keys %c) {
		if ($c{$_}>1) {
			push @dup,$_;
			++$mergeids{$_};	# make mark
			print "[!][$_] come with $c{$_} > 2 results.\n" if $c{$_}>2;
		}
	}
	if ($#dup == -1) {
		$sth3->execute($soapid,$fileid);
		++$itemsOutP0;
		print "[v]+>[$mergeid]\n" if $opt_v;
	} else {
		push @dup,$mergeid;
		$allDroped += $#dup;
		$maxDroped = $#dup if $maxDroped < $#dup;
# There do are cases that: "DBD::SQLite::db selectall_arrayref failed: Expression tree is too large (maximum depth 1000)"
		my ($csoapid,$cfileid,%candidates,$dsoapid,$dfileid);
		for $dmergeid (@dup) {
			($dfileid,$dsoapid)=split /_/,$dmergeid;
			$sth4->execute($dsoapid,$dfileid);	# len,hit,sumQasc,strand
			$qres = $sth4->fetchall_arrayref;	# ORDER BY strand ASC
			if ($#$qres != 1) {
				print "[!]Reads:[",1+$#$qres,"], not a PE soap for [$dmergeid] ! Ignored.(Maybe soap file with -t was merged with cat or cp)\n" if $#$qres > 1;	# since $mergeid has pass the check, always with at least 1 out.
				if ($#$qres == 0) {
					print "[!]Reads:[1], is a SE soap for [$dmergeid] at strand:$$qres[0][3] ! Ignored.(2)\n";	# EP ?
					#print "[!]Candidates are: [$dmergeid] of ";
					#print "[$_]\t" for (@dup);
					#print "\n";
				}
				# no more check for $#$qres == -1, we are smoking. Well, half smoking.
				next;	# also pass $#$qres = -1 or 0
			}
			$candidates{$dmergeid}=[ @$qres ];	# '+' -> [], '-' -> [len,hit,sumQasc]
		}
####		$dmergeid = &score_dups(\%candidates);	# Let's in_line ...
######	sub score_dups ($)
### Way1: for speed. Output should be the same as Way0.
	my (%len,%hit,%sumQasc,$len,$hit,$sumQasc,@t,$outT);
	for (keys %candidates) {
		$len = $candidates{$_}->[0]->[0] + $candidates{$_}->[1]->[0];
		push @{$len{$len}},$_;
	}
	@t=(sort {$b <=> $a} keys %len);	# sort numerically descending
	if ($#{$len{$t[0]}}==0) {
		if ($opt_v) {
			$outT=0;
			print "[v]->:";
			for $len (@t) {
				for (@{$len{$len}}) {
					printf '[v]drop%3u:',$outT if $outT > 0;
					print " [$_],len:$len\n";
					++$outT;
				}
			}
		}
		$dmergeid=${$len{$t[0]}}[0];
		goto RETURN;	# return $dmergeid;
	}
###
	for ( @{$len{$t[0]}} ) {
		$hit=$candidates{$_}->[0]->[1];
		push @{$hit{$hit}},$_;
	}
	@t=(sort {$a <=> $b} keys %hit);	# ASC
	if ($#{$hit{$t[0]}}==0) {
		if ($opt_v) {
			$outT=0;
			print "[v]->:";
			for $hit (@t) {
				for (@{$hit{$hit}}) {
					printf '[v]drop%3u:',$outT if $outT > 0;
					print " [$_],hit:$hit\n";
					++$outT;
				}
			}
		}
		$dmergeid=${$hit{$t[0]}}[0];
		goto RETURN;	# return $dmergeid;
	}
###
	for ( @{$hit{$t[0]}} ) {
		$sumQasc = $candidates{$_}->[0]->[2] + $candidates{$_}->[1]->[2];
		push @{$sumQasc{$sumQasc}},$_;
	}
	@t=(sort {$b <=> $a} keys %sumQasc);	# DESC
	if ($opt_v) {
		$outT=0;
		print "[v]->:";
		for $sumQasc (@t) {
			for (@{$sumQasc{$sumQasc}}) {
				printf '[v]drop%3u:',$outT if $outT > 0;
				print " [$_],sumQasc:$sumQasc\n";
				++$outT;
			}
		}
	}
	$dmergeid=${$sumQasc{$t[0]}}[0];
	goto RETURN;	# return $dmergeid;
### End Way1.
=pod
### Way0: for debug.
	my (%len,%hit,%sumQasc,$len,$hit,$sumQasc,@t,$outT);
	for (keys %candidates) {
		$len = $candidates{$_}->[0]->[0] + $candidates{$_}->[1]->[0];
		$hit = $candidates{$_}->[0]->[1];
		print "[!]$_ with different hits.\n" if $candidates{$_}->[0]->[1] != $candidates{$_}->[1]->[1];
		$sumQasc = $candidates{$_}->[0]->[2] + $candidates{$_}->[1]->[2];
		push @{$len{$len}},$_;
		push @{$hit{$hit}},$_;
		push @{$sumQasc{$sumQasc}},$_;
	}


### End Way0.
=cut
RETURN:
######	return $dmergeid;
		%candidates=();
		($cfileid,$csoapid)=split /_/,$dmergeid;
		$sth3->execute($csoapid,$cfileid);
		++$itemsOutP0;
	}
}
$dbh->commit;
#$dbh->do('CREATE INDEX IF NOT EXISTS c ON soap(chosen);'); if I just query once, no time benifit. http://stackoverflow.com/questions/1310306/will-i-save-any-time-on-a-index-that-selects-only-once
print STDERR 'D';
##FINISH deDup.

print STDERR "o\b";
my ($red,$redid,@aline);
#open OUT,'>',$outfile;	# opened
$sth5->execute;	# PE & SE	soapid,fileid,offset,'PE'/'SE',isTrim,pos
my ($pres,$len);
my ($Psoapid,$Pfileid,$Poffset,$PE,$Ppos);
while ( $pres=$sth5->fetchrow_arrayref ) {
	($Psoapid,$Pfileid,$Poffset,$PE,$isTrim,$Ppos)=@$pres;
	if ($opt_m) {
		$red=$PSE{$Pfileid}{$Poffset};
	} else {
		seek $FH{$Pfileid},$Poffset,0;
		$red=readline $FH{$Pfileid};
	}
	@aline=split(/\t/,$red);
	($redid,$len)=@aline[0,5];
	push @aline,$isTrim;
	$red=join "\t",@aline[0..9,-1];
	if ($PE eq 'PE') {
		$bpsOP += $len;	++$itemsOutP;
	} else {	# $PE eq 'SE'
		$bpsOS += $len;	++$itemsOutS;
	}
	print "[!]$Pfileid $PE [$Psoapid] <> [$redid], soap file changed.\n" if $Psoapid ne $redid;
	print OUT $Pfileid,'_',$red,"\n";	# rename soapid with prefix
}
close OUT;
print STDERR 'O';
##FINISH Output.
$the_time = [gettimeofday];
printf STDERR "\n-Finish Merge\033[32;1m @ %13.6f\033[0;0m sec.\n",tv_interval($start_time, $the_time);

close $_ for (values %FH);	# SOAP files closed.
$dbh->commit or warn $dbh->errstr;	# only 1 lib. No cycle, no DROP.
$dbh->disconnect;	# in memory db, no file to delete.
### End SQL
#$allDroped,$maxDroped,$itemsInP,$itemsOutP,$itemsInS,$bps{IP,OP,IS}
my ($itemsIn,$itemsOut,$bpsI,$bpsO);
$itemsIn=$itemsInP+$itemsInS;
$itemsOut=$itemsOutP+$itemsOutS;
$bpsI=$bpsIP+$bpsIS;
$bpsO=$bpsOP+$bpsOS;
my $report = "\033\\[32;1m\nIn Chr:[$opt_c]\t[PE,SE]\n
[${itemsInP},$itemsInS] ($itemsIn) lines from SOAP file(s). Trimed [${PEtrimed},${SEtrimed}].
[$o_f] PE reads filtered, causing [$o_fs] reads turn to SE and ignored.

PE pairs: [$itemsOutP0] chosen, [$allDroped] dropped, Max dropped per dupilcate is [$maxDroped].

[$itemsOutP,$itemsOutS] ($itemsOut) lines merged.

[$bpsIP,$bpsIS] ($bpsI) bps read in, [$bpsOP,$bpsOS] ($bpsO) bps write out.\033\\[0;0m
\n
!!!\t$opt_c\t$bpsO\t$bpsI\t$itemsIn\t$itemsOut\t$itemsInP\t$itemsInS\t$itemsOutP\t$bpsIP\t$bpsIS\t$bpsOP\t$maxDroped\t$allDroped\t!!!\n\n";
$report =~ s/\[/[\033[0;0m/g;
$report =~ s/\]/\033[32;1m]/g;
$report =~ s/\\\[\033\[0;0m/[/g;
$report =~ s/\\\033\[32;1m\]/]/g;
$report =~ s/\(/(\033[0;0m/g;
$report =~ s/\)/\033[32;1m)/g;
$report =~ s/\\\(\033\[0;0m/(/g;
$report =~ s/\\\033\[32;1m\)/)/g;
print STDERR $report;
#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
print "done !\n"
__END__
