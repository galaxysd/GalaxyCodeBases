#!/bin/env perl
use strict;
use warnings;
#use DBI;
use Galaxy::ShowHelp;

$main::VERSION=0.1.4;

our $opts='i:o:c:m:p:bvne';
our($opt_i, $opt_o, $opt_c, $opt_v, $opt_b, $opt_n, $opt_e, $opt_m, $opt_p);

our $desc='Genome creater for Single analyse';
our $help=<<EOH;
\t-c Consensus file with refbase
\t-p final SNP file (undef for directly from CNS)
\t-i filtered Indel file
\t-m merge.list for dividing mixed Chrs if any
\t-o Output file prefix for the sample, directories must exist
\t-n mask noncoverage bases to "n", otherwise little cased
\t-e skip hete indels, otherwise both homo and hete(little cased) used
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

warn "[x]Must specify -c xx.consensus.txt !\n" unless $opt_c;
warn "[x]Must specify -i xx.indel.txt.filter !\n" unless $opt_i;
warn "[x]Must specify -o ./xx !\n" unless $opt_o;
exit 1 unless $opt_i and $opt_o and $opt_c;
if ($opt_m) {
    die "[x]Invalid merge.list: [$opt_m] !\n" unless (-s $opt_m);
}
if ($opt_p) {
    die "[x]Invalid SNP file: [$opt_p] !\n" unless (-s $opt_p);
}

print STDERR "From [$opt_c][$opt_i] to [$opt_o]\n ";
print STDERR "merge.list is [$opt_m]\n " if $opt_m;
print STDERR "SNP file is [$opt_p]\n " if $opt_p;
print STDERR "noncoverage mask to n.   " if $opt_n;
print STDERR "hete indels will be skiped." if $opt_e;
warn "\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my (%Genome,%Depth,%Merged,%SNP,$snpChr);
#my %attr = (
#    RaiseError => 0,
#    PrintError => 1,
#    AutoCommit => 0
#);
#my $dbh = DBI->connect('dbi:SQLite:dbname=:memory:','','',\%attr) or die $DBI::errstr;
#my $sql=q/
#CREATE TABLE IF NOT EXISTS merge
#(  chrid TEXT,
#   scaffold TEXT,
#   start INTEGER,
#   end INTEGER  );
#/;
#$dbh->do($sql) or die $dbh->errstr;
#$dbh->commit;
#my $sthi = $dbh->prepare( 'INSERT INTO merge ( chrid,scaffold,start,end ) VALUES ( ?,?,?,? )' );
#my $stho = $dbh->prepare( 'SELECT DISTINCT scaffold,start,end FROM merge WHERE chrid=? AND ? BETWEEN start AND end' );

if ($opt_m) {
    open M,'<',$opt_m or die "Error opening $opt_m: $!\n";
    while (<M>) {
	chomp;
	my ($chrid,$scaffold,$start,$end)=split /\t/;
	#$sthi->execute($chrid,$scaffold,$start,$end);
	$Merged{$chrid}{$start}=[$scaffold,$end];
    }
    close M;
    warn "\n[!]merge.list done.\n";
}
if ($opt_p) {
    open M,'<',$opt_p or die "Error opening $opt_p: $!\n";
    while (<M>) {
	#chomp;
	my ($chrid,$pos,undef,$cons)=split /\t/;
	#$sthi->execute($chrid,$scaffold,$start,$end);
	$SNP{$chrid}{$pos}=$cons;
    }
    close M;
    warn "\n[!]SNP file done.\n";
}

open CNS,'<',$opt_c or die "Error opening $opt_c: $!\n";
my ($chr,$pos,$best,$depth,$ref,$indel,$bases,$homhet,$strand,$t);
while (<CNS>) {
    chomp;
    ($chr,$pos,$ref,$best,undef,undef,undef,undef,undef,undef,undef,undef,undef,$depth) = split /\t/;
    unless ($ref) {   # if file not completed.
	print STDERR '.';
	next;
    }
    if ($depth==0) {
	unless ($opt_n) {$best=lc $ref;}
	 else {$best=$ref='n';}
    } else {$Depth{$chr}->{$pos}=chr $depth;} # string is smaller than double
# WARING: Something might be wrong if $depth > 255 !
# However, people should not have such much money.
#$best = $ref; # DEBUG ONLY !!!
    if ($opt_p) {
    	if ($SNP{$chr}{$pos}) {$Genome{$chr}->[$pos]=$best;}
    	 else {$Genome{$chr}->[$pos]=$ref;}
    	#warn '.' if $best ne $SNP{$chr}{$pos};
    } else { $Genome{$chr}->[$pos]=$best; }
    #$Depth{$chr}->[$pos]=$depth;
}
close CNS;
warn "[!]CNS done.\n";

open INDEL,'<',$opt_i or die "Error opening $opt_i: $!\n";
while (<INDEL>) {
    ($chr,$pos,$indel,$bases,$strand,$homhet,undef,$depth,$t) = split /\t/;
    unless ($t) {   # if file not completed.
	print STDERR '\'';
	next;
    }
    next unless exists $Genome{$chr};
    if ($homhet eq 'hete') {
	unless ($opt_e) {$bases = lc $bases;}
	 else {next;}
    }
    if ($indel =~ /^([ID])(\d+)$/) {
	if ($1 eq 'I') {
	    $Genome{$chr}->[$pos] .= $bases;  # if length($Genome{$chr}->[$t]) == 1
	    $Depth{$chr}{$pos} .= chr($depth) x $2;
	} else {
	    my $i=$2-1;
	    #++$pos;    # deletion including it.
	    if ($opt_v) {
		print "[D]$indel:${1}_${2}@","$pos $strand\t";
		#print $Genome{$chr}->[$_] for ($pos..$pos+$i);
		print(($t=$Genome{$chr}->[$_])?$t:'.') for ($pos-5..$pos-1);
		print '<';
		print(($t=$Genome{$chr}->[$_])?$t:'.',' ') for ($pos..$pos+$i);
		print '>';
		print(($t=$Genome{$chr}->[$_])?$t:'.') for ($pos+$i+1..$pos+$i+5);
		print ' - ',$bases,"\n";
	    }
	    for $t ($pos..$pos+$i) {
		$Genome{$chr}->[$t]='';
		delete $Depth{$chr}{$pos};# if exists $Depth{$chr}{$pos};
	    }
	}
    } else {print STDERR '|';next;}
}
close INDEL;
warn "[!]INDEL done.\n";

my ($i,$dep,$scaf,$name);
for $chr (keys %Genome) {
    $Genome{$chr}->[0]='';  # no longer undef
    $t=0;
    $i=-1;
    unless (exists $Merged{$chr}) { # exists should faster than defined ?
	$name=$opt_o.'.'.$chr;
	open OUT,'>',$name.'.fa' or die "Error opening ${name}.fa: $!\n";
	print OUT ">$chr\n";
	open DEP,'>',$name.'.dep' or die "Error opening ${name}.dep: $!\n";
	print DEP ">$chr\n";
	#$t=0;
	#$i=-1;
	for (@{$Genome{$chr}}) {
	    ++$i;
	    next if $_ eq '';
	    ++$t;
	    print OUT $_;
	    $dep=$Depth{$chr}{$i} or $dep="\0";
	    $dep = join(' ',map ord,split //,$dep);
	    #print DEP '[',$dep,']';
	    if ($t > 79) {
		$t=0;
		print OUT "\n";
		print DEP "$dep\n";
	    } else {print DEP $dep,' ';}
	}
	close OUT;
	close DEP;
    } else {
	$t=$opt_o.'.'.$chr;
	system('mkdir','-p',$t);
	my ($end,$scafname,$toPrint)=(-1,'',0);	# -1 to fix 'print() on unopened filehandle OUT at ./conasm.pl line 178 & 179 of elsif ($i == $end).'
	$scaf=$Merged{$chr};
	for (@{$Genome{$chr}}) {
	    ++$i;
	    if (exists $$scaf{$i}) {	# Start point
		($scafname,$end)=@{$$scaf{$i}};
		$name=$opt_o.'.'.$chr.'/'.$scafname;
		open OUT,'>',$name.'.fa' or die "Error opening ${name}.fa: $!\n";
		print OUT ">$scafname\n";
		open DEP,'>',$name.'.dep' or die "Error opening ${name}.dep: $!\n";
		print DEP ">$scafname\n";
		$toPrint=1;	# when set to 0, no output for manmade junctions.
	    }
	    if ($i != $end) {
		next if ($_ eq '' or $toPrint==0);
		++$t;
		print OUT $_;
		$dep=$Depth{$chr}{$i} or $dep="\0";
		$dep = join(' ',map ord,split //,$dep);
		#print DEP '[',$dep,']';
		if ($t > 79) {
		    $t=0;
		    print OUT "\n";
		    print DEP "$dep\n";
		} else {print DEP $dep,' ';}
	    } else {	# End point
		if ($_ ne '') {
		    print OUT $_,"\n";
		    $dep=$Depth{$chr}{$i} or $dep="\0";
		    $dep = join(' ',map ord,split //,$dep);
		    print DEP "$dep\n";
		} else {
		    print OUT "\n";
		    print DEP "\n";
		}
		close OUT;
		close DEP;
		$t=$toPrint=0;
		next;
	    }
	}
    }
}
warn "[!]Output done.\n";

__END__
vf=11.9g for 1.3G Feb 20 07:23 consensus/QRS9.Gm03.txt, 47,781,076

find ./consensus/ -name '*.txt'|perl -lane '$m="";$vf="14g";$a=(split /\//)[-1];$b=(split /\./,$a)[0];open O,">./sh/$a.sh";if ($a=~/SGm/) {$m="-m soybean.merge.list ";$vf="4.7g"};print O "#\$ -N \"C_$a\"\n#\$ -cwd -r y -l vf=$vf,p=1\n#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH\n#\$ -o /dev/null -e /dev/null\n./conasm.pl ${m}-i ./Indel/${b}.indel.txt.filter -c $_ -bno ./out/$b 2>./log/$a.log";close O'

$ cat sh/QRS21.Gm13.txt.sh
#$ -N "C_QRS21.Gm13.txt"
#$ -cwd -r y -l vf=14g,s_core=1
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#$ -o /dev/null -e /dev/null
./conasm.pl -i ./Indel/QRS21.indel.txt.filter -c ./consensus/QRS21.Gm13.txt -bo ./out/QRS21 2>./log/QRS21.Gm13.txt.log

./conasmS.pl -c P142cns/P142.chromosome01.cns -p P142snp/P142.chromosome01.snp -i P142indel.list.filter -ne -o oP142j/m

cat chrorderNip | perl -lane '$p="142";$a=$_;open O,">./sh/do$p$a.sh";print O "#\$ -N \"C_$a\"\n#\$ -cwd -r y -l vf=12g,p=1\n#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH\n#\$ -o /dev/null -e /dev/null\n./conasmS.pl -c P${p}cns/P${p}.$a.cns -p P${p}snp/P${p}.$a.snp -i P${p}indel.list.filter -nbe -o oP${p}j/m 2>./log/${p}_$a.log";close O'

cat chrorder9311 | perl -lane '$p="143";$a=$_;open O,">./sh/do$p$a.sh";print O "#\$ -N \"C_$a\"\n#\$ -cwd -r y -l vf=12g,p=1\n#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH\n#\$ -o /dev/null -e /dev/null\n./conasmS.pl -c P${p}cns/P${p}.$a.cns -p P${p}snp/P${p}.$a.snp -i P${p}indel.list.filter -nbe -o oP${p}i/m 2>./log/${p}_$a.log";close O'

./conasm.pl -m soybean.merge.list -i ./Indel/QRS29.indel.txt.filter -c ./consensus/QRS29.SGm1.txt -bno ./out/QRS29 2>./log/QRS29.SGm1.txt.log
