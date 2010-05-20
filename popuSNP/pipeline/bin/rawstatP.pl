#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use File::Basename;
use Galaxy::ShowHelp;
use FindBin qw($Bin);

$main::VERSION=0.0.1;
my $SCRIPTS="$Bin/../scripts";

=pod
Note:
If a PE Lib come with only one .adapter.list file, all will use that one. If more than 2, only (sort)[0,-1] is used as _1 and _2.

./0rawfq , the path to fq file cannot contain more than 1 "_{$LibName}", since it is searched directly on fullpath.

$ perl -MTree::Suffix -e 'my $tree = Tree::Suffix->new(qw(stringxxxssx stringyx1xxssx axxxssxstring));my @lcs = $tree->lcs;print "[$_]\n" for @lcs;'
[string]

$ perl -MTree::Suffix -e 'my $tree = Tree::Suffix->new(qw(zzzzzzxxaaaaaastringaxxssx tyaaaaaastringzzzzzzyaxxssx zzzzzz1aaaaaaaaxxssxstring));my @lcs = $tree->lcs;print "[$_]\n" for @lcs;'
[zzzzzz]
[aaaaaa]
[axxssx]
[string]

=cut
our $opts='i:o:s:x:v:bqd';
our($opt_i, $opt_s, $opt_o, $opt_v, $opt_b, $opt_q, $opt_d, $opt_x);

our $desc='1.filter fq, 2.stats';
our $help=<<EOH;
\t-i FASTAQ file path (./0rawfq) with *.adapter.list and *.fq or *.fq.gz
\t-s Sample list (sample.lst) in format: /^Sample\\tLib\\tInsertSize\$/
\t-o Project output path (./1fqfilted), will mkdir if not exist
\t-q run qsub automatically
\t-x lib regex for Simulation Mode, undef for Normal Mode (undef)
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='./0rawfq' if ! $opt_i;
$opt_s='sample.lst' if ! $opt_s;
$opt_o='./1fqfilted' if ! $opt_o;
no warnings;
$opt_v=int $opt_v;
use warnings;

$opt_i=~s#/+$##;
$opt_o=~s#/+$##;
# `readlink -f` will be blank if target not exists.
system('mkdir','-p',$opt_o);

die "[x]-i $opt_i not exists !\n" unless -d $opt_i;

print STDERR "From [$opt_i] to [$opt_o] refer to [$opt_s]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (%SampleLib,%LibSample,%LibInsSize);
open SAMPLE,'<',$opt_s or die "Error opening $opt_s: $!\n";
while (<SAMPLE>) {
	chomp;
	s/\r//g;
	my ($sample,$lib,$min,$max)=split /\t/;
	unless (defined $max) {
		$max=int(0.5+$min*1.08);
		$min=int($min*0.92);
	} else { ($min,$max)=sort {$a <=> $b} ($min,$max); }
	push @{$SampleLib{$sample}},$lib;
	$LibInsSize{$lib}=[$min,$max];
	push @{$LibSample{$lib}},$sample;
}
close SAMPLE;
for (keys %LibSample) {
	my $t = scalar @{$LibSample{$_}};
	die "[!]Sample with Lib conflict($t) for Lib:[$_]. Check [$opt_s] !\n" if $t != 1;	# Check if a lib -> more samples
	$LibSample{$_}=${$LibSample{$_}}[0];
}
if ($opt_v) {
	print '-' x 80,"\n";
	for my $k (sort keys %SampleLib) {
		print scalar @{$SampleLib{$k}},"\tSample:[$k] => Lib:[",join('],[',@{$SampleLib{$k}}),"]\n"
	}
	print '-' x 80,"\n";
	for my $k (sort keys %LibInsSize) {
		print "Lib:[$k]\tInsize:[",join(',',@{$LibInsSize{$k}}),"]\n";
	}
	print '-' x 80,"\n";
}

my @files= qx/find $opt_i -iregex '.+\\.\\(fq\\(\\.gz\\)?\\|list\\)'/;	# find ./sprice/0raw/ -iregex '.+\.\(fq\(\.gz\)?\|list\)'
chomp @files;
#my @fq1 = `find $opt_i -name '*.fq'`;	# no need to sort
#my @fq2 = `find $opt_i -name '*.fq.gz'`;	# no need to sort
my @fqs=grep /\.fq/,@files;
my @Adapters=grep /\.list$/,@files;

my (%fqfile2rawfpe,%fq1,%fq2,%fqse,%fqpe);	# no ext.
#@fqfiles = map {[fileparse($_, qr/\.fq/)]} @fq;
#my %AdapterPath2Lists = map {(fileparse($_))[1,0]} @Adapters;
#print "[$_]\n" for map {(fileparse($_))[0,1]} @Adapters;
#print "[$_] => [$AdapterPath2Lists{$_}]\n" for keys %AdapterPath2Lists;

my (%AdapterPath2Lists,%fqname2adapter);
for (sort @Adapters) {	# sort so that _1 comes first
	my ($file, $path) = fileparse($_);
	push @{$AdapterPath2Lists{$path}},$file;
}
if ($opt_v) {
	print '[',join('],[',@{$AdapterPath2Lists{$_}}),'] <= [',$_,"]\n" for sort keys %AdapterPath2Lists;
	print '-' x 80,"\n";
}

for (@fqs) {
#print "$_\t";
	my ($file, $path, $ext) = fileparse($_, qr/\.fq(\.gz)?/);
	next if $file =~ /IndexPooling\d+(_[12])?$/i;
#print "$file, $path, $ext\n";
	$fqfile2rawfpe{$file}=[$path,$ext];
	my ($lib,$ab);
	unless ($opt_x) {
		$file =~ /_([^_]+)(_[12])?$/;	# $1, $2 is local, thus un-useable outside this `unless`.
		$lib = $1; $ab=$2;
	} else {
		$file =~ /($opt_x).*?(_[12])?$/;
		$lib = $1; $ab=$2;
	}
	unless ($ab) {
		push @{$fqse{$lib}},[$file];	# well, same API is better ?
#print "-[$lib] $file\t[$ab]\n" if $opt_v;
		next;
	}
	if ($ab eq '_1') {
		$fq1{$file}=$lib;
		$fqname2adapter{$file}=${$AdapterPath2Lists{$path}}[0] if defined $AdapterPath2Lists{$path};
	}
	if ($ab eq '_2') {
		$fq2{$file}=$lib;
		$fqname2adapter{$file}=${$AdapterPath2Lists{$path}}[-1] if defined $AdapterPath2Lists{$path};
	}
}
%AdapterPath2Lists=();	# useless now
if ($opt_v > 2) {
	print "[$_] => [$fqname2adapter{$_}][${$fqfile2rawfpe{$_}}[0]]\n" for sort keys %fqname2adapter;
	print '-' x 80,"\n";
}
for my $file (keys %fq1) {
	my $file2=$file;
	$file2=~s/_1$/_2/;
	my $lib=$fq1{$file};
	delete $fq1{$file};
	if (defined $fq2{$file2}) {
		push @{$fqpe{$lib}},[$file,$file2];
		delete $fq2{$file2};
	} else {
		warn "[!][$file${$fqfile2rawfpe{$file}}[1]] is not paired with _2 !\n";
		push @{$fqse{$lib}},[$file];
	}
}
for my $file (keys %fq2) {
	my $lib=$fq2{$file};
	delete $fq2{$file};
	warn "[!][$file${$fqfile2rawfpe{$file}}[1]] is not paired with _1 !\n";
	push @{$fqse{$lib}},[$file];
}
%fq1=%fq2=();	# useless now
my %fqbylib;#=%fqpe; If just copy, the value of hash will be the same pointer, thus conflict when pushing se in, which is, se would be pushed in to %fqpe
push @{$fqbylib{$_}},@{$fqpe{$_}} for (keys %fqpe);
push @{$fqbylib{$_}},@{$fqse{$_}} for (keys %fqse);
if ($opt_v) {
	for my $k (sort keys %fqbylib) {
		#print "\n[$k]\n[",join("]\n[",@{$fqbylib{$k}}),"]\n"
		print "[$k] => [";
		for (@{$fqbylib{$k}}) {
			print join(']-[',@$_);
		}
		print "]\n";
	}
	print '-' x 80,"\n";
}
### Got: %fqname2adapter, %fqpe, %fqse, %fqbylib, %SampleLib, %LibSample, %LibInsSize

### 1.filter fq
my ($opath,@cmdlines);
for my $lib (keys %fqbylib) {
	$opath=$opt_o.'/'.$lib.'='.$LibSample{$lib};
	system('mkdir','-p',$opath);
	for my $k (@{$fqbylib{$lib}}) {
		for my $name (@$k) {
			my ($path,$ext)=@{$fqfile2rawfpe{$name}};
			if ($fqname2adapter{$name}) {
				push @cmdlines,"perl $SCRIPTS/fqfilter.pl $path/$fqname2adapter{$name} $path/$name$ext > $opath/$name.fq 2>$opath/$name.nfo";
			} else {
				if ($ext eq '.fq') {
					system('ln','-s',`readlink -nf $path/$name$ext`,"$opath/$name.fq");
				}# else {push @cmdlines,"gzip -dc $path/$name$ext > $opath/$name.fq";}
				push @cmdlines,"perl $SCRIPTS/fqstat.pl $path/$name$ext > $opath/$name.fq 2>$opath/$name.nfo";
			}
		}
	}
}
if ($opt_v > 3) {
	print "[$_]\n" for @cmdlines;
	print '-' x 80,"\n";
}

open SH,'>',$opt_o.'/cmd.lst' or die "[x]Error $!\n";
print SH "$_\n" for @cmdlines;
close SH;
open SH,'>',$opt_o.'/job.sh' or die "[x]Error $!\n";
print SH "#!/bin/sh
#\$ -N \"filter\"
#\$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#\$ -cwd -r y -l vf=276M
#\$ -o /dev/null -e /dev/null
#\$ -S /bin/bash -t 1-",scalar @cmdlines,"
SEEDFILE=${opt_o}/cmd.lst
SEED=\$(sed -n -e \"\$SGE_TASK_ID p\" \$SEEDFILE)
eval \$SEED
";
close SH;
#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
