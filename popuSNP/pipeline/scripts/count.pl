#!/bin/env perl
#use lib '/share/raid010/resequencing/soft/lib';
#use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <in list> <Tag:Column(0):Precision(0),T:C:P> <out prefix>\n";
	exit;
}

my ($list,$arg,$out)=@ARGV;
# 3（depth），9（quality），10（copynum），15（RST_check）
$arg='Depth:3:0,Q:9:0,CopyNum:10:0.005,RST:15:0.0005' if $arg eq '0addcn';
my %Col;
for (split /,/,$arg) {
	s/\s//g;
	my @t=split /:/;
	# Better to change $t[2] to double here.
	$Col{$t[1]}=[$t[0],$t[2]];
	warn '[!] ',join("\t",$t[1],$t[0],$t[2]),"\n";
}

sub getVal($$) {
	my ($p,$d)=@_;
	return int($d/$p)*$p;
}

open L,'<',$list or die "[x]Error on $list: $!\n";
my (%Dat,%Count,$v);
while(my $file=<L>) {
	chomp $file;
	print STDERR "> $file\t";
	open I,'<',$file or die "[x]Error on $file: $!\n";
	while(<I>) {
		chomp;
		my @line=split /\t/;
		for my $k (keys %Col) {
			if ($Col{$k}->[1] == 0) { $v=$line[$k]; }	# faster
			 else { $v=&getVal($Col{$k}->[1],$line[$k]); }	# Well, C will inline while Perl won't.
			++$Dat{$k}->{$v};
			++$Count{$k};	# sum up for percents
		}
	}
	close I;
	warn "done.\n";
}
close L;
for my $k (keys %Dat) {
	open O,'>',join('',$out,$k,$Col{$k}->[0],'.dat') or die "[x]Error: $!\n";	# add $k so that always unique
	print O join("\t",$Col{$k}->[0],'Count','% of '.$Count{$k},'accumulation(%)'),"\n";
	my $suming=0;
	for my $x (sort {$a <=> $b} keys %{$Dat{$k}}) {
		$v=100*$Dat{$k}->{$x}/$Count{$k};
		$suming += $v;
		print O join("\t",$x,$Dat{$k}->{$x},$v,$suming),"\n";
	}
	close O;
}

__END__
# /nas/RD_09C/resequencing/soft/pipeline/population_snp/depth.pl
unless (@ARGV){
	print "<input dir> <output dir> <set nu>\n";
	exit;
}

my @files=`ls $ARGV[0]/*.add_cn`;
chomp @files;

my $input_dir = shift;
my $output_dir = shift;
my $nu_ar  = @ARGV;

for($nu_ar>0)
{
	my $line = shift @ARGV;
	print "deal: $line\n";
	count($output_dir,\@files,$line);
}

sub count
{
	my($out,$files,$setoff) = @_;
	my %depth=();
	my $sum;
open (OUT,">$out/$setoff.depth") or die;
foreach my $f (@$files)
{
	open (F,"$f") or die $!;
	while(my $l=<F> )
	{
		chomp ($l);
		my $dep=(split(/\t+/,$l))[$setoff];
		if (defined $depth{$dep}){$depth{$dep}++;}
		else{$depth{$dep}=1;}
	}
	close F;
}

foreach (sort {$a<=>$b} keys %depth){
	print OUT "$_\t$depth{$_}\n";
}

close OUT;
}
