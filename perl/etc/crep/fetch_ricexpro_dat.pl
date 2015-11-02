#!/usr/bin/perl -w
use strict;
use warnings;
#use IO::Handle;
#use Time::HiRes qw ( gettimeofday tv_interval );
#use Term::ANSIColor qw(:constants);
use LWP::UserAgent;
#use HTTP::Request::Common qw(POST);
use HTML::TreeBuilder::XPath;
#use LWP::Simple qw ( get );
#use HTTP::Cookies;
#use Text::CSV_XS;
use HTML::TableExtract;
use Data::Dump;

my @ArrayList = qw {RXP_1006 RXP_1009 RXP_1010 RXP_1008 RXP_1007 RXP_1012};
my $Order = '0 hr (Cy3)|0 hr (Cy5)|1 hr (Cy3)|1 hr (Cy5)|3 hr (Cy3)|3 hr (Cy5)|6 hr (Cy3)|6 hr (Cy5)|12 hr (Cy3)|12 hr (Cy5)';

open I,'<','resLOC.anno' or die $!;
open O,'>','resLOC.cydat' or die $!;
print O join("\t",'# FeatureNum',@ArrayList),"\n## $Order\n";

my %Acc2Loc;
while (<I>) {
	my ($tigLOC,undef,$FeatureNums) = split /\t/;
	#next unless defined $FeatureNums;
	my @FeatureNums = split /\|/,$FeatureNums;
	next unless @FeatureNums > 0;
	#print "$tigLOC,@FeatureNums\n";
	++$Acc2Loc{$_}{$tigLOC} for @FeatureNums;
}

my @AccNums = sort {$a <=> $b} keys %Acc2Loc;
#ddx \%Acc2Loc;
#print "@AccNums\n";
print join(',',@ArrayList),"\n";

my $ua = LWP::UserAgent->new;
$ua->agent("Mozilla/5.0");

for my $FeatureNum (@AccNums) {

#my $FeatureNum = 4865;
#my $Exp = 'RXP_1006';
my %ArrayDat;
for my $Exp (@ArrayList) {

my $req = HTTP::Request->new(GET => "http://ricexpro.dna.affrc.go.jp/$Exp/view-plot-data.php?featurenum=$FeatureNum");
my $res = $ua->request($req);
# check the outcome
if (! $res->is_success) {
		die "Error: Login as Guest failed !\n";
}

#ddx $res;
	my $tmp_html = $res->content;
	#print $tmp_html;
	$tmp_html = (grep /\<table\>\<tr\>\<th\>\<\/th\>\<th\>0 hr \(Cy3\)\<\/th\>/,(split /\n/,$tmp_html))[0];
#print $tmp_html,"\n";
=pod
		<table><tr><th></th><th>0 hr (Cy3)</th><th>0 hr (Cy5)</th><th>1 hr (Cy3)</th><th>1 hr (Cy5)</th><th>3 hr (Cy3)</th><th>3 hr (Cy5)</th><th>6 hr (Cy3)</th><th>6 hr (Cy5)</th><th>12 hr (Cy3)</th><th>12 hr (Cy5)</th></tr><tr><tr><th>rep1</th><td>1489.2</td><td>1468.1</td><td>2461.9</td><td>2467.7</td><td>2531.5</td><td>2546.4</td><td>3237.8</td><td>4325</td><td>2353</td><td>2589.5</td><tr><th>rep2</th><td>1584.1</td><td>1701.6</td><td>1875.8</td><td>2041.3</td><td>2742.3</td><td>2961.4</td><td>3275</td><td>3770</td><td>3569</td><td>4014</td><tr><th>mean</th><td>1536.6</td><td>1584.9</td><td>2168.9</td><td>2254.5</td><td>2636.9</td><td>2753.9</td><td>3256.4</td><td>4047.5</td><td>2961</td><td>3301.7</td><tr><th>median</th><td>1536.6</td><td>1584.9</td><td>2168.9</td><td>2254.5</td><td>2636.9</td><td>2753.9</td><td>3256.4</td><td>4047.5</td><td>2961</td><td>3301.7</td></tr></table>
=cut
	my $te = HTML::TableExtract->new();	#headers => [( undef, '0 hr (Cy3)', '0 hr (Cy5)', '1 hr (Cy3)', '1 hr (Cy5)', '3 hr (Cy3)', '3 hr (Cy5)', '6 hr (Cy3)', '6 hr (Cy5)', '12 hr (Cy3)', '12 hr (Cy5)' )] );
	$te->parse("$tmp_html");
	my $ts = $te->first_table_found();
	my @dat;
	unless ($ts) {
		warn "\n[!Table]$tmp_html\n";
	}
#ddx $ts;
#   "|0 hr (Cy3)|0 hr (Cy5)|1 hr (Cy3)|1 hr (Cy5)|3 hr (Cy3)|3 hr (Cy5)|6 hr (Cy3)|6 hr (Cy5)|12 hr (Cy3)|12 hr (Cy5)",
#   "rep1|1489.2|1468.1|2461.9|2467.7|2531.5|2546.4|3237.8|4325|2353|2589.5",
#   "rep2|1584.1|1701.6|1875.8|2041.3|2742.3|2961.4|3275|3770|3569|4014",
#   "mean|1536.6|1584.9|2168.9|2254.5|2636.9|2753.9|3256.4|4047.5|2961|3301.7",
#   "median|1536.6|1584.9|2168.9|2254.5|2636.9|2753.9|3256.4|4047.5|2961|3301.7",
	foreach my $row ($ts->rows) {
		my @rowdat = @{$row};
		my $title = shift @rowdat;
		unless (defined $title) {
			my $tmp = join("|",@rowdat);
			die if $tmp ne $Order;
			next;
		}
		push @dat,join("|",$title,@rowdat) if $title =~ /^(rep\d+|mean)$/;
	}
	$ArrayDat{$Exp}=\@dat;
	$req = undef;
#ddx \@dat;
}	# $Exp
print $FeatureNum," ---\n";
ddx \%ArrayDat;
	my $tmp;
	for my $Exp (@ArrayList) {
		$tmp .= "\t" . join(',',@{$ArrayDat{$Exp}} );
	}
	$tmp = $FeatureNum . $tmp;
	print O "$tmp\n";
}	# $FeatureNum
close I;
close O;

__END__
