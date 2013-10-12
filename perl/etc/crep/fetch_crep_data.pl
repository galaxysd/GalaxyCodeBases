#!/usr/bin/perl -w
use strict;
use warnings;
#use IO::Handle;
use Time::HiRes qw ( gettimeofday tv_interval );
use Getopt::Std;
use Term::ANSIColor qw(:constants);
use LWP::UserAgent;
use LWP::Simple qw ( get );
use HTTP::Cookies;
use Text::CSV_XS;
#use Data::Dumper;

$main::VERSION=1.3.2;

if (@ARGV == 0) {
	&VERSION_MESSAGE();
	&HELP_MESSAGE();
	die "\n";
}
our($opt_p, $opt_o, $opt_n, $opt_s);

$Getopt::Std::STANDARD_HELP_VERSION=1;
getopts('p:o:n:s:');
#print "$opt_i, $opt_o, $opt_e\n";

sub HELP_MESSAGE() {
	#my ($scr) = ($0 =~ m,([^/\\]+)$,);
	my $help=<<EOH;
\t-p NetAffx Rice Probe Set Information (.psi) file [./Rice.psi]
\t-o Output file (crep_all_tsv.txt)
\t-n Number of ProbeSet to retrieve at one time (50 <= 250)
\t-s To Start at . (for Rice 1..57381)
EOH
	print STDERR <<EOH;
\nUsage: \033[0;1m$0\033[0;0m [-OPTIONS [-MORE_OPTIONS]] [--] [PROGRAM_ARG1 ...]

The following single-character options are accepted:
\033[32;1m$help\033[0;0m
Options may be merged together.  -- stops processing of options.
Space is not required between options and their arguments.
EOH
}
sub VERSION_MESSAGE() {
	my $perlv = $];
	$perlv = sprintf "%vd", $^V if $] >= 5.006;
		my $ver = sprintf "%vd", $main::VERSION;
		my ($scr) = ($0 =~ m,([^/\\]+)$,);
	print STDERR <<EOH;
$scr version \033[0;1m$ver\033[0;0m, running under Perl version $perlv.
EOH
}

$opt_p='./Rice.psi' if ! defined $opt_p;
$opt_o='crep_all_tsv.txt' if ! defined $opt_o;
$opt_n=50 if ! defined $opt_n;
$opt_n=int($opt_n);
$opt_n=1 if $opt_n < 1;
$opt_n=500 if $opt_n > 500;
$opt_s=1 if ! defined $opt_s;
$opt_s=int($opt_s);
$opt_s=1 if $opt_s < 1;
&VERSION_MESSAGE();
warn "\nInput PSI: \033[32;1m$opt_p\033[0;0m\nOutput: \033[32;1m$opt_o\033[0;0m\nRetrieve per Turn: \033[32;1m$opt_n\033[0;0m\n\n";
warn "Start at \033[32;1m$opt_s\033[0;0m\n";
#warn 'Check fasta id matches: ',GREEN,BOLD,$opt_c?'Enabled':'Disabled',RESET,"\n";

my $start_time = [gettimeofday];
################### MAIN ####################
my (@psid,$tmp)=();
open PSI,'<',$opt_p or die "Error: $!\n";
while (<PSI>) {
	next if /^#/;
	push @psid,(split /\t/,$_)[1];
	#print "[$tmp_line]\n";
}
close PSI;
#for (@psid) {
#	print "[$_]\t";
#}
$tmp = $#psid + 1;
warn GREEN,BOLD,$tmp,RESET," ProbeSets loaded.\n";
$tmp = $opt_s-1;
if ($#psid >= $tmp-1) {
	splice(@psid,0,$tmp);
}
$tmp = $#psid + 1;
warn "Remains: \033[32;1m$tmp\033[0;0m.\n\n";
################### Login ####################
my $ua = LWP::UserAgent->new;
my ($req,$res,$count);
$ua->cookie_jar(HTTP::Cookies->new(file => "lwpcookies.txt", ignore_discard => 1,
								   autosave => 1));

$req = HTTP::Request->new(GET => 'http://crep.ncpgr.cn/crep-cgi/login2.pl?User_Name=guest&amp;Password=57510426775c5b0f');
# send request
$res = $ua->request($req);
# check the outcome
if (! $res->is_success) {
		die "Error: Login as Guest failed !\n";
}

################### Retrieve ####################
$req = HTTP::Request->new(POST => 'http://crep.ncpgr.cn/crep-cgi/element_multi_chronologer_view.pl');
$req->content_type('application/x-www-form-urlencoded');
my $ps_count = $#psid + 1;
my ($with_head,$to_retrieve,$csv_data_ref,$percent)=(1);
if ($opt_s==1) {
	open OUTPUT,'>',$opt_o or die "Error: $!\n";
} else {
	open OUTPUT,'>>',$opt_o or die "Error: $!\n";
	$tmp=tell OUTPUT;
	warn "$tmp";
	if ($tmp > 1) {
		$with_head=0;
	}
}
my $csv_subject = Text::CSV_XS->new;
print STDERR '==  0.00 %';
$tmp=-1;
while ($#psid >= $opt_n-1) {
=pod
	# login everytime to avoid session out.
	$req = HTTP::Request->new(GET => 'http://crep.ncpgr.cn/crep-cgi/login2.pl?User_Name=guest&amp;Password=57510426775c5b0f');
	$res = $ua->request($req);
	die "Error: Login as Guest failed !\n" if ! $res->is_success;

	$req = HTTP::Request->new(POST => 'http://crep.ncpgr.cn/crep-cgi/element_multi_chronologer_view.pl');
	$req->content_type('application/x-www-form-urlencoded');
=cut
	$to_retrieve = 'unique_id=' . join ' ',splice(@psid,0,$opt_n);
	#warn "[$to_retrieve]\n";
	$req->content($to_retrieve);
	$to_retrieve='';
	$res = $ua->request($req);
	($csv_data_ref,$count)=&crep_retrieve(\$res,$with_head);
	die "Error: $count <> $opt_n, Server Overflow !\n" if $count != $opt_n;
	$with_head=0;
	print OUTPUT $$csv_data_ref;
	$percent=sprintf "%5.2f",($ps_count-$#psid-1)*100/$ps_count;
	print STDERR "\b"x10, ($#psid % 2) ? ('>='):('=>')," $percent %";
	$tmp=1;
}
if ($#psid >= 0) {
=pod
	# login everytime to avoid session out.
	$req = HTTP::Request->new(GET => 'http://crep.ncpgr.cn/crep-cgi/login2.pl?User_Name=guest&amp;Password=57510426775c5b0f');
	$res = $ua->request($req);
	die "Error: Login as Guest failed !\n" if ! $res->is_success;
=cut
	$req = HTTP::Request->new(POST => 'http://crep.ncpgr.cn/crep-cgi/element_multi_chronologer_view.pl');
	$req->content_type('application/x-www-form-urlencoded');
	$to_retrieve = 'unique_id=' . join ' ', @psid;
	$req->content($to_retrieve);
	$to_retrieve='';
	$res = $ua->request($req);
	$with_head=0 if $tmp==1;
	($csv_data_ref,$count)=&crep_retrieve(\$res,$with_head);
	$tmp = $#psid + 1;
	die "Error: $count <> $tmp, Server Overflow !\n" if $count != $tmp;
	print OUTPUT $$csv_data_ref;
	print STDERR "\b"x10,">> 100 %   \n";
}
close OUTPUT;
#$req->content('unique_id=Os.54189.1.S1_at Os.42585.1.S1_at Os.12030.1.S1_at');
=pod
GET or POST is supported.(POST will be parsed over GET)
   DATA can be:
unique_id=Os.54189.1.S1_at Os.42585.1.S1_at Os.12030.1.S1_at
	OR
checkbox_Os.8271.1.S1_at=&checkbox_Os.8271.2.S1_x_at=&checkbox_Os.18205.1.S1_at=

since lc will be used before query, DATA is case-insensitive, still, 'ProbeSet' will be printed in the same case as those in DATA.
=cut
#$res = $ua->request($req);
sub crep_retrieve () {
	my ($res_ref,$with_head)=@_;
	my $tmp_html = $$res_ref->as_string;
#</table></center><br /><p /><br /><a href="/runblast/32494.csv">csv file</a><p /><br /></body></html>
	$tmp_html =~ m#</table></center><br /><p /><br /><a href="/runblast/(\d+).csv">csv file</a><p /><br /></body></html>$#s;
	my $csv_file=$1;
	if (! defined $csv_file) {
		print $tmp_html;
		die "Sever HTML wrong !";
	}
	#print "\nhttp://crep.ncpgr.cn/runblast/$csv_file.csv\n";
	my $csv = get "http://crep.ncpgr.cn/runblast/$csv_file.csv";
	#print $csv;
	my ($to_return,@row,$csvfh,$status)=('');
	#my $f = new File::Tabular($csvfh,{fieldSep=>','});
	#"ProbeSet","TIGR Locus","Minghui 63(endosperm, 14 days after pollination)","Shanyou 63(endosperm, 14 days after pollination)"
	#my $csv_subject = Text::CSV_XS->new;
	#my $status = $csv->parse ($line);
	my @csv=split /\n/,$csv;
	$csv='';
	shift @csv if $with_head == 0;
	for ( @csv ) {
	@row=();
		$status = $csv_subject->parse($_);
		for ($csv_subject->fields) {
				push @row,$_;
			}
	$to_return .= join "\t",@row;
	$to_return .= "\n";
	}
	my $item_count=$#csv + 1 - $with_head;
	@csv=();
	return (\$to_return,$item_count);
}

my $stop_time = [gettimeofday];
print STDERR "\n Time Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s)\n";

__END__
