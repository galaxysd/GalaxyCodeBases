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
use HTML::TableExtract;
use Data::Dump;

$main::VERSION=1.3.2;

if (@ARGV == 0) {
	&VERSION_MESSAGE();
	&HELP_MESSAGE();
	die "\n";
}
our($opt_p, $opt_o, $opt_n, $opt_s);

$Getopt::Std::STANDARD_HELP_VERSION=1;
getopts('p:o:s:');
#print "$opt_i, $opt_o, $opt_e\n";

sub HELP_MESSAGE() {
	#my ($scr) = ($0 =~ m,([^/\\]+)$,);
	my $help=<<EOH;
\t-p NetAffx Rice Probe Set Information (.psi) file [./Rice.psi]
\t-o Output file (crep_anno_tsv.txt)
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
$opt_o='crep_anno_tsv.txt' if ! defined $opt_o;
$opt_s=1 if ! defined $opt_s;
$opt_s=int($opt_s);
$opt_s=1 if $opt_s < 1;
&VERSION_MESSAGE();
warn "\nInput PSI: \033[32;1m$opt_p\033[0;0m\nOutput: \033[32;1m$opt_o\033[0;0m\n\n";
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
$ua->cookie_jar(HTTP::Cookies->new(file => "lwpcookies2.txt", ignore_discard => 1,
								   autosave => 1));

$req = HTTP::Request->new(GET => 'http://crep.ncpgr.cn/crep-cgi/login2.pl?User_Name=guest&amp;Password=57510426775c5b0f');
# send request
$res = $ua->request($req);
# check the outcome
if (! $res->is_success) {
		die "Error: Login as Guest failed !\n";
}

################### Retrieve ####################
# http://crep.ncpgr.cn/crep-cgi/elements_details.pl?unique_id=Os.24470.1.A1_at
$req = HTTP::Request->new(POST => 'http://crep.ncpgr.cn/crep-cgi/elements_details.pl');
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
while ($#psid >= 0) {
	$to_retrieve = shift(@psid);
	$req = HTTP::Request->new(GET => 'http://crep.ncpgr.cn/crep-cgi/elements_details.pl?unique_id='.$to_retrieve);
	#$req->content_type('application/x-www-form-urlencoded');

	$res = $ua->request($req);
	($csv_data_ref,$count)=&crep_retrieve(\$res,$with_head,$to_retrieve);
	die "Error: $count <> 1, Server Overflow !\n" if $count != 1;
	$with_head=0;
	print OUTPUT $$csv_data_ref;
#print "$$csv_data_ref";
	$percent=sprintf "%5.2f",($ps_count-$#psid-1)*100/$ps_count;
	print STDERR "\b"x10, ($#psid % 2) ? ('>='):('=>')," $percent %";
	$tmp=1;
}
if ($#psid >= 0) {
	die;
}
close OUTPUT;

sub crep_retrieve () {
	my ($res_ref,$with_head,$to_retrieve)=@_;
	my $tmp_html = $$res_ref->content;
=pod
<tr bgcolor="#FFFFFF"><th align="left" valign="center" width="20%">Associated TIGR Gene Locus <a class="descrp" onclick="help_window('/crep-cgi/help.pl?Element_details','Help','450','400')" href="javascript:void(0);"><img height="23" border="0" src="/images/question_s.gif" width="8" /></a></th><td align="left"><table border="0" cellspacing="1" cellpadding="5" width="95%"><tr><th align="center" valign="center" width="20%">TIGR Locus</th><th align="center" valign="center" width="20%">Matching Probes</th><th align="center" valign="center" width="60%">Description</th></tr><tr><td align="center" valign="center">NONE</td><td align="left" valign="center" width="20%">NONE</td><td align="left" valign="center" width="60%">NONE</td></tr></table></td></tr>

<tr bgcolor="#FFFFFF"><th align="left" valign="center" width="20%">Associated TIGR Gene Locus <a class="descrp" onclick="help_window('/crep-cgi/help.pl?Element_details','Help','450','400')" href="javascript:void(0);"><img height="23" border="0" src="/images/question_s.gif" width="8" /></a></th><td align="left"><table border="0" cellspacing="1" cellpadding="5" width="95%"><tr><th align="center" valign="center" width="20%">TIGR Locus</th><th align="center" valign="center" width="20%">Matching Probes</th><th align="center" valign="center" width="60%">Description</th></tr><tr><td align="center" valign="center" width="20%"><a href='http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?db=osa1r5&orf=LOC_Os12g17920' target='_blank'>LOC_Os12g17920</a></td><td align="center" valign="center" width="20%">1/11</td><td align="center" valign="center" width="60%">hypothetical protein</td></tr><tr><td align="center" valign="center" width="20%"><a href='http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?db=osa1r5&orf=LOC_Os11g17770' target='_blank'>LOC_Os11g17770</a></td><td align="center" valign="center" width="20%">1/11</td><td align="center" valign="center" width="60%">transposon protein, putative, Mutator sub-class</td></tr></table></td></tr>
=cut
	$tmp_html = (grep /Associated TIGR Gene Locus/,(split /\n/,$tmp_html))[0];
	#$tmp_html =~ m#^<tr bgcolor=\"\#FFFFFF\"><th[^>]*>Associated TIGR Gene Locus.+(<table.*?</table>)#s;
#print "\n\n[$tmp_html]\n";
	my $te = HTML::TableExtract->new( headers => [('TIGR Locus', 'Matching Probes', 'Description')] );
	$te->parse("$tmp_html");
	my $ts = $te->first_table_found();
	my ($to_return,$item_count,@dat) = ('',1);
	unless ($ts) {
		warn "\n[!Table]$tmp_html\n";
		$item_count = 0;
	}
#ddx $ts;
	foreach my $row ($ts->rows) {
		push @dat,join("|",@$row);
	}
	$to_return = join("\t",$to_retrieve,@dat)."\n";
#print "$to_return\n";
	return (\$to_return,$item_count);
}

my $stop_time = [gettimeofday];
print STDERR "\n Time Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s)\n";

__END__
