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
#use HTML::TableExtract;
use Data::Dump;

#open I,'<','resLOC.lst' or die $!;
#open O,'>','recLoc.anno' or die $!;

my $ua = LWP::UserAgent->new;
$ua->agent("Mozilla/5.0");

my ($req,$res,$count);

$req = HTTP::Request->new(POST => 'http://ricexpro.dna.affrc.go.jp/RXP_1006/gene-search.php');
$req->content_type('application/x-www-form-urlencoded');

my $keyword = 'LOC_Os01g01360';
#$keyword = 'LOC_Os01g01040';
$req->content("keyword=$keyword");
#$req = POST 'http://ricexpro.dna.affrc.go.jp/RXP_1006/gene-search.php', [ keyword => 'LOC_Os01g01030' ];
$res = $ua->request($req);
if ($res->is_success) {
	#print $res->decoded_content;
	my $tmp_html = $res->content;
	#print "$tmp_html";
=pod
      <div id="result">
        <table border="1" align="center" width="1000"><tr><th>Locus ID / Links</th><th>Locus<br>Select</th><th>FeatureNum<br>(Link to graph)</th><th>Feature<br>Select</th><th>Accession</th><th>Probe Sequence ID<br>(Link to SeqInfo)</th><th style="width:400px;">Description</th><th>MSU ID</th></tr>
<tr><td valign="top" rowspan="1"><span class="locus-link" tos17="1">Os01g0100500</span></td><td valign="top" rowspan="1"><input type="checkbox" class="locus-select" name="Os01g0100500"></td><td><a class="graph-link" barimg="images/barplot/RXP_1006-Os01g0100500-12943_bar.png" lineimg="images/lineplot/RXP_1006-Os01g0100500-12943_line.png" href="graph-view.php?featurenum=12943" target="_blank">12943</a></td>
<td><input type="checkbox" class="feature-select Os01g0100500-feature" name="12943"></td>
<td>AK067316</td><td><a href="probe-seq-info.php?seqid=S-10941" target="_blank">S-10941</a> (unique)</td><td><span class="desc descinfo">Similar to Pectinesterase-like protein.</span></td><td><a href="http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=LOC_Os01g01040" target="_blank">LOC_Os01g01040</a><br/><a href="http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=LOC_Os01g01030" target="_blank">LOC_Os01g01030</a><br/></td></tr>
</table>
      </div>
=cut
	my $tree= HTML::TreeBuilder::XPath->new;
	$tree->parse_content($tmp_html);
	my @toc = $tree->findnodes('//div[@id="result"]/table/tr/td');
	#ddx \@toc;
	#ddx $tree;
	my (@Locus,@FeatureNum,@Accession,@Desc,%Desc,@DescUniq);
	for my $el ( @toc ) {
		my $tmp = $el->as_HTML;
		#print $el->as_HTML," ---\n";
=pod
<td rowspan="2" valign="top"><span class="locus-link" tos17="0">Os01g0103100</span></td> ---
<td rowspan="2" valign="top"><input class="locus-select" name="Os01g0103100" type="checkbox" /></td> ---
<td><a barimg="images/barplot/RXP_1006-Os01g0103100-07015_bar.png" class="graph-link" href="graph-view.php?featurenum=7015" lineimg="images/lineplot/RXP_1006-Os01g0103100-07015_line.png" target="_blank">7015</a></td> ---
<td><input class="feature-select Os01g0103100-feature" name="7015" type="checkbox" /></td> ---
<td>AK070557</td> ---
<td><a href="probe-seq-info.php?seqid=S-5924" target="_blank">S-5924</a> (unique)</td> ---
<td><span class="desc descinfo">TGF-beta receptor, type I/II extracellular region family protein.</span></td> ---
<td><a href="http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=LOC_Os01g01360" target="_blank">LOC_Os01g01360</a><br /></td> ---
<td><a barimg="images/barplot/RXP_1006-Os01g0103100-36463_bar.png" class="graph-link" href="graph-view.php?featurenum=36463" lineimg="images/lineplot/RXP_1006-Os01g0103100-36463_line.png" target="_blank">36463</a></td> ---
<td><input class="feature-select Os01g0103100-feature" name="36463" type="checkbox" /></td> ---
<td>AK058723</td> ---
<td><a href="probe-seq-info.php?seqid=S-29658" target="_blank">S-29658</a> (unique)</td> ---
<td><span class="desc descinfo">TGF-beta receptor, type I/II extracellular region family protein.</span></td> ---
<td><a href="http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=LOC_Os01g01360" target="_blank">LOC_Os01g01360</a><br /></td> ---
=cut
		#print $el->as_trimmed_text,"\n";
		if ($tmp =~ /"locus-link".+\>(\w+)\<\//) {
			print "Locus=$1\n";
			push @Locus,$1;
		} elsif ($tmp =~ / barimg.+\>(\d+)\<\//) {
			print "FeatureNum=$1\n";
			push @FeatureNum,$1;
		} elsif ($tmp =~ /\<td\>(\w+)\<\/td\>/) {
			print "Accession=$1\n";
			push @Accession,$1;
		} elsif ($tmp =~ / class="desc descinfo".*\>([^<]+)\<\//) {
			print "Desc=$1\n";
			push @Desc,$1;
			++$Desc{$1};
		}
	}
	for (@Desc) {
		push @DescUniq,$_ if $Desc{$_} > 0;
		$Desc{$_} *= -1 if $Desc{$_} > 1;
	}
	my $tmp = join("\t",$keyword,join('|',@Locus),join('|',@FeatureNum),join('|',@Accession),join('|',@DescUniq));
	print $tmp,"\n";
	print '-' x 75,"\n";
	$tree->delete;
}
else {
	print "Error: " . $res->status_line . "\n";
}
my $tmp_html = $res->content;

#ddx $res;

__END__
awk '{print $1}' crep_all_tsv_new.txt.up2*.txt |sort|uniq > resLOC.lst
