#!/usr/bin/env perl
use strict;
use warnings;
use LWP;
use LWP::UserAgent;
use JSON qw( decode_json );
#use HTML::TreeBuilder 5 -weak;
use IO::Handle;
#use Data::Dumper;

my $ua = LWP::UserAgent->new;
$ua->agent("Mozilla/5.0 (Macintosh; Intel Mac OS X 10.12; rv:58.0) Gecko/20100101 Firefox/58.0");

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.xz$/) {
            open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
        open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.bz2$/) {
        open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}

sub getthings($) {
	my ($cnt) = @_;
	my ($tocSections,$jsonData)=('','');
	my @lines = split(/^/m,$cnt);
	for (my $i=0;$i<=$#lines;$i++) {
		if ($lines[$i] =~ /<script type="application\/ld\+json">/) {
			$jsonData = $lines[$i+1];
		}
		if ($lines[$i] =~ /Issue Section:/) {
			$lines[$i+1] =~ />([^<>]+)<\/a>/;
			$tocSections = $1;
			last;
		}
	}
	my $decoded_json = decode_json( $jsonData );
	#print ">>>$tocSections<<<\n";
	#print Dumper $decoded_json;
	return [$tocSections,$decoded_json];
}

sub fetchURL($) {
	my ($URL)=@_;
	my $req = HTTP::Request->new(GET => $URL);
	my $res = $ua->request($req);
	my $ret;
    if ($res->is_success) {
        #print $res->content;
		#my $tree = HTML::TreeBuilder->new_from_content($res->content);
		#$tree->dump; $tree->delete;
		$ret=getthings($res->content);
		return $ret;
    } else {
        print $res->status_line, "\n";
    }
}


my $fh = openfile('giga.tsv.bz2');
open O,'>','giga.authors.ini' || die("[x]Cannot Open Output File.");
binmode(O, ":utf8");
<$fh>;
while (<$fh>) {
	chomp;
	my @dat = split /\t/;
	print join(" | ",@dat[0,-1]),"\n";
	my $url = 'https://academic.oup.com/gigascience/article-lookup/doi/' . $dat[-1];
	my $ret=fetchURL($url);
	print O "[$dat[-1]]\nTitle=\"$dat[0]\"\nType=\"$ret->[0]\"\nAuthors={\n";
	#print Dumper $ret->[1];
	my $authors = ${$ret->[1]}{'author'};
	for (@$authors) {
		#print Dumper $_;
		print O join('"',"\t",$_->{'name'},"=",$_->{'affiliation'},"\n");
	}
	print O "}\n\n";
	O->flush();
}
close O;

__END__

eg:

[10.1186/s13742-015-0066-5]
Title="The ocean sampling day consortium"
Type="Commentary"
Authors={
	"Kopf, Anna"="1 Max Planck Institute for Marine Microbiology, Celsiusstrasse 1, D-28359Bremen, Germany   2 Jacobs University Bremen gGmbH, Campus Ring 1, D-28759 Bremen, Germany"
	"Bicak, Mesude"="3 University of Oxford, 7 Keble Road, OX1 3QG Oxford, Oxfordshire, UK"
	"Kottmann, Renzo"="1 Max Planck Institute for Marine Microbiology, Celsiusstrasse 1, D-28359Bremen, Germany"
	"Schnetzer, Julia"="1 Max Planck Institute for Marine Microbiology, Celsiusstrasse 1, D-28359Bremen, Germany   2 Jacobs University Bremen gGmbH, Campus Ring 1, D-28759 Bremen, Germany"
	"Øvreås, Lise"="26 Department of Biology, University of Bergen, Thormøhlensgate 53 B, 5020 Bergen, Norway"
	"Glöckner, Frank Oliver"="1 Max Planck Institute for Marine Microbiology, Celsiusstrasse 1, D-28359Bremen, Germany   2 Jacobs University Bremen gGmbH, Campus Ring 1, D-28759 Bremen, Germany"
}

