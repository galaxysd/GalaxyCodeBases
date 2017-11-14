#!/usr/bin/env perl
use strict;
use warnings;
use LWP;
use LWP::UserAgent;
use JSON qw( decode_json );
use IO::Handle;
use Data::Dumper;
use DBI;

my $dbh = DBI->connect("dbi:SQLite:dbname=giga.authors.sqlite","","",{RaiseError => 0,PrintError => 1,AutoCommit => 0}) or die $DBI::errstr;
$dbh->do("CREATE TABLE PubDat (DOI TEXT, Title TEXT, Type TEXT, Authors TEXT, RefList TEXT)") or die $dbh->errstr;
$dbh->commit;
my $sthi = $dbh->prepare( "INSERT INTO PubDat ( DOI,Title,Type,Authors,RefList ) VALUES ( ?,?,?,?,? )" );

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
	my ($tocSections,$reflist,$jsonData)=('','');
	my @lines = split(/^/m,$cnt);
	for (my $i=0;$i<=$#lines;$i++) {
		if ($lines[$i] =~ /<script type="application\/ld\+json">/) {
			$jsonData = $lines[$i+1];
		}
		if ($lines[$i] =~ /<div class="ref-list">/) {
			$reflist = $lines[$i];
			chomp($reflist);
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
	return [$tocSections,$decoded_json,$reflist];
}

sub fetchURL($$) {
	my ($URL,$times)=@_;
	my $req = HTTP::Request->new(GET => $URL);
	my $ret;
	for my $i (1 .. $times) {
		my $res = $ua->request($req);
	    if ($res->is_success) {
			$ret = getthings($res->content);
			return $ret;
	    } else {
	        print $res->status_line, " <<<--- $i of $times\n";
			$ret = ["\n",$res->status_line];
			sleep $i;
	    }
	}
	return $ret;
}


my $fh = openfile('giga.tsv.bz2');
open O,'>','giga.authors.ini' || die("[x]Cannot Open Output File.");
binmode(O, ":utf8");
<$fh>;
while (<$fh>) {
	chomp;
	my @dat = split /\t/;
	print join(" | ",@dat[0,16]),"\n";
	unless ($dat[16] =~ /\//) {
		print O "[$dat[16]]\nTitle=\"$dat[0]\"\nType=\"Missing or Misformatted DOI !\"\n\n";
		next;
	}
	my $url = 'https://academic.oup.com/gigascience/article-lookup/doi/' . $dat[16];
	my $ret=fetchURL($url,5);
	#my $ret=[''];
	if ($ret->[0] eq '') {
		print O "[$dat[16]]\nTitle=\"$dat[0]\"\nType=\"Wrong DOI !\"\n\n";
		next;
	} elsif ($ret->[0] eq "\n") {
		print O "[$dat[16]]\nTitle=\"$dat[0]\"\nType=\"Error: $ret->[1]\"\n\n";
		next;
	}
	print O "[$dat[16]]\nTitle=\"$dat[0]\"\nType=\"$ret->[0]\"\nAuthors={\n";
	#print Dumper $ret;
	my $authors = ${$ret->[1]}{'author'};
	my $AuthorStr = '';
	for (@$authors) {
		#print Dumper $_;
		print O join('"',"\t",$_->{'name'},"=",$_->{'affiliation'},"\n");
		$AuthorStr .= join('"','',$_->{'name'},"=",$_->{'affiliation'},"\n");
	}
	print O "}\nRefList=\{$ret->[2]\}\n\n";
	O->flush();
	$sthi->execute($dat[16],$dat[0],$ret->[0],$AuthorStr,$ret->[2]) or die $sthi->errstr;
	$dbh->commit;
}
close O;
$dbh->disconnect;

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

