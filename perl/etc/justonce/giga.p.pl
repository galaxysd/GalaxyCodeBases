#!/usr/bin/env perl
use strict;
use warnings;
use DBI;
use Data::Dumper;

my $dbh = DBI->connect("dbi:SQLite:dbname=giga.authors.sqlite","","",{RaiseError => 0,PrintError => 1,AutoCommit => 0}) or die $DBI::errstr;

my $sthi = $dbh->prepare("SELECT * FROM PubDat;");
$sthi->execute();
my $ret = $sthi->fetchall_arrayref;

if ($#$ret == -1) {
	die "[x]Empty !\n";
}

for (@$ret) {
	print "DOI: $_->[0]\n";
	print "Title: $_->[1]\n";
	print "Type: $_->[2]\n";
	print "Authors: { $_->[3] }\n";
	print "RefList: { $_->[4] }\n";
}

$dbh->rollback;
$dbh->disconnect;
