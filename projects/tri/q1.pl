#!/usr/bin/env perl
use strict;
use warnings;
use Bio::AlignIO;
use DBI;

use Data::Dump qw (ddx);
#use IO::Unread;

use Bio::EnsEMBL::Registry;

## Load the registry automatically
my $reg = "Bio::EnsEMBL::Registry";
#my $url = 'mysql://anonymous@ensembldb.ensembl.org';
my $url = 'mysql://root@localhost';
$reg->load_registry_from_url($url);

my $compara_db_adaptor = $reg->get_DBAdaptor('Multi', 'compara');
#my $genome_db_adaptor = $compara_db_adaptor->get_GenomeDBAdaptor();
my $dbh = $compara_db_adaptor->dbc->db_handle;

my $sql = qq{
SELECT DISTINCT ss.genome_db_id, gdb.taxon_id, gdb.name, gdb.assembly FROM species_set_tag sst
 JOIN species_set ss USING (species_set_id)
 JOIN genome_db gdb USING (genome_db_id)
 WHERE sst.value = 'mammals';
};

my $sth = $dbh->prepare($sql);
$sth->execute();

my (@IDsGnomeDB, %GnomeDBnfo);
while(my @retdat = $sth->fetchrow()) {
	#warn join("\t",@retdat),"\n";
	#$GnomeDBnfo{$retdat[0]} = [@retdat[1,2,3]];
	my $GnomedbID = shift @retdat;
	push @IDsGnomeDB,$GnomedbID;
	$GnomeDBnfo{$GnomedbID} = \@retdat;
}

ddx \%GnomeDBnfo;
