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
my $genome_db_adaptor = $compara_db_adaptor->get_GenomeDBAdaptor();
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

#my $MammaliaTaxID = 40674; # http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=40674&lvl=3&keep=1&srchmode=1&unlock
#my $MammaliaGDBs = $genome_db_adaptor->fetch_all_by_ancestral_taxon_id($MammaliaTaxID,1);
#ddx $MammaliaGDBs;

__END__

SELECT DISTINCT ss.species_set_id, ss1.genome_db_id, gdb1.taxon_id, gdb1.name, ss.genome_db_id, gdb.taxon_id, gdb.name, gdb.assembly FROM species_set ss
  JOIN species_set ss0 ON ss.genome_db_id=ss0.genome_db_id
  JOIN species_set ss1 ON ss.species_set_id=ss1.species_set_id
  JOIN species_set ss2 ON ss1.genome_db_id=ss2.genome_db_id
  JOIN species_set_tag sst0 ON sst0.species_set_id=ss0.species_set_id
  JOIN species_set_tag sst2 ON sst2.species_set_id=ss2.species_set_id
  JOIN genome_db gdb1 ON gdb1.genome_db_id=ss1.genome_db_id
  JOIN genome_db gdb ON gdb.genome_db_id=ss.genome_db_id
  WHERE sst0.value = 'mammals' AND sst2.value = 'mammals' AND ss1.genome_db_id != ss.genome_db_id;
