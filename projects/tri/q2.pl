#!/usr/bin/env perl
use strict;
use warnings;
use Bio::AlignIO;
use DBI;

use Data::Dump qw (ddx);
#use IO::Unread;

use Bio::Align::DNAStatistics;
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

$sql = qq{
	SELECT DISTINCT method_link_species_set_id FROM method_link_species_set
	  JOIN method_link USING (method_link_id)
	  JOIN species_set ss ON ss.species_set_id=method_link_species_set.species_set_id
	  JOIN species_set ss0 ON ss.genome_db_id=ss0.genome_db_id
	  JOIN species_set ss1 ON ss.species_set_id=ss1.species_set_id
	  JOIN species_set ss2 ON ss1.genome_db_id=ss2.genome_db_id
	  JOIN species_set_tag sst0 ON sst0.species_set_id=ss0.species_set_id
	  JOIN species_set_tag sst2 ON sst2.species_set_id=ss2.species_set_id
	  WHERE sst0.value = 'mammals' AND sst2.value = 'mammals' AND ss1.genome_db_id != ss.genome_db_id AND method_link.type='ENSEMBL_ORTHOLOGUES';
};
$sth = $dbh->prepare($sql);
$sth->execute();

my (@IDsMLSS);
while(my @retdat = $sth->fetchrow()) {
	push @IDsMLSS,$retdat[0];
}
warn 'MLSS_ID: ',join(', ',@IDsMLSS),"\n";

## The BioPerl alignment formatter
my $alignIO = Bio::AlignIO->newFh(-format => "clustalw");

my $HomologyAdaptor = $reg->get_adaptor("Multi", "compara", "Homology");
for my $mlss (@IDsMLSS) {
	my $homologies = $HomologyAdaptor->fetch_all_by_MethodLinkSpeciesSet($mlss);
	print '-----',$mlss,"-----\n";
    ## For each homology
    foreach my $this_homology (@{$homologies}) {
		next unless defined $this_homology->dn();
  	  my $description = $this_homology->description;
  	  # next unless ($description =~ /orth/);    # uncomment for orthologs only
        my ($a,$b) = @{$this_homology->gene_list};
        my $spa = $a->taxon->get_short_name;
        my $spb = $b->taxon->get_short_name;
        my $labela = $a->stable_id;
        $labela .= "(" . $a->display_label . ")" if $a->display_label;
        my $labelb = $b->stable_id;
        $labelb .= "(" . $b->display_label . ")" if $b->display_label;
		 print "$spa,$labela,$spb,$labelb,$description\n";
		#ddx $this_homology;
    }
    print "\n";
}

__END__
