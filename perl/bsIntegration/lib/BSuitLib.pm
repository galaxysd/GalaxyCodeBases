package main;
use BSuitInc;
use File::Path;
use File::Basename;
use Galaxy::IO;
use Galaxy::IO::FASTA;

sub do_pre() {
	#my $Config = $_[0];
	#ddx \$Config;
	my $RootPath = $Config->{'Output'}->{'WorkDir'};
	$RootPath =~ s/[\/\\]+$//g;
	warn "[!] WorkDir: [$RootPath]\n";
	File::Path::make_path("$RootPath/Ref",{verbose => 0,mode => 0755});
	my $HostRefName = basename($Config->{'RefFiles'}->{'HostRef'});
	my $VirusRefName = basename($Config->{'RefFiles'}->{'VirusRef'});
	my $RefFilesSHA = getFilesHash($HostRefName,$VirusRefName);
	my $Refprefix = getRef2char($HostRefName,$VirusRefName);
	my $Refile = "$RootPath/Ref/$RefFilesSHA/$Refprefix.fa";
warn "[$HostRefName,$VirusRefName] -> $Refprefix [$RefFilesSHA]\n";
	my $found = 0;
	if ( -f "$RootPath/Ref/Ref.ini" ) {
		my $RefConfig = Galaxy::IO::INI->new();
		$Config->read("$RootPath/Ref/Ref.ini");
	}
	if ($found==0) {
		File::Path::make_path("$RootPath/Ref/$RefFilesSHA",{verbose => 0,mode => 0755});
		my $Ref=openfile($Config->{'RefFiles'}->{'HostRef'});
		my $RefHash = readwholefa($Ref);
	}
}


1;
