package main;
use BSuitInc;
use File::Path;
use File::Basename;

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
warn "[$HostRefName,$VirusRefName] -> $Refprefix [$RefFilesSHA]\n";
	if ( -f "$RootPath/Ref/Ref.ini" ) {
		my $RefConfig = Galaxy::IO::INI->new();
		$Config->read("$RootPath/Ref/Ref.ini");
	}
}


1;
