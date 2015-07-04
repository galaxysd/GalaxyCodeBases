package main;
use File::Path;
use File::Basename;
use Digest::SHA; # qw(sha1_base64);

sub do_pre() {
	#my $Config = $_[0];
	ddx \$Config;
	my $RootPath = $Config->{'Output'}->{'WorkDir'};
	$RootPath =~ s/[\/\\]+$//g;
	warn "[!] WorkDir: [$RootPath]\n";
	File::Path::make_path("$RootPath/Ref",{verbose => 0,mode => 0755});
	my $HostRefName = basename($Config->{'RefFiles'}->{'HostRef'});
	my $VirusRefName = basename($Config->{'RefFiles'}->{'VirusRef'});
	my $RefFilesSHA = Digest::SHA::sha1_base64("$HostRefName,$VirusRefName");
warn "[$HostRefName,$VirusRefName] -> [$RefFilesSHA]\n";
	if ( -f "$RootPath/Ref/Ref.ini" ) {
		my $RefConfig = Galaxy::IO::INI->new();
		$Config->read("$RootPath/Ref/Ref.ini");
	}
}


1;
