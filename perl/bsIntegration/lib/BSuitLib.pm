package main;
use File::Path;

sub do_pre() {
	#my $Config = $_[0];
	ddx \$Config;
	my $RootPath = $Config->{'Output'}->{'WorkDir'};
	$RootPath =~ s/[\/\\]+$//g;
	warn "[!] WorkDir: [$RootPath]\n";
	File::Path::make_path("$RootPath/Ref",{verbose => 0,mode => 0755});
	if ( -f "$RootPath/Ref/Ref.ini" ) {
		my $RefConfig = Galaxy::IO::INI->new();
		$Config->read("$RootPath/Ref/Ref.ini");
	}
}


1;
