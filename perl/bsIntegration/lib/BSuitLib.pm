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
#warn "[$HostRefName,$VirusRefName] -> $Refprefix [$RefFilesSHA]\n";
	my $found = 0;
	my $RefConfig = Galaxy::IO::INI->new();
	if ( -f "$RootPath/Ref/Ref.ini" ) {
		$RefConfig->read("$RootPath/Ref/Ref.ini");
	}
	if ($found==0) {
		File::Path::make_path("$RootPath/Ref/$RefFilesSHA",{verbose => 0,mode => 0755});
		#$Config->{$RefFilesSHA};
		open O,'>',$Refile or die $!;
		my $FH=openfile($Config->{'RefFiles'}->{'HostRef'});
		my @ChrIDs;
		while (my $ret = FastaReadNext($FH)) {
			if ($$ret[0] =~ /(^chrEBV$)|(Un[_-])|(random$)/) {
				warn "[!]  skip Ref[$$ret[0]].\n";
				next;
			}
			my $len = length $$ret[1];
			warn "[!]  read Ref[$$ret[0]], $len bp.\n";
			print O ">$$ret[0]\n$$ret[1]\n";
			$RefConfig->{$RefFilesSHA}->{$$ret[0]} = $len;
			push @ChrIDs,$$ret[0];
		}
		$RefConfig->{$RefFilesSHA}->{RefChrIDs} = join(',',@ChrIDs);
		close $FH;
		@ChrIDs=();
		$FH=openfile($Config->{'RefFiles'}->{'VirusRef'});
		while (my $ret = FastaReadNext($FH)) {
			my $len = length $$ret[1];
			warn "[!]  read Virus[$$ret[0]], $len bp.\n";
			print O ">$$ret[0]\n$$ret[1]\n";
			$RefConfig->{$RefFilesSHA}->{$$ret[0]} = $len;
			push @ChrIDs,$$ret[0];
		}
		$RefConfig->{$RefFilesSHA}->{VirusChrIDs} = join(',',@ChrIDs);
		close $FH;
	}
	ddx $RefConfig;
}


1;
