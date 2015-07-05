package main;
use BSuitInc;
use File::Path;
use File::Basename;
use Galaxy::IO;
use Galaxy::IO::FASTA;

sub do_pre() {
	#my $Config = $_[0];
	ddx \$Config;
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
		$found=1 if exists $RefConfig->{$RefFilesSHA};
	}
	if ($found==0) {
		File::Path::make_path("$RootPath/Ref/$RefFilesSHA",{verbose => 0,mode => 0755});
		$RefConfig->{$RefFilesSHA} = { Refilename => $Refile, RefChrIDs => '', VirusChrIDs => '' };
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
		$RefConfig->write("$RootPath/Ref/Ref.ini");
	} else {
		warn "[!] Already Read References Pairs:[$HostRefName,$VirusRefName].\n";
	}
	#ddx \$RefConfig;
	warn "[!] Building index for [$Refile].\n";
	system("$RealBin/bin/bwameth.py",'index',$Refile);
	#warn "[!] Prepare done !\n";
}

sub do_aln() {
	my $RootPath = $Config->{'Output'}->{'WorkDir'};
	my $ProjectID = $Config->{'Output'}->{'ProjectID'};
	warn "[!] Working on: [$ProjectID] @ [$RootPath]\n";
	my $RefConfig = Galaxy::IO::INI->new();
	if ( -f "$RootPath/Ref/Ref.ini" ) {
		$RefConfig->read("$RootPath/Ref/Ref.ini");
	} else {die "[x] Prepare INI not found ! [$RootPath/Ref/Ref.ini]\n";}
	my $HostRefName = basename($Config->{'RefFiles'}->{'HostRef'});
	my $VirusRefName = basename($Config->{'RefFiles'}->{'VirusRef'});
	my $RefFilesSHA = getFilesHash($HostRefName,$VirusRefName);
	my $Refilename = $RefConfig->{$RefFilesSHA}->{'Refilename'};
	#warn "$Refilename\n";
	my (%tID);
	for (@{$Config->{'DataFiles'}->{'='}}) {
		/([^.]+)\.(\d)/ or die;
		$tID{$1}{$2} = $_;
	}
	#ddx \%tID;
	File::Path::make_path("$RootPath/${ProjectID}_aln",{verbose => 0,mode => 0755});
	open O,'>',"$RootPath/${ProjectID}_aln.sh" or die $!;
	print O "!#/bin/sh\n\n";
	for (keys %tID) {
		my $cmd = "$RealBin/bin/bwameth.py --reference $Refilename -t 24 -p ${_}_r --read-group $_ $Config->{'DataFiles'}->{$tID{$_}{1}} $Config->{'DataFiles'}->{$tID{$_}{2}}\n";
		print O $cmd;
	}
	close O;
	chmod 0755,"$RootPath/${ProjectID}_aln.sh";
	warn "[!] Please run [$RootPath/${ProjectID}_aln.sh] to do the aln.\n"
}

sub do_grep() {}

1;
