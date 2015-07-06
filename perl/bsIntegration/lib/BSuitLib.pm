package main;
use BSuitInc;
use File::Path;
use File::Basename;
use Galaxy::IO;
use Galaxy::IO::FASTA;

sub do_pre() {
	#my $Config = $_[0];
	#ddx \$Config;
	File::Path::make_path("$RootPath/Ref",{verbose => 0,mode => 0755});
	my $Refprefix = getRef2char($HostRefName,$VirusRefName);
	my $Refile = "$RootPath/Ref/$RefFilesSHA/$Refprefix.fa";
#warn "[$HostRefName,$VirusRefName] -> $Refprefix [$RefFilesSHA]\n";
	my $found = 0;
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
}

sub do_aln() {
	my $Refilename = warnFileExist($RefConfig->{$RefFilesSHA}->{'Refilename'});
	#warn "$Refilename\n";
	my (%tID);
	for (@{$Config->{'DataFiles'}->{'='}}) {
		/([^.]+)\.(\d)/ or die;
		$tID{$1}{$2} = $_;
	}
	#ddx \%tID;
	File::Path::make_path("$RootPath/${ProjectID}_aln",{verbose => 0,mode => 0755});
	open O,'>',"$RootPath/${ProjectID}_aln.sh" or die $!;
	print O "#!/bin/sh\n\n";
	for my $k (keys %tID) {
		my @FQ1c = split /\s*,\s*/,$Config->{'DataFiles'}->{$tID{$k}{1}};
		my @FQ2c = split /\s*,\s*/,$Config->{'DataFiles'}->{$tID{$k}{2}};
		die "[x]  DataFiles not paired ! [@FQ1c],[@FQ2c]\n" unless $#FQ1c == $#FQ1c;
		my $cmd;
		if (@FQ1c == 1) {
			$cmd = <<"CMD";
$RealBin/bin/bwameth.py --reference $Refilename -t 24 --read-group $k -p $RootPath/${ProjectID}_aln/$k @{[warnFileExist($FQ1c[0],$FQ2c[0])]} 2>$RootPath/${ProjectID}_aln/$k.log
CMD
			print O $cmd;
		} else {
			my @theBams;
			for my $i (0 .. $#FQ1c) {
				my $fID = basename($FQ1c[$i]);
				$fID =~ s/\.fq(\.gz)?$//i;
				$cmd = <<"CMD";
$RealBin/bin/bwameth.py --reference $Refilename -t 24 --read-group '\@RG\\tID:${k}_${i}_${fID}\\tSM:$k' -p $RootPath/${ProjectID}_aln/${k}_${i}_${fID} @{[warnFileExist($FQ1c[$i],$FQ2c[$i])]} 2>$RootPath/${ProjectID}_aln/${k}_${i}_${fID}.log
CMD
				push @theBams,"$RootPath/${ProjectID}_aln/${k}_${i}_${fID}.bam";
				print O $cmd;
			}
			my $theBamsJ = join(' ',@theBams);
			$cmd = <<"CMD";
samtools merge -l 9 $RootPath/${ProjectID}_aln/$k.bam $theBamsJ
samtools index $RootPath/${ProjectID}_aln/$k.bam
CMD
			print O $cmd;
		}
		$cmd = <<"CMD";
samtools sort -m 2415919104 -n $RootPath/${ProjectID}_aln/$k.bam -O bam -T $RootPath/${ProjectID}_aln/$k.sn >$RootPath/${ProjectID}_aln/$k.sn.bam 2>>$RootPath/${ProjectID}_aln/$k.log

CMD
		print O $cmd;
	}
	close O;
	chmod 0755,"$RootPath/${ProjectID}_aln.sh";
	warn "[!] Please run [$RootPath/${ProjectID}_aln.sh] to do the aln.\n"
}

sub do_grep() {	
	my @RefChrIDs = split(',',$RefConfig->{$RefFilesSHA}->{'RefChrIDs'});
	my @VirusChrIDs = split(',',$RefConfig->{$RefFilesSHA}->{'VirusChrIDs'});
	
	my (%tID);
	for (@{$Config->{'DataFiles'}->{'='}}) {
		/([^.]+)\.(\d)/ or die;
		$tID{$1}{$2} = $_;
	}
	#   "780_T" => { 1 => "780_T.1", 2 => "780_T.2" },
	#   "s01_P" => { 1 => "s01_P.1", 2 => "s01_P.2" },
	File::Path::make_path("$RootPath/${ProjectID}_grep",{verbose => 0,mode => 0755});
	for my $k (keys %tID) {
		my $myBamf = "$RootPath/${ProjectID}_aln/$k.sn.bam";
		warn $myBamf;
	}
}

1;
