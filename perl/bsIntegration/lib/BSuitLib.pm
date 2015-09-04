package main;
#use strict;
use warnings;
use BSuitInc;
use File::Path 2.08;	# http://search.cpan.org/~riche/File-Path-2.11/lib/File/Path.pm#API_CHANGES
use File::Basename;
use Galaxy::IO;
use Galaxy::IO::FASTA;
#use JSON;

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

sub do_grep($) {
	my $cfgfile = $_[0];
	my (%tID,%tFH);
	for (@{$Config->{'DataFiles'}->{'='}}) {
		/([^.]+)\.(\d)/ or die;
		$tID{$1}{$2} = $_;
	}
	#   "780_T" => { 1 => "780_T.1", 2 => "780_T.2" },
	#   "s01_P" => { 1 => "s01_P.1", 2 => "s01_P.2" },
	File::Path::make_path("$RootPath/${ProjectID}_grep",{verbose => 0,mode => 0755});
	my $GrepResult = Galaxy::IO::INI->new();
	my %ReadsIndex;
	for my $k (keys %tID) {
		my $myBamf = "$RootPath/${ProjectID}_aln/$k.bam";
		my $InsMean = $Config->{'InsertSizes'}->{$k} or die;
		my $InsSD = $Config->{'InsertSizes'}->{"$k.SD"} or die; # SD cannot be 0, so no need to test with defined.
		#warn "$myBamf,$InsMean,$InsSD";
		$GrepResult->{$k} = {
			'InBam' => $myBamf,
			InsMean => $InsMean,
			InsSD => $InsSD,
		};
		$GrepResult->{$k} = { DatFile => "$RootPath/${ProjectID}_grep/$k.sam" };
		open( IN,"-|","samtools view $myBamf") or die "Error opening $myBamf: $!\n";	# `-F768` later
		system( "samtools view -H $myBamf >".$GrepResult->{$k}{'DatFile'} );
		open GOUT,'>>',$GrepResult->{$k}{'DatFile'} or die "$!";
		print GOUT join("\t",'@PG','ID:bsuit',"CL:\"grep $cfgfile\""),"\n";
		print STDERR "[!] Reading [$myBamf] ...";
		while (my $line = <IN>) {
			my @Dat1 = split /\t/,$line;
			my $flag = 0;
			if ($Dat1[6] eq '=') {
				$flag |= 1 if abs(abs($Dat1[8])-$InsMean) > 3*$InsSD;
				$flag |= 2 if exists($VirusChrIDs{$Dat1[2]}) or exists($VirusChrIDs{$Dat1[6]});
			} else {
				$flag |= 4 if exists($VirusChrIDs{$Dat1[2]}) or exists($VirusChrIDs{$Dat1[6]});
			}
			next unless $flag;
			$flag |= 8 if $Dat1[5] !~ /^\d+M$/;
			my $id = join("\t",$k,$Dat1[0]);
			$ReadsIndex{$id} = [[0,$Dat1[2],$flag]];
		}
		close IN;
		open IN,'-|',"samtools view $myBamf" or die "Error opening $myBamf: $!\n";
		while (my $line = <IN>) {
			#my ($id, $flag, $ref, $pos, $mapq, $CIGAR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/,$line;
			#print "$id, $flag, $ref, $pos, $mapq, $CIGAR, $mref, $mpos, $isize\n";
			my @Dat1 = split /\t/,$line;
			my $r12R1;
			if ($Dat1[1] & 0x40) {
				$r12R1 = 1;
			} elsif ($Dat1[1] & 0x80) {
				$r12R1 = 2;
			} else {die $Dat1[1];}
			#--$r12R1;
			#my $line2 = <IN>;
			#die '[x]SAM/BAM file not paired !' unless defined($line2);
			#my @Dat2 = split /\t/,$line2;
			#next if $Dat1[4]<$CFGminMAPQ;
			my $flag = 0;
			if ($Dat1[6] eq '=') {
				$flag |= 1 if abs(abs($Dat1[8])-$InsMean) > 3*$InsSD;
				$flag |= 2 if exists($VirusChrIDs{$Dat1[2]}) or exists($VirusChrIDs{$Dat1[6]});
			} else {
				$flag |= 4 if exists($VirusChrIDs{$Dat1[2]}) or exists($VirusChrIDs{$Dat1[6]});
			}
			$flag |= 8 if $Dat1[5] !~ /^\d+M$/;
			my $curpos = tell(GOUT);
			my $id = join("\t",$k,$Dat1[0]);
			#warn "$flag $InsSD\t$curpos\n[${line}]\n" if $flag>1;
			next unless defined $ReadsIndex{$id};
			push @{$ReadsIndex{$id}},[$curpos,$r12R1,@Dat1[2,3,4,5],$flag];	# Chr, Pos, MapQ, CIGAR
			#print "($Dat1[2],$Dat1[3],$Dat1[5],$curpos)\n";
			print GOUT $line;
			#ddx $ReadsIndex{$id} if scalar @{$ReadsIndex{$id}} > 2;
			#last;
		}
		close IN;
		close GOUT;
		print STDERR "\b\b\bdone.\n";
		open $tFH{$k},'<',$GrepResult->{$k}{'DatFile'} or die "$!";
	}
	for my $tk (keys %ReadsIndex) {
		my ($flag,%tmp)=(0);
		for my $i (1 .. $#{$ReadsIndex{$tk}}) {
			my $x = join("\t",$ReadsIndex{$tk}->[$i][2],$ReadsIndex{$tk}->[$i][3]);
			$tmp{$x} = $i;
			$flag |= $ReadsIndex{$tk}->[$i][6];
		}
		my @s = sort sortChrPos(keys %tmp);
		$ReadsIndex{$tk}->[0][0] = $tmp{$s[0]};
		$ReadsIndex{$tk}->[0][1] = join("\t",$ReadsIndex{$tk}->[$tmp{$s[0]}][2],$ReadsIndex{$tk}->[$tmp{$s[0]}][3]);
		$ReadsIndex{$tk}->[0][2] = $flag;
#ddx $ReadsIndex{$tk};
# BSuitLib.pm:213: [
#   [1, "chr15\t20858350", 15],
#   [1477321, 1, "chr15", 20858350, 0, "48M42S", 12],
#   [33704587, 1, "gi|59585|emb|X04615.1|", 652, 0, "45S45M", 11],
#   [34264944, 2, "gi|59585|emb|X04615.1|", 717, 60, "56M34S", 12],
# ]
	}
	warn "[!] Bam Reading done.\n";
	#ddx \$GrepResult;
	#ddx (\%RefChrIDs,\%VirusChrIDs);
	open BOUT,'>',"$RootPath/${ProjectID}_grep/blocks.ini" or die "$!";
	my @IDsorted = sort { sortChrPos($ReadsIndex{$a}->[0][1],$ReadsIndex{$b}->[0][1]) } keys %ReadsIndex;
	my ($Cnt,$hChr,$hsPos,$hePos,%hChrRange,%vChrRange,@Store,@PureVirReads)=(0,"\t",0,0);
	for my $cid (@IDsorted) {
		my @minCPR = split /\t/,$ReadsIndex{$cid}->[0][1];
		push @minCPR,$ReadsIndex{$cid}->[0][0];
		if (exists $VirusChrIDs{$minCPR[0]}) {
			push @PureVirReads,$cid;
			next;
		}
		my ($thisehPos,$thisvChr,$thissvPos,$thisevPos);
		for my $i (1 .. $#{$ReadsIndex{$cid}}) {
			my $ret = mergeIn(\%hChrRange,$ReadsIndex{$cid}->[$i]);
		}
		my $rRead12 = $minCPR[2] ^ 3;
		if ($ReadsIndex{$cid}->[11] & 6) {
			my ($reflen,$readlen) = cigar2poses($ReadsIndex{$cid}->[$minCPR[2]+6]);
			$thisehPos = $minCPR[1] + $reflen;
			my ($Vreflen,$Vreadlen) = cigar2poses($ReadsIndex{$cid}->[$rRead12+6]);
			$thisvChr = $ReadsIndex{$cid}->[$rRead12+2];
			$thisvsPos = $ReadsIndex{$cid}->[$rRead12+4];
			$thisvePos = $thissvPos + $Vreflen;
		} else {
			my ($reflen,$readlen) = cigar2poses($ReadsIndex{$cid}->[$rRead12+6]);
			$thisehPos = $ReadsIndex{$cid}->[$rRead12+4] + $reflen;
		}
		if ($hChr eq "\t") {
			die "[x] The 1st sorted SAM record is Pure Virus PE.\n" if exists $VirusChrIDs{$minCPR[0]};
			$hChr = $minCPR[0];
			$hsPos = $minCPR[1];
			$hePos = $thisehPos;
			push @Store,$cid;
			next;
		}
		if ($hChr eq $minCPR[0] and $minCPR[1] <= $ehPos) {
			$hePos = $thisehPos;
			push @Store,$cid;
		} else {
			++$Cnt;
			print BOUT "[Block$Cnt]\nHostChrPoses=$hChr:$hsPos-$hePos\nVirusChrPoses=???\nSamFS=",
				join(',',map { my @t = split /\t/,$ReadsIndex{$_}->[0];$ReadsIndex{$_}->[$t[2]]; } @Store),"\n\n";
			@Store = ($cid);
			$hChr = $minCPR[0];
			$hsPos = $minCPR[1];
			$hePos = $thisehPos;
			%vChrRange = ();
		}
		my ($reflen,$readlen) = cigar2poses($ReadsIndex{$cid}->[$minCPR[2]+6]);
		my $start = $ReadsIndex{$cid}->[$minCPR[2]+6]
	}
	close BOUT;
	close $tFH{$_} for keys %tFH;
	warn $GrepResult->write_string;
	ddx $GrepResult;
	$GrepResult->write("$RootPath/${ProjectID}_grep.ini");
}

1;

__END__
samtools view -h /share/users/huxs/work/bsvir/bsI/SZ0010_aln/780_T.bam '*' | samtools bam2fq -O - | gzip -9 > /share/users/huxs/work/bsvir/bsI/SZ0010_aln/780_T.unmap.fq.gz
samtools view -h /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.bam '*'|samtools bam2fq -O -|gzip -9 > /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap.fq.gz &

./bin/bwameth.py --reference /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa -t 24 --read-group 780_T -p ~/work/bsvir/bsI/SZ0010_aln/780_T.unmap ~/work/bsvir/bsI/SZ0010_aln/780_T.unmap.fq.gz 2>~/work/bsvir/bsI/SZ0010_aln/780_T.unmap.log &    #/
./bin/bwameth.py --reference /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa -t 24 --read-group s01_P -p /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap.fq.gz 2>/share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap.log &

samtools view -h /share/users/huxs/work/bsvir/bsI/SZ0010_aln/780_T.bam 'gi|86261677|emb|AJ507799.2|'
