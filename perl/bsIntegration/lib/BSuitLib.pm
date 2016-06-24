package main;
use strict;
use warnings;
require BSuitInc;
use File::Path 2.08;	# http://search.cpan.org/~riche/File-Path-2.11/lib/File/Path.pm#API_CHANGES
use File::Basename;
use Galaxy::IO;
use Galaxy::IO::FASTA;
#use JSON;

sub do_pre() {
	#my $Config = $_[0];
	#ddx \$Config;
	File::Path::make_path("$main::RootPath/Ref",{verbose => 0,mode => 0755});
	my $Refprefix = getRef2char($main::HostRefName,$main::VirusRefName);
	my $Refile = "$main::RootPath/Ref/$main::RefFilesSHA/$Refprefix.fa";
#warn "[$main::HostRefName,$main::VirusRefName] -> $Refprefix [$main::RefFilesSHA]\n";
	my $found = 0;
	if ( -f "$main::RootPath/Ref/Ref.ini" ) {
		$main::RefConfig->read("$main::RootPath/Ref/Ref.ini");
		$found=1 if exists $main::RefConfig->{$main::RefFilesSHA};
	}
	if ($found==0) {
		File::Path::make_path("$main::RootPath/Ref/$main::RefFilesSHA",{verbose => 0,mode => 0755});
		$main::RefConfig->{$main::RefFilesSHA} = { Refilename => $Refile, RefChrIDs => '', VirusChrIDs => '' };
		open O,'>',$Refile or die $!;
		my $FH=openfile($main::Config->{'RefFiles'}->{'HostRef'});
		my @ChrIDs;
		while (my $ret = FastaReadNext($FH)) {
			if ($$ret[0] =~ /(^chrEBV$)|(Un[_-])|(random$)/) {
				warn "[!]  skip Ref[$$ret[0]].\n";
				next;
			}
			my $len = length $$ret[1];
			warn "[!]  read Ref[$$ret[0]], $len bp.\n";
			print O ">$$ret[0]\n$$ret[1]\n";
			$main::RefConfig->{$main::RefFilesSHA}->{$$ret[0]} = $len;
			push @ChrIDs,$$ret[0];
		}
		$main::RefConfig->{$main::RefFilesSHA}->{RefChrIDs} = join(',',@ChrIDs);
		close $FH;
		@ChrIDs=();
		$FH=openfile($main::Config->{'RefFiles'}->{'VirusRef'});
		while (my $ret = FastaReadNext($FH)) {
			my $len = length $$ret[1];
			warn "[!]  read Virus[$$ret[0]], $len bp.\n";
			print O ">$$ret[0]\n$$ret[1]\n";
			$main::RefConfig->{$main::RefFilesSHA}->{$$ret[0]} = $len;
			push @ChrIDs,$$ret[0];
		}
		$main::RefConfig->{$main::RefFilesSHA}->{VirusChrIDs} = join(',',@ChrIDs);
		close $FH;
		$main::RefConfig->write("$main::RootPath/Ref/Ref.ini");
	} else {
		warn "[!] Already Read References Pairs:[$main::HostRefName,$main::VirusRefName].\n";
	}
	#ddx \$main::RefConfig;
	unless ($main::DISABLE_REF_INDEX) {
		warn "[!] Building index for [$Refile].\n";
		system("$main::PathPrefix $RealBin/bin/bwameth.py index $Refile");	# unless -s $Refile.'.bwameth.c2t.sa';
	} else {
		warn "[!] Bwameth.py Index Building skipped for [$Refile].\n";
	}
	system("$main::PathPrefix samtools faidx $Refile") unless -s $Refile.'.fai';
}

sub do_aln() {
	my $Refilename = warnFileExist($main::RefConfig->{$main::RefFilesSHA}->{'Refilename'});
	#warn "$Refilename\n";
	my (%tID,%FQc,%maxReadNum);
	for (@{$main::Config->{'DataFiles'}->{'='}}) {
		/([^.]+)\.(\d)/ or die;
		$tID{$1}{$2} = $_;
		@{$FQc{$1}{$2}} = split /\s*,\s*/,$main::Config->{'DataFiles'}->{$_};
		$maxReadNum{$1} = 0 unless exists $maxReadNum{$1};
		$maxReadNum{$1} = $2 if $maxReadNum{$1} < $2;
	}
	#ddx \%tID;
	File::Path::make_path("$main::RootPath/${main::ProjectID}_aln",{verbose => 0,mode => 0755});
	open O,'>',"$main::RootPath/${main::ProjectID}_aln.sh" or die $!;
	print O "#!/bin/sh\n\nexport $main::PathPrefix\n\n";

	for my $k (keys %tID) {
		if ($maxReadNum{$k} == 1) {	# SE
			my @FQ2c;
			for my $f1 (@{$FQc{$k}{1}}) {
				#my $fID = basename($f1);
				my $f2 = " ";
				push @FQ2c,$f2;
			}
			@{$FQc{$k}{2}} = @FQ2c;
		}
	}
	#ddx \%FQc;
	for my $k (keys %tID) {
		#my @FQ1c = split /\s*,\s*/,$main::Config->{'DataFiles'}->{$tID{$k}{1}};
		#my @FQ2c = split /\s*,\s*/,$main::Config->{'DataFiles'}->{$tID{$k}{2}};
		my @FQ1c = @{$FQc{$k}{1}}; my @FQ2c = @{$FQc{$k}{2}};
		die "[x]  DataFiles not paired ! [@FQ1c],[@FQ2c]\n" unless $#FQ1c == $#FQ1c;
		my $cmd;
		if (@FQ1c == 1) {
			$cmd = <<"CMD";
$RealBin/bin/bwameth.py --reference $Refilename -t 24 --read-group $k -p $main::RootPath/${main::ProjectID}_aln/$k @{[warnFileExist($FQ1c[0],$FQ2c[0])]} 2>$main::RootPath/${main::ProjectID}_aln/$k.log
CMD
			print O $cmd;
		} else {
			my @theBams;
			for my $i (0 .. $#FQ1c) {
				my $fID = basename($FQ1c[$i]);
				$fID =~ s/\.fq(\.gz)?$//i;
				$cmd = <<"CMD";
$RealBin/bin/bwameth.py --reference $Refilename -t 24 --read-group '\@RG\\tID:${k}_${i}_${fID}\\tSM:$k' -p $main::RootPath/${main::ProjectID}_aln/${k}_${i}_${fID} @{[warnFileExist($FQ1c[$i],$FQ2c[$i])]} 2>$main::RootPath/${main::ProjectID}_aln/${k}_${i}_${fID}.log
CMD
				push @theBams,"$main::RootPath/${main::ProjectID}_aln/${k}_${i}_${fID}.bam";
				print O $cmd;
			}
			my $theBamsJ = join(' ',@theBams);
			$cmd = <<"CMD";
samtools merge -n -l 9 $main::RootPath/${main::ProjectID}_aln/$k.bam $theBamsJ
CMD
			print O $cmd;
		}
		$cmd = <<"CMD";
$RealBin/bin/bamfilter.pl 30 5 $main::RootPath/${main::ProjectID}_aln/$k.bam $main::RootPath/${main::ProjectID}_aln/S_${k}.bam
samtools sort -l 0 -m 2G $main::RootPath/${main::ProjectID}_aln/S_${k}.bam -T $main::RootPath/${main::ProjectID}_aln/P_${k} -o $main::RootPath/${main::ProjectID}_aln/P_${k}.bam
samtools index $main::RootPath/${main::ProjectID}_aln/P_${k}.bam
CMD
=pod
		$cmd = <<"CMD";
samtools sort -m 2415919104 -n $main::RootPath/${main::ProjectID}_aln/$k.bam -O bam -T $main::RootPath/${main::ProjectID}_aln/$k.sn >$main::RootPath/${main::ProjectID}_aln/$k.sn.bam 2>>$main::RootPath/${main::ProjectID}_aln/$k.log

samtools view -b -F256 $main::RootPath/${main::ProjectID}_aln/$k.sn.bam >$main::RootPath/${main::ProjectID}_aln/$k.snPstat.bam 2>>$main::RootPath/${main::ProjectID}_aln/$k.log

CMD
=cut
		print O $cmd;
	}
# Grep step0 Begin
	File::Path::make_path("$main::RootPath/${main::ProjectID}_grep",{verbose => 0,mode => 0755});
	my $WorkINI = Galaxy::IO::INI->new();
	$WorkINI->{'Output'} = $main::Config->{'Output'};
	$WorkINI->{'Ref'} = $main::RefConfig->{$main::RefFilesSHA};
	$WorkINI->{'InsertSizes'} = $main::Config->{'InsertSizes'};
	my %BamFiles;
	for my $k (keys %tID) {
		my $myBamf = "$main::RootPath/${main::ProjectID}_aln/P_$k.bam";
		$BamFiles{$k} = $myBamf;
	}
	$WorkINI->{'BamFiles'} = \%BamFiles;
	$WorkINI->write("$main::RootPath/${main::ProjectID}_grep/ToGrep.ini");
	my $cli = "$RealBin/bin/bsanalyser -p grep $main::RootPath/${main::ProjectID}_grep/ToGrep.ini";
	print O "\n$cli\n";
# Grep step0 End
	close O;
	chmod 0755,"$main::RootPath/${main::ProjectID}_aln.sh";
	warn "[!] Please run [$main::RootPath/${main::ProjectID}_aln.sh] to do the aln.\n"
}

sub do_grep($) {
	my $cfgfile = $_[0];
	my (%tID,%tFH);
	for (@{$main::Config->{'DataFiles'}->{'='}}) {
		/([^.]+)\.(\d)/ or die;
		$tID{$1}{$2} = $_;
	}
	#   "780_T" => { 1 => "780_T.1", 2 => "780_T.2" },
	#   "s01_P" => { 1 => "s01_P.1", 2 => "s01_P.2" },
	for my $k (keys %tID) {
		my $myBamf = "$main::RootPath/${main::ProjectID}_grep/$k.bam";
		print "[$myBamf]\n";
		open OUT,'>',"${myBamf}.grep" or die "Error opening ${myBamf}.grep: $!\n";
		open( IN,"-|","$main::PathPrefix samtools view $myBamf") or die "Error opening $myBamf: $!\n";
		my ($lastgid,@hReads,@vReads);
		my ($fhReads,$rhReads,$fvReads,$rvReads)=(0,0,0,0);	# /\bYD:Z:f\b/
		while (<IN>) {
			chomp;
			my @dat = split /\t/;
			/\tZc:i:(\d+)\b/ or die "[x]TAG:Zc:i not found.\n";
			my $thisGroup = $1;
#print "$lastgid <- $thisGroup\n";
			#print $thisGroup,"\t",join("][",@dat),"\n";
			if ($lastgid and ($lastgid != $thisGroup)) {
				my $skipflag = 0;
				if ($main::GrepMergeBetter) {
					$skipflag = 1 if ($fhReads < 1 or $rhReads < 1);
				} else {
					$skipflag = 1 if @hReads < 2;
				}
#print "$skipflag $lastgid <- $thisGroup\n";
				unless ($skipflag) {
					my $MergedHds = grepmerge(\@hReads);
					#ddx $MergedHds;
					print OUT join("\t",$hReads[0]->[2],$_,@{$MergedHds->{$_}}),"\n" for sort { $a <=> $b } keys %{$MergedHds};
					#die;
				}
				@hReads=();
				@vReads=();
				$lastgid = $thisGroup;
				($fhReads,$rhReads,$fvReads,$rvReads)=(0,0,0,0);
			} else {
				$lastgid = $thisGroup;
			}
			if (/\bZd:Z:H\b/) {
				push @hReads,\@dat;
				if (/\bYD:Z:f\b/) {++$fhReads}
				else {++$rhReads;}
			} elsif (/\bZd:Z:V\b/) {
				push @vReads,\@dat;
				if (/\bYD:Z:f\b/) {++$fvReads}
				else {++$rvReads;}
			}
		}
		close IN;
		close OUT;
	}
}

sub do_grep0($) {
	my $cfgfile = $_[0];
	my (%tID,%tFH);
	for (@{$main::Config->{'DataFiles'}->{'='}}) {
		/([^.]+)\.(\d)/ or die;
		$tID{$1}{$2} = $_;
	}
	#   "780_T" => { 1 => "780_T.1", 2 => "780_T.2" },
	#   "s01_P" => { 1 => "s01_P.1", 2 => "s01_P.2" },
	File::Path::make_path("$main::RootPath/${main::ProjectID}_grep",{verbose => 0,mode => 0755});
	my $GrepResult = Galaxy::IO::INI->new();
	my %ReadsIndex;
	for my $k (keys %tID) {
		my $myBamf = "$main::RootPath/${main::ProjectID}_aln/$k.bam";
		my $InsMean = $main::Config->{'InsertSizes'}->{$k} or die;
		my $InsSD = $main::Config->{'InsertSizes'}->{"$k.SD"} or die; # SD cannot be 0, so no need to test with defined.
		#warn "$myBamf,$InsMean,$InsSD";
		$GrepResult->{$k} = {
			'InBam' => $myBamf,
			InsMean => $InsMean,
			InsSD => $InsSD,
		};
		$GrepResult->{$k} = { DatFile => "$main::RootPath/${main::ProjectID}_grep/$k.sam" };
		open( IN,"-|","$main::PathPrefix samtools view $myBamf") or die "Error opening $myBamf: $!\n";	# `-F768` later
		system( "$main::PathPrefix samtools view -H $myBamf >".$GrepResult->{$k}{'DatFile'} );
		open GOUT,'>>',$GrepResult->{$k}{'DatFile'} or die "$!";
		print GOUT join("\t",'@PG','ID:bsuit',"CL:\"grep $cfgfile\""),"\n";
		print STDERR "[!] Reading [$myBamf] ...";
		while (my $line = <IN>) {
			my @Dat1 = split /\t/,$line;
			my $flag = 0;
			my $maxSC=0;
			while ($Dat1[5] =~ /(\d+)S/g) {
				$maxSC = $1 if $maxSC < $1;
			}
			$flag |= 1 if $maxSC > $main::minSoftClip;	# 加上SC的话，内存不乐观
			if ($Dat1[6] eq '=') {
				#$flag |= 1 if abs(abs($Dat1[8])-$InsMean) > 3*$InsSD;
				$flag |= 2 if exists($main::VirusChrIDs{$Dat1[2]}) or exists($main::VirusChrIDs{$Dat1[6]});
			} else {
				$flag |= 4 if exists($main::VirusChrIDs{$Dat1[2]}) or exists($main::VirusChrIDs{$Dat1[6]});
				$flag &= ~1 if $flag & 4;	# 人的PE，不同染色体的hit直接扔掉不管。
			}
			next unless $flag;
			$flag |= 8 if $Dat1[5] !~ /^\d+M$/;	# 数据只是占位置，所以可以去掉次行
			my $id = join("\t",$k,$Dat1[0]);
			$ReadsIndex{$id} = [[0,$Dat1[2],$flag]];	# 初始化，数据只是占位置。
		}
		close IN;
		open IN,'-|',"$main::PathPrefix samtools view $myBamf" or die "Error opening $myBamf: $!\n";	# `-F768` later
		while (my $line = <IN>) {
			#my ($id, $flag, $ref, $pos, $mapq, $CIGAR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/,$line;
			#print "$id, $flag, $ref, $pos, $mapq, $CIGAR, $mref, $mpos, $isize\n";
			my @Dat1 = split /\t/,$line;
			my $r12R1;
			if ($Dat1[1] & 0x40) {
				$r12R1 = 1;
			} elsif ($Dat1[1] & 0x80) {
				$r12R1 = 2;
			} else {$r12R1 = 1;}
			#--$r12R1;
			#my $line2 = <IN>;
			#die '[x]SAM/BAM file not paired !' unless defined($line2);
			#my @Dat2 = split /\t/,$line2;
			#next if $Dat1[4]<$CFGminMAPQ;
			my $flag = 0;
			my $maxSC=0;
			while ($Dat1[5] =~ /(\d+)S/g) {
				$maxSC = $1 if $maxSC < $1;
			}
			$flag |= 1 if $maxSC > $main::minSoftClip;
			if ($Dat1[6] eq '=') {
				#$flag |= 1 if abs(abs($Dat1[8])-$InsMean) > 3*$InsSD;
				$flag |= 2 if exists($main::VirusChrIDs{$Dat1[2]}) or exists($main::VirusChrIDs{$Dat1[6]});
			} else {
				$flag |= 4 if exists($main::VirusChrIDs{$Dat1[2]}) or exists($main::VirusChrIDs{$Dat1[6]});
				$flag &= ~1 if $flag & 4;
			}
			next unless $flag;
			$flag |= 8 if $Dat1[5] !~ /^\d+M$/;	# soft-clip
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
		if ($main::DEBUG) {
			open $tFH{$k},'<',$GrepResult->{$k}{'DatFile'} or die "$!";
		}
	}
	for my $tk (keys %ReadsIndex) {
		my ($flag,%tmpH,%tmpV)=(0);
		for my $i (1 .. $#{$ReadsIndex{$tk}}) {
			#my $x = join("\t",$ReadsIndex{$tk}->[$i][2],$ReadsIndex{$tk}->[$i][3]);
			my $x = $ReadsIndex{$tk}->[$i][4];
			if (exists $main::VirusChrIDs{$ReadsIndex{$tk}->[$i][2]}) {
				$tmpV{$x} = $i;	# 先认为最高质量的唯一。
			} else {
				$tmpH{$x} = $i
			}
			$flag |= $ReadsIndex{$tk}->[$i][6];
		}
		if ((keys %tmpH) > 0) {	# Select max MapQ as 1st & favors Human to Virus.
			my @s = sort { $b <=> $a } keys %tmpH;
			$ReadsIndex{$tk}->[0][0] = $tmpH{$s[0]};
			$ReadsIndex{$tk}->[0][1] = join("\t",$ReadsIndex{$tk}->[$tmpH{$s[0]}][2],$ReadsIndex{$tk}->[$tmpH{$s[0]}][3]);
		} else {
			my @s = sort { $b <=> $a } keys %tmpV;
			$ReadsIndex{$tk}->[0][0] = $tmpV{$s[0]};
			$ReadsIndex{$tk}->[0][1] = join("\t",$ReadsIndex{$tk}->[$tmpV{$s[0]}][2],$ReadsIndex{$tk}->[$tmpV{$s[0]}][3]);
		}
		#my @s = sort sortChrPos(keys %tmp);
		$ReadsIndex{$tk}->[0][2] = $flag;
#ddx $ReadsIndex{$tk};
# BSuitLib.pm:213: [
#   [1, "chr15\t20858350", 15],	# index, ChrPos, 总flag
#   [1477321, 1, "chr15", 20858350, 0, "48M42S", 12],
#   [33704587, 1, "gi|59585|emb|X04615.1|", 652, 0, "45S45M", 11],
#   [34264944, 2, "gi|59585|emb|X04615.1|", 717, 60, "56M34S", 12],
# ]
	}
	warn "[!] Bam Reading done.\n";
	#ddx \$GrepResult;
	#ddx (\%RefChrIDs,\%main::VirusChrIDs);
	open BOUT,'>',"$main::RootPath/${main::ProjectID}_grep/blocks.ini" or die "$!";
	if ($main::DEBUG) {
		open BOUTDBG,'>',"$main::RootPath/${main::ProjectID}_grep/blocks.txt" or die "$!";
	}
	my @IDsorted = sort { sortChrPos($ReadsIndex{$a}->[0][1],$ReadsIndex{$b}->[0][1]) } keys %ReadsIndex;
	my ($Cnt,%hChrRange,%vChrRange,@Store,@PureVirReads)=(0);
	for my $cid (@IDsorted) {
		my @minCPR = split /\t/,$ReadsIndex{$cid}->[0][1];
		push @minCPR,$ReadsIndex{$cid}->[0][0];
		if (exists $main::VirusChrIDs{$minCPR[0]}) {
			push @PureVirReads,$cid;
			next;
		}
		my ($thisehPos,$thisvChr,$thissvPos,$thisevPos);
		for my $i (1 .. $#{$ReadsIndex{$cid}}) {
			my ($DatRef,$isHost);
			if (exists $main::VirusChrIDs{$ReadsIndex{$cid}->[$i]->[2]}) {
				$DatRef = \%vChrRange;
				$isHost = 0;
			} else {
				$DatRef = \%hChrRange;
				$isHost = 1;
			}
			while ( mergeIn($isHost,$DatRef,$ReadsIndex{$cid}->[$i],\@Store,$cid,$i) ) {	# 目前只有桥墩，没有中间部分
				++$Cnt;
				print BOUT "[B$Cnt]\nHostRange=",formatChrRange(\%hChrRange),
					"\nVirusRange=",formatChrRange(\%vChrRange),"\nSamFS=",
					join(',',map { my @t = split /\n/;my @f=split /\t/,$t[0];join(':',$f[0],$ReadsIndex{$t[0]}->[$t[1]]->[0]); } @Store),"\n\n";
				if ($main::DEBUG) {
					print BOUTDBG "[B$Cnt]\nHostRange=",formatChrRange(\%hChrRange),
						"\nVirusRange=",formatChrRange(\%vChrRange),"\nSamFS=",
						join(',',map { my @t = split /\n/;my @f=split /\t/,$t[0];join(':',$f[0],$ReadsIndex{$t[0]}->[$t[1]]->[0]); } @Store),"\n";
					for (@Store) {
						my @t = split /\n/;
						my @f=split /\t/,$t[0];
						seek($tFH{$f[0]},$ReadsIndex{$t[0]}->[$t[1]]->[0],0);
						my $str = readline $tFH{$f[0]};
						print BOUTDBG $str;
					}
					print BOUTDBG "\n";
				}
				%hChrRange = %vChrRange = @Store = ();
			}
		}
	}
	close BOUT;
	if ($main::DEBUG) {
		close $tFH{$_} for keys %tFH;
		close BOUTDBG;
	}
	#warn $GrepResult->write_string;
	#ddx $GrepResult;
	$GrepResult->write("$main::RootPath/${main::ProjectID}_grep.ini");
	warn "[!] grep done with [$Cnt] items.\n";
}

sub do_analyse {
	my $Refilename = warnFileExist($main::RefConfig->{$main::RefFilesSHA}->{'Refilename'});
	my (%tID,%tFH);
	for (@{$main::Config->{'DataFiles'}->{'='}}) {
		/([^.]+)\.(\d)/ or die;
		$tID{$1}{$2} = $_;
	}
	my %ReadsIndex;
	for my $k (keys %tID) {
		my $GrepResult = "$main::RootPath/${main::ProjectID}_grep/$k.sam";
		my $InsMean = $main::Config->{'InsertSizes'}->{$k} or die;
		my $InsSD = $main::Config->{'InsertSizes'}->{"$k.SD"} or die; # SD cannot be 0, so no need to test with defined.
		open $tFH{$k},'<',$GrepResult or die "$!";
	}
	File::Path::make_path("$main::RootPath/${main::ProjectID}_analyse",{verbose => 0,mode => 0755});
	my $BlockINI = Galaxy::IO::INI->new();
	my $BlockINIFN = "$main::RootPath/${main::ProjectID}_grep/blocks.ini";
	if ( -f $BlockINIFN ) {
		$BlockINI->read($BlockINIFN);
	} else {die "[x] Grep INI not found ! [$BlockINIFN]\n";}
	warn "[!] minHostDepth = $main::minHostDepth\n";
	warn "[!] Running in DEVELOP(chr18 only) mode !\n" if ($main::DEVELOP);
	open OA,'>',"$main::RootPath/${main::ProjectID}_analyse.txt" or die;
	for my $Bid (@{$BlockINI->{']'}}) {
		my @HostRange = split /,/,$BlockINI->{$Bid}->{'HostRange'};
		my @VirusRange = split /,/,$BlockINI->{$Bid}->{'VirusRange'};
		my @SamFS = split /,/,$BlockINI->{$Bid}->{'SamFS'};
		next if (@HostRange<1 or @VirusRange<1);	# Need to do assembly
		my ($maxHostDepth,$maxItem) = (0,-1);
		for my $i (0 .. $#HostRange) {
			my ($chr,$range,$depth) = split /:/,$HostRange[$i];
			if ($maxHostDepth < $depth) {
				$maxHostDepth = $depth;
				$maxItem = $i;
			}
			$HostRange[$i] = [$chr,(split /-/,$range),$depth];
		}
		if ($main::DEVELOP) {
			next if $HostRange[$maxItem][0] ne 'chr18';
		}
		next if $maxHostDepth < $main::minHostDepth;
		#ddx $maxItem,\@HostRange;
		File::Path::make_path("$main::RootPath/${main::ProjectID}_analyse/idba/$Bid",{verbose => 0,mode => 0755});
		system("$main::PathPrefix samtools faidx $Refilename $HostRange[$maxItem][0]:$HostRange[$maxItem][1]-$HostRange[$maxItem][2] >$main::RootPath/${main::ProjectID}_analyse/idba/Ref.fa");
		my $FH;
		open $FH,'<',"$main::RootPath/${main::ProjectID}_analyse/idba/Ref.fa" or die $!;
		my $retHost = FastaReadNext($FH);
		close $FH;
		my @retVirus;
		$FH=openfile($main::Config->{'RefFiles'}->{'VirusRef'});
		while (my $ret = FastaReadNext($FH)) {
			push @retVirus,$ret;
		}
		close $FH;
		open FHo,'>',"$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/Refo.fa" or die $!;
		open FHf,'>',"$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/Reff.fa" or die $!;
		open FHr,'>',"$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/Refr.fa" or die $!;
		for ($retHost,@retVirus) {
			my ($id,$seq) = @$_;
			print FHo ">$id\n$seq\n";
			my $seqF = $seq;
			$seqF =~ tr /Cc/Tt/;
			print FHf ">${id}_F\n$seqF\n";
			my $seqR = $seq;
			$seqF =~ tr /Gg/Aa/;
			print FHr ">${id}_R\n$seqR\n";
		}
		close FHo; close FHf; close FHr;
		my %ReadsbyID;
		my %FHO = (
			o => [undef,0],
			f => [undef,0],
			r => [undef,0],
		);
		open $FHO{'o'}->[0],'>',"$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/Readso.fa" or die $!;
		open $FHO{'f'}->[0],'>',"$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/Readsf.fa" or die $!;
		open $FHO{'r'}->[0],'>',"$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/Readsr.fa" or die $!;
		my $maxReadLen = 0;
		for (@SamFS) {
			my ($fid,$pos) = split /:/,$_;
			seek($tFH{$fid},$pos,0);
			my $str = readline $tFH{$fid};
			chomp $str;
			my $strlen = length $str;
			$maxReadLen = $strlen if $maxReadLen < $strlen;
			my @dat = split /\t/,$str;
			next if $dat[1] & 256;	# FLAG, s  | 0x0100 | the alignment is not primary
			next if $dat[4] < 30;	# MAPQ
			my $seq = $dat[9];
			$seq = revcom($seq) if $dat[1] & 16;
			if ($dat[1] & 64) {
				$ReadsbyID{$dat[0]}->[1]=$seq;
			} elsif ($dat[1] & 128) {
				$ReadsbyID{$dat[0]}->[2]=$seq;
			}
			if (($dat[2] eq $HostRange[$maxItem][0]) or (exists $main::VirusChrIDs{$dat[2]})) {
				++$ReadsbyID{$dat[0]}->[0];
			}
		}
		for my $id (keys %ReadsbyID) {
			next unless defined($ReadsbyID{$id}->[1]) or defined($ReadsbyID{$id}->[2]);	# and->or to include SE.
			$ReadsbyID{$id}->[1] = '' unless defined($ReadsbyID{$id}->[1]);
			$ReadsbyID{$id}->[2] = '' unless defined($ReadsbyID{$id}->[2]);
			next unless defined($ReadsbyID{$id}->[0]);
			if ($main::FORCE_UNMETH) {
				if ($id =~ /^sf/) {
					$ReadsbyID{$id}->[1] =~ tr/gG/aA/;
					$ReadsbyID{$id}->[2] =~ tr/cC/tT/;
				} else {
					$ReadsbyID{$id}->[1] =~ tr/cC/tT/;
					$ReadsbyID{$id}->[2] =~ tr/gG/aA/;
				}
			}
			my $type = guessMethyl( $ReadsbyID{$id}->[1] . revcom($ReadsbyID{$id}->[2]) );
			my $str = join("\n",">$id/1",$ReadsbyID{$id}->[1],">$id/2",$ReadsbyID{$id}->[2],'');
			my $FH;
			if ($type eq '1CT') {
				$FH = $FHO{'f'}->[0];
				++$FHO{'f'}->[1];
			} elsif ($type eq '2GA') {
				$FH = $FHO{'r'}->[0];
				++$FHO{'r'}->[1];
			} else {
				$FH = $FHO{'o'}->[0];
				++$FHO{'o'}->[1];
			}
			print $FH $str;
		}
		#close $FHO{$_}->[0] for keys %FHO;
		my %Assem;
		my $idbaRead = '-r';
		$idbaRead = '-l' if $maxReadLen > 128;
		for my $fro (keys %FHO) {
			close $FHO{$fro}->[0];
			next if $FHO{$fro}->[1] <= 2;
			my $reff = "$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/Ref$fro.fa";
			my $readsf = "$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/Reads$fro.fa";
			next unless -s $readsf;
			my $outp = "$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/$fro";
			system("$RealBin/bin/idba_hybrid $main::idbacmd $idbaRead $readsf --reference $reff -o $outp >/dev/null");
			next unless -s "$outp/scaffold.fa";
			my ($tFH,@asm);
			open $tFH,'<',"$outp/scaffold.fa" or die $!;
			while (my $ret = FastaReadNext($tFH)) {
				push @asm,$ret;
			}
			$Assem{$fro} = \@asm;
		}
		my $ret = doAlign(\%Assem,[$retHost],\@retVirus);
		#die;
		print OA "[$Bid]\nRefCut: $ret->[0] $ret->[1]\nVirusCut: $ret->[2] $ret->[3] +$ret->[4]\n\n" if $ret->[1] > 0;
		if ($main::DEBUG) {
			open DBG,'>',"$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/idba.fa" or die $!;
			print DBG ">Host\n$$retHost[1]\n\n";
			#print DBG ">Virus$_->[0]\n$_->[1]\n" for @retVirus;
			print DBG "\n";
			for my $fro (sort keys %Assem) {
				print DBG ">Asm_${fro}$_ $Assem{$fro}->[$_]->[0]\n$Assem{$fro}->[$_]->[1]\n" for 0 .. $#{$Assem{$fro}};
			}
			print DBG "\nResult: ",join(",",@$ret),"\n";
		}
		warn "[$Bid]\n";
	}
	close OA;
	close $tFH{$_} for keys %tFH;
}

1;

__END__
samtools view -h /share/users/huxs/work/bsvir/bsI/SZ0010_aln/780_T.bam '*' | samtools bam2fq -O - | gzip -9 > /share/users/huxs/work/bsvir/bsI/SZ0010_aln/780_T.unmap.fq.gz
samtools view -h /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.bam '*'|samtools bam2fq -O -|gzip -9 > /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap.fq.gz &

./bin/bwameth.py --reference /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa -t 24 --read-group 780_T -p ~/work/bsvir/bsI/SZ0010_aln/780_T.unmap ~/work/bsvir/bsI/SZ0010_aln/780_T.unmap.fq.gz 2>~/work/bsvir/bsI/SZ0010_aln/780_T.unmap.log &    #/
./bin/bwameth.py --reference /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa -t 24 --read-group s01_P -p /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap.fq.gz 2>/share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap.log &

samtools view -h /share/users/huxs/work/bsvir/bsI/SZ0010_aln/780_T.bam 'gi|86261677|emb|AJ507799.2|'
