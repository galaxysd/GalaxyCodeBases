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
		$main::RefConfig->{RefFilesSHA}->{VirusChrIDs} = join(',',@ChrIDs);
		close $FH;
		$main::RefConfig->write("$main::RootPath/Ref/Ref.ini");
	} else {
		warn "[!] Already Read References Pairs:[$main::HostRefName,$main::VirusRefName].\n";
	}
	#ddx \$main::RefConfig;
	warn "[!] Building index for [$Refile].\n";
	system("$RealBin/bin/bwameth.py",'index',$Refile);
	system("samtools",'faidx',$Refile) unless -s $Refile.'.fai';
}

sub do_aln() {
	my $Refilename = warnFileExist($main::RefConfig->{$main::RefFilesSHA}->{'Refilename'});
	#warn "$Refilename\n";
	my (%tID);
	for (@{$main::Config->{'DataFiles'}->{'='}}) {
		/([^.]+)\.(\d)/ or die;
		$tID{$1}{$2} = $_;
	}
	#ddx \%tID;
	File::Path::make_path("$main::RootPath/${main::ProjectID}_aln",{verbose => 0,mode => 0755});
	open O,'>',"$main::RootPath/${main::ProjectID}_aln.sh" or die $!;
	print O "#!/bin/sh\n\n";
	for my $k (keys %tID) {
		my @FQ1c = split /\s*,\s*/,$main::Config->{'DataFiles'}->{$tID{$k}{1}};
		my @FQ2c = split /\s*,\s*/,$main::Config->{'DataFiles'}->{$tID{$k}{2}};
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
samtools merge -l 9 $main::RootPath/${main::ProjectID}_aln/$k.bam $theBamsJ
samtools index $main::RootPath/${main::ProjectID}_aln/$k.bam
CMD
			print O $cmd;
		}
		$cmd = <<"CMD";
samtools sort -m 2415919104 -n $main::RootPath/${main::ProjectID}_aln/$k.bam -O bam -T $main::RootPath/${main::ProjectID}_aln/$k.sn >$main::RootPath/${main::ProjectID}_aln/$k.sn.bam 2>>$main::RootPath/${main::ProjectID}_aln/$k.log

CMD
		print O $cmd;
	}
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
				$flag |= 2 if exists($main::VirusChrIDs{$Dat1[2]}) or exists($main::VirusChrIDs{$Dat1[6]});
			} else {
				$flag |= 4 if exists($main::VirusChrIDs{$Dat1[2]}) or exists($main::VirusChrIDs{$Dat1[6]});
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
				$flag |= 2 if exists($main::VirusChrIDs{$Dat1[2]}) or exists($main::VirusChrIDs{$Dat1[6]});
			} else {
				$flag |= 4 if exists($main::VirusChrIDs{$Dat1[2]}) or exists($main::VirusChrIDs{$Dat1[6]});
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
		if ($main::DEBUG) {
			open $tFH{$k},'<',$GrepResult->{$k}{'DatFile'} or die "$!";
		}
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
			while ( mergeIn($isHost,$DatRef,$ReadsIndex{$cid}->[$i],\@Store,$cid,$i) ) {
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
		system("samtools faidx $Refilename $HostRange[$maxItem][0]:$HostRange[$maxItem][1]-$HostRange[$maxItem][2] >$main::RootPath/${main::ProjectID}_analyse/idba/Ref.fa");
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
			next if $dat[1] & 256;
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
			next unless defined($ReadsbyID{$id}->[1]) and defined($ReadsbyID{$id}->[2]);
			next unless defined($ReadsbyID{$id}->[0]);
			my $type = guessMethyl( $ReadsbyID{$id}->[1] . revcom($ReadsbyID{$id}->[2]) );
			my $str = join("\n",">$id/1",$ReadsbyID{$id}->[1],">$id/2",$ReadsbyID{$id}->[2],'');
			if ($type eq 'CT') {
				print $FHO{'f'}->[0] $str;
				++$FHO{'f'}->[1];
			} elsif ($type eq 'GA') {
				print $FHO{'r'}->[0] $str;
				++$FHO{'r'}->[1];
			} else {
				print $FHO{'o'}->[0] $str;
				++$FHO{'o'}->[1];
			}
		}
		#close $FHO{$_}->[0] for keys %FHO;
		my %Assem;
		my $idbaRead = '-r';
		$idbaRead = '-l' if $maxReadLen > 128;
		for my $fro (keys %FHO) {
			close $FHO{$fro}->[0];
			next unless $FHO{$fro}->[1];
			my $reff = "$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/Ref$fro.fa";
			my $readsf = "$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/Reads$fro.fa";
			my $outp = "$main::RootPath/${main::ProjectID}_analyse/idba/$Bid/$fro";
			system("$RealBin/bin/idba_hybrid $main::idbacmd -r $readsf --reference $reff -o $outp");
			next unless -s "$outp/scaffold.fa";
			my ($tFH,@asm);
			open $tFH,'<',"$outp/scaffold.fa" or die $!;
			while (my $ret = FastaReadNext($tFH)) {
				push @asm,$ret;
			}
			$Assem{$fro} = \@asm;
		}
		die;
	}
	close $tFH{$_} for keys %tFH;
}

1;

__END__
samtools view -h /share/users/huxs/work/bsvir/bsI/SZ0010_aln/780_T.bam '*' | samtools bam2fq -O - | gzip -9 > /share/users/huxs/work/bsvir/bsI/SZ0010_aln/780_T.unmap.fq.gz
samtools view -h /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.bam '*'|samtools bam2fq -O -|gzip -9 > /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap.fq.gz &

./bin/bwameth.py --reference /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa -t 24 --read-group 780_T -p ~/work/bsvir/bsI/SZ0010_aln/780_T.unmap ~/work/bsvir/bsI/SZ0010_aln/780_T.unmap.fq.gz 2>~/work/bsvir/bsI/SZ0010_aln/780_T.unmap.log &    #/
./bin/bwameth.py --reference /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa -t 24 --read-group s01_P -p /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap /share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap.fq.gz 2>/share/users/huxs/work/bsvir/bsI/SZ0010_aln/s01_P.unmap.log &

samtools view -h /share/users/huxs/work/bsvir/bsI/SZ0010_aln/780_T.bam 'gi|86261677|emb|AJ507799.2|'
