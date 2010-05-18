#!/usr/bin/perl -w

use strict;
use Data::Dumper;
#use TableExtract;
require "/panfs/GAG/junli/raid010/pipeline/flow_script/TableExtract.pm";
use Getopt::Long;
use FindBin qw($Bin $Script);

my %opts = ();
GetOptions(\%opts,"fqDir:s","outDir:s","backDir:s","gsoap:i","sort:i","statistic:i","ISconfig:s","ref1:s","ref2:s","chrOrder:s","soapVersion:i","lib:s","mem:s","queue:s");

#by fangxd  Fri May 23 16:46:23 HKT 2008
#modify by lijun3 Tue Dec  1 15:38:15 CST 2009

unless(defined $opts{fqDir} && $opts{outDir} && $opts{ref1} && $opts{chrOrder}){
		print "\n$0 : auto run soap after solexa sequencing finished and provide basic statistics on mapping result,if want to bzip the reads and alignment result,please provide the directory to deposited them\n";
		print "perl $0\n";
		print "\t-fqDir\t\tthe directory contain FASTQ files,usually locate in /share/raid9/solexa-work/Run_seq_dir/Project_solexa_fq/PROJECT/LIBRARY\n";
		print "\t-outDir\t\tthe directory for soap alignment result files\n";
		print "\t-backDir\tthe directory to compress the fastq files,if not necessary,you can omit it \n";
		print "\t-gsoap\t\t switch for gsoap\n";
		print "\t-sort\t\t switch for sort the soap result\n";
		print "\t-statistic\t switch basic statistics for the alignment\n";
		print "\t-ISconfig\t configure file for Insert Size,useful for being without Eland alignment,the config file looks like\n";
		print "\t\t\t\t\tFastqName\tLibrary\tisPE\tRead1Length\tRead2Length\tInsertSize\tMinIS\tMaxIS\n";
		print "\t-lib\t configure file for Insert Size\n\t\t\tformat: libraryName insertSize [min max]\n";
		print "\t-ref1\t references in fasta format\n";
		print "\t\t/share/raid8/fangxd/data/personGenome/data/human.fa\n";
		print "\t-ref2\t references in binary format(soapVersion2 need)\n";
		print "\t\t /share/raid8/fangxd/data/personGenome/data/human.fa.index\n";
		print "\t-chrOrder\tchromomose sort order: /share/raid8/fangxd/data/personGenome/data/chroder \n";
		print "\t-soapVersion\t 1 for old version,2 for the new version\n\n";
		print "\t-mem\t\tMemory request in qusb,default 7g\n";
		print "\t-queue\t\tdesignation the queue: gag_snp.q (the all.q will be added automatically \n";
		exit 0;
}

my $projectDir = $opts{fqDir};
my $outputDir = $opts{outDir};
my $compressOutDir = $opts{backDir};
my $chrOrder = $opts{chrOrder};
my $refSeqFileFlat = $opts{ref1} if (defined $opts{ref1}); 
my $refSeqFileBinary = $opts{ref2} if (defined $opts{ref2}); 
my $queue=" ";
my $addque="";
 if(defined $opts{queue})
{


	$addque = $opts{queue};

	$queue=" -q all.q -q $addque";

	
}


my $soapVersion = 2;
if(defined $opts{soapVersion}){
		$soapVersion = $opts{soapVersion};
		if( $soapVersion == 2 && ! defined $opts{ref2}){
				die "SOAP2 need a pre-formated reference\n";
		}elsif($soapVersion == 1 && ! defined $opts{ref1}){
				die "SOAP1 need a reference in fasta format\n";
		}
}

print "fqDir\t$opts{fqDir}\n";
print "outDir\t$opts{outDir}\n";
my $seconds = 1;


######################################### globle variable


my %fqInformation = ();
mkdir $outputDir unless (-d $outputDir);
my $fqPathListFile = "$outputDir/fastq.list";
my $fqPathListFileOld = "$outputDir/fastq.list.old";
my $statisticsFile = "$outputDir/Statistics.txt";
my $statisticsFileRate = "$outputDir/Statistics.rate";
my @newFastq = ();
my %soapCommand = ();
my $soapBinOld = "$Bin/subBin/soap1 -w 5";	# only allow 5 equal best hit to save cpu time
my $soapBinNew = "$Bin/subBin/soap2.20 -p 4";
my $zipBin = `which bzip2`;
chomp $zipBin;
my $sortSoapBin = "python $Bin/subBin/SortSoap.py";
#	python  SortSoap.py soapFile chrOrder read1Length read2Length outputFile
my %libInfo = ();
my %fastqInformation = ();
my @fileToCompress = ();

######################################### globle variable


findNewFastq();
generateSoapShell();
submitJob();
#CompressFile(\@fileToCompress,$compressOutDir) if (defined $opts{backDir});

sub generateSoapShell{
		if(@newFastq > 0){
				foreach my $line (@newFastq){	# #FastqID        Library isPE    Read1Length     Read2Length     Median  MinInsertSize   MaxInsertSize   FastqFullPath
						my $refSeqFile = $refSeqFileFlat;
						my $soapComm = '';
						my @info = split /\s+/,$line;
						my $shellFile = '';
						my $readLength = $info[3];
						my $soapOutFile = '';
						my $fqFullPath = $info[-1];
						my $fqBasename = (split '/',$fqFullPath)[-1];
						my $unmapReadsOutDir = "$outputDir/unmapReads";
						if ($readLength>=60 && $soapVersion == 1){
								warn "soap version 1.0 can not handle read longer than 60 bp.. exit\n";
								exit 0;
						}
						my $largeLib = '';
						if($info[7]>1500){
								$largeLib = " -R ";
						}
						my $currentOutDir = "$outputDir/soap/";
						mkdir $currentOutDir unless (-d $currentOutDir);

						if($info[2]==1){ # PE
								$currentOutDir = "$outputDir/soap/PE";
								mkdir $currentOutDir unless (-d $currentOutDir);
								$shellFile = "$currentOutDir/soap_$fqBasename.sh";
								my $fqFullPath2 = $fqFullPath;
								$fqFullPath2 =~ s/\_1\.fq/\_2\.fq/;
								$soapOutFile = "$currentOutDir/$fqBasename.soap";
								if($soapVersion == 1){
										my $seed = $readLength>30?12:10;
										$soapComm = "$soapBinOld -a $fqFullPath -b $fqFullPath2 -d $refSeqFileFlat -o $soapOutFile  -2 $currentOutDir/$fqBasename.single -t -c 52 -p 8 -s $seed -m $info[6] -x $info[7]";
								}elsif($soapVersion == 2){
										my $seedLength = 32;
										my $maxMismatch = $readLength>70?3:1;
										my $minLength = 40;
										$soapComm = "$soapBinNew -a $fqFullPath -b $fqFullPath2 -D $refSeqFileBinary -o $soapOutFile  -2 $currentOutDir/$fqBasename.single -m $info[6] -x $info[7] $largeLib -t -s $minLength -l $seedLength -v $maxMismatch"; 
					#this modify is just for no trim					
										#$soapComm = "$soapBinNew -a $fqFullPath -b $fqFullPath2 -D $refSeqFileBinary -o $soapOutFile  -2 $currentOutDir/$fqBasename.single -m $info[6] -x $info[7] $largeLib -l $seedLength -v $maxMismatch";
										if(defined($opts{gsoap})){
												mkdir $unmapReadsOutDir unless(-d $unmapReadsOutDir);
												$unmapReadsOutDir = "$unmapReadsOutDir/PE";
												mkdir $unmapReadsOutDir unless(-d $unmapReadsOutDir);
												my $unmapFile = "$unmapReadsOutDir/$fqBasename.unmap";
												#$soapComm .= " -u $unmapFile";
												$soapComm .= " -g 5 -u $unmapFile";
										}
								}
								push @fileToCompress,$fqFullPath,$fqFullPath2,$soapOutFile,"$currentOutDir/$fqBasename.single";
						}else{

								my $currentOutDir = "$outputDir/soap/SE";
								mkdir $currentOutDir unless (-d $currentOutDir);

								$shellFile = "$currentOutDir/soap_$fqBasename.sh";
								$soapOutFile = "$currentOutDir/$fqBasename.soap";

								if($soapVersion == 1){
										my $seed = $readLength>30?12:10;
										$soapComm = "$soapBinOld -a $fqFullPath -d $refSeqFileFlat -o $soapOutFile -t -p 8 -s $seed";
								}elsif($soapVersion == 2){
										my $seedLength = 32;
										my $maxMismatch = $readLength>70?3:1;
										my $minLength = 40;
										$soapComm = "$soapBinNew -a $fqFullPath -D $refSeqFileBinary -o $soapOutFile -t -s $minLength -l $seedLength -v $maxMismatch"; 
				#this modify if just for been
										#$soapComm = "$soapBinNew -a $fqFullPath -D $refSeqFileBinary -o $soapOutFile -l $seedLength -v $maxMismatch";
										if(defined($opts{gsoap})){
												mkdir $unmapReadsOutDir unless(-d $unmapReadsOutDir);
												$unmapReadsOutDir = "$unmapReadsOutDir/SE";
												mkdir $unmapReadsOutDir unless(-d $unmapReadsOutDir);
												my $unmapFile = "$unmapReadsOutDir/$fqBasename.unmap";
												#$soapComm .= " -u $unmapFile";
												$soapComm .= " -g 5 -u $unmapFile";
										}
								}
								push @fileToCompress,$fqFullPath,$soapOutFile;
						}
						my $memeryRequest = '';
						#if(defined $opts{mem}){
						#	$memeryRequest = $opts{mem};
						#	$memeryRequest =~ s/g//i;
						#}
						$memeryRequest = $soapVersion==1?"14g":"7g";

						if(defined $opts{mem}){
                                                        $memeryRequest = $opts{mem};
                                                        #$memeryRequest =~ s/g//i;
                                                }

						$soapCommand{$shellFile} = join "\t",$memeryRequest,$soapComm;
						$fastqInformation{$shellFile} = $line . "\t" . $soapOutFile;
				}
		}
}
sub submitJob{

		my %jobIDs = ();
		my $psudoJobID = 0;
		foreach my $shellFile (keys %soapCommand){
				my $soapOutFile = (split "\t",$fastqInformation{$shellFile})[-1];
				if(-f $soapOutFile){
						$psudoJobID++;
						$jobIDs{"psudoJobID.$psudoJobID"} = $fastqInformation{$shellFile};
				}else{
						print "$shellFile\n";
						my ($memeryRequest,$soapComm) = split "\t",$soapCommand{$shellFile};
						open OUT,">$shellFile" or die "Error in writing [ $shellFile ]:$!\n";
						print OUT "$soapComm\n";
						close OUT;
						#next;
						my $qsub = "qsub$queue -cwd -l vf=$memeryRequest $shellFile";
						my $jobID = `$qsub`;
						print "$qsub jobID $jobID\n";
						$jobID = (split /\s+/,$jobID)[2];
						while($jobID !~ /\d+/){
								$jobID = `$qsub`;
								$jobID = (split /\s+/,$jobID)[2];
								print "\tresub $qsub\n";
								sleep 10;
						}
						$jobIDs{$jobID} = $fastqInformation{$shellFile};
						sleep 10;
				}
		}
		waitJobsDoneSortSoap2(%jobIDs);
}

sub waitJobsDoneSortSoap2{
		my (%jobsID) = @_;
		while(keys %jobsID){

		foreach my $jobid (keys %jobsID){
			my $qstat = `qstat -j $jobid`;
			chomp $qstat;

               		my @qstat_lines = split /\n/, $qstat;


				if(defined($qstat_lines[1]))
				{
		      #  print "test : $qstat_lines[1]\n";       
	        	        my $job = (split /\s+/,$qstat_lines[1])[1];
                
        	      #  print "job: $job\n";
                        	}

				else
				{

				print "Finish jobid $jobid\n";
        	                sortSoap($jobsID{$jobid}) if(defined $opts{sort});
                	       mapStatics($jobsID{$jobid}) if (defined $opts{statistic});
	                       delete $jobsID{$jobid};
		


				}
			
	

			

			}
	
				sleep $seconds;
		}

	
}






sub waitJobsDoneSortSoap{

		my (%jobsID) = @_;
		while( keys %jobsID > 0){
				my $qstat = `qstat`;
				chomp $qstat;
				my @qstatLines = split /\n/,$qstat;
				if(@qstatLines < 2){
						sleep 10;
						next;
				}
				my %runningJobs = ();
				foreach my $line (@qstatLines){
						my $job = (split /\s+/,$line)[0];
						next if (! $job);
						if($line =~ /Eqw/){
								my $qresub = `qresub $job`;
								chomp $qresub;
								my $qresubJobID = (split /\s+/,$qresub)[2];
								$jobsID{$qresubJobID} = $jobsID{$job};
								delete $jobsID{$job};
								`qdel $job`;
								$runningJobs{$qresubJobID} = 1;
						}else{
								if(exists $jobsID{$job}){
										$runningJobs{$job} = 1;
								}
						}
				}
				foreach my $jobID(keys %jobsID){
						if(! exists $runningJobs{$jobID}){
								print "Finish jobid $jobID\n";
								if ($jobID !~ /^\d+$/ && $jobID !~ /psudoJobID/){
										delete $jobsID{$jobID};
										next;
								}
								sortSoap($jobsID{$jobID}) if(defined $opts{sort});
								mapStatics($jobsID{$jobID}) if (defined $opts{statistic});
								delete $jobsID{$jobID};
						}
				}
				sleep $seconds;
		}
}






sub mapStatics{

#die "exit in mapStatics\n";
		my ($line) = @_;
		my @info = split /\s+/,$line;
		my $fastqID = $info[0];
		my $library = $info[1];
		my $isPE = $info[2];
		my %finishedStatistics = ();
		if(-f $statisticsFile){
				open IN,"$statisticsFile" or die "Error in reading [ $statisticsFile ] $!\n";
				while(my $line = <IN>){
						if ($line !~ /^\d+/){ 
								next;
						}
						chomp $line;
						$finishedStatistics{(split /\s+/,$line)[0]} = 1;
				}
				close IN;
		}
		return 0  if (exists $finishedStatistics{$fastqID});

		print "$fastqID\tStatistics...\n";
		my ($fastqFile,$soapOutFile) = @info[-2,-1];
		my @soapOutFiles = ();
		push @soapOutFiles,$soapOutFile;
		if($isPE){
				my $soapSingleFile = $soapOutFile;
				$soapSingleFile =~ s/\.soap/\.single/;
				push @soapOutFiles,$soapSingleFile;
		}

		my %readLength = ();
		$readLength{a} = $info[3];
		$readLength{b} = $info[4];

		my $tmpCount = `wc -l $fastqFile`;
		chomp $tmpCount;
		my $readCount = (split /\s+/,$tmpCount)[0];
		$readCount /= 4;
		my $baseCount = $readCount * $readLength{a};
		if($isPE){
				$baseCount = $readCount * ($readLength{a}+$readLength{b});
				$readCount *= 2; 
		}

		my $mapReadsCount = 0;
		my $mapBasesCount = 0;
		my $uniqMapReadsCount = 0;
		my $uniqMapBasesCount = 0;
		my $trimReadsCount = 0;
		my $trimBasesCount = 0;
		my %mismatchCount = ();
		my $dupRate = "NA";
		$mismatchCount{0} = 0;
		$mismatchCount{1} = 0;
		$mismatchCount{2} = 0;

		my $peMapReadsCount = 0;
		my $peMapReadsRate = 0;

		my $peMapBasesCount = 0;

		foreach my $soapOutFile(@soapOutFiles){
				print "\treading $soapOutFile\n";
				open IN,"$soapOutFile" or die "Error in reading file [ $soapOutFile ] $!\n";
				my @lines = ();
				my $isPESoap = 0;
				if($isPE && $soapOutFile =~ /\.soap/){
						$isPESoap = 1;
				}
				my @checkDuplicate = ();
				while(my $line = <IN>){
						@lines = split /\s+/,$line;
						next if (@info <10 || @info > 12);
						$mismatchCount{$lines[9]}++;
						$mapReadsCount++;
						$mapBasesCount += $lines[5];
						if($isPESoap){
								$peMapBasesCount += $lines[5];
								$peMapReadsCount++; 
								if($.<=2000000){
										push @checkDuplicate,"$info[0],$info[7],$info[8]";
								}
						}
						if($lines[3] == 1){
								$uniqMapReadsCount++; 
								$uniqMapBasesCount+=$lines[5];
						}
						if($lines[5] != $readLength{$lines[4]}){
								$trimReadsCount++;
								$trimBasesCount+=$readLength{$lines[4]} - $lines[5];
						}
				}
				if($isPESoap){
						my $dup = 0;
						my $total = 0;
						my %pcrDup = ();
						foreach my $chr(keys %pcrDup){
								my %freq = ();
								for(my $i=0;$i<@checkDuplicate;$i+=2){
										my($id1,$chr1,$pos1) = split ',',$checkDuplicate[$i];
										my($id2,$chr2,$pos2) = split ',',$checkDuplicate[$i+1];
										if($id1 eq $id2 && $chr1 eq $chr2){
												my $key = '';
												if($pos1>$pos2){
														$key = "$pos2,$pos1";
												}else{
														$key = "$pos1,$pos2";
												}
												$pcrDup{$chr1}{$key}++;
										}
								}
								@checkDuplicate = ();	# release memory
										foreach my $chr(keys %pcrDup){
												foreach my $pair(keys %{$pcrDup{$chr}}){
														$total++;
														if( $pcrDup{$chr}{$pair} > 1 ){
																$dup += $pcrDup{$chr}{$pair} - 1;
														}
												}
										}
						}
						if($total == 0){
								$dupRate = "NA";
						}else{
								$dupRate = sprintf("%.4f",$dup / $total * 100);
						}
				}
		}
		if(-f $statisticsFile){
				open STAT,">>$statisticsFile" or die "Error in writting [ $statisticsFile ] $!\n";
				open RATE,">>$statisticsFileRate" or die "Error in writting [ $statisticsFileRate ] $!\n";
		}else{
				open STAT,">$statisticsFile" or die "Error in writting [ $statisticsFile ] $!\n";
				open RATE,">$statisticsFileRate" or die "Error in writting [ $statisticsFileRate ] $!\n";
				my $title1 = "Fastq\tLibrary\tIsPE\tLength1\tLength2\t#Reads\t#Bases\t#MapableReads\t#MapableBases\t#UniqMapReads\t#PEMapableBases\t#UniqMapableBases\t#TrimedBases\t#U0Reads\t#U1Reads\t#U2Reads\n";
				print STAT "$title1\n";
				my $title2 = "Fastq\tLibrary\tIsPE\tLength1\tLength2\t#Reads\t#Bases\tMapableReads(%)\tMapableBases(%)\tpeMapBasesRate(%)\tUniqMapRead\tUniqMapableBases(%)\tTrimedBases(%)\tU0Reads(%)\tU1Reads(%)\tU2Reads(%)\tPCRdupRate\n";
				print RATE "$title2\n";
		}

		my $mapReadsRate = sprintf("%.2f",$mapReadsCount / $readCount * 100);
		my $mapBasesRate = sprintf("%.2f",$mapBasesCount / $baseCount * 100);
		my $uniqBasesRate = sprintf("%.2f",$uniqMapBasesCount / $mapBasesCount * 100);
		my $trimBasesRate = sprintf("%.2f",$trimBasesCount / $mapBasesCount * 100);
		my $peMapBasesRate = sprintf("%.2f",$peMapBasesCount / $mapBasesCount * 100);
		my $u0Rate = sprintf("%.2f",$mismatchCount{0} / $readCount * 100);
		my $u1Rate = sprintf("%.2f",$mismatchCount{1} / $readCount * 100);
		my $u2Rate = sprintf("%.2f",$mismatchCount{2} / $readCount * 100);
		my $uniqMapReadRate = sprintf("%.2f",$uniqMapReadsCount/$mapReadsCount);

		print STAT "$fastqID\t$library\t$isPE\t$readLength{a}\t$readLength{b}\t$readCount\t$baseCount\t$mapReadsCount\t$mapBasesCount\t$peMapBasesCount\t$uniqMapReadsCount\t$uniqMapBasesCount\t$trimBasesCount\t$mismatchCount{0}\t$mismatchCount{1}\t$mismatchCount{2}\n";
		print RATE "$fastqID\t$library\t$isPE\t$readLength{a}\t$readLength{b}\t$readCount\t$baseCount\t$mapReadsRate\t$mapBasesRate\t$peMapBasesRate\t$uniqMapReadRate\t$uniqBasesRate\t$trimBasesRate\t$u0Rate\t$u1Rate\t$u2Rate\t$dupRate\n";

		close STAT;
		close RATE;
}

sub sortSoap{
		my ($info) = @_;
# 081103_I328_FC30JLNAAXX_L5_HSKnomeR4AADCABPE_1.fq	HSKnomeR4AADCABPE	1	75	75	222	198	249	/share/raid9/solexa-work/Run_seq_dir/Project_solexa_fq/SZC08002_HSKnomeR4_resequencing_DNA/HSKnomeR4AADCABPE/081103_I328_0002_FC30JLNAAXX/081103_I328_FC30JLNAAXX_L5_HSKnomeR4AADCABPE/081103_I328_FC30JLNAAXX_L5_HSKnomeR4AADCABPE_1.fq	/share/raid11/resequencing/Knome4/NS12665//PE/081103_I328_FC30JLNAAXX_L5_HSKnomeR4AADCABPE_1.fq.soap
		my @info = split /\s+/,$info;
		my $isPE = $info[2];
		my ($read1Len,$read2Len) = @info[3,4];
		my $soapOutFile = $info[-1];
		my $soapOutFileBasename = (split '/',$soapOutFile)[-1];
		my $type = $isPE>0?"PE":"SE";
		my $sortOutDir = "$outputDir/SoapSort";
		mkdir $sortOutDir unless (-d $sortOutDir);
		$sortOutDir = "$outputDir/SoapSort/$type";
		mkdir $sortOutDir unless (-d $sortOutDir);
		my $shellFile = "$sortOutDir/sort_$soapOutFileBasename.sh";
		my $soapSortOutFile = "$sortOutDir/$soapOutFileBasename.sort";
#my $sortComm1 = "$sortSoapBin -in $soapOutFile -out $soapSortOutFile -order $chrOrder -a $read1Len -b $read2Len";
		my $sortComm1 = "$sortSoapBin $soapOutFile $chrOrder $read1Len $read2Len $soapSortOutFile";
		my $sortComm2 = "";
		my @commands = ();
		push @commands,$sortComm1 if(! -f $soapSortOutFile);
		if($isPE){
				my $soapSortOutFile2 = $soapSortOutFile;
				my $soapOutFile2 = $soapOutFile;
				$soapSortOutFile2 =~ s/\.soap/\.single/g;
				$soapOutFile2 =~ s/\.soap/\.single/g;
				$sortComm2 = "$sortSoapBin $soapOutFile2 $chrOrder $read1Len $read2Len $soapSortOutFile2";
				if(! -f $soapSortOutFile2){
						push @commands,$sortComm2;
				}
		}
		if(@commands > 0){
				open OUT,">$shellFile" or die "Error in writting file: [ $shellFile ] $!\n";
				print OUT "$sortComm1\n";
				print OUT "$sortComm2\n";
				close OUT;
				my $memeryRequest = "5g";
				if($read1Len > 60){
						$memeryRequest = "6g";
				}
				my $qsub = "qsub$queue -cwd -l vf=$memeryRequest  $shellFile";
				print "Sorting SoapOut\n$qsub\n";
				`$qsub`;
		}
}

sub findNewFastq{

		my %oldFq = ();
		my %newFq = ();
		my %oldFqListFileLibCount = ();
		my %newFqListFileLibCount = ();

		my @newFqFiles = `find $projectDir -name "*.fq"`;
		chomp @newFqFiles;
		print "#########################################\n";
		print "Fastqs\n";
		my $fqIndex = 0;
		map{if($_!~/_2\.fq/){$fqIndex++;print "\t$fqIndex\t$_\n"}}@newFqFiles;
		print "#########################################\n\n";
		foreach my $fqFile(@newFqFiles){
#next if ($fqFile =~ /\/maq/);	# tmep add
#next if ($fqFile =~ /FC20B3BAAXX/);	# tmep add
				my $fqBasename = (split '/',$fqFile)[-1];
				$newFq{$fqBasename} = $fqFile;
				my $libName = findLibName($fqBasename);
				$newFqListFileLibCount{$libName}++; 
		}

		if(-f $fqPathListFile){
				open IN,"$fqPathListFile" or die "Error in Reading [ $fqPathListFile ] $! $.\n";
				while(my $line = <IN>){
						chomp $line;
						my @info = split /\s+/,$line;
						$oldFq{$info[0]} = $info[1];
						my $libName = findLibName($info[0]);
						$oldFqListFileLibCount{$libName}++;
				}
				close IN;
		}
		my $old = keys %oldFq;
		my $new = keys %newFq;
		if(1){

				system("cp $fqPathListFile $fqPathListFile.old") if (-f $fqPathListFile);
				open FQLIST,">>$fqPathListFile" or die "Error in appending [ $fqPathListFile ] $!\n";
				foreach my $fq(sort keys %newFq){

						my $fqBasename = (split '/',$fq)[-1];
						next if ($fqBasename =~ /\_2.fq/);	
						if( ! exists $oldFq{$fqBasename}){

								if(! exists $fqInformation{$fqBasename}){
										getFqInformation($newFq{$fq});
								}

								if(! exists $fqInformation{$fqBasename}){
										warn "\t!!!!!Lacking information of [ $fqBasename ] !!!!!!\n";
										next;
								}else{
										$fqInformation{$fqBasename} = $fqInformation{$fqBasename} . "\t". "$newFq{$fq}"; 
										my @elementCount = split /\s+/,$fqInformation{$fqBasename};
										if(@elementCount < 8){
												warn "\t!!!!!Lacking information of $fqBasename,please check it !!!!!\n";
												next;
										}
										print FQLIST "$fqBasename\t$fqInformation{$fqBasename}\n";
										push @newFastq,"$fqBasename\t$fqInformation{$fqBasename}";	# find new Read to run soap
								}
						}
				}
				close FQLIST;
		}
}

sub getFqInformation{
		my ($fqFile) = @_;
		my $isPE = 0;
		if($fqFile =~ /\_1.fq/){
				$isPE = 1;
		}
		my $info = '';
		my $fqBasename = (split '/',$fqFile)[-1];
		#my $library = (split '/',$fqFile)[-4];
		my $library = (split /\_/, $fqFile)[-2];
		if($isPE == 0){
				open IN,"$fqFile" or die "Error in reading [ $fqFile ] $!\n";
				<IN>;
				my $seq = <IN>;
				chomp $seq;
				my $readLen = length $seq;
				$info = join "\t",$library,$isPE,$readLen,0,0,0,0;
				$fqInformation{$fqBasename} = $info;
		}else{
				my @path = (split '/',$fqFile);
				pop @path;
				my $reportPath = join '/',@path;
				pop @path;
				my $path = join '/',@path;
				my @html = `find $path -name Summary.htm`;
				chomp @html;
				my $htmlFile = $html[0];
				my %laneInfo = ();
				my @reportFile = glob("$reportPath/*.report");
				my $reportFile = '';
				if(@reportFile>0){
						$reportFile = $reportFile[0];
				}
				if(defined $opts{ISconfig} && -f $opts{ISconfig}){
						readConfig();
				}
				if(defined $opts{lib} && -f $opts{lib}){
						print "read lib\n";
						readLibInfo();
				}
				if(! exists $fqInformation{$fqBasename}){
						if(defined $opts{lib} && !exists $libInfo{$library}){
								print "no $library\n";
								return 0;
						}
						if(exists $libInfo{$library}){
	#modify 091019							open T,"$fqFile" or die "$!\n";
	
								open T,"$fqFile" or print "$!\n";
								my $tmp = <T>;
								$tmp = <T>;
								close T;
								chomp $tmp;
								my $readLen = length $tmp;
								$fqInformation{$fqBasename} = join "\t",$library,$isPE,$readLen,$readLen,$libInfo{$library};
								print "$fqBasename\t$fqInformation{$fqBasename}\n"; 
						}elsif(-f $reportFile){
								print "ReportFile $reportFile\n";
								open REPORT,"$reportFile" or die "Error in reading [ $reportFile ] $!\n";
								my $laneIndex = '';
								my $insertSize = 0;
								my @infoArray = ();
								while(my $line = <REPORT>){
										chomp $line;
										if($line =~ /Lane\s+(\d+)/){
												$laneIndex = $1;
										}
										if($line =~ /Length\s+(\d+)\;(\d+)/){
												my $readLength1 = $1;
												my $readLength2 = $2;
												push @infoArray,$readLength1,$readLength2;
										}
										if($line =~ /InsertSize\b/){
												if($line =~ /unknown/i){
														if (defined $opts{ISconfig} && -f $opts{ISconfig}){
																readConfig();
																print "Find Insert information in $opts{ISconfig}\n";
														}
														if(! exists $fqInformation{$fqBasename}){
																warn "\n!!!!  No insert size information in $reportFile !!!!\n\n";
														}
														last;
												}

												$insertSize = (split /\s+/,$line)[-1];
												push @infoArray,$insertSize;
										}
										if($line =~ /InsertSizeSD/){
												my $sd = (split /\s+/,$line)[-1];
												my ($blow,$above) = ($sd =~ /\-(\d+)\/\+(\d+)/);

												
												my $lowThresh = $insertSize - $blow * 3;
				if($lowThresh<0){ $lowThresh= 3*$blow - $insertSize;}
												my $highThresh = $insertSize + $above * 3;
												push @infoArray,$lowThresh,$highThresh;
										}
								}
								close REPORT;
								$laneInfo{$laneIndex} = join "\t",@infoArray;
								$fqInformation{$fqBasename} = join "\t",$library,$isPE,$laneInfo{$laneIndex};	# without fqFullpath

						}elsif(-f $htmlFile){
								print "html $htmlFile\n";
								extractInfoFromHtml($htmlFile,\%laneInfo);
								my ($currentLaneIndex) = ($fqBasename =~ /_L(\d+)_/);
								foreach my $laneIndex(sort{$a<=>$b} keys %laneInfo){
										my $fqNameCopy = $fqBasename;
										$fqNameCopy =~ s/_L${currentLaneIndex}_/_L${laneIndex}_/;
										$fqInformation{$fqNameCopy} = join "\t",$library,$isPE,$laneInfo{$laneIndex};	# without fqFullpath
								}
						}else{
								print "No Summary file found in [ $path ]\n";
						}
				}
		}
}
sub extractInfoFromHtml{

		my ($htmlFile,$hLaneInfo) = @_;
		my $te = HTML::TableExtract->new();
		$te->parse_file($htmlFile);
		open OUT,">te.out" or die "$!\n";
		foreach my $ts ($te->tables) {
				print OUT "Table (", join(',', $ts->coords), "):\n";
				foreach my $row ($ts->rows) {
						print OUT join(',', @$row), "\n";
				}
		}
		close OUT;
		open IN,"te.out" or die "$!\n";
		my $laneIndex = 0;
		my %laneReadLength = ();
		my @targetLane = ();	# deal with some lane without ELAND alignment
				while(my $line = <IN>){
						chomp $line;

#Table (0,2):
						if($line =~ /Table\s\(0,2\)\:/){
								<IN>;
								foreach my $index (1..8){
										my $info = <IN>;
#7,unknown,HUMAN_refseq_eland,ELAND_PAIR,26, 26,'((CHASTITY>=0.6))',4,Lane 7
										my @info = split ',',$info;
										if($info =~ /ELAND_PAIR/){
												$laneReadLength{$info[0]} = join "\t",@info[4,5];
												push @targetLane,$info[0];
										}
								}

						}elsif($line =~ /Below-median/){
								my $info = <IN>;
								chomp $info;
								my $laneNum = $targetLane[$laneIndex];
								my ($median,$lowThresh,$highThresh) = (split ',',$info)[0,3,4];
								$$hLaneInfo{$laneNum} = join "\t",$laneReadLength{$laneNum},$median,$lowThresh,$highThresh;
								$laneIndex++;
						}
				}
		close IN;
		unlink "te.out";
}
sub readConfig{
		open IN,"$opts{ISconfig}" or die "Error in reading $opts{ISconfig}\n";
		while(my $line = <IN>){
				chomp $line;
				my($fqBasename,$library,$isPE,$read1Length,$read2Length,$insertSize,$min,$max) = split /\s+/,$line;
				my $lowThresh;
				my $highThresh;
				

				if(defined($min) and defined($max))
				{
				$lowThresh = $min;		#$insertSize - $insertSize * 0.11;
				$highThresh= $max;		#$insertSize + $insertSize * 0.11;
				}
				else
				{

				$lowThresh = $insertSize - $insertSize * 0.11;
				$highThresh= $insertSize + $insertSize * 0.11;

				}
				
				$fqInformation{$fqBasename} = join "\t",$library,$isPE,$read1Length,$read2Length,$insertSize,$lowThresh,$highThresh;
		}
		close IN;
}
sub findLibName{
		my ($fqFile) = @_;
		my $libName = '';
		if($fqFile =~ /_1\.fq/){
				$libName = (split '_',$fqFile)[-2];
		}else{
				$libName = (split '_',$fqFile)[-1];
				$libName =~ s/\.fq//;
		}
		return $libName;
}
sub readLibInfo{
		open IN,"$opts{lib}" or die "$!\n";
		while(my $line = <IN>){
				chomp $line;
				#my ($lib,$ins) = split /\s+/,$line;
				#my $min = $ins - $ins * 0.08 ;
				#my $max = $ins + $ins * 0.08 ;
				my ($lib, $ins, $min, $max) = split /\s+/, $line;
				unless (defined $max) {
					$min = $ins - $ins * 0.08 ;
					$max = $ins + $ins * 0.08 ;
				}
				$libInfo{$lib} = "$ins\t$min\t$max";
		}
		close IN;
}
sub CompressFile{
		my ($a_file,$dir) = @_;
		my @files = @$a_file;
		my @bzipComm = ();
		if(@files>0){
				mkdir $dir unless(-f $dir);
				foreach my $file(@files){
						if( ! -f $file){
								warn "No $file \n";
								next;
						}
						my $fileBasename = (split '/',$file)[-1];
						my $outFileName = '';
						my $outDir = '';

						if($file =~ /\.fq$/){
								$outDir = "$dir/Fastq";
								mkdir $outDir unless(-f $outDir);
								if($file =~ /\_1\.fq/ || $file =~ /\_2\.fq/){
										$outDir .= '/PE';
								}else{
										$outDir .= '/SE';
								}
						}elsif($file =~ /soap$/ || $file =~ /single$/){
								$outDir = "$dir/Alignment";
								mkdir $outDir unless(-f $outDir);
								if($file =~ /\_1\.fq\.soap/ || $file =~ /single/){
										$outDir .= '/PE';
								}else{
										$outDir .= '/SE';
								}
						}

						mkdir $outDir unless(-f $outDir);
						$outFileName = "$outDir/$fileBasename.bz2";
						next if ( -f $outFileName);
						my $comm = "$zipBin -c $file > $outFileName";
						push @bzipComm,$comm;
				}
		}
		if(@bzipComm>0){
				foreach my $comm(@bzipComm){
						print "$comm\n";
						`$comm`;
				}
		}
}
