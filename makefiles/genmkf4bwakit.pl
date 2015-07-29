#!/usr/bin/env perl
use strict;
#use warnings;
#use IO::Unread qw(unread);
use Data::Dump qw(ddx);

#die "Usage: $0 <fq.gz path>\n" if @ARGV < 1;
my ($inp)=@ARGV;
$inp='fq';

my $Ref = 'ref/hs37d5.fa';
my $BwaKit = "./bwa.kit/run-bwamem -sad -t12";
my $SamTools = './bwa.kit/samtools';
my $SamtoolsMerge = "$SamTools merge -l 9";
my $readsCounter = '/share/users/huxs/git/toGit/c_cpp/faststater/readsCounter';
my $GATK = 'java -Xmx16g -jar /opt/jar/GenomeAnalysisTK.jar';

opendir my($dh), $inp or die "Couldn't open dir '$inp': $!";
my @files = readdir $dh;
closedir $dh;
# http://perlmeme.org/faqs/file_io/directory_listing.html

@files = grep(/\.f(ast|)q\b/i, @files);

#ddx \@files;

my (%Pairs,@AllTargets);
for ( @files ) {
	m~^([^/]+?)([\W_.]*?)([12])?(\.clean)?(\.f(ast|)q(\.gz)?)$~i or die $_;
	my ($m,$rp,$r,$ext) = ($1,$2,$3,$5);
	#$m =~ s/\W//g and die $m;	# 减号也是\W, >_<
	print "$1, $2, $3, $4, $5, $6, $7\t$m,$rp,$r,$ext, $_\n";
die unless ($r == 1) or ($r==2);	# no SE now
	push @{$Pairs{$m}},[$_,$r];
}

my (%SamplesFN,%Samples);
for (keys %Pairs) {
	my ($sid) = (split /_/,$_)[-2];
	push @{$SamplesFN{$sid}},$_;
}
@AllTargets = sort keys %SamplesFN;

ddx \%Pairs;
ddx \%SamplesFN;
print join(',',@AllTargets),"\n";

open S,'>','Samples.list.pre';
print S "$_\t\n" for @AllTargets;
close S;

if (-s 'Samples.list') {
	open S,'<','Samples.list';
	while (<S>) {
		chomp;
		my ($id,$sample)=split /\t/;
		$Samples{$sample} = $SamplesFN{$id};
	}
} else {
	$Samples{$_} = $SamplesFN{$_} for @AllTargets;
}
ddx \%Samples;
print join(',',@AllTargets),"\n";
@AllTargets = sort keys %Samples;
print join(',',@AllTargets),"\n";

open M,'>','Makefile' or die $!;
print M 'all: ',join(' ',@AllTargets,'_AdditionalTG_'),"\n";
my (@AdditionalTG,@PHONY);

mkdir 'bam';
mkdir 'stat';
mkdir 'merged';

for my $sid ( @AllTargets ) {
	my @Runs =  @{$Samples{$sid}};
	for my $Target (@Runs) {
		my @Files = @{$Pairs{$Target}};
		die if @Files != 2;	# no SE now
		my (@fqFiles,@newFiles);
		for (@Files) {
			my ($fqname,$read12) = @$_;
			$newFiles[$read12-1] = "${Target}_$read12";
			$fqFiles[$read12-1] = $fqname;
			print M "
stat/${Target}_$read12.fqstat: $inp/$fqname
\t$readsCounter -o stat/${Target}_$read12.fqstat $inp/$fqname\n";
			push @AdditionalTG, "stat/${Target}_$read12.fqstat";
		}
	print M "
bam/${Target}.aln.bam: ",join(' ',(map { "$inp/$_" } @fqFiles)),"
\t$BwaKit -R'\@RG\\tID:$Target\\tLB:${sid}-1\\tPL:ILLUMINA\\tSM:$sid\\tPI:350' -o bam/${Target} $Ref ",join(' ',(map { "$inp/$_" } @fqFiles) )," | sh
";
	}
	print M "
merged/$sid.bam: ",join(' ',(map { "bam/$_.aln.bam" } @Runs)),"
\t$SamtoolsMerge merged/$sid.bam ",join(' ',(map { "bam/$_.aln.bam" } @Runs)),"
\t$SamTools index merged/$sid.bam
";
	print M "
merged/$sid.bam.bai: merged/$sid.bam
\t$SamTools index merged/$sid.bam
";
	push @PHONY,$sid;
	print M "$sid: merged/$sid.bam.bai\n";
}

print M "\n";
print M join(' ','_AdditionalTG_:',@AdditionalTG),"\n";
print M ".PHONY: all clean",join(' ',@PHONY),"\n";
print M "clean:
\t-rm bam/*.aln.*.bam
";
close M;

__END__
make -j4

$ ./bwa.kit/run-bwamem -sad -t24 -R'@RG\tID:FCAP086\tLB:FCAP086\tPL:ILLUMINA\tSM:FCAP086' -o bam/FCAP086_H2LGFCCXX_L3 ref/Felis_catus80_chr.fa fq/FCAP086_H2LGFCCXX_L3_1.fq.gz fq/FCAP086_H2LGFCCXX_L3_2.fq.gz
./bwa.kit/seqtk mergepe fq/FCAP086_H2LGFCCXX_L3_1.fq.gz fq/FCAP086_H2LGFCCXX_L3_2.fq.gz \
  | ./bwa.kit/trimadap 2> bam/FCAP086_H2LGFCCXX_L3.log.trim \
  | ./bwa.kit/bwa mem -p -t24 -R'@RG\tID:FCAP086\tLB:FCAP086\tPL:ILLUMINA\tSM:FCAP086' ref/Felis_catus80_chr.fa - 2> bam/FCAP086_H2LGFCCXX_L3.log.bwamem \
  | ./bwa.kit/samblaster 2> bam/FCAP086_H2LGFCCXX_L3.log.dedup \
  | ./bwa.kit/samtools sort -@ 4 -m1G - bam/FCAP086_H2LGFCCXX_L3.aln;

ARR=`echo {10814..10823}| tr ' ' ,` && pidstat -h -r -u  -p $ARR 1



@RG Read group. Unordered multiple @RG lines are allowed.
	ID* Read group identifier. Each @RG line must have a unique ID. The value of ID is used in the RG tags of alignment records. Must be unique among all read groups in header section. Read group IDs may be modified when merging SAM files in order to handle collisions.
	CN Name of sequencing center producing the read.
	DS Description.
	DT Date the run was produced (ISO8601 date or date/time).
	FO Flow order. The array of nucleotide bases that correspond to the nucleotides used for each flow of each read. Multi-base flows are encoded in IUPAC format, and non-nucleotide flows by various other characters. Format: /\*|[ACMGRSVTWYHKDBN]+/
	KS The array of nucleotide bases that correspond to the key sequence of each read.
	LB Library.
	PG Programs used for processing the read group.
	PI Predicted median insert size.
	PL Platform/technology used to produce the reads. Valid values: CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT, ONT, and PACBIO.
	PM Platform model. Free-form text providing further details of the platform/technology used.
	PU Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier.
	SM Sample. Use pool name where a pool is being sequenced

http://gatkforums.broadinstitute.org/discussion/1317/collected-faqs-about-bam-files
Dad's data:
@RG	ID:FLOWCELL1.LANE1	PL:ILLUMINA	LB:LIB-DAD-1 SM:DAD	PI:200
@RG	ID:FLOWCELL1.LANE2	PL:ILLUMINA	LB:LIB-DAD-1 SM:DAD	PI:200
@RG	ID:FLOWCELL1.LANE3	PL:ILLUMINA	LB:LIB-DAD-2 SM:DAD	PI:400
@RG	ID:FLOWCELL1.LANE4	PL:ILLUMINA	LB:LIB-DAD-2 SM:DAD	PI:400

Mom's data:
@RG	ID:FLOWCELL1.LANE5	PL:ILLUMINA	LB:LIB-MOM-1 SM:MOM	PI:200
@RG	ID:FLOWCELL1.LANE6	PL:ILLUMINA	LB:LIB-MOM-1 SM:MOM	PI:200
@RG	ID:FLOWCELL1.LANE7	PL:ILLUMINA	LB:LIB-MOM-2 SM:MOM	PI:400
@RG	ID:FLOWCELL1.LANE8	PL:ILLUMINA	LB:LIB-MOM-2 SM:MOM	PI:400

