#!/usr/bin/perl
use strict;
use warnings;
use Math::Random qw(random_uniform);
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.1.5;

our $opts='i:o:n:c:r:bva';
our($opt_i, $opt_o, $opt_n, $opt_c, $opt_r, $opt_v, $opt_b, $opt_a);

our $help=<<EOH;
\t-i SNP PWM file (./snp.pwm)
\t-n Number of Diploid Individue (100)
\t-c Chromosome ID (undef for ALL)
\t-r Reference FASTA file
\t-o Output prefix with existed path (./out/Din)
\t-a Append ChrID to output path
\t  Out will be ./out/ChrID/Din.popsnp & ./out/ChrID/Din/ind_{NUM}_[AB].fa
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH
# \t-d dbSNP SQLite data file (_tdbSNP.sqlite)

ShowHelp();

$opt_i='./snp.pwm' if ! defined $opt_i;
$opt_o='./out/Din' if ! $opt_o;
$opt_n='100' if ! $opt_n;
$opt_c='ALL' if ! $opt_c;
my $out=$opt_o;
$out =~ s#(.*)/([^/]+)#$1/$opt_c/$2# if $opt_a;

print STDERR "From [$opt_i].[$opt_c] with [$opt_r] to [$out]/ind_{NUM}_[AB].fa of [$opt_n]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

system("mkdir -p $out");

#@a=random_uniform($n)
#http://doc.bioperl.org/bioperl-live/Bio/Tools/IUPAC.html
my %REV_IUB = (AA	=> 'A',
				TT	=> 'T',
				CC	=> 'C',
				GG 	=> 'G',
				AC	=> 'M',
				AG	=> 'R',
				AT	=> 'W',
				CG	=> 'S',
				CT	=> 'Y',
				'GT'	=> 'K',
				ACG	=> 'V',
				ACT	=> 'H',
				AGT	=> 'D',
				CGT	=> 'B',
				ACGT=> 'N',
			);
# chromosome01	203	T	59,54,5,0,0	T:0.635593,C:0.364407
my $start_time = [gettimeofday];

my %PWM;	# -> {chr}{pos} -> [] -> [base,SigmaP]
open P,'<',$opt_i or die "Error opening $opt_i: $!\n";
while (<P>) {
	my ($chr,$pos,$ref,$stat,$pwm)=split /\t/;
	if ($opt_c) {
		next if $chr ne $opt_c;
	};
	--$pos;	# substr starts from 0;
	$PWM{$chr}{$pos}=[];
	for (split /\,/,$pwm) {
		push @{$PWM{$chr}{$pos}},[split /:/];
	}
	my $i=0;
	for (@{$PWM{$chr}{$pos}}) {
		$i+=$$_[1];
		$$_[1]=$i;
	}
	#warn "[!]The sum is [$i].\n" if $i != 1;
}
close P;
warn "[!]PWM loaded.\n";

open SNP,'>',$out.'.popsnp' or die "Error opening ${out}.popsnp: $!\n";
my @FH;
for (my $i = $opt_n; $i > 0; $i--) {
	my ($fha,$fhb);
	my ($nameA,$nameB)=($out.'/ind_'.$opt_c.'_'.$i.'_A.fa',$out.'/ind_'.$opt_c.'_'.$i.'_B.fa');
	open $fha,'>',$nameA or die "Error opening $nameA: $!\n";
	open $fhb,'>',$nameB or die "Error opening $nameB: $!\n";
#print $fha ">Din_${i}_A\n"; print $fhb ">Din_${i}_B\n";
	push @FH,[$fha,$fhb];
}
warn "[!]Output files created.\n";

my %Genome;
open G,'<',$opt_r or die "Error opening $opt_i: $!\n";
while (<G>) {
	s/^>//;
	my $title = $_;
	my $seqname = $1 if($title =~ /^(\S+)/);
#	if ($opt_c) {
#		if
#	}
print STDERR "[!]Loading Seq. >$seqname ...\n";
	$/=">";
	my $seq=<G>;
	chomp $seq;
	$seq=~s/\s//g;
	$/="\n";
	$Genome{$seqname}=\$seq;
	#$seq='';
}
close G;

sub PrintArray ($$$) {
	my ($FHr,$strr,$useAB)=@_;
	if ($useAB) {
		for (@$FHr) {
			my ($fha,$fhb)=@$_;
			print $fha $$strr,"_A"; print $fhb $$strr,"_B";	# no \n here, turn to the `&PrintArray() unless $pos%80` below.
		}
	} elsif (ref($strr) eq 'SCALAR') {
		for (@$FHr) {
			my ($fha,$fhb)=@$_;
			print $fha $$strr; print $fhb $$strr;
		}
	} elsif (ref($strr) eq 'ARRAY') {
		for (@$FHr) {
			my ($fha,$fhb)=@$_;
			print $fha shift(@$strr);
			print $fhb shift(@$strr);
		}
	}
}

my (@theBases,$pwm,@SNP,@iuba);
for my $chr (keys %Genome) {
	my $GenomeR=$Genome{$chr};
	my $length=length $$GenomeR;
	print STDERR "[!]>$chr: $length .\b";
	my $t=">Din_${chr}";
	&PrintArray(\@FH,\$t,1);
	warn "\n";
	for (my $pos = 0; $pos < $length; $pos++) {
		&PrintArray(\@FH,\$/) unless $pos%80;	# first is 0, which % 80.
		my $refbase=substr $$GenomeR,$pos,1;
		if ($PWM{$chr}{$pos}) {
			my $i=$pos+1;	# I was planing to use $pos=0 and $PWM{$chr}{$pos-1} and then ($pos-1)+1. It seems that I had written those 3 lines in different time thus I forgot the last trick.
			#print SNP "$chr\t$i\t$refbase\t";	# refbase should be tabbed
			#$refbase='-'; &PrintArray(\@FH,\$refbase);
			$pwm=$PWM{$chr}{$pos};	# []->[base,SigmaP]
			my (%pwmBases,@pwmSvalues,$pwmSvalue);
			for (@$pwm) {
#warn @$_;
				$pwmBases{$$_[1]}=$$_[0];
				push @pwmSvalues,$$_[1];
			}
#warn @pwmSvalues;
			@theBases=();
			for my $uniform ( random_uniform(2*$opt_n) ) {
				for (@pwmSvalues) {
					$pwmSvalue=$_;
					last if $pwmSvalue >= $uniform;
				}
#warn $pwmSvalue;
				push @theBases,$pwmBases{$pwmSvalue};
			}
			my $flag=0;
			my @t=@theBases;
			@SNP=();
			while (@t) {
				@iuba=();
				push @iuba,shift @t;
				push @iuba,shift @t;
				my $iub=join '',(sort @iuba);
				$flag = 1 if $iub ne $refbase;
				push @SNP,$REV_IUB{$iub};
			}
			if ($flag == 1) {	# maybe useful ?
				print SNP "$chr\t$i\t$refbase\t",join(" ",@SNP),"\n";
			} else {
				warn "[!]Empty PSNP appears:$chr\t$i\t$refbase\n";
			}
#warn "$chr\t$pos\t$refbase ",join(" ",@SNP);

			&PrintArray(\@FH,\@theBases);
		} else {
			&PrintArray(\@FH,\$refbase);
		}
	}
}

close SNP;
for (@FH) {
	my ($fha,$fhb)=@$_;
	close $fha; close $fhb;
}

__END__
cat chrorder | while read a;do echo "#$ -N \"${a}_si\"" >./shell/${a}_si100.sh;echo "#$ -cwd -r y -l vf=2G,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o ./out/${a}_si100.log -e ./out/${a}_si100.err" >> ./shell/${a}_si100.sh; echo ./simindu.pl -n 100 -ba -c $a -r ./Reference/$a -o ./out/D${a} >> ./shell/${a}_si100.sh; done
