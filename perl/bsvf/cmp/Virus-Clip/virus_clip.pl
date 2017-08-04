#usage: perl virus_clip.pl input.txt blastn blastnDB annovar annovarDB

my $input = $ARGV[0];
my $blastn = $ARGV[1];
my $blastdb = $ARGV[2];
my $annovar = $ARGV[3];
my $annovardb = $ARGV[4];

my %left_count;
my %left_length;
my %left_seq;
my %right_count;
my %right_length;
my %right_seq;
my %virus_left_length;
my %virus_left_seq;
my %virus_right_length;
my %virus_right_seq;

my @temp = ();

open(INPUT, $input);

while (my $line = <INPUT>) {
	chomp($line);

	my @column = split("\t",$line);

	my $flag = $column[0];
	my $pos = $column[1];
	my $mapq = $column[2];
	my $cigar = $column[3];
	my $seq = $column[4];
	my $dir = "";
	my $length = length($seq);

	if ($flag & 16) {
		$dir = "REV";
	}
	else {
		$dir = "FWD"
	}

	my $clip1 = "";
	my $clip2 = "";

    $cigar =~ /^(\d+)(\D).*$/;	#left softclip

	if ($2 eq "S") {
		$clip1 = $1;
	}
	
	$cigar =~ /.*\D(\d+)(\D)$/;	#right softclip
	
	if ($2 eq "S") {
		$clip2 = $1;
	}	


	#count the number of bases that aligned to virus genome
	my $size = 0;
	my @cigar_flag = ($cigar =~ /(\d+\D)/g);

	foreach my $item (@cigar_flag) {
		if ($item =~ /(\d+)[DM]/) {
			$size += $1;
		}
	}
	

	if ($clip1 ne "") {
		my $breakpoint = $pos;				#pos denotes the leftmost matched bases i.e. the breakpoint on HBV
		$left_count{$breakpoint}++;
		if ($clip1 > $left_length{$breakpoint}) {
			$left_length{$breakpoint} = $clip1;
			$left_sequence{$breakpoint} = substr($seq, 0, $clip1);
		}

		my $v_length = $length - $clip1 - $clip2;
		if ($v_length > $virus_right_length{$breakpoint}) {
			$virus_right_length{$breakpoint} = $v_length;
			$virus_right_sequence{$breakpoint} = substr($seq, $clip1, $v_length);
		}
	}
	
	if ($clip2 ne "") {
		my $breakpoint = $pos + $size - 1;					#breakpoint denotes the leftmost matched bases i.e. the breakpoint on HBV
		$right_count{$breakpoint}++;
		if ($clip2 > $right_length{$breakpoint}) {
			$right_length{$breakpoint} = $clip2;
			$right_sequence{$breakpoint} = substr($seq, $length - $clip2, $clip2);
		}

		my $v_length = $length - $clip1 - $clip2;
		if ($v_length > $virus_left_length{$breakpoint}) {
			$virus_left_length{$breakpoint} = $v_length;
			$virus_left_sequence{$breakpoint} = substr($seq, $clip1, $v_length);
		}
	}

}

close INPUT;


open(ANNO_TEMP, ">virus_clip.annovar");

my $breakpoint = "";

foreach $breakpoint (keys %left_count) {
	open(TEMP, ">temp.txt");
	print TEMP "$left_sequence{$breakpoint}\n";
	close TEMP;

	system("$blastn -db $blastdb -query temp.txt -out blastn.out");

	my $human_chr = "";
	my $human_pos = "";

	open(BLASTN, "blastn.out");
	while (my $line = <BLASTN>) {
		chomp $line;

		if ($line =~ />/) {				#extract matched human chr
			my @column = split(" ",$line);
			$human_chr = $column[1];
		}
		elsif ($line =~ /^Sbjct/) {		#extract matched human position
			my @column = split("  ",$line);
			$human_pos = $column[1] + length($left_sequence{$breakpoint}) - 1;
			last;
		}
	}
	
	if ($human_chr ne "" && $human_pos ne "") {
		push (@temp, "Human\t$human_chr\t$human_pos\t$left_sequence{$breakpoint}\tVirus\tchrVirus\t$breakpoint\t$virus_right_sequence{$breakpoint}\t$left_count{$breakpoint}");
		print ANNO_TEMP "$human_chr\t$human_pos\t$human_pos\t0\t0\n";
	}
}


foreach $breakpoint (keys %right_count) {
	open(TEMP, ">temp.txt");
	print TEMP "$right_sequence{$breakpoint}\n";
	close TEMP;

	system("$blastn -db $blastdb -query temp.txt -out blastn.out");

	my $human_chr = "";
	my $human_pos = "";

	open(BLASTN, "blastn.out");
	while (my $line = <BLASTN>) {
		chomp $line;

		if ($line =~ />/) {				#extract matched human chr
			my @column = split(" ",$line);
			$human_chr = $column[1];
		}
		elsif ($line =~ /^Sbjct/) {		#extract matched human position
			my @column = split("  ",$line);
			$human_pos = $column[1];
			last;
		}
	}	

	if ($human_chr ne "" && $human_pos ne "") {
		push (@temp, "Virus\tchrVirus\t$breakpoint\t$virus_left_sequence{$breakpoint}\tHuman\t$human_chr\t$human_pos\t$right_sequence{$breakpoint}\t$right_count{$breakpoint}");
		print ANNO_TEMP "$human_chr\t$human_pos\t$human_pos\t0\t0\n";
	}
}


system("$annovar -buildver hg19 -out integration.gene -neargene 10000 virus_clip.annovar $annovardb");

open(OUTPUT, ">virus_clip.out");

print OUTPUT "Left_element\tLeft_chr\tLeft_pos\tLeft_seq\tRight_element\tRight_chr\tRight_pos\tRight_seq\tread_count\tGene_region\tHuman_gene\n";	#Header

open(ANNO, "integration.gene.variant_function");
my $i = 0;

while (my $line = <ANNO>) {
	chomp($line);
	
	my @column = split("\t",$line);

	print OUTPUT "$temp[$i]\t$column[0]\t$column[1]\n";
	$i++;
}