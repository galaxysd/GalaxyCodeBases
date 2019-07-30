use strict;
use warnings;

my $qc_dir = shift;
my $output = shift;
open OUT,">$output" or die($!);

$qc_dir =~ s/\/$//g;
my @path = split /\//,$qc_dir;
my $now = $path[scalar @path - 2];
my @useless = splice @path,scalar @path - 2;
my $total_dir = join "/",@path;
my @content = readpipe("ls -F $total_dir");

my %chips;
for (@content){
	chomp($_);
	if ($_ =~ /^\d+CG\S+\/$/ or $_ =~ /^\d+CG\S+\@$/){
		$_ =~ s/\/$//;
		$_ =~ s/\@$//;
		my ($date) = $_ =~ /^(\d+)CG/;
		$chips{$_} = $date;
	}
}

my @dirs = ($now);
if ($now =~ /^\d+CG/){
	my ($nowtime) = $now =~ /^(\d+)CG/;
	foreach my $dir (sort {$chips{$b}<=>$chips{$a}} keys %chips){
		next if ($dir eq $now);
		next unless (-e "$total_dir/$dir/6record");
		if ($chips{$dir} eq $nowtime){
			push @dirs,$dir;
		}elsif ($chips{$dir} < $nowtime && scalar @dirs < 3){
			push @dirs,$dir;
		}
	}
}

my (%nFs,%aFs,%store);
foreach my $take_dir (@dirs){
	my $list = "$total_dir/$take_dir/family.lst";
	open LI,"<$list" or die($!);
	while (my $info = <LI>){
		chomp($info);
		next unless ($info =~ /\S+/);
		my ($M,$F,$C) = split /\s+/,$info;
		$aFs{$M} = "p$C.M";
		$aFs{$F} = "p$C.F";
		if ($take_dir eq $now){
			$nFs{$M} = "p$C.M";
			$nFs{$F} = "p$C.F";
		}
		$store{$M} = "$total_dir/$take_dir/4tsv";
		$store{$F} = "$total_dir/$take_dir/4tsv";
	}
	close LI;
}

my @now_fathers = sort keys %nFs;
my @all_fathers = sort keys %aFs;
print OUT "\t";
print OUT join("\t",@all_fathers),"\n";
foreach my $F1 (@now_fathers){
	print OUT "$F1\t";
	my @value;
	foreach my $F2 (@all_fathers){
		my (%one,%two);
		&get_geno(\%one,$F1);
		&get_geno(\%two,$F2);
		my ($err,$total) = (0,0);
		foreach my $locus (keys %one){
			if (defined $two{$locus}){
				$total++;
				unless ($one{$locus} eq $two{$locus}){
					$err++;
				}
			}
		}
		if ($total > 0){
			push @value,$err / $total;
		}else{
			push @value,"NA";
		}
	}
	print OUT join("\t",@value),"\n";
}
close OUT;

##########################################
sub get_geno {
	my ($hash,$tag) = @_;
	my $file = "$store{$tag}/$aFs{$tag}.tsv";
	open TE,"<$file" or die($!);
	while (<TE>){
		chomp;
		my @data = split /\t/,$_;
		next unless ($data[3] > 100);
		my @geno = splice @data,4;

		my %check;
		for (@geno){
			my @info = split /;/,$_;
			my @alleles = split /\//,$info[0];
			my @sort = sort @alleles;
			my $allele = join "/",@sort;
			$check{$allele}++;
		}
		foreach my $allele (keys %check){
			if ($check{$allele} == scalar @geno){
				$hash->{$data[0]} = $allele;
			}
		}
	}
	close TE;
}
