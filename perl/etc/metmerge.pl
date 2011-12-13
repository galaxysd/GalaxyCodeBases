#!/bin/env perl
use strict;
use warnings;

use Data::Dumper;
#
die "Usage: $0 <mode(pe,se,total)> <input> <output>\n" if @ARGV != 3 or $ARGV[0] !~ /pe|se|total/;

my ($Mode,$Infile,$Outfile)=@ARGV;
warn "Mode=[$Mode], [$Infile]->[$Outfile]\n";

if ($Infile=~/.bz2$/) {
	open( IN,"-|","bzip2 -dc $Infile") or die "Error opening $Infile: $!\n";
} elsif ($Infile=~/.gz$/) {
	open( IN,"-|","gzip -dc $Infile") or die "Error opening $Infile: $!\n";
} else {open( IN,"<",$Infile) or die "Error opening $Infile: $!\n";}

sub splitLine($$) {
	my ($line,$isPE)=@_;
	my @lineItem=split /\t/,$line;
	my @Head=split /:/,$lineItem[0];
	my $CW=pop @Head;
	pop @Head; pop @Head;	# >FCC01P5ACXX:8:1101:10000:187253#CGATGTAT/1:rep1:+:len68:W	
	my $rep=pop @Head;
	$rep =~ s/^rep// or die "[x]ID format error.\n";
	my $ID=join ':',@Head;
	if ($isPE) {
		$ID =~ s#/[12]$##;
	}
	return [$ID,$CW,$rep,\@lineItem];
}

sub mergePEitems($) {
	my @Items=@{$_[0]};
	die "[x]PE lines missing.\n" if scalar @Items % 2;
	my $Paires = (scalar @Items) / 2;
	my (@OutItem,@a,@b);
	for (my $i=0;$i<$Paires;$i++) {
		@a=@{$Items[$i]};
		@b=@{$Items[$i + $Paires]};
		die "[x]PE data error.\n" if ($a[1] ne $b[1]) or ($a[2] ne $b[2]);
		push @OutItem,[@a,$b[3]];
	}
	return \@OutItem;
}

my $lastLine;
sub ReadItems($) {
	my ($isPE)=@_;
	my ($lastID,$id,$line,@Items,@dat)=('','');

	if (defined $lastLine) {
		$line=$lastLine;
		$lastLine=undef;
	} else {
		$line=<IN>;
		return [] unless defined $line;
		chomp $line;
	}
	($id,@dat)=@{&splitLine($line,$isPE)};
	$lastID=$id;
	push @Items,[$id,@dat];

	while ($lastID eq $id) {
		$line=<IN>;
		last unless defined $line;
		chomp $line;
		($id,@dat)=@{&splitLine($line,$isPE)};
		if ($lastID eq $id) {
			push @Items,[$id,@dat];
		} else {
			$lastLine=$line;
		}
	}
	if ($isPE) {
		my $ref=mergePEitems(\@Items);
		return $ref;	# [] of [$ID,$CW,$rep,\@lineItem1,\@lineItem2]
	} else {
		return \@Items;	# [] of [$ID,$CW,$rep,\@lineItem]
	}
}

sub splitTLine($) {
	my ($line)=@_;
	my @lineItem=split /\t/,$line;
	my $ID=$lineItem[0];	# >FCC01P5ACXX:8:1101:10000:15797#CGATGTAT/ 
	my $t=@lineItem;
	if ($t == 19) {	# PE
		return [$ID,1,\@lineItem];
	} elsif ($t == 12) {	# SE
		return [$ID,0,\@lineItem];
	} else {
		die "[x]Line format error.\n";
	}
}
sub LoadTItems() {
	my ($lastID,$id,$line,@Items,@dat)=('','');

	if (defined $lastLine) {
		$line=$lastLine;
		$lastLine=undef;
	} else {
		$line=<IN>;
		return [] unless defined $line;
		chomp $line;
	}
	($id,@dat)=@{&splitTLine($line)};
	$lastID=$id;
	push @Items,[$id,@dat];

	while ($lastID eq $id) {
		$line=<IN>;
		last unless defined $line;
		chomp $line;
		($id,@dat)=@{&splitTLine($line)};
		if ($lastID eq $id) {
			push @Items,[$id,@dat];
		} else {
			$lastLine=$line;
		}
	}
	return \@Items;	# [] of [$ID,$isPE,\@lineItem]
}

open OUT,'>',$Outfile or die "Error opening $Outfile: $!\n";

sub main_se() {
	my @dat;
	my $count=0;
	while(@dat=@{&ReadItems(0)}) {
		++$count;
		print Dumper(\@dat),'-' x 75,"$count\n";
	}
	warn "\nTotal Groups: $count\n";
}
sub main_pe() {
	my @dat;
	my $count=0;
	while(@dat=@{&ReadItems(1)}) {
		++$count;
		print Dumper(\@dat),'-' x 75,"$count\n"
	}
	warn "\nTotal Groups: $count\n";
}
sub main_total() {
	my @dat;
	my $count=0;
	while(@dat=@{&LoadTItems()}) {
		++$count;
		print Dumper(\@dat),'-' x 75,"$count\n";
		if (@dat==1) {
			print OUT join("\t",@{$dat[0]->[2]}),"\n";
		} else {
			for (@dat) {
				if ($$_[1] == 1) {
					print OUT join("\t",@{$$_[2]}),"\n";
					last;
				}
			}
		}
	}
	warn "\nTotal Groups: $count\n";
}

if ($Mode eq 'pe') {
	main_pe();
} elsif ($Mode eq 'se') {
	main_se();
} elsif ($Mode eq 'total') {
	main_total();
}

close IN;
close OUT;
