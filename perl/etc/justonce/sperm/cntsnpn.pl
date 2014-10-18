#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump;

#die "Usage: $0 <blood vcf.gz> <sperm vcf.gz>\n" if @ARGV < 2;
die "Usage: $0 <vcf1> <vcf2> [vcf3.gz ..]\n" if @ARGV < 2 or @ARGV > 32;
my @INS=@ARGV;
my $Length = scalar @INS;
my (@FH,%SNP);

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.xz$/) {
	    open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}

sub readfile($) {
	my $id = $_[0];
	my $fh = $FH[ $id ];
	print STDERR "Reading [$id] ...";
	my ($cntID,$totalID) = (0,0);
	while (<$fh>) {
		next if /^#/;
		my ($chr,$pos,undef,$ref,$alt,$QUAL,undef,$INFO,undef,$data) = split /\t/;
		next if $chr =~ /[XYM]/i;
		next if $INFO =~ /INDEL;/;
		++$totalID;
		if (defined $data) {
			next if $QUAL < 20;
			my $GQ = (split /:/,$data)[-1];
			next if $GQ < 20;
		}
		$chr =~ s/^chr//i;
		next unless $chr =~ /^\d+/;
		#my $GT;
		#($chr,$pos,$GT) = @{getGT($_)};
		my $GT = 1;
		$SNP{$chr}{$pos}{$id} = $GT;
		++$cntID;
	}
	print STDERR "\b\b\b\b, done ! [$totalID -> $cntID]\n";
	print OUT "[$INS[$id]]:[$totalID -> $cntID]\n";
}

for (@INS) {
	my $t;
	$t = openfile($_);
	push @FH,$t;
}
open OUT,'>','stat.txt';
print OUT "[@INS] -> [stat.txt]\n";

warn "[@INS] -> [stat.txt]\n";

readfile($_) for (0 .. $#FH);
close $_ for @FH;

#ddx \%SNP;
#              141050953 => { 1 => "TT", 2 => "TT" },

my %Stat;
for my $chr (keys %SNP) {
	for my $pos (keys %{$SNP{$chr}}) {
		my %Dat = %{$SNP{$chr}{$pos}};
		my $flag = 0;
		for my $k (keys %Dat) {
			$flag |= 2 ** $k;
		}
		++$Stat{$flag};
#ddx $SNP{$chr}{$pos}; print "$flag\n\n";
	}
}

for (sort {$a<=>$b} keys %Stat) {
	my $v = sprintf("%0${Length}b",$_);
	my $str = join("\t",$v,$Stat{$_},$_);
	print "$str\n";
	print OUT "$str\n" 
}
close OUT;

__END__
perl ./cntsnp3.pl sperm23.5cf.filter sperm24.5cf.filter sperm28.5cf.filter
2	243901
4	268409
6	38623
8	253689
10	33956
12	35961
14	7207

perl ./cntsnp3.pl sperms01.5cf.filter sperms02.5cf.filter sperms03.5cf.filter
2	676931
4	635059
6	60638
8	623031
10	65199
12	58892
14	36960
