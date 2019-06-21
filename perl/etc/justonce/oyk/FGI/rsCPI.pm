package main;

use strict;
#no strict "subs";
use warnings;
#use Data::Dump qw(ddx);
=head1 NAME && VERSION
   cpi.pm
   Version:  V1.1
=head1 DESCRIPTION
	chosen snp
=head1 AUTHOR && CONTACT
    Author :  gaoshengjie huxuesong 2018/05/06
    modified at 2019/06/18
    Contact:  gaoshengjie@genomics.org.cn
=head1 UPDATE LOG
    For missing record, MAF = 0.45, Ref/Alt follow VCF.
=cut

our (%Markers,%MarkerAF);
# %Markers
#   {
#     rs1000381   => ["chr5", 56885878, 0.124304011506483, 0.187151034246758],
#     rs10005624  => ["chr4", 125929227, 0.122668415174707, 0.186323232412646],
#     rs10011941  => ["chr4", 156398437, 0.116374475912992, 0.183033072043504],
#   }
# %MarkerAF
#   {
#     rs9995391   => { A => 0.5391, T => 0.4609 },
#     rs9929923   => { A => 0.4429, C => 0.5522, G => "0.0050" },
#     rs9996524   => { A => "0.5450", G => "0.4550" },
#   }

#print $CHR{"SNP5538"}{"ref"}."\n";
#SNP1110 501     G,A     999     G/G;9996,7      A/A;8,9955      G/A;5895,61;2.58311161233576e-19;G/G;5895,61
our $err=0;
our $special=0;
our $pass=0;
sub getcpiD(@) {
	my @a = @_;
	if(@a<6){
		$pass++;
		return [1,0,$err,$special,$pass];
	}
	my @alleles = qw{A C G T};
	my @geno=split /\,/,$a[2];
	my ($pe);
	if (exists $MarkerAF{$a[0]}){
		$pe = $Markers{$a[0]}->[2];
		foreach my $allele (@alleles){
			unless (defined $MarkerAF{$a[0]}{$allele}){
				$MarkerAF{$a[0]}{$allele} = 0.01;
			}
		}
	}else{
		$MarkerAF{$a[0]}{$geno[0]} = 0.55;
		for my $i (1..scalar @geno - 1){
			$MarkerAF{$a[0]}{$geno[$i]} = 0.45;
		}
		$pe = 0.125;
	}
	my @dad=split /\;/,$a[5];
	my @ch=split /\;/,$a[6];
	my @genodad=split /\//,$dad[0];
	my @genoch=split /\//,$ch[0];
	if ($genoch[0] ne $genodad[0] && $genoch[0] ne $genodad[1] && $genoch[1] ne $genodad[0] && $genoch[1] ne $genodad[1]){
		$err++;
		return [1e-4,$pe,$err,$special,$pass];
	} elsif ($genoch[0] eq $genoch[1]){
		if ($genodad[0] eq $genodad[1]){
			return [1/$MarkerAF{$a[0]}{$genoch[0]},$pe,$err,$special,$pass];
		}elsif ($genodad[0] ne $genodad[1]){
			return [0.5/$MarkerAF{$a[0]}{$genoch[0]},$pe,$err,$special,$pass];
		}
	} elsif ($genoch[0] ne $genoch[1]){
		my @child = sort @ch;
		my @Gdad = sort @dad;
		if ($genodad[0] eq $genodad[1]){
			return [0.5/$MarkerAF{$a[0]}{$genodad[0]},$pe,$err,$special,$pass];
		} elsif ($genodad[0] ne $genodad[1]){
			my @child = sort @ch;
			my @Gdad = sort @dad;
			if ($child[0] eq $Gdad[0] && $child[1] eq $Gdad[1]){
				my $pi = ($MarkerAF{$a[0]}{$child[0]} + $MarkerAF{$a[0]}{$child[1]}) / (4 * $MarkerAF{$a[0]}{$child[0]} * $MarkerAF{$a[0]}{$child[1]});
				return [$pi,$pe,$err,$special,$pass];
			}else{
				push @child,@Gdad;
				my %count;
				my @need = grep { ++$count{$_} > 1 } @child;
				return [0.25/$MarkerAF{$a[0]}{$need[0]},$pe,$err,$special,$pass];
			}
		}
	}
	return [999,$pe,-111,-222,-333];
}
sub getcpiT(@) {
	no warnings 'uninitialized';
	my @a = @_;
	if(@a<6){
		$pass++;
		return [1,0,$err,$special,$pass];
	}
	my @geno=split /\,/,$a[2];
	my $pe;
	my @alleles = qw{A C G T};
	if (exists $MarkerAF{$a[0]}) {
		$pe = $Markers{$a[0]}->[3];
		foreach my $allele (@alleles){
			unless (defined $MarkerAF{$a[0]}{$allele}){
				$MarkerAF{$a[0]}{$allele} = 0.01;
			}
		}
	} else {
		$MarkerAF{$a[0]}{$geno[0]} = 0.55;
		for my $i (1..scalar @geno - 1){
			$MarkerAF{$a[0]}{$geno[$i]} = 0.45;
		}
		$pe = 0.1875;
	}
	my @mm=split /\;/,$a[4];
	my @dad=split /\;/,$a[5];
	my @ch=split /\;/,$a[6];
	my @genomm=split /\//,$mm[0];
	my @genodad=split /\//,$dad[0];
	my @genoch=split /\//,$ch[0];
	if($genoch[0] ne $genomm[0] && $genoch[1] ne $genomm[1]){
		$pass++;
		return [1,0,$pe,$err,$special,$pass];
	} elsif ($genoch[0] eq $genoch[1]){
		if ($genoch[0] ne $genodad[0] && $genoch[0] ne $genodad[1]){
			$err++;
			return [1e-4,$pe,$err,$special,$pass];
		} elsif ($genodad[0] eq $genodad[1]){
			return [1/$MarkerAF{$a[0]}{$genoch[0]},$pe,$err,$special,$pass];
		} elsif ($genodad[0] ne $genodad[1]){
			return [0.5/$MarkerAF{$a[0]}{$genoch[0]},$pe,$err,$special,$pass];
		}
	} elsif ($genoch[0] ne $genoch[1]){
		my $fromF = $genoch[0] eq $genomm[0] ? $genoch[1] : $genoch[0];
		if ($fromF ne $genodad[0] && $fromF ne $genodad[1]){
			$err++;
			return [1e-4,$pe,$err,$special,$pass];
		} elsif ($genodad[0] eq $genodad[1]){
			return [1/$fromF,$pe,$err,$special,$pass];
		} elsif ($genodad[0] ne $genodad[1]){
			return [0.5/$fromF,$pe,$err,$special,$pass];
		}
	}
	return [999,$pe,-111,-222,-333];
}

1;
