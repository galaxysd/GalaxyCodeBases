#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <fq.lst> <fq.nfo> <fq.stat>\n";
	exit;
}

my ($fqlst,$fqnfo,$statout) = @ARGV;
my (%DATrbrf,%Librmx,%FQnfo,$min,$max);
open LST,'<',$fqlst or die "[x]Error opening $fqlst: $!\n";
while (<LST>) {
	chomp;
	my ($PESE,$sample,$lib,$FL,$ins,$ext,$path,@fqs)=split /\t/;
	($min,$max)=split ',',$ins;
	$Librmx{$sample}{$lib}{$FL}=[-1,$min,$max];	# readlen,minIS,maxIS
	$FQnfo{$sample}{$lib}{$FL}=[$ext,$path,@fqs];
	$DATrbrf{Lane}{$sample}{$lib}{$FL}=[0,0,0,0,0];	# rawReads,rawBP,filteredReads,filteredBP,FileCount
	$DATrbrf{Lib}{$sample}{$lib}=[0,0,0,0,0] unless exists $DATrbrf{Lib}{$sample}{$lib};	# may be faster ?
	$DATrbrf{Sample}{$sample}=[0,0,0,0,0] unless exists $DATrbrf{Sample}{$sample};
}

sub sumup ($$) {
	my ($arrayr,$hashr)=@_;
	$$arrayr[0] += $$hashr{InReads}||0;
	$$arrayr[1] += $$hashr{InBPs}||0;
	$$arrayr[2] += $$hashr{OutReads}||0;
	$$arrayr[3] += $$hashr{OutBP}||0;
	++$$arrayr[4];
}

my (@NFO,%NFO,$tmp);
open O,'>',$fqnfo or die "[x]Error opening $fqnfo: $!\n";
for my $sample (sort keys %FQnfo) {
	for my $lib (keys %{$FQnfo{$sample}}) {
		for my $FL (keys %{$FQnfo{$sample}{$lib}}) {
			my ($ext,$path,@fqs)=@{$FQnfo{$sample}{$lib}{$FL}};
			for (@fqs) {
				#print "[$_]\n";
				open NFO,'<',"$path$_.nfo" or (warn "[!]Error opening $path$_.nfo: $!\n" and next);
				@NFO=<NFO>;
				chomp @NFO;
				%NFO = map {split /\t/} grep /\w\t\d/,@NFO;
				#print "$_ => $NFO{$_}\n" for keys %NFO;
				close NFO;
				unless (grep /^All done/,@NFO) {
					warn "[!] $path$_$ext Not Finished.\n";
					next;
				}
				unless (exists $NFO{MaxReadLen} and exists $NFO{InReads} and exists $NFO{InBPs} and exists $NFO{OutReads} and exists $NFO{OutBP}) {
					print STDERR "[!] $path$_.nfo ERROR. ";
					if (exists $NFO{MaxReadLen}) { warn "MaxReadLen:[$NFO{MaxReadLen}] found and used.\n";}
					 else {warn "MaxReadLen NOT found, file skipped.\n";next;}
				}
				$Librmx{$sample}{$lib}{$FL}->[0]=$NFO{MaxReadLen} if $Librmx{$sample}{$lib}{$FL}->[0] < $NFO{MaxReadLen};
				&sumup($DATrbrf{Sample}{$sample},\%NFO);
				&sumup($DATrbrf{Lib}{$sample}{$lib},\%NFO);
				&sumup($DATrbrf{Lane}{$sample}{$lib}{$FL},\%NFO);
			}
			#next if $Librmx{$sample}{$lib}{$FL}->[0]==-1;
			$tmp=join(',',@{$Librmx{$sample}{$lib}{$FL}});
			if ($#fqs == 1) {	# PE
				print O join("\t",'PE',$sample,$lib,$FL,$tmp,@{$FQnfo{$sample}{$lib}{$FL}}),"\n";
			} else {
				print O join("\t",'SE',$sample,$lib,$FL,$tmp,@{$FQnfo{$sample}{$lib}{$FL}}),"\n";
			}
		}
	}
}
close O;

open O,'>',$statout or die "[x]Error opening $statout: $!\n";
print O "# SubItemOrder: rawReads,rawBP,filteredReads,filteredBP,FileCount\n";
my ($Rsample,$Rlib,$RFL);
for my $sample (sort keys %FQnfo) {
	$Rsample=join ',',@{$DATrbrf{Sample}{$sample}};
	for my $lib (keys %{$FQnfo{$sample}}) {
		$Rlib=join ',',@{$DATrbrf{Lib}{$sample}{$lib}};
		for my $FL (keys %{$FQnfo{$sample}{$lib}}) {
			$RFL=join ',',@{$DATrbrf{Lane}{$sample}{$lib}{$FL}};
			print O join("\t",$sample,$Rsample,$lib,$Rlib,$FL,$RFL),"\n";
		}
	}
}
close O;

__END__
对./1fqfilted/fqs.stat：
# SubItemOrder: rawReads,rawBP,filteredReads,filteredBP,FileCount
按个体统计：grep -v '#' ./1fqfilted/fqs.stat | awk '{print $1"\t"$2}' |sort|uniq
按文库统计：grep -v '#' ./1fqfilted/fqs.stat | awk '{print $3"\t"$4}' |sort|uniq
按lane统计：grep -v '#' ./1fqfilted/fqs.stat | awk '{print $5"\t"$6}'
