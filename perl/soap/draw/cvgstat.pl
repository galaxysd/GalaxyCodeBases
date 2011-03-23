#!/bin/env perl
use strict;
use warnings;

unless (@ARGV>0){
	print "perl $0 depsingle\n";
	exit;
}
my ($depf)=@ARGV;
my $statf=$depf;
$statf =~ s/\.[^.\/\\]+$//;
$statf .= '.depstat';

my ($ChrID,%Stat,%StatA)=('Empty');
open DEPTH,'<',$depf or die "[x]Error opening $depf: $!\n";
open STAT,'>',$statf or die "[x]Error opening $statf: $!\n";
open STATA,'>',$statf.'all' or die "[x]Error opening ${statf}all: $!\n";
while (<DEPTH>) {
	chomp;
	if (/^>/) {
		if (exists $Stat{$ChrID}) {
			print STAT "[Chr$ChrID]\n";	# for Reading in R, SectionName must be string
			for my $v (sort {$a <=> $b} keys %{$Stat{$ChrID}}) {
				print STAT 'd',$v,"=",$Stat{$ChrID}{$v},"\n";	# for Reading in R, Key must be string
				#print $v,"=",$Stat{$ChrID}{$v},"\n";
			}
			print STAT "\n";
			print STDERR "\b\b\bdone.\n";
			delete $Stat{$ChrID};
		}
		print STDERR "$_ ...";
		s/^>//;
		$ChrID=$_;
		next;
	} else {
		for (split /\s+/) {
			++$Stat{$ChrID}{$_};
			++$StatA{$_};
		}
	}
}
close DEPTH;

		# for the last Chr.
		if (exists $Stat{$ChrID}) {
			print STAT "[$ChrID]\n";
			for my $v (sort {$a <=> $b} keys %{$Stat{$ChrID}}) {
				print STAT $v,"=",$Stat{$ChrID}{$v},"\n";
			}
			print STAT "\n";
			print STDERR "\b\b\bo\n";
			delete $Stat{$ChrID};
		}

close STAT;

for my $v (sort {$a <=> $b} keys %StatA) {
	print STATA $v,"\t",$StatA{$v},"\n";
}

close STATA;

__END__
> x=seq(0,30,1)
> plot(x,dpois(x,12),type='p')
> lines(x,dpois(x,12))
> lines(c(12,12),c(-1,1),col='blue')

HUMfcmRABDMAAPE
Lib='HUMfcmRAADDAAPE'
a=read.delim(paste('g',Lib,'.depstatall',sep=''),F)
a[1,2]=a[1,2]-239850802
EffLen=sum(as.numeric(a$V2))
maxX=20
maxY=0.15
x=seq(0,maxX+1,1)
c=sum(as.numeric(a$V1*a$V2))/EffLen
pdf(paste('g',Lib,'.pdf',sep=''),width=8,height=6,paper="special")
plot(a$V1,a$V2/EffLen,xlim=c(0,maxX),ylim=c(0,maxY),type='b',xlab='depth',ylab='density',main=paste('g',Lib,sep=''),
	sub=paste('depth=',round(c,4),', ','cvg=',round(100*(1-a[1,2]/EffLen),4),' % < ',round(100*(1-exp(-c)),4),' %',sep='')
	)
#lines(a$V1,a$V2/EffLen)
lines(x,dpois(x,c),col='blue')
lines(c(c,c),c(-1,1),col='red')
axis(1,at=c)
dev.off()

https://stat.ethz.ch/pipermail/r-help/2007-June/134110.html

Since in my problem the structure of the INI sections is almost static and
always present, I extended your example to create an in-memory list of
everything in the INI file with this function:

# Prototype of how to read INI files to process olfactometer data
# efg, 13 June 2007
# Thanks to Gabor Grothendieck for helpful suggestions in the R-Help
# mailing list on how to parse the INI file.
Parse.INI <- function(INI.filename)
{
  connection <- file(INI.filename)
  Lines  <- readLines(connection)
  close(connection)

  Lines <- chartr("[]", "==", Lines)  # change section headers

  connection <- textConnection(Lines)
  d <- read.table(connection, as.is = TRUE, sep = "=", fill = TRUE)
  close(connection)

  L <- d$V1 == ""                    # location of section breaks
  d <- subset(transform(d, V3 = V2[which(L)[cumsum(L)]])[1:3],
                           V1 != "")

  ToParse  <- paste("INI.list$", d$V3, "$",  d$V1, " <- '",
                    d$V2, "'", sep="")

  INI.list <- list()
  eval(parse(text=ToParse))

  return(INI.list)
}


Here's an example of using the above function (I'll put the sample input
file below):

INI1 <- Parse.INI("sample.ini")

# Explore INI contents
summary(INI1)

INI1$SystemSetup$OlfactometerCode
INI1$DefaultLevels
unlist(INI1$DefaultLevels)
INI1$Map

INI1$Map$port1
as.integer( unlist( strsplit(INI1$Map$port1, ",") ) )

= = = = =
Sample output:

> INI1 <- Parse.INI("sample.ini")
>
> # Explore INI contents
> summary(INI1)
              Length Class  Mode
SystemSetup   1      -none- list
Files         8      -none- list
DefaultLevels 4      -none- list
OdorNames     2      -none- list
Map           3      -none- list
>
> INI1$SystemSetup$OlfactometerCode
[1] "3"
> INI1$DefaultLevels
$FC00
[1] "50"

$FC01
[1] "100"

$FC02
[1] "50"

$FC10
[1] "50"

> unlist(INI1$DefaultLevels)
 FC00  FC01  FC02  FC10
 "50" "100"  "50"  "50"
> INI1$Map
$port0
[1] "0,0,0,0,0,0,0,0,0,0,0,0"

$port1
[1] "0,0,0,0,0,0,0,0,0,0,0,0"

$port2
[1] "0,0,0,0,0,0,0,0,0,0,0,0"

>
> INI1$Map$port1
[1] "0,0,0,0,0,0,0,0,0,0,0,0"
> as.integer( unlist( strsplit(INI1$Map$port1, ",") ) )
 [1] 0 0 0 0 0 0 0 0 0 0 0 0

= = = = =
Sample input file, sample.ini:

[SystemSetup]
OlfactometerCode=3
[Files]
prelog0=Part0.txt
date0=2:06:27.461 PM 6/9/2007
note0=group1-1
name0=group1
prelog1=Part1.txt
date1=2:09:16.809 PM 6/9/2007
note1=group1-1
name1=group1-1
[DefaultLevels]
FC00=50
FC01=100
FC02=50
FC10=50
[OdorNames]
port0=None
port1=None
[Map]
port0=0,0,0,0,0,0,0,0,0,0,0,0
port1=0,0,0,0,0,0,0,0,0,0,0,0
port2=0,0,0,0,0,0,0,0,0,0,0,0

plot(as.matrix(INI1$Chr11),type='l')
