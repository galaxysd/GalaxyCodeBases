#!/usr/local/bin/perl -w
use strict;

my ($i, $rp, $count, $name, $oname
	);
 
$name = "";
$oname = "";
my (%seq) = ();
my ($mindiff) = 0;
my ($mod) = 0;
open (FH, $ARGV[0]) || die "cannot open"; 
while (<FH>)
{
   my ($line) = $_;
   chomp $line;
   if(substr($line,0,1) =~ />/)
   {
     $mod=0;
   }
   else
   {
     $mod=1;
   }
   last;
}
close FH;

if($mod==0)
{

  foreach $name (@ARGV) {
	open F,$name;
	my ($tag) = "";
	my ($tail) = "";
	my ($length) = 0;
	my ($pos) = -1;
	while (<F>) {
		chomp;
		my ($line) = $_;
		if ($line =~ /^\>(\S+)\s*(.*)/) {
			if($pos >= 0) {
				print "$tag $pos $length $name $tail\n";
				$length = 0;
			}
			$tag = $1;
			$tail = "";
			$tail = $2 if defined $2;
			$pos = tell(F) - length($line) - 1;
			next;
		}
		if ($line =~ /^[acgtACGTnN]/) {
			$length += length($line);
			next;
		}
	}
	close F;
	print "$tag $pos $length $name $tail\n";
  }
}
else
{
  foreach $name (@ARGV) {
        if($name =~ /(\S+)\.gz$/ or $name =~ /(\S+)\.Z$/) {
                open F,"gunzip -c $name |";
                $name = $1;
        } else {
                open F,$name;
        }
        my ($tag) = "";
        my ($tail) = "";
        my ($length) = 0;
        my ($pos) = -1;
        my ($tellF) = 0.;
        while (<F>) {
                $tellF += length($_);
                chomp;
                my ($line) = $_;
                if ($line =~ /^\@(\S+)\s*(.*)/) {
                        if($pos >= 0) {
                                printf "$tag %.0f $length $name $tail\n",$pos;
                                $length = 0;
                        }
                        $tag = $1;
                        $tail = "";
                        $tail = $2 if defined $2;
                        $pos = $tellF - length($line) - 1;
                        next;
                }
                if ($line =~ /^[acgtACGTnN]/) {
                        $length += length($line);
                        next;
                }
                if ($line =~ /^\+$tag/) {
                        $line = <F>;
                        $tellF += length($line);
                        next;
                }
        }
        close F;
        printf "$tag %.0f $length $name $tail\n",$pos;
  }
}
