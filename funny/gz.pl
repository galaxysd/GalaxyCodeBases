#!/usr/bin/perl -w
use strict;
use warnings;
#use IO::Handle;
use IPC::Open2;
use POSIX ":sys_wait_h";

=pod
http://tools.ietf.org/html/rfc1951#section-3.2
5 Bits: HLIT, # of Literal/Length codes - 257 (257 - 286)
5 Bits: HDIST, # of Distance codes - 1        (1 - 32)
4 Bits: HCLEN, # of Code Length codes - 4     (4 - 19)
...
e3 02 00 for "\n";

? bytes  compressed data
4 bytes  crc32
4 bytes  uncompressed input size modulo 2^32
=cut
my $minGZ=10;	# (5+5+4)b + 8 byte
sub Compress($$) {
	my ($Str,$Bin)=@_;
	#$Str=v0 if $Str eq '';
	if (length $Str <= $minGZ) {
		$_[1] = $Str;
		return 0;
	}
	my ($pid,$chld_out, $chld_in);
	my $sleep_count = 0;
	do {
		$pid = open2($chld_out, $chld_in, 'gzip','-6nc');
		unless (defined $pid) {
			warn "[!]Cannot fork gzip: $!";
			die "[x]Bailing out" if $sleep_count++ > 6;
			sleep 10;
		}
	} until defined $pid;
	print $chld_in $Str;
	close $chld_in;
	local $/;	# enable "slurp" mode
	$Bin=<$chld_out>;
	close $chld_out;
	$_[1] = substr $Bin,10;
	#waitpid $pid, 0;
	return $pid;
}
my $GZHeader = v31.139.8.0.0.0.0.0.0.3;	# -9 => 2.3; -2~8 => 0.3
#cat t.pl|while read a;do echo "$a"|gzip -6nc|hexdump -C;done
#00000000  1f 8b 08 00 00 00 00 00  00 03 e3 02 00 93 06 d7
sub deCompress($$) {
	my ($Bin,$Str)=@_;
	if (length $Bin <= $minGZ) {
		$_[1] = $Bin;
		return 0;
	}
	my ($pid,$chld_out, $chld_in);
	my $sleep_count = 0;
	do {
		$pid = open2($chld_out, $chld_in, 'gzip','-dc');
		unless (defined $pid) {
			warn "[!]Cannot fork gzip: $!";
			die "[x]Bailing out" if $sleep_count++ > 6;
			sleep 10;
		}
	} until defined $pid;
	print $chld_in $GZHeader,$Bin;
	close $chld_in;
	local $/;	# enable "slurp" mode
	$Str=<$chld_out>;
	close $chld_out;
	$_[1] = $Str;
	#waitpid $pid, 0;
	return $pid;
}
sub waitGout() {
	my $kid;
	do {
		$kid = waitpid(-1, WNOHANG);
		#warn $kid;
	} while $kid > 0;
}
=pod
$ man waitpid

pid_t waitpid(pid_t pid, int *status, int options);
       The value of pid can be:

       < -1   meaning wait for any child process whose  process  group  ID  is
              equal to the absolute value of pid.

       -1     meaning wait for any child process.

       0      meaning  wait  for  any  child process whose process group ID is
              equal to that of the calling process.

       > 0    meaning wait for the child whose process  ID  is  equal  to  the
              value of pid.

       The  value  of  options  is an OR of zero or more of the following con©\
       stants:

       WNOHANG     return immediately if no child has exited.

       WUNTRACED   also return if a child has  stopped  (but  not  traced  via
                   ptrace(2)).   Status for traced children which have stopped
                   is provided even if this option is not specified.

       WCONTINUED (since Linux 2.6.10)
                   also return if a stopped child has been resumed by delivery
                   of SIGCONT.
=cut

my ($ttt,$aaa);
my($chld_out, $chld_in);
#$|=1;
my ($a,$b,$c,$sumA,$sumB);
open IN,'<',$ARGV[0] or die "$!";
#my $i=1;
while (<IN>) {
	chomp;
	$a=length($_);
	Compress($_,$ttt);
	$b=length $ttt;
	deCompress($ttt,$aaa);
	$c=length $aaa;
	chomp $aaa;
	print int(1000*$b/$a)/10," %\t$a\t$b\t$c,[$aaa]\n" if $a > 0;
	$sumA += $a; $sumB += $b;
	#waitpid(-1,0);
	#waitGout() unless $i % 100;
	#++$i;
	waitGout();
	#waitpid(-1,0);
	#waitpid(-1,0);
}
close IN;
print '-'x75,"\n",int(10000*$sumB/$sumA)/100," %\t$sumA\t$sumB\n";

__END__

