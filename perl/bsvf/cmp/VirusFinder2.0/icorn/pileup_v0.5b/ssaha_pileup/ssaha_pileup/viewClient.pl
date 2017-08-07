#!/usr/bin/perl
use strict;
use warnings;
use Socket;
use Term::ANSIColor qw(:constants);
# colours available ####################################
# BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE
########################################################
$Term::ANSIColor::AUTORESET = 1;

my %params;
my $server;
my $recvline;
my $numSNP;
my $numBlock;
my @SNPs;
my @blocks;
my $block;
my $SNP;
my $i;

$SIG{'INT'}  = 'sigint_handler';
$SIG{'PIPE'} = 'sigpipe_handler';
$SIG{'TERM'} = 'sigterm_handler';

default_params();
input_params(@ARGV);

socket(SH,AF_INET,SOCK_STREAM,getprotobyname('tcp')) or die "Socket err: $!\n";
$server = sockaddr_in ($params{PORT},inet_aton($params{SERVER}));
connect (SH, $server) or die "Can't connect: $!\n";

defined(send(SH, param_line(),0)) or die "Can't send : $!\n";
defined ($recvline=<SH>) or die "Can't recieve : $!\n";

system("clear");
print BOLD "YABA: Yet Another Browser Application.\n\n";

while(defined ($recvline=<SH>)) {
    defined (send(SH, $recvline,0)) or die "Can't send : $!\n";
    if($params{VIEW_MOD} == 2) {
	if($recvline =~ m$<[acgtACGTnN-]+>$) {
	    @SNPs = $recvline =~ m$<([acgtACGTnN-]+)>$g;
	    @blocks = split /<[acgtACGTnN-]+>/, $recvline;
	    $numSNP = scalar(@SNPs);
	    $numBlock = scalar(@blocks);
	    die "Corrupt output\n" if($numBlock != $numSNP + 1);

	    for($i=0;$i<$numSNP;$i++) {
		print $blocks[$i];
		if($SNPs[$i] =~ m/[acgtACGT]/) {
		    print BOLD RED $SNPs[$i];
		}
		elsif($SNPs[$i] =~ m/-/) {
		    print BOLD GREEN $SNPs[$i];
		}
		elsif($SNPs[$i] =~ m/[nN]/) {
		    print BOLD BLUE $SNPs[$i];
		}
		else {
		    die "Error parsing output\n";
		}
	    }
	    print $blocks[$i];
	}
	else {
	    print $recvline;
	}
    }
    else {
	print $recvline;
    }
}

close(SH);

sub default_params {
  %params=(
  PORT           => '1172',
  SERVER         => 'pingu',
  CONTIG         => 1,
  LOCUS          => 1,
  QUAL           => 23,
  VIEW_MOD       => 1,
  SSAHA2_VERSION => '1.0.1',
  GAPS           => 1);
}

sub input_params {
  print_help() if(scalar(@_) < 2);
  for(my $i=0; $i<(@_); $i++) {
    $params{PORT}       = $_[$i+1] if($_[$i] eq "-port");
    $params{SERVER}     = $_[$i+1] if($_[$i] eq "-server");
    $params{CONTIG} = $_[$i+1] if($_[$i] eq "-contig");
    $params{LOCUS}  = $_[$i+1] if($_[$i] eq "-pos");
    $params{QUAL} = $_[$i+1] if($_[$i] eq "-qual");
    $params{VIEW_MOD}  = $_[$i+1] if($_[$i] eq "-view");
    $params{SSAHA2_VERSION} = $_[$i+1] if($_[$i] eq "-version");
    $params{GAPS} = $_[$i+1] if($_[$i] eq "-gaps");
  }
  $params{CONTIG} -= 1;
  if($params{CONTIG} < 0) {
      die "Warning '-contig' now starts counting at 1.\n";
  }
}

sub param_line {
  "COMMAND_INPUT $params{QUAL} $params{LOCUS} $params{CONTIG} $params{VIEW_MOD} $params{SSAHA2_VERSION} $params{GAPS}\n"
}

sub print_help {
    printf("Usage ./viewClient.pl <-server $params{SERVER}> <-port $params{PORT}> <-contig $params{CONTIG}> <-pos $params{LOCUS}> <-view $params{VIEW_MOD}> <-qual $params{23}> <-version 1.0.2> <-gaps 1>\n");
    exit(0);
}

sub sigint_handler {
  # NB sig_pipe is sent to server automatically by system before TCP
  # socket is broken.
    print "\n\n***************************************",
          "\n\n         Client was interrupted        ",
          "\n\n***************************************\n";
    die $!;
}

sub sigterm_handler {
  # NB sig_pipe is sent to server automatically by system before TCP
  # socket is broken.
    print "\n\n***************************************",
          "\n\n           Client was killed           ",
          "\n\n***************************************\n";
    die $!;
}

sub sigpipe_handler {
    print "\n\n***************************************",
          "\n\n               Broken Pipe             ",
          "\n\n***************************************\n";
    die $!;
}
