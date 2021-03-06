#!/usr/bin/env perl
use strict;
use warnings;
use threads;

die "Usage: $0 <TS file> <Out>.mkv [Split=15]\n" if @ARGV <2;

my ($fin,$fout,$nsplit)=@ARGV;
$nsplit = 15 unless defined $nsplit;

my $FFMPEG_BIN = 'ffmpeg-10bit';
$FFMPEG_BIN = './ffmpeg';

my $TcpPort0 = 12000;
#my @jobs;
my ($ffarmIDs,$PreOut,@ffarms,$foutIn);

my $fps1001o = 60000;
my $fps1001n = $fps1001o/$nsplit;

for my $j (1 .. $nsplit) {
	my ($Port1,$Port2,$OutTag) = ($TcpPort0+$j,$TcpPort0+$nsplit+$j,"O$j");
	#push @jobs,[$Port1,$Port2,$OutTag];
	$ffarmIDs .= "[$OutTag]";
	$PreOut .= " -map \'[$OutTag]\' -f yuv4mpegpipe tcp://localhost:$Port1?listen";
	push @ffarms,"$FFMPEG_BIN -i tcp://localhost:$Port1 -vf owdenoise=8:0.309:0,fps=\'$fps1001n/1001\' -f yuv4mpegpipe tcp://localhost:$Port2?listen >$fout.p$j.log 2>&1";
	$foutIn .= " -i tcp://localhost:$Port2";
}

sub doffarm() {
	#system("AV_LOG_FORCE_COLOR=1 $_");
	print "$_ &\n";
}

my $cmd = "$FFMPEG_BIN -threads 16 -i $fin -lavfi \"bwdif=send_field:tff,mcdeint=fast:tff:10,scale=w=iw:h=ih/2,drawtext=\'fontfile=arial.ttf:text=%{n}_%{pict_type}_%{pts}:fontcolor=Aqua:fontsize=32\',select=n=4:e=\'mod(n,$nsplit)+1\'$ffarmIDs\"$PreOut >$fout.pre.log 2>&1";
print "$cmd &\nsleep 3\n";

my %thread;
foreach (@ffarms) {
	$thread{$_} = threads->new(\&doffarm, $_);
}

sleep 1;

#$cmd = "$FFMPEG_BIN $foutIn -lavfi [0:v][1:v][2:v][3:v]interleave=4,fps=\'60000/1001\' -pix_fmt yuv420p10le -vcodec libx264 -preset placebo -x264opts \'crf=23:vbv-maxrate=15000:vbv-bufsize=15000:threads=32:colormatrix=bt709:colorprim=bt709:transfer=bt709:psnr=1:ssim=1\' -vframes 1600 -f matroska -y $fout.mkv >$fout.cmp.log 2>&1";
$cmd = "$FFMPEG_BIN $foutIn -lavfi [0:v][1:v][2:v][3:v]interleave=4,fps=\'60000/1001\' -pix_fmt yuv420p10le -strict -1 -y -f yuv4mpegpipe - | mpv -";
#system("AV_LOG_FORCE_COLOR=1 $cmd");
print "sleep 2\n$cmd\n";

my %count;
foreach (sort keys %thread) {
	$count{$_} = $thread{$_}->join;
}

#print "[m]done !\n";

