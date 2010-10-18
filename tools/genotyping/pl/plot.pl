#!/usr/bin/perl -w
use strict;
use SVG;
unless(@ARGV)
{
	die "perl $0 <chrlen file><input file><out file><chr><width><proportion>\n";
}
my $chrlen;
open IN,"$ARGV[0]" or die "The chrlen file can't open:$!\n";
while(<IN>)
{
	chomp;
	my @temp=split;
	if($temp[0]=~/$ARGV[3]/i)
	{
		$chrlen=$temp[1];
	}
}
close IN;
open CIN,"$ARGV[1]" or die "The input file can't open:$!\n";

my $temp='#';
while ($temp =~ /^#/) {
	$temp=<CIN>;
}
chomp $temp;
#chomp(my $temp=<CIN>);

my @name=split /\s+/,$temp;
#shift @name;
my $samnum=@name;
#print "@name\n$samnum\n";
#ast;
my $side=100;
my $width=$chrlen/$ARGV[5]+$side*3;
my $height=$samnum*$ARGV[4]+200;
my $svg=SVG->new('width',$width,'height',$height);
$svg->rect('x',0,'y',0,'width',$width,'height',$height,'stroke','white','fill','white');

my %ReFormatGT=(0=>1, 1=>2, 0.5=>3,'NA'=>0, 'N'=>-1);
my %FormatGT=(1=>0, 2=>1, 3=>0.5, 0=>'NA', -1=>'N');
#my @color=("#2A2F7F","#68BE4A","#8C63A4","#B6292B","#FCDAD5");
my %Color=(
	1  => '#CC0000',
	3  => '#1100CC',
	2  => '#008800',
	0  => '#CCCCCC',
	-1 => '#000022',
);
#my @label=("9311","pair64","no_decide","het","no_cov");
my @Label=( [1,'9311'],[2,'PA64'],[3,'Hete'],[0,'N/A'],[-1,'N Zone'] );
pop @Label;
my $i=0;
foreach(@Label) {
	my ($id,$label)=@$_;
	$svg->line('x1',$side,'y1',$side/2+$ARGV[4]*$i*2,'x2',$side+20,'y2',$side/2+$ARGV[4]*$i*2,'stroke',$Color{$id},'stroke-width',$ARGV[4]);
	$svg->text('x',$side+25,'y',$side/2+$ARGV[4]*$i*2+3,'font-family','Arial','stroke-width',0.3,'font-size',12,'-cdata',$label);
	++$i;
}
#plot left
$svg->line('x1',$side,'y1',$side,'x2',$side,'y2',$height-100,'stroke','black','stroke-width',1);
foreach(1..int(($height-200)/($ARGV[4]*20)))
{
	$svg->line('x1',$side,'y1',$height-100-20*$_*$ARGV[4],'x2',$side-5,'y2',$height-100-20*$_*$ARGV[4],'stroke','black','stroke-width',0.7);
        $svg->text('x',$side-30,'y',$height-100-20*$_*$ARGV[4]+5,'font-family','Arial','stroke-width',0.3,'font-size',14,'-cdata',20*$_);
}
my $sit=($height-200)/2+100+30;
$svg->text('x',50,'y',$sit,'font-family','Arial','stroke-width',0.5,'font-size',16,'-cdata','Individuals','transform',"rotate(270,50,$sit)");
#plot down
$svg->line('x1',$side,'y1',$height-100,'x2',$width-$side*2,'y2',$height-100,'stroke','black','stroke-width',1);
my $par=int($chrlen/5000000);
my $len=int($chrlen/($ARGV[5]*$par));
foreach(0..$par)
{
	$svg->line('x1',$side+$_*$len,'y1',$height-100,'x2',$side+$_*$len,'y2',$height-100+5,'stroke','black','stroke-width',0.7);
	$svg->text('x',$side+$_*$len-5,'y',$height-100+20,'font-family','Arial','stroke-width',0.3,'font-size',14,'-cdata',5*$_);
}
$svg->text('x',($width-$side*3)/2+$side,'y',$height-100+50,'font-family','Arial','stroke-width',0.5,'font-size',16,'-cdata','Position(Mb)');
#plot right
$svg->line('x1',$width-200,'y1',$side,'x2',$width-200,'y2',$height-100,'stroke','black','stroke-width',1);
my $m=0;
foreach my $sam(1..$samnum)
{
	$svg->line('x1',$width-2*$side,'y1',$height-$side-$ARGV[4]*$sam+$ARGV[4]/2,'x2',$width-2*$side+$ARGV[4]*4,'y2',$height-$side-$ARGV[4]*$sam-$ARGV[4]*2,'stroke','black','stroke-width',0.7);
	if($sam%2==1)
	{
	 $svg->line('x1',$width-2*$side+$ARGV[4]*4,'y1',$height-$side-$ARGV[4]*$sam-$ARGV[4]*2,'x2',$width-2*$side+$ARGV[4]*14,'y2',$height-$side-$ARGV[4]*$sam-$ARGV[4]*2,'stroke','black','stroke-width',0.7);
	$svg->text('x',$width-2*$side+$ARGV[4]*15,'y',$height-$side-$ARGV[4]*$sam-$ARGV[4]*2+$ARGV[4]/2+2,'font-family','Arial','stroke-width',0.5,'font-size',10,'-cdata',$name[$sam-1]);
	}
	else
	{
		$svg->line('x1',$width-2*$side+$ARGV[4]*4,'y1',$height-$side-$ARGV[4]*$sam-$ARGV[4]*2,'x2',$width-2*$side+$ARGV[4]*6,'y2',$height-$side-$ARGV[4]*$sam-$ARGV[4]*2,'stroke','black','stroke-width',0.7);
        $svg->text('x',$width-2*$side+$ARGV[4]*6,'y',$height-$side-$ARGV[4]*$sam-$ARGV[4]*2+$ARGV[4]/2+2,'font-family','Arial','stroke-width',0.5,'font-size',10,'-cdata',$name[$sam-1]);
	}
}
#plot up
$svg->line('x1',$side,'y1',$side,'x2',$width-200,'y2',$side,'stroke','black','stroke-width',1);
#plot data
#chomp(my $line=<CIN>);
#my @last=split /\s+/,$line;
#my @color=("#2A2F7F","#00A06B");
my $mark=0;
my @new;
my @last;
my @start;
while(<CIN>)
{
	next if /^#/;
	chomp;
	my @now=split;
	if($mark==0)
	{
		@last=@now;
		foreach my $sam(1..$samnum)
        	{
			$start[$sam]=1;
		}
		$start[0]="seat";
	}
	foreach my $sam(1..$samnum)
	{
		if($last[$sam] ne $now[$sam])
		{
			my $n=$last[$sam];
			my $x1=$side+$start[$sam]/$ARGV[5];
			my $y1=$height-$side-$ARGV[4]*$sam+$ARGV[4]/2;
			my $x2=$now[0]/$ARGV[5]+$side;
			my $y2=$y1;
			my $colr=$Color{$ReFormatGT{$n}};
#		print "$x1\t$y1\n";
=pod
			if($n eq 0.5)
			{
				$colr=$color[3];
			}
			elsif($n eq "NA")
			{
				$colr=$color[4];
			}
			else
			{
				$colr=$color[$n];
			}
=cut
			$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke',$colr,'stroke-width',$ARGV[4]);# unless $n =~ /^N/;
			$last[$sam]=$now[$sam];
			$start[$sam]=$now[0];
		}
	}
#	print "@last\n";
#	last;
	$mark=1;
	@new=@now;
}
foreach my $sam(1..$samnum)
{
	if($last[$sam] eq $new[$sam])
        {
        	my $n=$last[$sam];
                my $x1=$side+$start[$sam]/$ARGV[5];
                my $y1=$height-$side-$ARGV[4]*$sam+$ARGV[4]/2;
                my $x2=$new[0]/$ARGV[5]+$side;
                my $y2=$y1;
#               print "$x1\t$y1\n";
		my $colr=$Color{$ReFormatGT{$n}};
                $svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke',$colr,'stroke-width',$ARGV[4]);# unless $n =~ /^N/;
	}
}

open(OUT,">$ARGV[2]")or die;
print OUT $svg->xmlify();
close OUT;

__END__
/nas/RD_09C/resequencing/user/zhangxm/soft/zijiao/new_selfing.pl
perl /nas/RD_09C/resequencing/user/zhangxm/soft/zijiao/new_selfing.pl /nas/RD_09C/resequencing/Genome/Rice/Genome_9311/chr.nfo /panfs/POPULATION/PROJECT/Rice_RIL/chongqiong_RIL/data/Chr01.out.new /panfs/POPULATION/PROJECT/Rice_RIL/chongqiong_RIL/pic/Chr01.out.new.svg Chr01 5 40000
cat ../9311/chrorder |while read a; do ./plot.pl ../9311/chr.nfo ril.$a.rgt ril.$a.r.svg $a 6 40000; done
find . -name '*.?1.svg'|perl -lane '$a=$_;s/\.svg$/\.png/;system "convert $a $_"'
