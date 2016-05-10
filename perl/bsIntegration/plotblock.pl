#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <block file>\n" if @ARGV<1;

my ($in) = @ARGV;
mkdir '/tmp/plots';
mkdir './plots';

$in = './example18/simVir4_grep/blocks.txt.gz';
my $chrOnly = 'chr18';
my $viruschr = 'gi|59585|emb|X04615.1|';

my $IMG_Width = 1200;
my $IMG_Height = 900;
my $IMG_Boder = 50;

my $Paint_Width = $IMG_Width - 2*$IMG_Boder;
my $Paint_Height = $IMG_Height - 2*$IMG_Boder;
my $Stroke_Width = $Paint_Width/240;

sub openFH($) {
	my $inf = $_[0];
	my $FH;
	if ($inf =~ /\.gz$/i) {
		open $FH,'-|',"gzip -dc $inf" or die "Error opening $inf: $!\n";
	} elsif ($inf =~ /\.xz$/i) {
		open $FH,'-|',"xz -dc $inf" or die "Error opening $inf: $!\n";
	} elsif ($inf =~ /\.bz2$/i) {
		open $FH,'-|',"bzip2 -dc $inf" or die "Error opening $inf: $!\n";
	} else {
		open $FH,'<',$inf or die "Error opening $inf: $!\n";
	}
	return $FH;
}

my $fhin = openFH($in);
my @cache=();
my ($id,$chr,$range,$cnt,$vchr,$vrange,$vcnt);

sub getmaxkey($) {
	my $inh = $_[0];
	my @keys = sort { $inh->{$b} <=> $inh->{$a} } keys %{$inh};
	return $keys[0];
}

sub CIAGRdepth($$$$) {
	my ($CIAGR,$leftPos,$DatHash,$SimInfo) = @_;
	my $curPos = $leftPos;
	my (%LR,$ret);
	while ($CIAGR =~ /(\d+)(\D)/g) {
		if ($2 eq 'M') {
			my $end = $curPos + $1 -1;
			for ($curPos .. $end) {
				++$DatHash->{$_};
				if (defined $SimInfo) {
					if ($_ >= $$SimInfo[0] and $_ < $$SimInfo[1]) {
						++$LR{'L'};
					} elsif ($_ >= $$SimInfo[1] and $_ <= $$SimInfo[2]) {
						++$LR{'R'};
					}
				}
			}
			$curPos = $end +1;
		} elsif ($2 eq 'S') {
			$curPos += $1;
		} elsif ($2 eq 'D') {
			$curPos += $1;
		} else {
			die "[x]CIAGR[$1 $2] not allowed.\n" if $2 ne 'I';
		}
	}
	if (%LR) {
		$ret = getmaxkey(\%LR);
	}
	return $ret;
}

sub getdepth($) {
	my ($readsRef)=@_;
	my (%HumDepth,%VirusDepth,%sEdges,%pEdges,%cpItemCache);
	%pEdges = ('LL'=>0, 'RR'=>0);
	my (%HumDep,%VirDep,%HumCut,%HumLeft,%HumRight,%VirStrand,%VirLeft,%VirRight);
	my $doflag = 0;
	for (@$readsRef) {
		my ($qname,$flag,$rname,$pos,$mapq,$CIAGR,undef,undef,undef,undef,undef,@opts) = split /\t/;
		$qname =~ /_Ref_(\d+)_(\d+)_(\d+)_Vir_([+-])_(\d+)_(\d+)_/ or die;
		++$HumLeft{$1};
		++$HumCut{$2};
		++$HumRight{$3};
		++$VirStrand{$4};
		++$VirLeft{$5};
		++$VirRight{$6};
		my @XAs = grep(/^[SX]A:Z:/,@opts);
		#print "+++ @XAs\n";
		my @t = grep(/\b[fr](\Q$chrOnly\E|\Q$viruschr\E),/,@XAs);
		if (scalar @t>=0) {
			$doflag = 1;
			#print "--- $doflag, @t\n";
			$cpItemCache{"$qname\t$rname\t$pos"} = $CIAGR;
		}
	}
	return unless $doflag; # in fact, all items have $doflag==1;
	my ($humLeft,$humCut,$humRight,$virStrand,$virLeft,$virRight);
	$humLeft = getmaxkey(\%HumLeft);
	$humCut = getmaxkey(\%HumCut);
	$humRight = getmaxkey(\%HumRight);
	$virStrand = getmaxkey(\%VirStrand);
	$virLeft = getmaxkey(\%VirLeft);
	$virRight = getmaxkey(\%VirRight);
	my $SimInfo = [$humLeft,$humCut,$humRight,$virStrand,$virLeft,$virRight];
	print ">>> $humLeft,$humCut,$humRight,$virStrand,$virLeft,$virRight $doflag:\n";
	for (@$readsRef) {
		my ($qname,$flag,$rname,$pos,$mapq,$CIAGR,$MRNM,$MPOS,undef,undef,undef,@opts) = split /\t/;
		my ($LR1,$LR2,@itemsEdges,@itempEdges);
		$MRNM = $rname if $MRNM eq '=';
		if ($rname eq $chrOnly) {
			$LR1 = CIAGRdepth($CIAGR,$pos,\%HumDepth,$SimInfo);
			next unless defined $LR1;
		} else {
			$LR1 = 'V';
			CIAGRdepth($CIAGR,$pos,\%VirusDepth,undef);
		}
		if (exists $cpItemCache{"$qname\t$MRNM\t$MPOS"}) {
			my $MCIAGR = $cpItemCache{"$qname\t$MRNM\t$MPOS"};
			if ($MRNM eq $chrOnly) {
				$LR2 = CIAGRdepth($MCIAGR,$MPOS,\%HumDepth,$SimInfo);
			} else {
				$LR2 = 'V';
				CIAGRdepth($MCIAGR,$MPOS,\%VirusDepth,undef);
			}
			if (defined $LR2) {
				my @t = sort ($LR1,$LR2);
				++$pEdges{join('',@t)};
			}
		}
		my @XAHit;
		for (@opts) {
			if (/^[SX]A:Z:([\w\,\;\|\.\+\-]+)$/) {
				push @XAHit,split(/;/,$1);
			}
		}
		if (@XAHit) {
			#print "$qname\n";
			for (@XAHit) {
				my ($XAchr,$XApos,$XAcigar,$XAnm) = split /,/;
				if ($XAchr =~ /^[fr](\Q$chrOnly\E|\Q$viruschr\E)$/) {
					#print "$1\t$XAchr,$XApos,$XAcigar,$XAnm\n";
					my $LR3;
					if ($1 eq $chrOnly) {
						$LR3 = CIAGRdepth($XAcigar,abs($XApos),\%HumDepth,$SimInfo) if $LR1 eq 'V';
					} elsif ($1 eq $viruschr and $LR1 ne 'V') {
						$LR3 = 'V';
						CIAGRdepth($XAcigar,abs($XApos),\%VirusDepth,undef);
					}
					if (defined $LR3) {
						my @t = sort ($LR1,$LR3);
						push @itemsEdges,join('',@t);
					}
				}
			}
			if (@itemsEdges) {
				my $weight = 1 / @itemsEdges;
				$sEdges{$_} += $weight for @itemsEdges;
			}
		}
	}
	if (exists $sEdges{'LV'} and exists $sEdges{'RV'}) {
		ddx [\%sEdges,\%pEdges];
		return (\%HumDepth,\%VirusDepth,\%sEdges,\%pEdges,$SimInfo);
	} else {
		warn ".\n";
		return;
	}
}

sub doplot($$$$$$) {
	my ($fh,$HumDepth,$VirusDepth,$sEdges,$pEdges,$SimInfo) = @_;
	my ($MaxDepth,$HumMin,$HumMax,$VirMin,$VirMax)=(0,@$SimInfo[0,2,4,5]);
	my @HumCovered = sort {$a <=> $b} keys %{$HumDepth};
	my @VirCovered = sort {$a <=> $b} keys %{$VirusDepth};
	for (@HumCovered) {
		$MaxDepth = $HumDepth->{$_} if $MaxDepth < $HumDepth->{$_};
	}
	for (@VirCovered) {
		$MaxDepth = $VirusDepth->{$_} if $MaxDepth < $VirusDepth->{$_};
	}
	$HumMin = $HumCovered[0] if $HumMin > $HumCovered[0];
	$HumMax = $HumCovered[-1] if $HumMax < $HumCovered[-1];
	$VirMin = $VirCovered[0] if $VirMin > $VirCovered[0];
	$VirMax = $VirCovered[-1] if $VirMax < $VirCovered[-1];
	my $HumRange = $HumMax - $HumMin +1;
	my $VirRange = $VirMax - $VirMin +1;
	my $RatioH = $Paint_Width / $HumRange;
	my $RatioV = $Paint_Width / $VirRange;
	my $XRatio = ($RatioH < $RatioV)?$RatioH:$RatioV;
	my @Cuts;
	for (@$SimInfo[0,1,2]) {
		push @Cuts,($_ - $HumMin);
	}
	for (@$SimInfo[4,5]) {
		push @Cuts,($_ - $VirMin);
	}
	print "Hum:$HumCovered[0] - $HumCovered[-1]\nVir:$VirCovered[0] - $VirCovered[-1]\n@$SimInfo\t@Cuts\n";
	$_ *= $XRatio for @Cuts;
	my @BarHeight = ($Paint_Height*3/8, $Paint_Height*5/8, $Paint_Height*7/16);
	my $YRatio = ($Paint_Height*3/8)/$MaxDepth;
	my $HumLeft = ($Paint_Width - $HumRange*$XRatio)/2;
	my $VirLeft = ($Paint_Width - $VirRange*$XRatio)/2;
	print $fh '<line stroke="gray" stroke-width="',$Stroke_Width,'" x1="',$VirLeft+$IMG_Boder,'" y1="',$IMG_Boder+$BarHeight[0]+$Stroke_Width/2,'" x2="',$VirLeft+$IMG_Boder+$VirRange*$XRatio,'" y2="',$IMG_Boder+$BarHeight[0]+$Stroke_Width/2,'"/>',"\n";
	print $fh '<line stroke="gray" stroke-width="',$Stroke_Width,'" x1="',$HumLeft+$IMG_Boder,'" y1="',$IMG_Boder+$BarHeight[1]-$Stroke_Width/2,'" x2="',$HumLeft+$IMG_Boder+$HumRange*$XRatio,'" y2="',$IMG_Boder+$BarHeight[1]-$Stroke_Width/2,'"/>',"\n";

	print $fh '<line stroke="firebrick" stroke-width="',$Stroke_Width,'" x1="',$Cuts[3]+$VirLeft+$IMG_Boder,'" y1="',$IMG_Boder+$BarHeight[0]+$Stroke_Width,'" x2="',$Cuts[4]+$VirLeft+$IMG_Boder,'" y2="',$IMG_Boder+$BarHeight[0]+$Stroke_Width,'"/>',"\n";
	print $fh '<line stroke="forestgreen" stroke-width="',$Stroke_Width,'" x1="',$Cuts[0]+$HumLeft+$IMG_Boder,'" y1="',$IMG_Boder+$BarHeight[1]-$Stroke_Width,'" x2="',$Cuts[1]+$HumLeft+$IMG_Boder,'" y2="',$IMG_Boder+$BarHeight[1]-$Stroke_Width,'"/>',"\n";
	print $fh '<line stroke="mediumblue" stroke-width="',$Stroke_Width,'" x1="',$Cuts[1]+$HumLeft+$IMG_Boder,'" y1="',$IMG_Boder+$BarHeight[1]-$Stroke_Width,'" x2="',$Cuts[2]+$HumLeft+$IMG_Boder,'" y2="',$IMG_Boder+$BarHeight[1]-$Stroke_Width,'"/>',"\n";

	print $fh "<defs>\n";
	print $fh '<path id="pLV" d="M',$Cuts[0]+$HumLeft+$IMG_Boder+$Stroke_Width/2,',',$IMG_Boder+$BarHeight[1]-2*$Stroke_Width,' Q',$HumLeft+$IMG_Boder,',',$IMG_Boder+$BarHeight[2],' ',$Cuts[3]+$VirLeft+$IMG_Boder,',',$IMG_Boder+$BarHeight[0]+2*$Stroke_Width,'"
fill="none" stroke="limegreen" stroke-width="',$Stroke_Width,'" />',"\n";
	print $fh '<path id="sLV" d="M',$VirLeft+$IMG_Boder+($Cuts[4]+15*$Cuts[3])/16,',',$IMG_Boder+$BarHeight[0]+2*$Stroke_Width,' Q',$HumLeft+$IMG_Boder+$Cuts[1],',',$IMG_Boder+$BarHeight[2],' ',$Cuts[1]+$HumLeft+$IMG_Boder-$Stroke_Width/2,',',$IMG_Boder+$BarHeight[1]-2*$Stroke_Width,'"
fill="none" stroke="limegreen" stroke-width="',$Stroke_Width,'" />',"\n";

	print $fh '<path id="sRV" d="M',$Cuts[1]+$HumLeft+$IMG_Boder+$Stroke_Width/2,',',$IMG_Boder+$BarHeight[1]-2*$Stroke_Width,' Q',$HumLeft+$IMG_Boder+$Cuts[1],',',$IMG_Boder+$BarHeight[2],' ',$VirLeft+$IMG_Boder+$Cuts[3]+($Cuts[4]-$Cuts[3])*15/16,',',$IMG_Boder+$BarHeight[0]+2*$Stroke_Width,'"
fill="none" stroke="royalblue" stroke-width="',$Stroke_Width,'" />',"\n";
	print $fh '<path id="pRV" d="M',$Cuts[4]+$VirLeft+$IMG_Boder,',',$IMG_Boder+$BarHeight[0]+2*$Stroke_Width,' Q',$HumLeft+$IMG_Boder+$Cuts[2],',',$IMG_Boder+$BarHeight[2],' ',$Cuts[2]+$HumLeft+$IMG_Boder-$Stroke_Width/2,',',$IMG_Boder+$BarHeight[1]-2*$Stroke_Width,'"
fill="none" stroke="royalblue" stroke-width="',$Stroke_Width,'" />',"\n";
	print $fh "</defs>\n<use xlink:href=\"#pLV\"/>\n<use xlink:href=\"#sLV\"/>\n<use xlink:href=\"#pRV\"/>\n<use xlink:href=\"#sRV\"/>\n";
	print $fh '<text text-anchor="middle" dominant-baseline="text-before-edge" font-family="Verdana" font-size="20">',"\n";
	print $fh "<textPath startOffset=\"50%\" xlink:href=\"#pLV\">pLV: ${$pEdges}{'LV'}</textPath>\n";
	print $fh "<textPath startOffset=\"50%\" xlink:href=\"#sLV\">sLV: ${$sEdges}{'LV'}</textPath>\n";
	print $fh "<textPath startOffset=\"50%\" xlink:href=\"#pRV\">pRV: ${$pEdges}{'RV'}</textPath>\n";
	print $fh "<textPath startOffset=\"50%\" xlink:href=\"#sRV\">sRV: ${$sEdges}{'RV'}</textPath>\n";
	print $fh '</text>',"\n";
	print $fh '<text x="',0.5*($Cuts[0]+$Cuts[1])+$HumLeft+$IMG_Boder,'" y="',$IMG_Boder+$BarHeight[1]-2*$Stroke_Width,'" font-family="Verdana" font-size="20" text-anchor="middle" alignment-baseline="text-after-edge">pLL: ',${$pEdges}{'LL'},'</text>',"\n";
	print $fh '<text x="',0.5*($Cuts[1]+$Cuts[2])+$HumLeft+$IMG_Boder,'" y="',$IMG_Boder+$BarHeight[1]-2*$Stroke_Width,'" font-family="Verdana" font-size="20" text-anchor="middle" alignment-baseline="text-after-edge">pRR: ',${$pEdges}{'RR'},'</text>',"\n";
	print $fh '<text x="',0.5*($Cuts[3]+$Cuts[4])+$VirLeft+$IMG_Boder,'" y="',$IMG_Boder+$BarHeight[0]+$Stroke_Width,'" font-family="Verdana" font-size="20" text-anchor="middle" alignment-baseline="text-before-edge">pVV: ',${$pEdges}{'VV'},'</text>',"\n";

	for my $bp (@VirCovered) {
		my $depth = $VirusDepth->{$bp};
		my $x = ($bp - $VirMin +0.4) * $XRatio +$VirLeft+$IMG_Boder;
		#print "$bp:$depth -> $x\n";
		print $fh '<line stroke="red" stroke-width="',$XRatio*0.8,'" x1="',$x,'" y1="',$IMG_Boder+$BarHeight[0],'" x2="',$x,'" y2="',$IMG_Boder+$BarHeight[0]-$depth*$YRatio,'"/>'," $bp\n";
	}
	#print "@VirCovered\n";
	for my $bp (@HumCovered) {
		my $depth = $HumDepth->{$bp};
		my $x = ($bp - $HumMin +0.4) * $XRatio +$HumLeft+$IMG_Boder;
		#print "$bp:$depth -> $x\n";
		print $fh '<line stroke="blue" stroke-width="',$XRatio*0.8,'" x1="',$x,'" y1="',$IMG_Boder+$BarHeight[1],'" x2="',$x,'" y2="',$IMG_Boder+$BarHeight[1]+$depth*$YRatio,'"/>'," $bp\n";
	}
}

while (<$fhin>) {
	chomp;
	if (/^\[(B\d+)\]/) {
		my ($HumDepth,$VirusDepth,$sEdges,$pEdges,$SimInfo);
		if (defined $id and defined $vchr) {
			#print "[$id]\t$chr,$range,$cnt\t$vchr,$vrange,$vcnt\n",join("\n",@cache),"\n\n\n";
			warn "> $id - $1\n";
			($HumDepth,$VirusDepth,$sEdges,$pEdges,$SimInfo) = getdepth(\@cache);
			if (defined $HumDepth) {
				my $fh;
				open $fh,'>',"./plots/$id.svg";
				print $fh '<?xml version="1.0"?>
<svg width="',$IMG_Width,'" height="',$IMG_Height,'" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
';
				doplot($fh,$HumDepth,$VirusDepth,$sEdges,$pEdges,$SimInfo);
				print $fh "</svg>\n";
				close $fh;
			}
if ($id eq 'B1553') {
	#ddx $VirusDepth;
	#die "$id\n";
}
		}
		$id = $1;
		@cache = ();
		$chr='';
		$vchr=undef;
	} elsif (/^HostRange=([\w\:\-\,\|\.]+)$/) {
		my $t = $1;
		$id = undef if $t =~ /,/;
		($chr,$range,$cnt) = split /:/,$t;
		$id = undef if $chr ne $chrOnly;
	} elsif (/^VirusRange=([\w\:\-\,\|\.]+)$/) {
		my $t = $1;
		$id = undef if $t =~ /,/;
		($vchr,$vrange,$vcnt) = split /:/,$t;
		$id = undef if $vchr ne $viruschr;
	} elsif (/^SamFS=/) {
		next;
	} else {
		push @cache,$_ unless /^$/;
	}
}
