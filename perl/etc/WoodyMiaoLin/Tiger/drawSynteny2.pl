#!/usr/bin/perl
use strict;
use warnings;

die "Usage: $0 <gene> <SNP> <output>\n" if @ARGV<3;
my $in1=shift;
my $in2=shift;
my $out=shift;

open I1, "<", $in1 or die "Error opening $in1 : $!\n";
open I2, "<", $in2 or die "Error opening $in2 : $!\n";
open O, ">", $out or die "Error opening $out : $!\n";


my $nSpeci = 4; # 物种数
my $nGene = 22; # 基因数
my $yScale = 0.0002; # y轴方向的缩放系数，等于基因的高度除以DNA碱基数
my $genWid = 10; # 基因矩形宽度的一半
my $spaChr = 200; # 染色体间距

# 读取染色体与基因信息
my @yCoo;
for my $a (0..$nSpeci) {
  for my $b (0..$nGene) {
    my $c = <I1>;
    chomp $c;
    @{$yCoo[$a][$b]} = split /\t/, $c;
  }
}
close I1;

# 读取SNP信息
my %SNPcoo;
while (<I2>) {
  chomp;
  my @a = split /\t/;
  $SNPcoo{$a[0]} = $a[1];
}
close I2;

print O "<?xml version=\"1.0\"?>\n\n<svg xmlns=\"http://www.w3.org/2000/svg\">\n\n";

# 画染色体 Scaffold75长5658748
for (2..$nSpeci) {
  my $x = $spaChr*$_;
  my $y1 = 0;
  my $y2 = 9302904*$yScale;
  print O "<line x1=\"$x\" y1=\"$y1\" x2=\"$x\" y2=\"$y2\" stroke=\"black\" stroke-width=\"5\"/>\n";
}
print O "<line x1=\"$spaChr\" y1=\"0\" x2=\"$spaChr\" y2=\"", 5658748*$yScale, "\" stroke=\"black\" stroke-width=\"5\"/>\n";
print O "<line x1=\"$spaChr\" y1=\"", 5758748*$yScale, "\" x2=\"$spaChr\" y2=\"", 9302904*$yScale, "\" stroke=\"black\" stroke-width=\"5\"/>\n";

# 画基因
for my $a (1..$nSpeci) {
  for my $b (1..$nGene) {
    my $x = $spaChr*$a - $genWid;
    my $y = $yScale*$yCoo[$a][$b][0];
    my $w = 2*$genWid;
    my $h = $yScale*($yCoo[$a][$b][1]-$yCoo[$a][$b][0]);
    print O "<rect x=\"$x\" y=\"$y\" width=\"$w\" height=\"$h\" fill=\"black\"/>\n";
  }
}

# 画基因连线
for my $b (1..$nGene) {
  for my $a (1..$nSpeci-1) {
    for my $c (0..1) {
      my $x1 = $spaChr*$a + $genWid;
      my $y1 = $yScale*$yCoo[$a][$b][$c];
      my $x2 = $spaChr*($a+1) - $genWid;
      my $y2 = $yScale*$yCoo[$a+1][$b][$c];
      print O "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" stroke=\"black\" stroke-width=\"0.5\"/>\n";
    }
  }
}

# 打印基因名
for (1..$nGene) {
  my $x = $nSpeci*$spaChr + 2*$genWid;
  my $y = 0.5*$yScale*($yCoo[$nSpeci][$_][0] + $yCoo[$nSpeci][$_][1]) + 2;
  my $gene = $yCoo[0][$_][0];
  print O "<text x=\"$x\" y=\"$y\" font-family=\"Courier\" font-size=\"8\" fill=\"black\">$gene</text>\n";
}

# 打印染色体名
for (1..$nSpeci) {
  my $x = ($_-0.4)*$spaChr;
  my $scaffold = $yCoo[$_][0][0];
  print O "<text x=\"$x\" y=\"-10\" font-family=\"Courier\" font-size=\"12\" fill=\"black\">$scaffold</text>\n";
}

# 标记SNP
for (sort {$a <=> $b} keys %SNPcoo) {
  my $x1 = $spaChr - 3*$genWid;
  my $x2 = $spaChr;
  my $y = $yScale*$_;
  print O "<line x1=\"$x1\" y1=\"$y\" x2=\"$x2\" y2=\"$y\" stroke=\"black\" stroke-width=\"1\"/>\n";
  my $SNP = "$_($SNPcoo{$_})";
  my $ty = $yScale*$_ + 2;
  print O "<text x=\"120\" y=\"$ty\" font-family=\"Courier\" font-size=\"8\" fill=\"black\">$SNP</text>\n";
}

print O "\n</svg>\n";
close O;
