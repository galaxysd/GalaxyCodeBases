#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <section_len> <f/r patten> <nodes_overlap> [collapse_repeats] [out_prefix]\n" if @ARGV < 3;
my ($SecLen,$FRpatten,$OverlapNodes,$Collapse,$out)=@ARGV;
$out='olc' unless $out;

#my @Seq=qw/1 2 3 4 A B C D 5 6 7 8 A B C D 9 10 11 12 D C B A 13 14 15 16/;
#my @Seq=qw/1 2 3 A B C 4 5 6 A B C 7 8 9 C B A 10 11 12/;
my @Seq;
die "[x] section_len is [2,26].\n" if $SecLen<2 or $SecLen>26;
$FRpatten=lc $FRpatten;
my $u=0;
push @Seq,++$u for (1..$SecLen);
for (split //,$FRpatten) {
    if ($_ eq 'f') {
        push @Seq,chr(64+$_) for (1..$SecLen);
    } elsif ($_ eq 'r') {
        push @Seq,chr(65+$SecLen-$_) for (1..$SecLen);
    }
    push @Seq,++$u for (1..$SecLen);
}
my $Seq=join('-',@Seq);
my $filename="$out.$FRpatten${SecLen}o$OverlapNodes";
$filename .= 'c' if $Collapse;
print "Seq: $Seq\nOut: $filename.{gv,png}\n";

my @U=grep(/^\d+$/,@Seq);
my @R=grep(/[^\d]/,@Seq);
my %t;
++$t{$_} for @U;
@U=sort keys %t;
%t=();
++$t{$_} for @R;
@R=sort keys %t;
my $RepeatCount=(values %t)[0];
for (values %t) {
    $RepeatCount=$_ if $RepeatCount>$_;
}
%t=();
#print "[$_] " for @seq;

sub getName($) {
    my $tag=$_[0];
    return [$tag] if $Collapse;
    if ($tag =~ /^\d+$/) {
        return [$tag];
        #return ["u$tag"];
    } else {
        my @t;
        #push @t,"r${tag}$_" for (1..$RepeatCount);
        push @t,"${tag}$_" for (1..$RepeatCount);
        return \@t;
    }
}

open O,'>',$filename.'.gv' or die "$!";
print O <<HEAD;
graph "OLC" {
\trankdir=LR;
\tgraph [ fontname = "Arial", label = "OLC plot of $RepeatCount Repeats with Overlap=$OverlapNodes\\nSeq: $Seq" ];
HEAD
print O "\tnode [shape = ellipse]; ";
print O " @{&getName($_)}" for (@R);
print O ";\n\tnode [shape = box];\n";
unless ($Collapse) {
#print O "\t{rank=same; ";
    print O "\t{rank=same; @{&getName($_)} ;}\n" for (@R);
#print O ";}\n";
}
=pod
for (my $i=0;$i<=$#Seq-$OverlapNodes;$i++) {
    #print "[$Seq[$i]],",join(",",@{&getName($Seq[$i])}),"\n";
    for my $a (@{&getName($Seq[$i])}) {
        my @t;
        for my $j (1..$OverlapNodes) {
            push @t,@{&getName($Seq[$i+$j])};
        }
        print O "\t$a -> { ",join('; ',@t)," };\n";
    }
}
=cut
my %Edges;
for (my $i=0;$i<$#Seq;$i++) {
    #print "[$Seq[$i]],",join(",",@{&getName($Seq[$i])}),"\n";
    for my $a (@{&getName($Seq[$i])}) {
        for my $j (1..$OverlapNodes) {
            push @{$Edges{$a}},@{&getName($Seq[$i+$j])} if $i+$j<=$#Seq;
        }
    }
}

sub getStr($$) {
    my ($a,$b)=sort @_;
    my $v=0.1;
    if ("$a$b" =~ /^\d+$/) {
        $v=10;
    } elsif ($a !~ /^\d+$/ and $b !~ /^\d+$/) {
        $v=1;
    }
    #print "$a -- $b ,$v\n";
    return ["$a -- $b",$v];
}
my %str;
for my $a (sort keys %Edges) {
    %t=();
    ++$t{$_} for @{$Edges{$a}};
    #$Edges{$a} = [keys %t];
    #print O "\t$a -> { ",join('; ',keys %t)," };\n";    # dot
    #print O "\t$a -- $_ ;\n" for sort keys %t;
    for my $b (keys %t) {
        my ($tstr,$v)=@{&getStr($a,$b)};
        $str{$tstr} += $v;
    }
}
print O "\t$_ [weight=$str{$_}];\n" for sort keys %str;
%t=();
print O "}\n";
close O;
system('dot','-Tpng',"-o$filename.png","$filename.gv");

