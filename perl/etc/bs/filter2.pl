#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <snp file> <Cancer_Watson,Cancer_Crick,Normal_Watson,Normal_Crick> >out.txt\n" if @ARGV < 2;
my @maxDepth = split /,/,$ARGV[1];
warn "[!]MaxDepth:[",join(',',@maxDepth),"]\n";

no warnings 'qw';
my @thePOS = qw(chrom position);
my @LastEight = qw(Number_of_watson[A,T,C,G]_Normal Number_of_crick[A,T,C,G]_Normal Mean_Quality_of_Watson[A,T,C,G]_Normal Mean_Quality_of_Crick[A,T,C,G]_Normall Number_of_watson[A,T,C,G]_Cancer Number_of_crick[A,T,C,G]_Cancer Mean_Quality_of_Watson[A,T,C,G]_Cancer Mean_Quality_of_Crick[A,T,C,G]_Cancer);
my @SELECTED = (qw(ref var),@LastEight[4,5,0,1],'somatic_status',@LastEight[6,7,2,3]);
my %gOrder;
@gOrder{qw(A T C G)} = qw(0 1 2 3);

sub readnext($) {
	my $in = $_[0];
	if ( ! eof($in->[0]) ) {
		my $record = readline($in->[0]);
		die "readline failed: $!" unless defined $record;
		chomp($record);
		my %hash;
		@hash{@{$in->[1]}} = split /\t/, $record;
		@{$in}[2,3] = @hash{@thePOS};
		$in->[4] = [@hash{@SELECTED}];
		$in->[5] = $record;
		#ddx \%hash;
		return 1;
	} else {
		@{$in}[2,3] = qw(_EOF_ 0);
		$in->[4] = [qw(0 0 NA NA 0 0 NA NA NA)];
		$in->[5] = '';
		return 0;
	}
}

sub initfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.zst$/) {
			open( $infile,"-|","zstd -qdc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	chomp(my $t = <$infile>);
	#print "$t\n";
	my @tt = split("\t",$t);
	map { s/_read(\d)$/_reads$1/ } @tt;
	push @tt,qw(eW eC eA eB);
	#print join("|",@tt),"\n";
	@tt[-8 .. -1] = @LastEight;
	print join("\t",@tt),"\n";
	my $ret = [$infile,\@tt,undef,-1,[],[]];
	#readnext($ret) or die "[x]File [$filename] is empty. $!\n";
	return $ret;
}

my $FH = initfile($ARGV[0]);
#ddx $FH;

while($FH->[3]) {
	readnext($FH) or last;
	next unless defined $FH->[3];	# skip empty lines
	if (length($FH->[4][1])>1) {
		print $FH->[5],"\n";
		next;
	}
	my $t1 = $gOrder{$FH->[4][0]};	# ref
	my $t2 = $gOrder{$FH->[4][1]};	# var
	my $GT = join('',(sort {$a cmp $b} ($FH->[4][0],$FH->[4][1])));
	my @d1 = (split(',',$FH->[4][2]))[$t1,$t2,0,1,2,3];	# Watson_Cancer的Number, $d1[1]是Ref，$d1[2]是Var，$d1[3]是A，$d1[4]是T，$d1[5]是C，$d1[6]是G
	my @d2 = (split(',',$FH->[4][3]))[$t1,$t2,0,1,2,3];	# Crick_Cancer的Number
	my @d3 = (split(',',$FH->[4][4]))[$t1,$t2,0,1,2,3];	# Watson_Normal的Number
	my @d4 = (split(',',$FH->[4][5]))[$t1,$t2,0,1,2,3];	# Crick_Normal的Number
	my $s1 = $d1[0] + $d1[1];
	my $s2 = $d2[0] + $d2[1];
	my $sa = $s1+$s2;
	my $s3 = $d3[0] + $d3[1];
	my $s4 = $d4[0] + $d4[1];
	my @q1 = (split(',',$FH->[4][7]))[$t1,$t2,0,1,2,3];	# Watson_Cancer的Mean_Quality, $q1[1]是Ref，$q1[2]是Var，$q1[3]是A，$q1[4]是T，$q1[5]是C，$q1[6]是G
	my @q2 = (split(',',$FH->[4][8]))[$t1,$t2,0,1,2,3];	# Crick_Cancer的Mean_Quality
	my @q3 = (split(',',$FH->[4][9]))[$t1,$t2,0,1,2,3];	# Watson_Normal的Mean_Quality
	my @q4 = (split(',',$FH->[4][10]))[$t1,$t2,0,1,2,3];	# Crick_Normal的Mean_Quality
	if ($FH->[4][6] eq 'Germline') {
		;
	} elsif ($FH->[4][6] eq 'Somatic') {
		next if $sa<15;
		if ($s1>$maxDepth[0] or $s2>$maxDepth[1] or $s3>$maxDepth[2] or $s4>$maxDepth[3]) {
			next;
		}
		if ($GT eq 'AA') {
			;
		} elsif ($GT eq 'TT') {
			;
		} elsif ($GT eq 'AT') {
			# 1.1 突变碱基正负链大等于2，1.2 突变频率大于0.2，深度大于15。1.3 突变碱基平均质量值大于20 1.4 除了AT，其他大都是0
			next if $d1[1]<2 or $d2[1]<2;
			next if (($d1[1]+$d2[1])/$sa)<0.2;
			next if ($q1[1]+$q2[1])<40;
			next if ($d1[2]+$d1[3]+$d1[4]+$d1[5])>($d1[0]+$d1[1]);
		} elsif ($GT eq 'CC') {
			;
		} elsif ($GT eq 'GG') {
			;
		} elsif ($GT eq 'CT') {
			;
		} elsif ($GT eq 'AG') {
			;
		} elsif ($GT eq 'GT') {
			;
		} elsif ($GT eq 'CG') {
			;
		} elsif ($GT eq 'AC') {
			;
		} else {die;}
	} elsif ($FH->[4][6] eq 'LOH') {
		;
	} else {
		;
	}

	print $FH->[5],"\n";
	#ddx $FH;
}

__END__

先分类型：
第一， germline过滤条件：（原则：原则尽量宽松）
只考虑两个文件都是正负链都有支持突变类型。就通过。

第二，somatic过滤条件（原则：非常严格）
（一）癌症组（EA，EB，EC，ED），
1. 如果genotype是AT，AA，TT
1.1 突变碱基正负链大等于2，
1.2 突变频率大于0.2，深度大于15。
1.3 突变碱基平均质量值大于20

1. 如果genotype是AT，
1.1 突变碱基正负链大等于2，
1.2 突变频率大于0.2，深度大于15。
1.3 突变碱基平均质量值大于20
1.4 除了AT，其他大都是0

2.  如果是AA，
2.1 支持A正负链都大于5
2.2 总深度大于15
2.3 碱基质量值大于20
2.4 除了A，其他都是0

3. 如果是TT
3.1 支持T正负链都大于5
3.2 总深度大于15
3.3 碱基质量值大于20
3.4 除了T，其他都是0



3. genotype是CC
正链(C+T)>5,G=0,A=0
负链C>5,A=0,T=0,G=0
depth>15

6. GG
正链：G>5, 其他都是0。
负链：G+A>5,其他都是0。
depth >15

4. genotype是CT
正链：C+T>5,其他都是0
负链：C>2,T>2，其他都是0，C/（C+T）>0.2, T/(C+T)>0.2
detph>15

7. AG
正链：A>2,G>2,其他都是0, G/（G+A）>0.2, A/(G+A)>0.2
负链：A+G>5，其他都是0.
detph>15

2. 如果genotype是AC，
总depth>15
正链A>2,(C+T)>2,G=0
负链A>2,C>2.T=0,G=0
全部A  质量值>30
负链C  质量值>30
# MAF>0.2 =>
(C+T)/total >0.2 and A/total >0.2

5，genotype是TG,
正链：T>2, G>2, 其他都是0, T/(T+G)> 0.2, G/(T+G)>0.2
负链：T>2,(G+A)>2,C=0, T/(T+G+A)>0.2, (G+A)/(T+G+A) > 0.2
depth >15


8. CG
正链：（C+T）>2; G>2, A=0, G/(T+G+C)> 0.2, (C+T)/(C+T+G)>0.2
负链：C>2,(G+A)>2, T=0, C/(C+A+G)>0.2; (G+A)/(C+A+G)>0.2
depth >15


（二）正常组（ Number_of_watson[A,T,C,G]       Number_of_crick[A,T,C,G]        Mean_Quality_of_Watson[A,T,C,G] Mean_Quality_of_Crick[A,T,C,G]）
1. 如果ref是A或者T，
1.1 除了ref碱基，其他碱基都是0.
1.2 ref碱基正负链都大于5，
1.3 ref总深度大于15

2.如果ref是C
2.1 除了C碱基和正链的T碱基，其他碱基都是0
2.2 负链C碱基质量值大于30，并且个数大于5，正链（C+T大于5）

3. 如果ref是G
3.1 除了G碱基和负链A碱基，其他碱基个数都是0
3.2 正链G碱基质量值大于30，并且个数大于5，负链（A+G大于5个）

第三，LOH的过滤条件（这个条件更麻烦，后面单独讨论）
