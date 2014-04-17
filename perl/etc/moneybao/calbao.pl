#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dump;

die "Usage: $0 <Money> <年化收益>%\n" if @ARGV < 2;
my $total = shift;
my $annelrate = shift;

=pod
=((1+B1/10000)*(1+B2/10000)*(1+B3/10000)*(1+B4/10000)*(1+B5/10000)*(1+B6/10000)*(1+B7/10000))^(365/7)-1
=(((G1+1)^(7/365))^(1/7)-1)*10000
华夏财富宝货币
http://www.chinaamc.com/product/fundLishijingzhi.do?fundcode=000343&isQuery=&pageIndex=1&begindate=2013-04-01&enddate=

天弘增利宝货币
http://www.thfund.com.cn/website/funds/fundnet.jsp?fundcode_select=000198&fundcode=null&startdate=2014-03-1&enddate=2014-04-30&x=31&y=8
http://www.thfund.com.cn/website/hd/zlb/newzlbrev2.jsp
=cut
my $fitness = ((($annelrate/100 +1) ** (7/365)) ** (1/7) -1) * 10000;

print "$annelrate % => $fitness\n";

my $money0 = $total;
for my $i (1..31) {
	my $add = int(($total*$fitness/10000)*100 +1)/100;
	print "$i:\t",join(', ',$total,$add,$total+$add),"\n";
	$total += $add;
}

print '-' x 75,"\n$total - $money0 = ",int(($total - $money0)*100 +1)/100,"\n";
