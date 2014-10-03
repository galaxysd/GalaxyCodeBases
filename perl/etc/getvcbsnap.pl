#! /usr/bin/perl -w

my $baseURL='http://www.svwfl.com/cdn';
my $project=shift @ARGV;
for my $n (@ARGV) {
	system "wget $baseURL/$project/$n.jpg";
	system "wget $baseURL/$project/${n}v.jpg";
}

__END__
getvcbsnap.pl code 1543 7166 12620 18688 53750 60806 68689 77932

