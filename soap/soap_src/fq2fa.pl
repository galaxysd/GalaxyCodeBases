perl -e 'my $count=0; while(<>){if($count%4==0) {if(!/^@(\S+)/) {die "unknown string $_ at line $count\n";} print ">$1\n";} elsif($count%4==1) {print $_;} $count++;}'

