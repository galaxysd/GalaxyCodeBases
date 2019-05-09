package FGI::OYK;

use strict;
use warnings;
use Data::Dump qw(ddx);

require Exporter;

my $SgeQueue = 'fgi.q';
my $SgeProject = 'fgi';
my $SgeJobPrefix = 'nip'.(time() % 999);

#----- systemic variables -----
our (@ISA, @EXPORT, $VERSION, @EXPORT_OK, %EXPORT_TAGS);
@ISA = qw(Exporter);
@EXPORT = qw(SHcutadapt);
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw(cpi getcpi)],
                 OTHER   => [qw(mendnumber)]);
$VERSION = "1.0";

our ($SHcutadapt);

$SHcutadapt = <<"END_SH";
#!/bin/sh
#$ -S /bin/bash
#$ -q $SgeQueue -P $SgeProject
#$ -N 0${SgeJobPrefix}cut
#$ -l vf=600M,num_proc=2
#$ -binding linear:3
#$ -cwd -r y
#$ -v PERL5LIB,PATH,LD_LIBRARY_PATH


END_SH

1;
