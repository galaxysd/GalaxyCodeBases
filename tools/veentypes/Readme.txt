For an all default run:	./vnstat.pl -b

Help:	./vnstat.pl

Full options:	./vnstat.pl -i gff.lst -c chr.len -l 10 -o stat.txt -d details.lst -t /var/tmp/swap.tmp -b

Usage: ./vnstat.pl [-OPTIONS [-MORE_OPTIONS]] [--] [PROGRAM_ARG1 ...]

The following single-character options are accepted:
        -i GFF files list (gff.lst) [SampleName\tPath_to_file\n](8 lines max !)
        -c Chromosome length list (chr.len) [ChrID\s+Len\n]
        -l Minimal overlap length (10)
        -o Output Stat (stat.txt)
        -d Details dump to (details.lst)
        -t tmpfile for swap, will be removed (./._tmp_)
        -v show verbose info to STDOUT
        -b No pause for batch runs

Options may be merged together.  -- stops processing of options.
Space is not required between options and their arguments.
