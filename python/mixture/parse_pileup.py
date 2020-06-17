#!/usr/bin/env python3
# Modified by Daniel Gomez-Sanchez: adding opening bgzipped pileups
# <https://github.com/magicDGS/PoolHMM/blob/master/Modules/parse_pileup.py>
# <https://github.com/bihealth/vcfpy/blob/master/vcfpy/parser.py>
import sys
import getopt
import re
import os.path
import gzip

import zstandard
import io

def openPileup(pileup_file, log=False):
    """
    Open pileup in text format or gzipped
    pileup_prefix is the prefix of .pileup.gz or .pileup
    mode is the way of opening the file
    log=True prints if compressed or not compressed are used
    return a filehandler
    exit if IOError found
    """
    #pileup file loading
    try:
        if pileup_file.endswith('.zstd'):
            fh = open(pileup_file, 'rb')
            dctx = zstandard.ZstdDecompressor()
            stream_reader = dctx.stream_reader(fh)
            pileup = io.TextIOWrapper(stream_reader, encoding='utf-8')
        elif pileup_file.endswith('.gz'):
            pileup = gzip.open(pileup_file, 'rt')
        else:
            pileup = open(pileup_file, 'rt')
        if log:
            print ("Input file: %s" % pileup_file)
        ## return the fh
        return pileup
    except IOError:
        print ("Could not open input file %s" %  pileup_file)
        sys.exit()

def mpileup_parser(line):
    if line == "\n":
        print ("pileup parser empty line provided")
        sys.exit(8)
    #...

    #line is splited to get ch, pos, rc, cov, nucs & qual infos
    try:
        line_items=line.split()
    except ValueError:
        print ("could not parse pileup entry %s" % (line))
        sys.exit(9)
    #...

    LineSplited = len(line_items)
    if (LineSplited % 3) == 0 and LineSplited >= 6:
        SampleCnt = (LineSplited-3) / 3
        (ch, pos, rc) = line_items
    SampleRows=['cov', 'nucs', 'qual']
    SampleNO = range(1,SampleCnt+1)
    SampleID = [[str(i)+str(j) for j in SampleRows] for i in SampleNO]
        #(ch, pos, rc, cov, nucs, qual) = line_items
    else:
        print ("wrong number of columns in pileup line (SampleCnt=%d): %s" % (SampleCnt,line))
        sys.exit()

    #nucs is filtered
    #step 1
    s1 = re.findall(r"[-+](\d+)", nucs)

    for n in s1:
        l1 = ['[-+]',n,'[ACGTNacgtn]{',n, '}']
        pattern = r''.join(l1)
        nucs = re.sub(pattern, '', nucs)
    #...

    #step2
    nucs = re.sub('\^.', '', nucs)

    #step3
    nucs = re.sub('\$', '', nucs)
    
    #step 4
    nucs = re.sub('\.', rc.upper(), nucs)
    
    #step 5
    nucs = re.sub(',', rc.lower(), nucs)
    #...

    #exception size of sequence != size of quality
    if len(nucs) != len(qual):
        print ("Size of sequence does not equal size of quality: %s, %s" % (nucs,line))
        sys.exit(10)
    #...
    
    #filter the pileup file by quality 
    i, ac, tc, cc, gc, nco, dell, co = 0, 0, 0, 0, 0, 0, 0, 0

    nucs_filtered = []
    qual_filtered = []
    for qc in qual:
        nc = nucs[i]
        quality = self.encode(qc)
        if quality >= self.minQual:
            co += 1
            if nc == "A" or nc == "a":
                ac += 1
            elif nc == 'T' or nc == 't':
                tc += 1
            elif nc == "C" or nc == "c":
                cc += 1
            elif nc == "G" or nc == "g":
                gc += 1
            elif nc == "N" or nc == "n":
                nco += 1
            elif nc == "*":
                dell += 1
            else: 
                print ("Could not parse pileup; Unknown allele : %s in %s" % (a,line))
                sys.exit(1)
            #...
            if nc != "*" and nc != "n" and nc != "N":
                nucs_filtered.append(nc)
                qual_filtered.append(quality)

        #...
        i += 1
    #...

    alar = [{'a':'A', 'c':ac}, {'a':'T', 'c':tc}, {'a':'C', 'c':cc}, {'a':'G', 'c':gc}]
    eucov = ac + tc + cc + gc
    
    if len(nucs_filtered) != len(qual_filtered):
        print ("Error: Length of nucleotides does not agree with lenght of quality string!")
        sys.exit(11)
    #...
    
    if len(nucs_filtered) != eucov:
        print ("Error : Coverage does not agree with length of nucleotides : %d n: %d " % (eucov, len(nucs_filtered)))
        sys.exit(12)
    #...
    
    # pos chr refc nucs qual totcov eucov alleles A T C G N del valid valid_alleles derived_allele allele_code
    entry={
        'pos':pos,
        'ch':ch,
        'refc':rc,
        'nucs':nucs_filtered,
        'qual':qual_filtered,
        'totcov':co,
        'eucov':eucov,
        'alleles':alar,
        'del':dell,
        'N':nco
    };

    alleles = 0
    if ac >= self.minCount:
        alleles += 1
    #...

    if tc >=  self.minCount:
        alleles += 1
    #...

    if cc >=  self.minCount:
        alleles += 1
    #...

    if gc >=  self.minCount:
        alleles += 1
    #...
    
    allele_code = "na"
    if alleles == 1:
        allele_code = "M"
    elif alleles == 2:
        allele_code = "S"
    elif alleles >2:
        allele_code = "T"
    entry['allele_code'] = allele_code        
    #..

    valid = 1
    if  entry['del'] >0 or entry['eucov'] < self.minCoverage or entry['eucov'] > self.maxCoverage:
        valid = 0
    entry["valid"] = valid

    entry['ancestral_allele']="N"
    entry['derived_allele']="N"
    entry['removed_alleles']= 0
    entry['unfolded']=1
    entry['refc']=entry['refc'].upper()

    return entry
#...


if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    import sys
    with openPileup(sys.argv[1]) as TextIn:
        for line in TextIn:
            record = mpileup_parser(line)
            print(record)
    exit()
