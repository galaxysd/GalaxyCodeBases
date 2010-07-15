#!/usr/bin/python
import sys
import struct
import math
import time

def read_header(infile):
    try:
        assert 'glf' == infile.read(3)
    except AssertionError:
        print >>sys.stderr, "File is not GLF"
        sys.exit(1)
    (major_ver, minor_ver) = struct.unpack("1i1i", infile.read(8))
    #print >>sys.stderr, "GLF v%d.%d" %(major_ver, minor_ver)
    while(True):
        (lineTypeNameLen,) = struct.unpack("1i", infile.read(4))
        (lineTypeName,) = struct.unpack("%ds" %(lineTypeNameLen), infile.read(lineTypeNameLen))
        if "CHROMOSOMES" == lineTypeName[0:11]:
            break
        else:
            (dataLen,) = struct.unpack("1i", infile.read(4))
            dataLen = int(dataLen)
            (data,) = struct.unpack("%ds" %(dataLen), infile.read(dataLen))
        print >>sys.stderr, lineTypeName, data
    (chrNum,) = struct.unpack("1i", infile.read(4))
    return chrNum

def read_chr(infile):
    (chrNameLen,) = struct.unpack("1i", infile.read(4))
    (chrName,) = struct.unpack("%ds" %(chrNameLen), infile.read(chrNameLen))
    (chrLen,) = struct.unpack("1i", infile.read(4))
    print >>sys.stderr, chrNameLen, chrName, chrLen
    return (chrName, chrLen)

def read_base(infile):
    global pos;
    (tmp1, tmp2) = struct.unpack("2c", infile.read(2))
    ref = (ord(tmp1)>>4)&0xF
    copyNum =  (ord(tmp2)&0xF)
    depth = ((ord(tmp1)&0xF)<<4)|((ord(tmp2)>>4)&0xF)
    #print >>sys.stderr, "%c\t%d\t%d\t" %(decode[ref], copyNum, depth),
    LLR = map(ord, struct.unpack("10c", infile.read(10)) )
    #print >>sys.stderr, "\t".join(map(str, LLR))
    #sys.exit(1)
    return (decode[ref],copyNum, depth, LLR)

Files = []
if len(sys.argv)!=4:
    print >>sys.stderr, "python", sys.argv[0],"SNPinfo FILELIST OUTPUT"
    sys.exit(1)
else:
    SNPinfo = open(sys.argv[1],"r")
    GLFlist = open(sys.argv[2],"r")
    Files = []
    for line in GLFlist:
        Files.append(open(line.rstrip(),"rb"))
        assert(Files[-1])
    N_indi = len(Files)
    output = open(sys.argv[3],"w")

decode = ('N','A','C','N','G','N','N','N','T','N','N','N','N','N','N','N')
abbv = ('A','M','W','R','M','C','Y','S','W','Y','T','K','R','S','K','G','N')
code = (0,5,15,10,1,3,2,7,6,11)
rev_code = (0,4,6,5,4,1,8,7,6,8,3,9,5,7,9,2)

# Read the SNPs
snp_hash = {}
for line in SNPinfo:
    tabs = line.split()
    #chr9    121463515       0       1       0       0.01    0.01    0       0       rs236
    if not snp_hash.has_key(tabs[0]):
        snp_hash[tabs[0]] = {}
    snp = 0
    for i in xrange(5,9):
        if float(tabs[i]) > 0:
            snp = ((snp<<2)|(i-5))
        else:
            pass
    snp = (snp&0xF)
    assert snp < 16
    snp_hash[tabs[0]][int(tabs[1])-1] = snp

for chr_name in snp_hash.iterkeys():
    snps = snp_hash[chr_name].items()
    snps.sort()
    snp_hash[chr_name] = snps

# Get number of chromosomes
LLR_mat = []
for i in xrange(N_indi):
    LLR_mat.append([0.0,]*10)
chrNum = -1
for file in Files:
    if -1 == chrNum:
        chrNum = read_header(file)
    else:
        ##### Emergency
        try:
            assert chrNum == read_header(file)
        except AssertionError:
            print >>sys.stderr, "Inconsistent number of chromosomes. Only the first chromosome is called."
            chrNum =1

print >>sys.stderr, "Total chromosomes %d" %(chrNum)
for chr_count in xrange(chrNum):
    chrName = ""; chrLen = -1
    for file in Files:
        if -1==chrLen:
            (chrName, chrLen) = read_chr(file)
        else:
            assert (chrName, chrLen) == read_chr(file)
        #file.seek(12*(start), 1)
    prev_pos = 0
    chrName = chrName[:-1]
    n_pos = len(snp_hash[chrName])
    print >>sys.stderr, "%s length=%d building consensus for all individuals in %d SNP sites" %(chrName, chrLen, n_pos)
    finished = 0
    for (pos, type) in snp_hash[chrName]:
        if pos >= chrLen:
            print >>sys.stderr, "Warning!!! Your SNPs are located outside the chromosome!!! I am stopping!"
            break
        #count = [0.0,]*4
        total_depth = 0
        copyNum = 0
        print >>output, "%s\t%d\t%c" %(chrName, pos+1, abbv[type]),
        ref = 'N'
        for i in xrange(N_indi):
            Files[i].seek(12*(pos-prev_pos),1)
            (new_ref,indi_CN, indi_dep, LLR) = read_base(Files[i])
            best_genotype = 16
            second_genotype = 16
            for glf_type in xrange(10):
                genotype = code[glf_type]
                if (genotype^type) == 0: # HET SNP
                    pass
                elif ((genotype^type)&3 == 0) or (((genotype^type)>>2)&3 == 0): # HOM SNP or one allele correct
                    if (genotype&3) == ((genotype>>2)&3): # HOM SNP
                        LLR[glf_type] += 0
                    else: # One allele mutation
                        LLR[glf_type] += 27
                else: # Two allele mutation
                    if (genotype&3) == ((genotype>>2)&3): # HOM mutation
                        LLR[glf_type] += 30
                    else: # HET two allele mutation
                        LLR[glf_type] += 54
                if best_genotype == 16 or LLR[glf_type] < LLR[rev_code[best_genotype]] :
                    second_genotype = best_genotype
                    best_genotype = genotype
                elif second_genotype == 16 or LLR[glf_type] < LLR[rev_code[second_genotype]]:
                    second_genotype = genotype
                else:
                    pass
            print >>output, "\t%c\t%c\t%d" %(abbv[best_genotype], abbv[second_genotype], LLR[rev_code[second_genotype]]-LLR[rev_code[best_genotype]]),
            total_depth += indi_dep
            #copyNum += (indi_dep*(indi_CN+0.5))
            if 'N' == ref:
                ref = new_ref
            else:
                if 'N' != new_ref:
                    assert ref == new_ref
            #LLR_sum = 0.0
        prev_pos = pos+1
        print >>output, "\t%d" %(total_depth)
        finished += 1
        if 0 == finished % 10000:
            print >>sys.stderr, "%d sites finished" %(finished)
    print >>sys.stderr, "Finished", chrName
print >>sys.stderr, "Mission Complete!"