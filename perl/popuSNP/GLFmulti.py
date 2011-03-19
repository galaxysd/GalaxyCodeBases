#!/usr/bin/python
import sys
import struct
import math
import time
decode = ('N','A','C','N','G','N','N','N','T','N','N','N','N','N','N','N')
code = (0,5,15,10,1,3,2,7,6,11)
rev_code = (0,4,6,5,4,1,8,7,6,8,3,9,5,7,9,2)

# potential useful
def factorial(N_hap):
    log10_fact = [0.0]*(N_hap+1)
    for i in xrange(1,N_hap+1):
        log10_fact[i] = log10_fact[i-1]+math.log10(i)
    return log10_fact

def choose(N_hap, log10_fact):
    log10_choose = [0.0]*(N_hap+1)
    for i in xrange(1,N_hap):
        log10_choose[i] = log10_fact[N_hap] - log10_fact[i] - log10_fact[N_hap-i]
    return log10_choose

def prior_mat_gen(N_hap):
    log10_prob = []
    log10_fact = factorial(2)
    log10_choose = choose(2,log10_fact)
    for minor_count in xrange(N_hap+1):
        log10_prob.append([0.0,]*3)
        if minor_count == 0:
            freq = 0.1/N_hap
        elif minor_count == N_hap:
            freq = 1-0.1/N_hap
        else:
            freq = minor_count/float(N_hap)
        for i in xrange(3):
            log10_prob[minor_count][i] = log10_choose[i] + i*math.log10(freq)+(2-i)*math.log10(1.0-freq)
    return log10_prob

def log2normal():
    power = [0.0,]*256
    for i in xrange(256):
        power[i] = math.pow(10,-i/10.0)
    return power

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
    return (decode[ref],copyNum, depth, LLR)
Files = []
if len(sys.argv)!=5:
    print >>sys.stderr, "python", sys.argv[0],"START END FILELIST OUTPUT"
    sys.exit(1)
else:
    start = int(sys.argv[1])-1; end = int(sys.argv[2])
    GLFlist = open(sys.argv[3],"r")
    Files = []
    for line in GLFlist:
        Files.append(open(line.rstrip(),"rb"))
        assert(Files[-1])
    output = open(sys.argv[4],"w")

N_indi = len(Files)
N_hap = 2*N_indi
minLLR = int(20+10*math.log10(N_hap))
log10_prob = prior_mat_gen(N_hap)
power = log2normal()
lo = [0,]*4 # Searching space of the 4 bases, the lower cut
hi = [0,]*4
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
        file.seek(12*(start), 1)
    chrName = chrName[:-1]
    print >>sys.stderr, "%s %d start: %d end: %d" %(chrName, chrLen, start+1, end)

    for pos in xrange(start, end):
        if pos >= chrLen:
            break
        count = [0.0,]*4
        total_depth = 0
        copyNum = 0

        ref = 'N'
        for i in xrange(N_indi):
            (new_ref,indi_CN, indi_dep, LLR) = read_base(Files[i])
            total_depth += indi_dep
            copyNum += (indi_dep*(indi_CN+0.5))
            if 'N' == ref:
                ref = new_ref
            else:
                if 'N' != new_ref:
                    assert ref == new_ref
            #LLR_sum = 0.0
            for type in xrange(10):
                LLR_mat[i][type] = power[LLR[type]]
                if type < 4: # homo
                    count[code[type]&3] += LLR_mat[i][type]
                else:
                    count[code[type]&3] += LLR_mat[i][type]
                    count[(code[type]>>2)&3] += LLR_mat[i][type]
        # select the best 2 bases
        base1 = -1; base2 = -1
        for base in xrange(4):
            if -1==base1 or count[base] > count[base1]:
                base2 = base1
                base1 = base
            elif -1==base2 or count[base] > count[base2]:
                base2 = base
            else:
                pass

        G_LL = [0.0,]*(N_hap+1)
        besttype = [rev_code[(base1<<2)|base1], rev_code[(base1<<2)|base2],rev_code[(base2<<2)|base2]]
        for minor_count in xrange(N_hap+1):
            if minor_count >= N_hap:
                break
            else:
                for indi in xrange(N_indi):
                    indi_G_LL = 0.0
                    for i in xrange(3):
                        indi_G_LL += (math.pow(10,log10_prob[minor_count][i])*LLR_mat[indi][besttype[i]])
                    G_LL[minor_count] += math.log10(indi_G_LL)
        #print count
        #print G_LL
        est_G = -1
        for minor_count in xrange(N_hap+1):
            if minor_count >= N_hap-2:
                break
            else:
                if G_LL[minor_count] > G_LL[minor_count+1] > G_LL[minor_count+2]:
                    est_G = minor_count
                    break
                else:
                    pass
        if est_G == -1:
            est_G = 0
        #print "%s\t%d\t%c\t%c\t%c\t%d\t%d\t%f\t%f" %(chrName, pos, ref, "ACTG"[base1], "ACTG"[base2], N_hap-est_G, est_G, count[base1], count[base2])
        if total_depth == 0:
            copyNum = 15
        else:
            copyNum /= total_depth
        try:
            print >>output, "%s\t%d\t%c\t%d\t%f\t%c\t%c\t%d\t%d\t%d" %(chrName, pos+1, ref, total_depth, copyNum, "ACTG"[base1], "ACTG"[base2], N_hap-est_G, est_G, int(10*(G_LL[est_G]-G_LL[0])+0.5))
        except ValueError:
            print >>sys.stderr, G_LL[0],G_LL[est_G]
            sys.exit(1)
    print >>sys.stderr, "Finished", chrName
print >>sys.stderr, "Mission Complete!"

