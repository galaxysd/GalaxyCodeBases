#!/opt/rocks/bin/python
import sys

try:
    snp = open(sys.argv[1],"r")
    min_soap_qual = int(sys.argv[2])
    max_soap_rep = float(sys.argv[3])
    min_dist = int(sys.argv[4])
    min_depth = int(sys.argv[5])
    max_soap_depth = int(sys.argv[6])
except IndexError:
    print >>sys.stderr, "python",sys.argv[0],"Snp MinQual Max_soap_rep MinDist MinDepth MaxDepth"
    sys.exit(1)

decode = ["A", "M", "W", "R", "M", "C", "Y", "S", "W", "Y", "T", "K", "R", "S", "K", "G"]
code = {}
for i in xrange(16):
    code[decode[i]]=i

for line in snp:
    tabs = line.split()
    count = len(tabs)
    if count< 14:
        continue
    ori = tabs[2]; soap_base = tabs[3]; qual = int(tabs[4])
    base1Depth = tabs[8]
    base2Depth = tabs[12]
    depth = int(tabs[13]); rank_sum = float(tabs[14])
    rep = float(tabs[15]); isSNP = (1==int(tabs[16]) or 5==int(tabs[16]));dist = int(tabs[17])
    if qual>= min_soap_qual and depth <= max_soap_depth and rep <= max_soap_rep and (dist >= min_dist or isSNP):
        if code[soap_base] % 5 == 0:
            # HOM
            if depth >= min_depth:
                print line.strip()
            else:
                pass
        else:
            # HET
            if  base1Depth >= min_depth and base2Depth >= min_depth:
                print line.strip()
            else:
                pass
    else:
        pass

