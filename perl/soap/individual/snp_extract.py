#!/usr/bin/python
import sys
import re

try:
    soap_cns = open(sys.argv[1], "r")
except IndexError:
    print >>sys.stderr, "python", sys.argv[0], "SoapCNS"
    sys.exit(1)

regexp = re.compile("([ATGCXN]):(\d+)\*(\d+)")
prev_pos = -1000000;
prev_info = ""
discreptancy = []
for line in soap_cns:
    tabs = line.split()
    if tabs[3] != tabs[2] and tabs[2] != "N":
        dist = int(tabs[1])-prev_pos
        if dist < 0: # New chr
            dist = 1000000
        try:
            if discreptancy[5] > dist and discreptancy[0] == tabs[0]:
                discreptancy[5] = dist
            else:
                pass
            print >>sys.stdout, prev_info+"\t"+str(discreptancy[5])
            discreptancy = []
        except IndexError:
            if len(discreptancy) == 0:
                pass
            else:
                print >>sys.stderr, "Index Error in discreptancies. Exit 255"
                sys.exit(255)
        discreptancy=[tabs[0],tabs[1],tabs[2],tabs[3], tabs[12], dist]
        prev_pos = int(tabs[1])
        prev_info = line.rstrip()
    else:
        pass
soap_cns.close()
print >>sys.stdout, prev_info+"\t"+str(discreptancy[5])
