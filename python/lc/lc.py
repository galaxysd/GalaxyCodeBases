#!/usr/bin/env python3
import site,os; site.addsitedir('/'.join((os.path.expanduser("~"),'git/toGit/python/Gsite')))
#import sys; print('\n'.join(sys.path))
import Galaxy.asizeof

import csv,sys,math
from collections import OrderedDict
from collections import defaultdict

naStr = 'NA'
ZoneNum = 50.0

import pprint
from xopen import xopen

class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

print(sys.argv)

MethRate = Vividict()
MethRateCnt = 0
with xopen('meth.cg.rate.bz2','rt') as f:
#with open('meth.cg.rate','rt') as f:
    tsvin = csv.reader(f, delimiter='\t')
    for row in tsvin:
        if row[0].startswith('#'):
            continue
        #print(row)
        #Rate = (int(row[3])+int(row[6]))/(int(row[4])+int(row[7]))
        MethRate[row[0]][int(row[1])] = float(row[9])
        MethRateCnt += 1
        #pprint.pprint(MethRate)
print(Galaxy.asizeof.asizeof(MethRate))

#CDSdat = OrderedDict()
CDSdatCnt = 0
CDSdat = defaultdict(set)

with xopen('CDS.bed.gz','rt') as f:
#with open('CDS.bed','rt') as f:
    tsvin = csv.reader(f, delimiter='\t')
    for row in tsvin:
        ZoneLength = 1 + int(row[2]) - int(row[1])
        #print((ZoneLength,row))
        theKey = '\t'.join(row[0:3])
        #CDSdat[theKey]=','.join(filter(None,(CDSdat.setdefault(theKey),row[3])))
        #CDSdat[theKey] |= {row[3]}
        CDSdat[theKey].add(row[3])
        CDSdatCnt += 1
        #pprint.pprint(CDSdat)
print(Galaxy.asizeof.asizeof(CDSdat))
print((MethRateCnt,CDSdatCnt,len(MethRate),len(CDSdat)))

with open('meth.zones','wt') as outf:
    for k, v in CDSdat.items():
        v = ','.join(sorted(v))
        (Chrid,pLeft,pRight) = k.split('\t')
        pLeft = int(pLeft)
        pRight = int(pRight)
        #print((k,v,CDSdat[k],Chrid,pLeft,pRight))
        ZoneLen = 1 + pRight - pLeft
        ZoneStep = math.ceil(float(ZoneLen)/ZoneNum)
        ZoneCnt = math.ceil(ZoneLen/ZoneStep)
        ZoneValues = []
        ZoneValueCntV = ZoneCnt
        startPos = pLeft
        while startPos <= pRight:
            endPos = startPos + ZoneStep
            if endPos > pRight: endPos = pRight
            ZoneSum = 0
            ZoneHits = 0
            for p in range(startPos,endPos):
                if p in MethRate[Chrid]:
                    ZoneSum += MethRate[Chrid][p]
                    ZoneHits += 1
            if ZoneHits == 0:
                ZoneValues.append(naStr)
                ZoneValueCntV -= 1
            else:
                ZoneValue1 = float(ZoneSum/ZoneHits)
                ZoneValues.append("%.8f" % ZoneValue1)
            startPos += ZoneStep
        ZoneValues = [v,k,','.join(map(str,(ZoneStep,ZoneCnt,ZoneValueCntV)))] + ZoneValues
        if ZoneValueCntV > 0:
            print('\t'.join(ZoneValues),file=outf,flush=True)

#for theKey in CDSdat.keys():
#    print(theKey,CDSdat[theKey])
'''
scp t.py blc4:/ldfssz1/FGI/users/huxuesong/t/lc/

f=magic.from_file('CDS.bed.gz')
print(f)
f=magic.from_file('meth.cg.rate.bz2')
print(f)
f=magic.from_file('adapter.fa')
print(f)

15615200
81835576
(164859, 202393, 25, 198788)
'''
