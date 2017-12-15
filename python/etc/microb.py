#!/usr/bin/env python3

import vcf
import gzip,csv
import pprint

SNPdata = {}
SNPchrID = set()

ChrIDs ={
"DS362249":"DS362249.1", "DS362248":"DS362248.1", "DS362247":"DS362247.1", "DS362246":"DS362246.1", "DS362245":"DS362245.1", "DS362244":"DS362244.1", "DS362243":"DS362243.1", "DS362242":"DS362242.1", "DS362241":"DS362241.1", "DS362240":"DS362240.1", "DS362239":"DS362239.1", "DS362238":"DS362238.1", "DS362237":"DS362237.1", "DS362236":"DS362236.1", "DS362235":"DS362235.1", "DS362234":"DS362234.1", "DS362233":"DS362233.1", "DS362232":"DS362232.1", "DS362231":"DS362231.1", "DS362230":"DS362230.1", "DS362229":"DS362229.1", "DS362228":"DS362228.1", "DS362227":"DS362227.1", "DS362226":"DS362226.1", "DS362225":"DS362225.1", "DS362224":"DS362224.1", "DS362223":"DS362223.1", "DS362222":"DS362222.1", "DS362221":"DS362221.1", "DS362220":"DS362220.1", "DS362219":"DS362219.1", "DS362218":"DS362218.1", "DS362217":"DS362217.1", "DS499581":"DS499581.1", "DS499580":"DS499580.1", "DS499579":"DS499579.1", "DS499578":"DS499578.1", "DS499577":"DS499577.1", "DS499576":"DS499576.1", "DS499575":"DS499575.1", "DS499574":"DS499574.1", "DS499573":"DS499573.1", "DS499572":"DS499572.1", "DS499571":"DS499571.1", "DS499570":"DS499570.1", "CP001107":"CP001107.1", "CP001104":"CP001104.1", "CP001105":"CP001105.1", "CP001106":"CP001106.1"
}

SNPcnt = 0
f = gzip.open('snpcandidatforpcr.out.gz', 'rt')
tsvin = csv.reader(f, delimiter='\t')
for row in tsvin:
    if len(row)==0: continue
    #print(row)
    ChrID = row[0].split('.')[1]
    thePos = row[2]
    #theRef = row[3]
    Alleles = row[4].split(',')
    Allele1 = Alleles[0].split('|')
    Allele2 = Alleles[1].split('|')
    AlleleP = { Allele1[0]:sum(map(int,Allele1[2:])), Allele2[0]:sum(map(int,Allele2[2:])) }
    if ChrID in ChrIDs:
        theKey = '\t'.join((ChrIDs[ChrID],thePos))
        SNPdata[theKey] = AlleleP
        SNPchrID.add(ChrIDs[ChrID])
        #print((ChrID,thePos,theRef,theKey,AlleleP))
        #print(Allele1)
        #print(Allele2)
        SNPcnt += 1
f.close()
#pprint.pprint(SNPchrID)
print('[!]SNP read:%d.' % (SNPcnt))

SNPcnt = 0
DepthStat = {}
GTStat = {}
vcf_reader = vcf.Reader(filename='microb.vcf.gz',compressed=True)
for theChrid in sorted(SNPchrID):
    for record in vcf_reader.fetch(theChrid):
        theKey = '\t'.join((record.CHROM,str(record.POS)))
        #print((theKey))
        if theKey in SNPdata:
            SNPcnt += 1
            #print(record)
            theBases = [record.REF]
            for abp in record.ALT:
                theBases.append(str(abp))
            SNPdatastr = []
            for abp in SNPdata[theKey]:
                SNPdatastr.append( ''.join([abp,':',str(SNPdata[theKey][abp]) ]) )
            print('%4d   %s\t%s\t%s'%(SNPcnt, theKey.replace('\t',':'), ','.join(SNPdatastr),','.join(theBases)))
            #print(SNPdata[theKey])
            #print(record.genotype('ERR589860'))
            SampleGTs = set()
            Depths = []
            for sample in record.samples:
                #print(sample)
                theADstr = ','.join(map(str, sample['AD'] if isinstance(sample['AD'], list) else list(str(sample['AD']))))
                """
                if type(sample['AD']) is list:
                    theADstr = ','.join(map(str,sample['AD']))
                elif type(sample['AD']) is int:
                    theADstr = str(sample['AD'])
                else:
                    theADstr = 'ERROR'
                """
                print('%s:%s\t%4d %s=%s'%(sample.sample,sample.gt_bases,sample['DP'],sample['GT'],theADstr))
                Depths.append(sample['DP'])
                if sample.gt_bases != None:
                    SampleGTs.add(sample.gt_bases)
            if len(Depths):
                AvgDepth = sum(Depths)/len(Depths)
                theDepth = int(round(AvgDepth,-1))
                DepthStat[theDepth] = 1 + DepthStat.setdefault(theDepth,0)
                GTStat[len(SampleGTs)] = 1 + GTStat.setdefault(len(SampleGTs),0)
            print((len(SampleGTs),AvgDepth,theDepth))
            print((GTStat,DepthStat))
