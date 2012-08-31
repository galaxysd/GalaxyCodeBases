'''
Created on 2012-7-12

@author: LiJinsen
'''

import sys, os

#read data file
filename1 = []
def mydir(arg, dirname, names):
    files = [os.path.normpath(os.path.join(dirname, file)) for file in names]
    for filename in files:
        if filename.find("FaMo")!=-1:
            filename1.append(filename)

if len(sys.argv)==1:
    path=os.getcwd()
else:
    path = sys.argv[1]
os.path.walk(path, mydir, 0)
data = {}

for filename in filename1:
    print "Loading File:", filename
    content = open(filename,'r').read()
    tempdata = {}
    '''
    #Chr        Pos            FaHap    MoHap    SpNum    S01    S02    S03    S04    S05    S06    S07    S08    S09    S10    S11    S12    S13    S14    S15    S16    S17    S18    S19    S20    S21    S22    S23    S24    S25    S26    S27    S28    S29    S30    S31    S32    S33    S34    S35    S36    S37    S38    S39    S40    S41    S42    S43    S44    S45    S46    S47    S48    S49    S50    S51    S52    S53    S54    S55    S56    S57    S58    S59    S60    S61    S62    S63    S64    S65    S66    S67    S68    S69    S70    S71    S72    S73    S74    S75    S76    S77    S78    S79    S80    S81    S82    S83    S84    S85    S86    S87    S88    S89    S90    S91    S92    S93    S94    S95    S96    S97    S98    S99
    #chr13    19020145    T        G        9        N    T    N    N    N    N    N    N    N    N    N    N    N    G    N    N    N    N    N    N    N    G    N    G    N    N    N    N    N    T    N    N    N    N    N    N    N    N    N    N    N    N    N    N    G    N    N    N    N    T    N    N    N    G    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    N    G    N    N    N    N    N    N
    #chr13    19020627    G        A        39        N    G    N    A    G    N    N    N    G    N    A    G    A    A    A    N    N    G    N    N    N    A    N    A    G    G    N    G    N    G    G    G    N    N    N    G    N    G    N    G    N    N    A    N    N    N    N    N    N    G    A    A    N    A    N    N    G    N    N    A    N    N    N    G    N    N    A    G    G    N    N    A    N    N    N    G    N    N    N    N    N    N    A    N    N    N    N    A    N    N    N    N    A    G    N    G    N    N    N
    #chr13    19020776    G        T        24        N    G    N    T    G    N    N    N    N    N    T    G    N    T    T    T    N    N    N    G    N    T    N    T    G    N    N    N    G    G    N    N    N    N    N    G    N    N    N    N    N    N    N    N    N    N    T    N    N    N    N    N    N    N    N    N    G    N    N    T    N    N    N    N    N    N    T    G    N    T    N    N    T    N    N    N    N    N    N    N    N    G    N    N    N    N    N    T    N    N    N    N    N    N    N    N    N    N    N
    '''
    itemkey = []
    chrName = ''
    for line in content.split('\n'):
        if line=='':
            continue
        tmp = line.split('\t')
        if tmp[0]=='#Chr':
            itemkey = tmp[:4]
            continue
        if chrName=='':
            chrName = tmp[0]
        tempdata.update({tmp[1]:tmp[2]+tmp[3]})
            
    data.update({chrName:tempdata})

#read result file
filename1 = []
def mydir2(arg, dirname, names):
    files = [os.path.normpath(os.path.join(dirname, file)) for file in names]
    for filename in files:
        if filename.find("result")!=-1:
            filename1.append(filename)

path = 'data'
allSNP = 0
mayErrorSNP = 0
unsolvedSNP = 0
os.path.walk(path, mydir2, 0)
w = open('Result.all.txt.TRASH','w')
for filename in filename1:
    print 'Reading File:', filename
    chrName = filename.split('.')[2]
    content = open(filename,'r').read().split('\n')
    '''
    #Sperm    Chr    Pos    Observation    Guess
    S86    chrX    62304    1    F
    S86    chrX    173892    1    F
    S86    chrX    178135    1    F
    '''
    count = {'F':0,'M':0, 'U':0}
    last = ''
    change = []
    for line in content:
        item = line.split('\t')
        if item[0]=='' or item[0][0]=='#':
            continue
        if (item[4]=='F' and item[3]=='0') or (item[4]=='M' and item[3]=='1'):
            mayErrorSNP+=1
        if item[4]=='U':
            unsolvedSNP+=1
        allSNP+=1
        if last=='':
            last = item
            count[last[4]]=1
        else:
            count[item[4]]+=1
            if item[4]=='U' and last[4]!='U':
                change = last
            if last[4]=='U' and item[4]!='U':
                if change!=[] and change[4]!=item[4]:
                    w.write(change[0])
                    w.write('\t')
                    w.write(change[1])
                    w.write('\t')
                    w.write(str((int(change[2])+int(item[2]))/2))
                    w.write('\t')
                    w.write(str( int(item[2])-int(change[2]) ))
                    w.write('\t')
                    w.write(change[2]+','+data[chrName][change[2]])
                    w.write('\t')
                    w.write(item[2]+','+data[chrName][item[2]])
                    w.write('\t1\t')
                    w.write(str(count['U']))
                    w.write('\t\n')
                count = {'F':0,'M':0, 'U':0}
                change = []
            if item[4]!=last[4] and item[4]!='U' and last[4]!='U':
                w.write(last[0])
                w.write('\t')
                w.write(last[1])
                w.write('\t')
                w.write(str( (int(last[2])+int(item[2]))/2 ))
                w.write('\t')
                w.write(str( int(item[2])-int(last[2]) ))
                w.write('\t')
                w.write(last[2]+','+data[chrName][last[2]])
                w.write('\t')
                w.write(item[2]+','+data[chrName][item[2]])
                w.write('\t1\t0\t\n')
                count = {'F':0,'M':0, 'U':0}
            last = item
            

print "All SNP", allSNP
print 'Error SNP', mayErrorSNP
print 'Illegible SNP', unsolvedSNP
    
        