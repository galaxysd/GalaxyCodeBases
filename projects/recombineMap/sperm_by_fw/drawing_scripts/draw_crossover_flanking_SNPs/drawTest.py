'''
Created on 2012-7-8

@author: LiJinsen
'''
import sys, os, cPickle
from drawSVG import *;
'''
mySVG = SVG()
mySVG.addChildElement('rect',
                     {'x':20,
                      'y':40,
                      'width':80,
                      'height':50,
                      'fill': 'blue'})
mySVG.addChildElement('text',
                     {'x':22, 'y':12,},
                      "Test text")
mySVG.write('test.svg')
'''

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

if os.path.isfile(os.getcwd()+'tmp.pickle'):
    print "Loading File: tmp.pickle"
    data = cPickle.load(open("tmp.pickle", "r"))
else:
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
        for line in content.split('\n'):
            if line=='':
                continue
            tmp = line.split('\t')
            if tmp[0]=='#Chr':
                itemkey = tmp
                for item in itemkey:
                    tempdata.update({item:[]})
                continue
            for i in range(len(itemkey)):
                tempdata[itemkey[i]].append(tmp[i])
        data.update({tempdata[itemkey[0]][0]:tempdata})
    #cPickle.dump(data, open("tmp.pickle", "w"))
        
#read result file
filename = sys.argv[2]
print "Loading File:", filename
file = open(filename,'r')
content = file.read()
result = []
'''
S01    chr1    67236662    88414    "67192455,TC"    "67280869,GA"    1
S01    chr2    1120458    143548    "1048684,CA"    "1192232,TC"    1
S01    chr2    147107306    39677    "147087468,GA"    "147127145,TC"    1
S01    chr3    12810562    17438    "12801843,TC"    "12819281,CT"    1
'''
for line in content.split('\n'):
    if len(line)<1:
        continue
    result.append(line.split('\t'))

#find changed position and draw it
def findPosList(chrName, spermName, upPos, downPos):
    deepData = data[chrName]
    spermData = deepData[spermName]
    FaHap = deepData['FaHap']
    MoHap = deepData['MoHap']
    Pos = deepData['Pos']
    upIndex = Pos.index(upPos)
    downIndex = Pos.index(downPos)
    counti = 0
    countj = 0
    back_data = []
    forward_data = []
    i = upIndex
    j = downIndex
    while(counti<100 and i>=0):
        if spermData[i]!='N':
            counti+=1;
            if spermData[i]==FaHap[i]:
                k = 0
            else:
                k = 1
            back_data.append((Pos[i],k,FaHap[i]+MoHap[i]))
        i-=1
    while(countj<100 and j<len(spermData)):
        if spermData[j]!='N':
            countj+=1;
            if spermData[j]==FaHap[j]:
                k = 0
            else:
                k = 1
            forward_data.append((Pos[j],k,FaHap[j]+MoHap[j]))
        j+=1
    back_data.reverse()
    
    return_data1 = back_data
    return_data2 = []
    for tmp in range(upIndex+1, downIndex):
        if spermData[tmp]!='N':
            if spermData[tmp] == FaHap[tmp]:
                k = 0
            else:
                k = 1
            return_data2.append((Pos[tmp],k,FaHap[tmp]+MoHap[tmp]))

    return_data3 = forward_data
    
    return (return_data1, return_data2, return_data3)


w = open("Calibaration.txt",'w')
wFileWithCalib = open(sys.argv[2]+'.trust','w')
startx = 0
starty = 30

def drawingrect(posdata):
    global startx, starty
    for kitem in posdata:
        w.write(str(kitem[1]))
        col = 'blue'
        if kitem[1]==0:
            col = 'red'
        mySVG.addChildElement('rect', {'x': startx,
                                        'y':starty-6,
                                        'width':2,
                                        'height':10,
                                        'fill':col})
        startx += 3
    w.write(' ')
    startx += 3

def countFalse(posdata):
    countf = 0
    countm = 0
    for item in posdata:
        if item[1]==1:
            countf+=1
        else:
            countm+=1
    if countf+countm==0:
        return str(1)
    ratio = max(countf, countm)*1.0/(countf+countm)
    return str(ratio)

    
SVGfilename = ''



for item in result:
    if item[6]=='0':
        continue
    if SVGfilename=='':
        mySVG = SVG({'width':2000, 'height':5000})
        SVGfilename = item[0]
        
    if SVGfilename!=item[0]:
        print "Writing SVG...", SVGfilename+'.svg'
        mySVG.attributes['viewBox'] = '0 0 '+str(startx+5)+' '+str(starty+5)
        mySVG.write(SVGfilename+'.svg')
        mySVG = SVG({'width':2000, "Height":3000})
        starty = 30
        SVGfilename = item[0]

    print 'Analyzing', item[0], item[1]
    upPos = item[4].replace('"', '').split(',')[0]
    downPos = item[5].replace('"', '').split(',')[0]
    interval = int(item[3])
    if interval:
        beforedata, posdata, afterdata = findPosList(item[1], item[0], upPos, downPos)
        #if posdata==[]:
        #    continue
        startx = 10
        tmpstr = item[0]+' '+item[1]+' '+upPos+' '+item[3]
        for ii in range(len(item)-1):
            wFileWithCalib.write(item[ii])
            wFileWithCalib.write('\t')
        wFileWithCalib.write(countFalse(beforedata)+'\t')
        wFileWithCalib.write(countFalse(posdata)+'\t')
        wFileWithCalib.write(countFalse(afterdata)+'\t')
        wFileWithCalib.write('\n')
        mySVG.addChildElement('text', {'x':startx, 'y':starty}, tmpstr)
        startx = 250
        w.write(item[0]+'\t'+item[1]+'\t'+upPos+'\t'+item[3]+'\n')
        drawingrect(beforedata)
        boxstartx = startx
        drawingrect(posdata)
        mySVG.addChildElement('line', {'x1': boxstartx-3,
                                        'y1':starty-8,
                                        'x2':boxstartx-3,
                                        'y2':starty+7,
                                        'stroke':'black',
                                        'stroke-width':1})
        mySVG.addChildElement('line', {'x1': boxstartx-3,
                                        'y1':starty-8,
                                        'x2':startx-2,
                                        'y2':starty-8,
                                        'stroke':'black',
                                        'stroke-width':1})
        mySVG.addChildElement('line', {'x1': startx-2,
                                        'y1':starty-8,
                                        'x2':startx-2,
                                        'y2':starty+7,
                                        'stroke':'black',
                                        'stroke-width':1})
        mySVG.addChildElement('line', {'x1': boxstartx-3,
                                        'y1':starty+7,
                                        'x2':startx-2,
                                        'y2':starty+7,
                                        'stroke':'black',
                                        'stroke-width':1})
        drawingrect(afterdata)
        w.write('\n')
        starty += 20
w.close()
mySVG.write('test.svg')
    
            