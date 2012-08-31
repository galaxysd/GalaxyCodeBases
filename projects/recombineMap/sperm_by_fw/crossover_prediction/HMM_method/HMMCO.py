'''
Created on 2012-7-11

@author: LiJinsen
'''
import sys, os
from mpmath import mpf

##Baum-Welch Part
def Baum_Welch(obs, states, start_p, trans_p, emit_p):
    T = {}
    for state in states:
        T[state] = start_p[state]
    
    
    pass

##viterbi part
def forward_viterbi(obs, states, start_p, trans_p, emit_p):
    T = {}
    for state in states:
        ##    prob.    V.path    V.prob.
        T[state] = (start_p[state], [state], start_p[state])
    for output in obs:
        U = {}
        for next_state in states:
            total = mpf(0)
            argmax = None
            valmax = mpf(0)
            for source_state in states:
                (prob, v_path, v_prob) = T[source_state]
                prob = mpf(prob)
                v_prob = mpf(v_prob)
                tmp1 = emit_p[source_state][output]
                tmp2 = trans_p[source_state][next_state]
                p = tmp1 * tmp2
                prob *= p
                v_prob *= p
                total += prob
                if v_prob > valmax:
                    argmax = v_path + [next_state]
                    valmax = v_prob
            U[next_state] = (total, argmax, valmax)
        T = U
    ## apply sum/max to the final states:
    total = 0
    argmax = None
    valmax = 0
    for state in states:
        (prob, v_path, v_prob) = T[state]
        total += prob
        if v_prob > valmax:
            argmax = v_path
            valmax = v_prob
    return (total, argmax, valmax)

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
#diffFile = open('difference.txt','r').read().split('\n')
#needReRun = []
#for line in diffFile:
#    tmp = line.split('\t')
#    needReRun.append(tmp)
for filename in filename1:
    print "Loading File:", filename
    content = open(filename,'r').read()
    tempdata = {}
    posData = {}
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
                posData.update({item:[]})
            continue
        for i in range(len(itemkey)):
            if i>=5:
                if tmp[i]==tmp[2]:
                    tmp[i]='1'#father is 1
                if tmp[i]==tmp[3]:
                    tmp[i]='-1'#mother is -1
                if tmp[i]=='N':
                    continue
            tempdata[itemkey[i]].append(tmp[i])
            posData[itemkey[i]].append(tmp[1])
    chrName = tempdata['#Chr'][0]
    transition_probability = {'F' : {'F': 0.9999998,'U': 0.00000009,'M': 0.00000002},
                              'M' : {'M': 0.9999998,'U': 0.00000009,'F': 0.00000002},
                              'U' : {'U': 0.9, 'F': 0.05, 'M': 0.05}}
    emission_probability = {'F' : {'1': 0.9, '-1': 0.1},
                            'M' : {'-1': 0.9, '1': 0.1},
                            'U' : {'1':0.5, '-1':0.5}}
    start_probability = {'F' : 0.3,
                         'M' : 0.3,
                         'U' : 0.4}
    for item in itemkey:
        if item[0]=='S' and item[1]!='p':
            #if [item, chrName, 'More'] not in needReRun:
            #    continue
            states = ('F','U','M')
            observations = tempdata[item]
            print 'Calculating',item
            tmpdt = forward_viterbi(observations, states, start_probability, transition_probability, emission_probability)
            print 'Calculating',item, 'Finished!'
            print 'Writing File',"result."+item+'.'+chrName
            w = open("result."+item+'.'+chrName,'w')
            w.write('#Sperm\tChr\tPos\tObservation\tGuess\n')
            posTemp = posData[item]
            for i in range(len(tmpdt[1])-1):
                w.write(item)
                w.write('\t')
                w.write(chrName)
                w.write('\t')
                w.write(posTemp[i])
                w.write('\t')
                w.write(observations[i])
                w.write('\t')
                w.write(tmpdt[1][i])
                w.write('\n')
            w.close()
            print 'Writing Completed.'
            