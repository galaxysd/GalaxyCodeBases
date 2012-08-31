'''
Created on 2012-7-12

@author: LiJinsen
'''
import sys, os

filename = sys.argv[1]
filename2 = sys.argv[2]

#raw data
file = open(filename,'r')
content1 = file.read()
result1 = {}
for line in content1.split('\n'):
    if len(line)<1:
        continue
    item = line.split('\t')
    if not result1.has_key(item[0]):
        result1.update({item[0]:{}})
    if not result1[item[0]].has_key(item[1]):
        result1[item[0]].update({item[1]:[0]})
    result1[item[0]][item[1]][0]+=1
    result1[item[0]][item[1]].append(item)
    
#HMM data
file = open(filename2,'r')
content2 = file.read()
result2 = {}
for line in content2.split('\n'):
    if len(line)<1:
        continue
    item = line.split('\t')
    if not result2.has_key(item[0]):
        result2.update({item[0]:{}})
    if not result2[item[0]].has_key(item[1]):
        result2[item[0]].update({item[1]:[0]})
    result2[item[0]][item[1]][0]+=1
    result2[item[0]][item[1]].append(item)

#print result1
    
w = open('difference.txt','w')
w1 = open('differenceDetail.raw.txt','w')
w2 = open('differenceDetail.new.txt','w')
def writeData(item1, item2):
    for obj in item1:
        flag = False
        for tmpobj in item2:
            if obj[2]==tmpobj[2]:
                flag = True
                break
        if flag==True:
            continue
        for c in obj:
            w2.write(c)
            w2.write('\t')
        w2.write('\n')
    
for itemkey1 in result1.keys():
    for itemkey2 in result1[itemkey1].keys():
        if not result2[itemkey1].has_key(itemkey2):
            w.write(itemkey1)
            w.write('\t')
            w.write(itemkey2)
            w.write('\tNot In\n')
            continue
        if result1[itemkey1][itemkey2][0] != result2[itemkey1][itemkey2][0]:
            w.write(itemkey1)
            w.write('\t')
            w.write(itemkey2)
            w.write('\t')
            if result1[itemkey1][itemkey2] < result2[itemkey1][itemkey2]:
                w.write('More')
            else:
                w.write('Less')
            w.write('\n')
            writeData(result2[itemkey1][itemkey2][1:], result1[itemkey1][itemkey2][1:])
        