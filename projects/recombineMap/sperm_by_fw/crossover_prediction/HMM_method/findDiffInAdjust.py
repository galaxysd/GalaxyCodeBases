'''
Created on 2012-7-13

@author: LiJinsen
'''
content1 = open('differenceDetail.new.txt.adjust','r').read().split('\n')
content2 = open('differenceDetail.new.txt','r').read().split('\n')

data1 = []
for line in content1:
    if line!='' and line[-1]!='\t':
        line+='\t'
    if line!='' and line[-2]=='\t':
        line = line[:-2]
    data1.append(line)
    
w = open('bad.result.txt','w')
for line in content2:
    if line not in data1:
        w.write(line)
        w.write('\n')