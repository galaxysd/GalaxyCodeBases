import sys, os

filename = sys.argv[1]

file = open(filename,'r')
content = file.read()

result = []
firstline = ''
for line in content.split('\n'):
	if len(line)<1:
		continue
	result.append(line.split('\t'))
	
def mycmp(a,b):
	if a[0]==b[0]:
		if a[1]==b[1]:
			return cmp(int(a[2]),int(b[2]))
		else:
			if a[1]=='chrX':
				return 1
			if b[1]=='chrX':
				return -1
			if len(a[1])==len(b[1]):
				return cmp(a[1],b[1])
			else:
				return cmp(len(a[1]),len(b[1]))
	else:
		return cmp(a[0],b[0])

result.sort(mycmp)

file.close()
w = open(filename,'w')
for line in result:
	print line
	for item in line:
		if item=='':
			continue
		w.write(item)
		w.write('\t')
	w.write('\n')
	
