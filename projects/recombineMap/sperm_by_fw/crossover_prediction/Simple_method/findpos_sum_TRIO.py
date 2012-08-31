import sys, os

filename = sys.argv[1]

content = open(filename,'r').read()


data = []
#data[i][j] is the 
	#j#		#j#			#j#
#i##Chr		Pos			FaHap	MoHap	SpNum	S01	S02	S03	S04	S05	S06	S07	S08	S09	S10	S11	S12	S13	S14	S15	S16	S17	S18	S19	S20	S21	S22	S23	S24	S25	S26	S27	S28	S29	S30	S31	S32	S33	S34	S35	S36	S37	S38	S39	S40	S41	S42	S43	S44	S45	S46	S47	S48	S49	S50	S51	S52	S53	S54	S55	S56	S57	S58	S59	S60	S61	S62	S63	S64	S65	S66	S67	S68	S69	S70	S71	S72	S73	S74	S75	S76	S77	S78	S79	S80	S81	S82	S83	S84	S85	S86	S87	S88	S89	S90	S91	S92	S93	S94	S95	S96	S97	S98	S99
#i##chr13	19020145	T		G		9		N	T	N	N	N	N	N	N	N	N	N	N	N	G	N	N	N	N	N	N	N	G	N	G	N	N	N	N	N	T	N	N	N	N	N	N	N	N	N	N	N	N	N	N	G	N	N	N	N	T	N	N	N	G	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	N	G	N	N	N	N	N	N
#i##chr13	19020627	G		A		39		N	G	N	A	G	N	N	N	G	N	A	G	A	A	A	N	N	G	N	N	N	A	N	A	G	G	N	G	N	G	G	G	N	N	N	G	N	G	N	G	N	N	A	N	N	N	N	N	N	G	A	A	N	A	N	N	G	N	N	A	N	N	N	G	N	N	A	G	G	N	N	A	N	N	N	G	N	N	N	N	N	N	A	N	N	N	N	A	N	N	N	N	A	G	N	G	N	N	N
#i##chr13	19020776	G		T		24		N	G	N	T	G	N	N	N	N	N	T	G	N	T	T	T	N	N	N	G	N	T	N	T	G	N	N	N	G	G	N	N	N	N	N	G	N	N	N	N	N	N	N	N	N	N	T	N	N	N	N	N	N	N	N	N	G	N	N	T	N	N	N	N	N	N	T	G	N	T	N	N	T	N	N	N	N	N	N	N	N	G	N	N	N	N	N	T	N	N	N	N	N	N	N	N	N	N	N

for line in content.split('\n'):
	data.append(line.split('\t'))
	
last = [-1]*99
now = [-1]*99
pos = []
ma = []
macursor = []

for i in range(99):
	pos.append([])
	ma.append([])
	macursor.append([])



for i in range(1,len(data)):
	for j in range(5,len(data[i])):
		temp = data[i][j]
		if temp=='N':
			continue
		fahap = data[i][2]
		mohap = data[i][3]
		if len(macursor[j-5])==10:
			del macursor[j-5][0]
		if temp==fahap:
			macursor[j-5].append(100)
			now[j-5] = 0
		if temp==mohap:
			macursor[j-5].append(0)
			now[j-5] = 1
		ma[j-5].append(sum(macursor[j-5])/len(macursor[j-5]))
		pos[j-5].append(data[i][1])
		'''if last[j-5]==-1 and now[j-5]!=-1:
			last[j-5] = now[j-5]
			continue;
		if now[j-5]==-1:
			continue;
		if last[j-5]!=now[j-5]:
			recombinationpos[j-5].append(data[i][1])
			last[j-5]=now[j-5]'''

w = open(data[1][0]+'.changepos.sum.trio','w')
j = 0 #sperm number.
recombinepos = []

for sperm in ma:
	start = sperm[0]
	#w.write("S"+str(j+1)+"\t")
	count = 0
	total = len(sperm)
	for l in range(1,len(sperm)):
		if sperm[l]==100-start:
			count += 1
			if count<1000 and l<=total-20:
				continue
			k = l
			while(sperm[k]!=start):
				k-=1
			recombinepos.append([int(pos[j][k]),j+1])
			#w.write(pos[j][k])
			#w.write('\t')
			start = 100 - start
			count=0
	#w.write('\n')
	j+=1
	
def mycmp(E1, E2):
	return cmp(E1[0],E2[0])
recombinepos.sort(mycmp)
w.write("window\t[pos,sperm]")
interval = 10000
startpos = -interval
for l in recombinepos:
	if l[1]==27 or l[1]==68 or l[1]==87:
		continue
	if l[1]==41:
		if data[1][0]=='chr13' or data[1][0]=='chr14' or data[1][0]=='chr16':
			continue
	if l[1]==63:
		if data[1][0] =='chr14':
			continue
	if l[1]==65 and data[1][0]=='chr6':
		continue
	if l[1]==97:
		if data[1][0]=='chr10' or data[1][0]=='chr16' or data[1][0]=='chr8' or data[1][0]=='chr1':
			continue
		
	if l[0]>startpos and l[0]<startpos+interval:
		w.write(str(l[0]))
		w.write(',')
		w.write("S"+str(l[1]))
		w.write('\t')
	else:
		w.write('\n')
		startpos=l[0]/interval*interval
		w.write(str(startpos)+'-'+str(startpos+interval)+'\t')
		w.write(str(l[0]))
		w.write(',')
		w.write("S"+str(l[1]))
		w.write('\t')


	
w.close()

'''
for sperm in recombinationpos:
	print "Sperm S",i,
	print sperm
	print '\n'
	i+=1
'''	