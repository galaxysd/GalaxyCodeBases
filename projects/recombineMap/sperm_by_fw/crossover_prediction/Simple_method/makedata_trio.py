import sys, os

filename1 = []
def mydir(arg, dirname, names):
	files = [os.path.normpath(os.path.join(dirname, file)) for file in names]
	for filename in files:
		if filename.find("score20")!=-1:
			filename1.append(filename)


	
if len(sys.argv)==1:
	path=os.getcwd()
else:
	path = sys.argv[1]
os.path.walk(path, mydir, 0)
w = open('all.sum.changepos.standard.trio','w')
result = []
w.write("#SpermID\tChr\tPos\tInterval\tUpSNP\tDownSNP\tIsCurated\n");
for filename in filename1:
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
	snp = []

	for i in range(99):
		pos.append([])
		ma.append([])
		macursor.append([])
		snp.append({})

	for i in range(1,len(data)):
		for j in range(5,len(data[i])):
			temp = data[i][j]
			snp[j-5].update({data[i][1]:data[i][2]+data[i][3]})
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

	
	j = 0 #sperm number.
	start = 0#start finding pos
	for sperm in ma:
		start = sperm[0]
		count = 0
		total = len(sperm)
		eachsperm = []
		thechangepos = 0
		for l in range(1,len(sperm)):
			if sperm[l]==start and count==0:
				thechangepos = l
			if sperm[l]==100-start:
				count += 1
				if count<50 and l<=total-20:
					continue
				#k = l
				#while(sperm[k]!=start):
				#	k-=1
				k = thechangepos
				if j<9:
					eachsperm.append("S0"+str(j+1))
				else:
					eachsperm.append("S"+str(j+1))
				eachsperm.append(data[1][0])
				eachsperm.append(str( (int(pos[j][k])+int(pos[j][k+1]))/2 ))
				eachsperm.append(str( int(pos[j][k+1])-int(pos[j][k])) )
				eachsperm.append(pos[j][k]+','+snp[j][pos[j][k]])
				eachsperm.append(pos[j][k+1]+','+snp[j][pos[j][k+1]])
				'''for position in range(start, len(snp[j])):
					item = snp[j][position]
					if item[0]==pos[j][k-1]:
						eachsperm.append(pos[j][k-1]+','+item[1])
						continue
					if item[0]==pos[j][k]:
						eachsperm.append(pos[j][k]+','+item[1])
						start = position
						break'''
				start = 100 - start
				count=0
				print eachsperm
				for item in eachsperm:
					w.write(item)
					w.write("\t")
				w.write('\n')
				eachsperm = []
	
		#result.append(eachsperm)
		j+=1
'''
def mycmp(E1, E2):
	if E1[0]>E2[0]:
		return 1
	elif E1[0]==E2[0]:
		if E1[1]==E2[1]:
			return cmp(int(E1[2]),int(E2[2]))
	return -1
#result.sort(mycmp)

for line in result:
	for item in line:
		w.write(item)
		w.write("\t")
	w.write("\n")
'''