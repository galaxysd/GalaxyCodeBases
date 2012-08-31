import sys, os

changeposfile = []
def mydir(arg, dirname, names):
	files = [os.path.normpath(os.path.join(dirname, file)) for file in names]
	for filename in files:
		if filename.find('changepos.trio')!=-1:
			changeposfile.append(open(filename,"r"))


	
if len(sys.argv)==1:
	path=os.getcwd()
else:
	path = sys.argv[1]
	
os.path.walk(path, mydir, 0)

w = open("sum.changepos.trio","w")
sum = {}

for file in changeposfile:
	for line in file:
		item = line.split('\t')
		if not sum.has_key(item[0]):
			sum.update({item[0]:0})
		print item, len(item)
		sum[item[0]]+=len(item)-2
items = sum.items()
items.sort()
for item in items:
	w.write(str(item[0]))
	w.write('\t')
	w.write(str(item[1]))
	w.write('\n')
	
print items