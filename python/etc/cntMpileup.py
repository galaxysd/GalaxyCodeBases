#!/usr/bin/env pypy3

import lz4.frame
import os
import csv
import itertools
import sys
import collections

CoreKeys = ["{}{}".format(x,y) for x,y in itertools.product(['fwd','rev'], ['A','C','G','T','N'])]
CounterKeys = CoreKeys + ['cntIns','cntDel'];
theSamples = ['B7','D3','Normal'];

def main():
	if len(sys.argv) < 3 :
		print('Usage:',sys.argv[0],'<samtools.mpileup.lz4> <BED file>',file=sys.stderr,flush=True);
		exit(0);
	inDepthFile = sys.argv[1]
	inBEDFile = sys.argv[2]

	theZone = collections.defaultdict(dict)
	with open(inBEDFile, newline='') as csvfile:
		tsvin = csv.DictReader(csvfile, delimiter='\t', fieldnames=['ChrID','Start','End'])
		for row in tsvin:
			theZone[row['ChrID']][row['Start']] = 1
	#print(theZone)
	#print("# Samples: [{}]".format(','.join(theSamples)))
	print('# Format: {}'.format(','.join(CoreKeys)))
	print("\t".join(['ChrID','Pos']+theSamples))
	with lz4.frame.open(inDepthFile, newline='', mode='rt') as fp:
		tsvin = csv.DictReader(fp, delimiter='\t', fieldnames=['ChrID','Pos','Ref'] + 
			["{}_{}".format(x,y) for x,y in itertools.product(theSamples, ['Depth','Bases','Qs'])]
		)
		sbc = {key:collections.Counter() for key in theSamples}
		bc = {key:{} for key in theSamples}
		try:
			for row in tsvin:
				#print(row)
				if row['ChrID'] in theZone:
					if row['Pos'] in theZone[row['ChrID']]:
						#print(row)
						print("\t".join((row['ChrID'],row['Pos'])),end='\t')
						for k in theSamples:
							bc[k] = Base_Counter(row['{}_Depth'.format(k)],row['{}_Bases'.format(k)])
							sbc[k] += {x: bc[k][x] for x in CounterKeys}
							#print(bc[k])
							print(",".join( [str(bc[k][x]) for x in CoreKeys] ),end='\t')
						print(flush=True);
		except KeyboardInterrupt:
			print('\n#[!]Ctrl+C pressed.',file=sys.stderr,flush=True)
			sys.stdout.flush()
			pass
		#print(sbc);

# Base Counting Subroutine *[Completed]
# https://github.com/photonchang/allelecount/blob/master/allelecount.py
def Base_Counter(inDepth,inBases):
	# Cleaning up Base String + Indel Counting
	CleanString = ''
	countIn = 0
	countDel = 0
	IndelHolder = []
	IndelDeterminant = 0
	CleanBool = False
	
	for currentIndex, Strholder in enumerate(inBases):
		# Skipping of '^' Signage
		if CleanBool == True:
			CleanBool = False
			continue
		
		if Strholder == '^':
			CleanBool = True
			continue
		
		# Skipping Indel
		if IndelDeterminant > 0:
			IndelDeterminant -= 1
			continue
		
		if Strholder == '+':
			countIn += 1
			
			# Determining Length of Indel
			# Since Illumina NGS has an upper limit of less than 999bp read length
			IndelDeterminant = 0
			
			if (currentIndex + 4) <= len(inBases):
				if inBases[currentIndex + 1: currentIndex + 4].isnumeric() == True:
					IndelDeterminant = int(inBases[currentIndex + 1: currentIndex + 4]) + 3
				elif inBases[currentIndex + 1: currentIndex + 3]. isnumeric() == True:
					IndelDeterminant = int(inBases[currentIndex + 1: currentIndex + 3]) + 2
				elif inBases[currentIndex + 1: currentIndex + 2].isnumeric() == True:
					IndelDeterminant = int(inBases[currentIndex + 1: currentIndex + 2]) + 1
			elif (currentIndex + 3) <= len(inBases):
				if inBases[currentIndex + 1: currentIndex + 3]. isnumeric() == True:
					IndelDeterminant = int(inBases[currentIndex + 1: currentIndex + 3]) + 2
				elif inBases[currentIndex + 1: currentIndex + 2].isnumeric() == True:
					IndelDeterminant = int(inBases[currentIndex + 1: currentIndex + 2]) + 1
			
			IndelHolder.append(inBases[currentIndex:currentIndex + IndelDeterminant + 1])
			continue
		
		if Strholder == '-':
			countDel += 1
			
			# Determining Length of Indel
			# Since Illumina NGS has an upper limit of less than 999bp read length
			IndelDeterminant = 0
			
			if (currentIndex + 4) <= len(inBases):
				if inBases[currentIndex + 1: currentIndex + 4].isnumeric() == True:
					IndelDeterminant = int(inBases[currentIndex + 1: currentIndex + 4]) + 3
				elif inBases[currentIndex + 1: currentIndex + 3]. isnumeric() == True:
					IndelDeterminant = int(inBases[currentIndex + 1: currentIndex + 3]) + 2
				elif inBases[currentIndex + 1: currentIndex + 2].isnumeric() == True:
					IndelDeterminant = int(inBases[currentIndex + 1: currentIndex + 2]) + 1
			elif (currentIndex + 3) <= len(inBases):
				if inBases[currentIndex + 1: currentIndex + 3]. isnumeric() == True:
					IndelDeterminant = int(inBases[currentIndex + 1: currentIndex + 3]) + 2
				elif inBases[currentIndex + 1: currentIndex + 2].isnumeric() == True:
					IndelDeterminant = int(inBases[currentIndex + 1: currentIndex + 2]) + 1
			
			IndelHolder.append(inBases[currentIndex: currentIndex + len(str(IndelDeterminant)) + 1])
			continue
		
		CleanString += Strholder
	
	else:
	# Transferring Back Cleaned String
		inBases = CleanString
		
		# '$' Signage Stripping
		inBases = inBases.replace('$', '')
		
	# Base Count Var Initialization
	bigA = 0
	bigC = 0
	bigG = 0
	bigT = 0
	smallA = 0
	smallC = 0
	smallG = 0
	smallT = 0
	delBase = 0
	
	# Base Counting
	bigA = inBases.count('A')
	bigC = inBases.count('C')
	bigG = inBases.count('G')
	bigT = inBases.count('T')
	bigN = inBases.count('N')
	
	smallA = inBases.count('a')
	smallC = inBases.count('c')
	smallG = inBases.count('g')
	smallT = inBases.count('t')
	smallN = inBases.count('n')
	
	delBase = inBases.count('*')
	
	# Internal Check - Throws Out Error (NOT STD-IN/OUT Compatible: Should break pipeline)
	InternalCounter = 0
	InternalCounter = bigA + bigC + bigG + bigT + bigN + smallA + smallC + smallG + smallT + smallN + delBase
	if InternalCounter != int(inDepth):
		print('Error:')
		print('Reported count: ' + inDepth)
		print('Internal counter: ' + str(InternalCounter))
		print('Internal sum: ' + str(len(inBases)))
		print('Number of insertions: ' + str(countIn))
		print('Number of deleions: ' + str(countDel))
		print('Post-processed bases: ' + inBases)
		#sys.exit()
	
	# Indel Compilation
	IndelSetDict = set(IndelHolder)
	tmpIndelString = ''
	FinalIndelHolder = []
	
	for EveryIndel in IndelSetDict:
		tmpIndelString = ''
		tmpIndelString = str(IndelHolder.count(EveryIndel)) + ":" + EveryIndel
		FinalIndelHolder.append(tmpIndelString)
	
	# Return Output
	FinalOutput = collections.OrderedDict.fromkeys( CounterKeys + ['Indels'] )
	FinalOutput['fwdA'] = bigA; FinalOutput['fwdC'] = bigC; FinalOutput['fwdG'] = bigG; FinalOutput['fwdT'] = bigT; FinalOutput['fwdN'] = bigN;
	FinalOutput['revA'] = smallA; FinalOutput['revC'] = smallC; FinalOutput['revG'] = smallG; FinalOutput['revT'] = smallT; FinalOutput['revN'] = smallN;
	FinalOutput['cntIns'] = countIn; FinalOutput['cntDel'] = countDel; FinalOutput['Indels'] = ';'.join(FinalIndelHolder); 	
	return FinalOutput

if __name__ == "__main__":
	main()  # time python3 ./samdepthplot.py t.tsv.gz 1

'''
./cntMpileup.py z.mp9k.lz4 z.bed
'''
