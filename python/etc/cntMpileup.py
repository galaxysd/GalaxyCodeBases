#!/usr/bin/env pypy3

import lz4.frame
import os
import csv
import itertools
import sys
import collections

CounterKeys = ["{}{}".format(x,y) for x,y in itertools.product(['fwd','rev'], ['A','C','G','T','N'])] + ['cntIns','cntDel'];
theSamples = ['B7','D3','Normal'];

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

with lz4.frame.open('z.mp9k.lz4', newline='', mode='rt') as fp:
	tsvin = csv.DictReader(fp, delimiter='\t', fieldnames=['ChrID','Pos','Ref'] + 
		["{}_{}".format(x,y) for x,y in itertools.product(theSamples, ['Depth','Bases','Qs'])]
	)
	sbc = {key:collections.Counter() for key in theSamples}
	bc = {key:{} for key in theSamples}
	sbcB7 = {key:0 for key in CounterKeys};
	sbcD3 = {key:0 for key in CounterKeys};
	sbcNm = {key:0 for key in CounterKeys};
	try:
		for row in tsvin:
			#print(row)
			bcB7 = Base_Counter(row['B7_Depth'],row['B7_Bases'])
			bcD3 = Base_Counter(row['D3_Depth'],row['D3_Bases'])
			bcNm = Base_Counter(row['Normal_Depth'],row['Normal_Bases'])
			for k in CounterKeys:
				sbcB7[k] += bcB7[k];
				sbcD3[k] += bcD3[k];
				sbcNm[k] += bcNm[k];
			for k in theSamples:
				bc[k] = Base_Counter(row['{}_Depth'.format(k)],row['{}_Bases'.format(k)])
				sbc[k] += {x: bc[k][x] for x in CounterKeys}
	except KeyboardInterrupt:
		print('\n[!]Ctrl+C pressed.',file=sys.stderr,flush=True)
		pass
	#print('[!]Lines Read:[{}], MaxDepth is [{}].'.format(RecordCnt,MaxDepth),file=sys.stderr,flush=True)
	print(sbcB7);
	print(sbcD3);
	print(sbcNm);
	print(sbc);
