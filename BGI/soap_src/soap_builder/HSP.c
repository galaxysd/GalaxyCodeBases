/*

   HSP.c		BWTBlastn functions

   This module contains miscellaneous BWTBlastn functions.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "TextConverter.h"
#include "MiscUtilities.h"
#include "Socket.h"
#include "r250.h"
#include "HSP.h"
#include "HSPstatistic.h"


extern  double stat_expectationValue;

void HSPFillCharMap(unsigned char charMap[255]) {

	int i;

	for (i=0; i<255; i++) {
		charMap[i] = nonMatchDnaCharIndex;
	}
	for (i=0; i<16; i++) {
		charMap[(int)dnaChar[i]] = (unsigned char)i;
		charMap[(int)dnaChar[i] - 'A' + 'a'] = (unsigned char)i;
	}

}

void HSPFillComplementMap(unsigned char complementMap[255]) {

	int i;

	for (i=0; i<255; i++) {
		complementMap[i] = nonMatchDnaCharIndex;
	}
	for (i=0; i<16; i++) {
		complementMap[(int)dnaComplement[i]] = (unsigned int)i;
		complementMap[(int)dnaComplement[i] - 'A' + 'a'] = (unsigned int)i;
	}

}


HSP *HSPLoad(MMPool *mmPool, const char *PackedDNAFileName, const char *AnnotationFileName, const char *AmbiguityFileName) {

	HSP *hsp;

	FILE *annotationFile = NULL ;

	hsp = MMPoolDispatch(mmPool, sizeof(HSP));

	// Load packed DNA
	if (PackedDNAFileName != NULL && PackedDNAFileName[0] != '\0' && PackedDNAFileName[0] != '-') {
		hsp->packedDNA = DNALoadPacked(PackedDNAFileName, &hsp->dnaLength, TRUE);
	} else {
		hsp->packedDNA = NULL;
		hsp->dnaLength = 0;
	}
	return hsp;

}


void HSPFree(MMPool *mmPool, HSP *hsp) {

	if (hsp->packedDNA != NULL) {
		DNAFreePacked(hsp->packedDNA, hsp->dnaLength);
	}
	MMUnitFree(hsp->seqOffset, (hsp->numOfSeq+1) * sizeof(SeqOffset));
	MMUnitFree(hsp->annotation, (hsp->numOfSeq+1) * sizeof(Annotation));
	MMUnitFree(hsp->ambiguity, (hsp->numOfAmbiguity+2) * sizeof(Ambiguity));

	MMPoolReturn(mmPool, hsp, sizeof(hsp));

}

unsigned int HSPParseFASTAToPacked(const char * FASTAFileName, const char * annotationFileName, const char * packedDNAFileName, const char* ambiguityFileName, const char * repeatFileName, const unsigned int FASTARandomSeed, const int maskLowerCase) {

	FILE *FASTAFile, *annotationFile, *packedDNAFile, *ambiguityFile, *repeatFile;

	NewAnnotation *chrAnnotation;
	int chrAnnAllocated = 256;
	int blockAllocated = 256;

	char c, *ch;
	int chrNum, blockNum;
	unsigned int i, l;
	int nCount;
	unsigned int chrLen, usefulCharNum, numCharInBuffer, totalNumChar;
	unsigned int numRepeat = 0;
	unsigned char charMap[255];
	char *chrSeq, *p;
	unsigned int chrAllocated = 65536;
	unsigned char buffer[PACKED_BUFFER_SIZE];
	unsigned char packedBuffer[PACKED_BUFFER_SIZE / 4];
	chrLen = usefulCharNum = numCharInBuffer = totalNumChar = chrNum = blockNum = i = l = nCount = 0;
	FASTAFile = (FILE*)fopen64(FASTAFileName, "r");
	if (FASTAFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open FASTAFileName!\n");
		exit(1);
	}

	annotationFile = (FILE*)fopen64(annotationFileName, "w");
	if (annotationFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open annotationFileName!\n");
		exit(1);
	}

	packedDNAFile = (FILE*)fopen64(packedDNAFileName, "wb");
	if (packedDNAFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open packedDNAFileName!\n");
		exit(1);
	}

	ambiguityFile = (FILE*)fopen64(ambiguityFileName, "w");
	if (ambiguityFile == NULL) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open ambiguityFileName!\n");
		exit(1);
	}
/*
	repeatFile = (FILE *)fopen64(repeatFileName, "w");
	if (maskLowerCase && repeatFile == NULL ) {
		fprintf(stderr, "ParseFASTToPacked() : Cannot open repeatFileName!\n");
		exit(1);
	}
	//*/

	HSPFillCharMap(charMap);

	c = (char)getc(FASTAFile);
	if (c != '>') {
		fprintf(stderr, "ParseFASTToPacked() : FASTA file does not begin with '>'!\n");
		exit(1);
	}
	chrAnnotation = (NewAnnotation *)malloc(sizeof(NewAnnotation)*chrAnnAllocated);
	chrSeq = (char*)malloc(sizeof(char)*chrAllocated);
	chrNum = blockNum = usefulCharNum = numCharInBuffer = 0;
	while(!feof(FASTAFile)){
		if (feof(FASTAFile)) break;
		if (chrNum == chrAnnAllocated){
			chrAnnAllocated <<= 1;
			chrAnnotation = (NewAnnotation *)realloc(chrAnnotation, sizeof(NewAnnotation)*chrAnnAllocated);
//			printf("%d\n", chrNum);
		}

		l=0;
		c = (char)getc(FASTAFile);
		while(!feof(FASTAFile) && c!='\t' && c!=' ' && c!='\n' && l<MAX_SEQ_NAME_LENGTH){
			chrAnnotation[chrNum].chrName[l]=c;
			l++;
			c=(char)getc(FASTAFile);
		}
		chrAnnotation[chrNum].chrName[l]='\0';
		while(c!='\n'){
			c=(char)getc(FASTAFile);
		}
		chrLen = 0;
		while(c!='>' && !feof(FASTAFile)){
			if (c!='\n'){
				if (!maskLowerCase && c>='a' && c<='z'){
					c += 'A'-'a';
				}
				if (chrLen >= chrAllocated){
					chrAllocated <<= 1;
					chrSeq = (char*)realloc(chrSeq, sizeof(char)*chrAllocated);
				}
				*(chrSeq+chrLen) = c;
				chrLen += 1;
			}
			c=(char)getc(FASTAFile);
		}
		if (maskLowerCase) {
			p = chrSeq; 
			unsigned int repeat_beg , repeat_end = 0;
			repeat_end = repeat_beg = 0;
			i = 0;
			while (p != chrSeq + chrLen) {
				if (*p >= 'a' && *p <= 'z') {
					repeat_beg = i;
					while (*p >= 'a' && *p <= 'z' && p != chrSeq + chrLen) {
						*p += 'A' - 'a';
						i ++; p ++;
					}
					repeat_end = i;
				}
				i++; p++;
				numRepeat++;
				fprintf(repeatFile, "%d\t%u\t%u\n", chrNum, repeat_beg, repeat_end);
			}
		}
		if (chrLen <= 75) {
			fprintf(stderr, "ParseFASTToPacked(): reference <= 75 bp is filtered. Continue ...\n");
			continue;
		}
		//*
		i=0;
		p=chrSeq;
		while (ambiguityCount[charMap[(int)*p]] == 1 && i++ != chrLen) p++;
		if (i == chrLen) {
			blockNum = 1;
			chrAnnotation[chrNum].blockInChr = (ChrBlock *)malloc(sizeof(ChrBlock)*blockNum);
			chrAnnotation[chrNum].chrStart = usefulCharNum;
			chrAnnotation[chrNum].blockNum = blockNum;
			chrAnnotation[chrNum].blockInChr[0].blockStart = usefulCharNum;
			chrAnnotation[chrNum].blockInChr[0].ori = 0;
			usefulCharNum += chrLen;
			chrAnnotation[chrNum].chrEnd = usefulCharNum-1;
			chrAnnotation[chrNum].blockInChr[0].blockEnd = usefulCharNum-1;
			i=0;
			while(i<chrLen){
				if (numCharInBuffer >= PACKED_BUFFER_SIZE) {
					ConvertTextToBytePacked(buffer, packedBuffer, charMap, 4, PACKED_BUFFER_SIZE);
					fwrite(packedBuffer, 1, PACKED_BUFFER_SIZE / 4, packedDNAFile);
					numCharInBuffer = 0;
				}
				buffer[numCharInBuffer++] = chrSeq[i++];
			}
		} else {
			i=0;
			p = chrSeq;
			while (ambiguityCount[charMap[(int)*p]]!=1 && i++!=chrLen) p++;
			if (i<10) i=0;
			blockNum = 1;
			chrAnnotation[chrNum].blockInChr = (ChrBlock *)malloc(sizeof(ChrBlock)*blockAllocated);
			chrAnnotation[chrNum].chrStart = usefulCharNum;
			chrAnnotation[chrNum].blockInChr[blockNum-1].ori = i;
			chrAnnotation[chrNum].blockInChr[blockNum-1].blockStart = usefulCharNum;
			int len=0;
			while (i<chrLen) {
				if(ambiguityCount[charMap[(int)*p]] == 1){
					if (numCharInBuffer >= PACKED_BUFFER_SIZE) {
						ConvertTextToBytePacked(buffer, packedBuffer, charMap, 4, PACKED_BUFFER_SIZE);
						fwrite(packedBuffer, 1, PACKED_BUFFER_SIZE / 4, packedDNAFile);
						numCharInBuffer = 0;
					}
					buffer[numCharInBuffer++] = *p++;
					i++;
					usefulCharNum++;
					len++;
				}else{
					nCount = 0;
					while((ambiguityCount[charMap[(int)*p]]!=1) && i<chrLen){
						nCount++;
						i++;
						p++;
					}
					if (nCount<10) {
						do {
							if (numCharInBuffer >= PACKED_BUFFER_SIZE) {
								ConvertTextToBytePacked(buffer, packedBuffer, charMap, 4, PACKED_BUFFER_SIZE);
								fwrite(packedBuffer, 1, PACKED_BUFFER_SIZE / 4, packedDNAFile);
								numCharInBuffer = 0;
							}
							buffer[numCharInBuffer++] = 'G';
							usefulCharNum++;
							len++;
						} while(--nCount>0);
					} else {
						if (i<chrLen) {
							chrAnnotation[chrNum].blockInChr[blockNum-1].blockEnd = usefulCharNum -1;
							chrAnnotation[chrNum].blockInChr[blockNum-1].ori = i-nCount-len;
							if (blockNum == blockAllocated){
								blockAllocated <<= 1;
								chrAnnotation[chrNum].blockInChr = (ChrBlock *)realloc(chrAnnotation[chrNum].blockInChr, sizeof(ChrBlock)*blockAllocated);
							}
							blockNum++;
							len=0;
							chrAnnotation[chrNum].blockInChr[blockNum-1].blockStart = usefulCharNum;
						} else {
							i-=nCount;
							break;
						}
					}
				}
			}
			chrAnnotation[chrNum].blockInChr[blockNum-1].blockEnd = usefulCharNum-1;
			chrAnnotation[chrNum].blockInChr[blockNum-1].ori = i-len;
			chrAnnotation[chrNum].blockNum = blockNum;
			chrAnnotation[chrNum].chrEnd = usefulCharNum-1;
		}
//*/
		chrNum++;
		totalNumChar += chrLen;
	}
	if (numCharInBuffer > 0) {
		ConvertTextToBytePacked(buffer, packedBuffer, charMap, 4, numCharInBuffer);
		fwrite(packedBuffer, 1, (numCharInBuffer + 3) / 4, packedDNAFile);
		numCharInBuffer = 0;
	}
	if (totalNumChar % 4 == 0) {
		c = 0;
		fwrite(&c, 1, 1, packedDNAFile);
	}
	c = (char)(totalNumChar % 4);
	fwrite(&c, 1, 1, packedDNAFile);
	fclose(packedDNAFile);
	fprintf(annotationFile, "%u\t%d\t%d\n", totalNumChar, chrNum, FASTARandomSeed);
	int j=0;
	int total = 0;
	for (i=0;i<chrNum;i++) {
		fprintf(annotationFile, "%d\t%s\n", (int)strlen(chrAnnotation[i].chrName), chrAnnotation[i].chrName);
		total += chrAnnotation[i].blockNum;
	}
	fprintf(annotationFile, "%d\n", total);
	for(i=0;i<chrNum;i++){
		for(j=0;j<chrAnnotation[i].blockNum;j++){
			fprintf(annotationFile,"%d\t%u\t%u\t%u\n",i, chrAnnotation[i].blockInChr[j].blockStart, chrAnnotation[i].blockInChr[j].blockEnd, chrAnnotation[i].blockInChr[j].ori);
		}
		free(chrAnnotation[i].blockInChr);
	}

	if (maskLowerCase) {
		fprintf(repeatFile, "%u\n", numRepeat);
		fclose(repeatFile);
	}

	free(chrAnnotation);
	fclose(annotationFile);
	return chrNum;
}

NewAnnotation *AnnLoad(const char *annotationFileName, int *chrNum){
	FILE *annotationFile;
	if ((annotationFile = fopen(annotationFileName, "r")) == NULL){
		fprintf(stderr, "Load Annotation File error\n");
		exit(1);
	}
	NewAnnotation *chrAnn;
	int tmp, blockNum, totalNumChar, FASTARandomSeed;
	fscanf(annotationFile, "%u\t%d\t%d\n", &totalNumChar, &tmp, &FASTARandomSeed);
	chrAnn = (NewAnnotation *)malloc(sizeof(NewAnnotation)*tmp);
	int i,j;
	for(i=0;i<tmp;i++){
		fscanf(annotationFile, "%s\t%u\t%u\t%d\n", chrAnn[i].chrName, &chrAnn[i].chrStart, &chrAnn[i].chrEnd, &chrAnn[i].blockNum);
		chrAnn[i].blockInChr = (ChrBlock*) malloc(sizeof(ChrBlock)*chrAnn[i].blockNum);
		for(j=0;j<chrAnn[i].blockNum;j++){
			fscanf(annotationFile,"%u\t%u\t%u\n", &chrAnn[i].blockInChr[j].blockStart, &chrAnn[i].blockInChr[j].blockEnd, &chrAnn[i].blockInChr[j].ori);
		}
	}
	fclose(annotationFile);
	*chrNum=tmp;
	return chrAnn;
}



unsigned int HSPPackedToFASTA(const char* FASTAFileName, const char* annotationFileName, const char* packedDNAFileName, const char* ambiguityFileName) {

	HSP *hsp;
	FILE *FASTAFile;

	hsp = HSPLoad(NULL, packedDNAFileName, annotationFileName, ambiguityFileName);

	// Generate FASTA from packed
	FASTAFile = (FILE*)fopen64(FASTAFileName, "w");
	if (FASTAFile == NULL) {
		fprintf(stderr, "Cannot open FASTA file!\n");
		exit(1);
	}

	// to be done...
	fprintf(stderr, "HSPPackedToFASTA(): Function not complete!\n");

	fclose(FASTAFile);

	HSPFree(NULL, hsp);

	return 0;

}

