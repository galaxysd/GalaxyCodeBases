/*

   BWTFormatdb.c		Build index for FASTA database

   This program builds index for FASTA database for use of BWTBlastn.

   Copyright (C) 2006, Wong Chi Kwong.

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
#include "TypeNLimit.h"
#include "BWTConstruct.h"
#include "MiscUtilities.h"
#include "DNACount.h"
#include "TextConverter.h"
#include "MemManager.h"
#include "iniparser.h"
#include "HSP.h"
#include "Timing.h"
#include "lookupBuilder.h"
#include "highOccBuilder.h"

#define MAX_FILENAME_LEN 1024

// Database and ini
dictionary *ParseInput(int argc, char** argv);
void ParseIniFile(char *iniFileName);
void ProcessIni();
void ValidateIni();
void PrintIni();
void PrintShortDesc();
void PrintHelp();

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName);
int buildLookupTable(char * PackedDNAFileName, char * LookupTableFileName, int lookupTableSize);

	// Parameters
	char IniFileName[MAX_FILENAME_LEN+1];
	int Confirmation;
	
	// BuildTasks parameters
	int ParseFASTA = TRUE;
	int BuildBWT = TRUE;
	int BuildSaValue = TRUE;
	int BuildSaIndex = TRUE; 
	int BuildLookUp = TRUE;
	int BuildHOT = TRUE;

	// Memory parameters
	unsigned int PoolSize = 2097152;				// 2M  - fixed; not configurable through ini

	// Display parameters
	int ShowProgress = FALSE;

	// Database parameters
	char FASTAFileName[MAX_FILENAME_LEN+1] = "";
	char DatabaseName[MAX_FILENAME_LEN+1] = "";
	char AnnotationFileName[MAX_FILENAME_LEN+1] = "*.index.ann";
	char AmbiguityFileName[MAX_FILENAME_LEN+1] = "*.index.amb";
	char RepeatFileName[MAX_FILENAME_LEN+1] = "*.index.rep";
	char PackedDNAFileName[MAX_FILENAME_LEN+1] = "*.index.pac";
	char BWTCodeFileName[MAX_FILENAME_LEN+1] = "*.index.bwt";
	char BWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.index.fmv";
	char SaValueFileName[MAX_FILENAME_LEN+1] = "*.index.sa";
	char SaIndexFileName[MAX_FILENAME_LEN+1] = "*.index.sai";
	
	char RevPackedDNAFileName[MAX_FILENAME_LEN+1] = "*.index.rev.pac";
	char RevBWTCodeFileName[MAX_FILENAME_LEN+1] = "*.index.rev.bwt";
	char RevBWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.index.rev.fmv";
	
	char LookupTableFileName[MAX_FILENAME_LEN+1] = "*.index.lkt";
	char RevLookupTableFileName[MAX_FILENAME_LEN+1] = "*.index.rev.lkt";

    char HighOccHashTableFileName[MAX_FILENAME_LEN+1] = "*.index.hot";
	// Parse FASTA parameters
	unsigned int FASTARandomSeed = 0;
	int MaskLowerCase = FALSE;

	// Build BWT parameters
	unsigned int OccValueFreq = 256;
	float TargetNBit = 2.5;
	unsigned int InitialMaxBuildSize = 10000000;
	unsigned int IncMaxBuildSize = 10000000;

	// Build SA value parameters
	unsigned int SaValueFreq = 8;

	// Build SA index parameters
	unsigned int SaIndexNumOfChar = 12;
	
	//Look Up Table parameters
	unsigned int LookUpTableSize=13;
    unsigned int ReversedLookUpTableSize=13;
	
	//High Occurrences Pattern Hash Table parameters
	unsigned int HashPatternLength=35;
	unsigned int HashOccThreshold=4;

void ExtractionHighOccPattern(const BWT *bwt, int index, unsigned int l, unsigned int r) {
    if (l<=r && r-l+1>=HashOccThreshold) {
        if (index>0) {
            //If the packing is not finished
            unsigned int new_l,new_r;
            
            unsigned int occCount_pstart[4];
            unsigned int occCount_pend[4];
            BWTAllOccValue(bwt,l,occCount_pstart);
            BWTAllOccValue(bwt,r + 1,occCount_pend);
            
            unsigned char ec;
            for (ec=0;ec<4;ec++) {
        		new_l = bwt->cumulativeFreq[ec] + occCount_pstart[ec] + 1;
        		new_r = bwt->cumulativeFreq[ec] + occCount_pend[ec];

                ExtractionHighOccPattern(bwt,index-1,new_l,new_r);
            }
        } else if (index==0) {
            //If the packing pattern reaches the maximum length 
			//l,r
			ho_append(l,r);
        }
    }
}

void ho_writetofile(const char * fileName,BWT * bwt) {
    //fprintf(stderr,"Writing Hash Table to File..\n");
     unsigned int acc_index=0;
     unsigned int acc_occIndex=0;
     unsigned int i;
     
    FILE * outFile;
    if(!(outFile = fopen(fileName, "w"))) return;
    
    fwrite(&ho_hashTableSize,sizeof(unsigned int),1,outFile);
    fwrite(&ho_hash_a,sizeof(unsigned int),1,outFile);
    fwrite(&ho_hash_b,sizeof(unsigned int),1,outFile);
    fwrite(&ho_prime,sizeof(unsigned int),1,outFile);
    fwrite(&ho_ttlItem,sizeof(unsigned int),1,outFile);
    fwrite(&ho_ttlOccurrence,sizeof(unsigned int),1,outFile);
    
    //fprintf(stderr,"Writing Hash Table to File..Phrase One\n");
    for (i=0;i<ho_hashTableSize;i++) {
        if (ho_hashtable[i].count==0) {
            char c = 0;
            fwrite(&c,1,1,outFile);
        } else {
            char c = 1;
            fwrite(&c,1,1,outFile);
            fwrite(&(acc_index),sizeof(unsigned int),1,outFile);
            fwrite(&(ho_hashtable[i].count),sizeof(unsigned int),1,outFile);
            acc_index+=ho_hashtable[i].count;
        }
    }
    //fprintf(stderr,"Writing Hash Table to File..Phrase Two\n");
    for (i=0;i<ho_hashTableSize;i++) {
       hashItemPtr node = ho_hashtable[i].item;
       while (node!=NULL) {
            fwrite(&node->sa_l,sizeof(unsigned int),1,outFile);
            fwrite(&node->sa_r,sizeof(unsigned int),1,outFile);
            fwrite(&(acc_occIndex),sizeof(unsigned int),1,outFile);
            acc_occIndex+=(node->sa_r)-(node->sa_l)+1;
            node=node->next;
       }
        
    }
    //fprintf(stderr,"Writing Hash Table to File..Phrase Three\n");
    //unsigned int * buffer = malloc(sizeof(unsigned int)*1024);
   // unsigned int bufferIndex =0 ;
    for (i=0;i<ho_hashTableSize;i++) {
       hashItemPtr node = ho_hashtable[i].item;
       while (node!=NULL) {
             unsigned int j;
             for (j=node->sa_l;j<=node->sa_r;j++) {
                unsigned int tp = BWTSaValue(bwt,j);
		          fwrite(&tp,sizeof(unsigned int),1,outFile);
                // buffer[bufferIndex++]=tp;
                // if (bufferIndex>=1024) {
                     //fwrite(buffer,sizeof(unsigned int)*1024,1,outFile);
                //     bufferIndex=0;
                // }
             }
            node=node->next;
       }
    }
    //free(buffer);
    fclose(outFile);//*/
    //fprintf(stderr,"Finished Writing Hash Table to File.\n");
}

void ho_buildandpopulate(BWT * bwt) {
    
    //fprintf(stderr,"Build and Populate Hash Table.\n");
    //Initialize the Hash Table of size doubled the number of pattern
     ho_hashTableSize=ho_ttlItem*2;
     ho_hashtable=malloc(sizeof(struct hashTyp)*ho_hashTableSize);
     unsigned int i;
     for (i=0;i<ho_hashTableSize;i++) {
         ho_hashtable[i].count=0;
         ho_hashtable[i].item=NULL;
     }
     
     ho_cellSpacePtr node = ho_root;
     while (node!=NULL) {
           for (i=0;i<node->count;i++) {
               unsigned int index = ho_hash(node->saL[i]);
               hashItemPtr newItem = malloc(sizeof(struct hashItem));
               newItem->sa_l=node->saL[i];
               newItem->sa_r=node->saR[i];
               newItem->next=ho_hashtable[index].item;
               ho_hashtable[index].item=newItem;
               
               ho_hashtable[index].count++;
           }
           node=node->next;
     }
     
     ho_writetofile(HighOccHashTableFileName,bwt);
}

void buildReversePacked(const char *inputFileName, unsigned int *textLength, const unsigned int convertToWordPacked) {

	FILE *inputFile;
	FILE *outputFile;
	unsigned char tempChar[4];
	unsigned int writeBuffer=0;
	unsigned int writeBufferIndex=0;
	unsigned int *packedText;
	unsigned int packedFileLen;
	unsigned char lastByteLength;
	unsigned int wordToProcess;
	long long i;
	int j,k;

	inputFile = (FILE*)fopen(inputFileName, "rb");
	outputFile = (FILE*)fopen(RevPackedDNAFileName, "wb");

	if (inputFile == NULL) {
		fprintf(stderr, "buildReversePacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(inputFile, -1, SEEK_END);
	packedFileLen = ftell(inputFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "buildReversePacked(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, inputFile);

	*textLength = (packedFileLen - 1) * 4 + lastByteLength;

	wordToProcess = (*textLength + 16 - 1) / 16;
	packedText = MMUnitAllocate((wordToProcess + 1) * sizeof(unsigned int));	// allocate 1 more word at end
	if (packedText == NULL) {
		fprintf(stderr, "buildReversePacked(): Cannot MMUnitAllocate packedText\n");
		exit(0);
	}
	packedText[wordToProcess - 1] = 0;
	packedText[wordToProcess] = 0;

	fseek(inputFile, 0L, SEEK_SET);
	fread(packedText, 1, packedFileLen, inputFile);
	fclose(inputFile);
	
	
    //printf("lastByteLength = %u\n",lastByteLength);
    long long currentLocation = (wordToProcess)*16;
    //printf("currentLocation = %u\n",currentLocation);
    //printf("*textLength = %u\n",*textLength);
	for (i=wordToProcess-1; i>=0; i--) {
		*(unsigned int*)tempChar = packedText[i];
		packedText[i] = (tempChar[0] << 24) | (tempChar[1] << 16) | (tempChar[2] << 8) | tempChar[3];

		unsigned int mk = packedText[i];
		for (j=0;j<16;j++) {
			if (writeBufferIndex>=16) {
				*(unsigned int*)tempChar = writeBuffer;
				writeBuffer = (tempChar[0] << 24) | (tempChar[1] << 16) | (tempChar[2] << 8) | tempChar[3];
				fwrite(&writeBuffer,sizeof(unsigned int),1,outputFile);
				writeBufferIndex=0;
			}
            writeBuffer<<=2;
            unsigned char c = mk & 3;
            mk>>=2;
            currentLocation--;
            if (currentLocation>=*textLength) {continue;}
            writeBuffer|=c;writeBufferIndex++;
            
        }
    }
    //*/
    
	
    
	MMUnitFree(packedText, ((wordToProcess + 1) * sizeof(unsigned int)));
	//printf("Finished main loop..\n");
	if (writeBufferIndex>0) {
        printf("Wiping..%d\n",writeBufferIndex);
        for (k=writeBufferIndex;k<16;k++) writeBuffer<<=2;
        *(unsigned int*)tempChar = writeBuffer;
        if (writeBufferIndex<4) {
        	fwrite(&tempChar[3],sizeof(unsigned char),1,outputFile);
        } else if (writeBufferIndex<8) {
        	fwrite(&tempChar[3],sizeof(unsigned char),1,outputFile);
        	fwrite(&tempChar[2],sizeof(unsigned char),1,outputFile);
        } else if (writeBufferIndex<12) {
        	fwrite(&tempChar[3],sizeof(unsigned char),1,outputFile);
        	fwrite(&tempChar[2],sizeof(unsigned char),1,outputFile);
        	fwrite(&tempChar[1],sizeof(unsigned char),1,outputFile);
        } else if (writeBufferIndex<16) {
        	fwrite(&tempChar[3],sizeof(unsigned char),1,outputFile);
        	fwrite(&tempChar[2],sizeof(unsigned char),1,outputFile);
        	fwrite(&tempChar[1],sizeof(unsigned char),1,outputFile);
        	fwrite(&tempChar[0],sizeof(unsigned char),1,outputFile);
        }
    }
    fwrite(&lastByteLength,sizeof(unsigned char),1,outputFile);
    fclose(outputFile);
    //free(packedText);
	return;
}

int main(int argc, char** argv) {

	char c;
	MMPool *mmPool;
	dictionary *programInput;
	double startTime;
	double elapsedTime = 0, totalElapsedTime = 0;

	char filename[MAX_FILENAME_LEN+1];
	BWT *bwt = NULL;
	unsigned int textLength = 0;
	unsigned int numSeq;

	BWTInc *bwtInc = NULL;

	// Program input
	programInput = ParseInput(argc, argv);
	PrintShortDesc();

	// Ini
	if (strcmp(argv[0] + strlen(argv[0]) - 4, ".exe") == 0) {
		*(argv[0] + strlen(argv[0]) - 4) = '\0';
	}
	sprintf(filename, "%s.ini", argv[0]);
	ParseIniFile(filename);
	//printf("\n");
	ProcessIni();
	ValidateIni();
	PrintIni();

	if (Confirmation == TRUE) {
		printf("Press Y to go or N to cancel. ");
        c = (char)getchar();
		while (c != 'y' && c != 'Y' && c != 'n' && c!= 'N') {
			c = (char)getchar();
		}
		if (c == 'n' || c == 'N') {
			exit(0);
		}
	}

	startTime = setStartTime();

	MMMasterInitialize(1, 0, FALSE, NULL);
	mmPool = MMPoolCreate(PoolSize);

	// Parse FASTA file to produce packed DNA and annotation file
	if (ParseFASTA == TRUE) {

		printf("Parsing FASTA file..\n");

		numSeq = HSPParseFASTAToPacked(FASTAFileName, AnnotationFileName, PackedDNAFileName, AmbiguityFileName, RepeatFileName, FASTARandomSeed, MaskLowerCase);
		
		
        //Parse packed DNA to construct the packed reversed DNA
        
		unsigned int textLen;
        buildReversePacked(PackedDNAFileName,&textLen,TRUE);
	    //printf("Reversed Packed DNA generated..\n");
		printf("Finished. Parsed %u sequences.\n", numSeq);


		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

	}
	
	if (BuildLookUp == TRUE) {
    	//Construct Look-Up Table
    	printf("Building Look-Up..\n");
    	buildLookupTable(PackedDNAFileName,LookupTableFileName,LookUpTableSize);
    	//printf("Look-Up Table is built.\n");
    	buildLookupTable(RevPackedDNAFileName,RevLookupTableFileName,ReversedLookUpTableSize);
    	//printf("Reversed Look-Up Table is built.\n");
    	
		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Finished.\nElapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");
    }

	// Construct BWTInc from text
	if (BuildBWT == TRUE) {

		printf("Building BWT..\n");

		bwtInc = BWTIncConstructFromPacked(mmPool, PackedDNAFileName, ShowProgress, 
										   TargetNBit, InitialMaxBuildSize, IncMaxBuildSize);

		printf("Finished constructing BWT in %u iterations.  ", bwtInc->numberOfIterationDone);
		
		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		printf("Saving BWT..\n");
		BWTSaveBwtCodeAndOcc(bwtInc->bwt, BWTCodeFileName, BWTOccValueFileName);
		printf("Finished saving BWT.  ");
		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		textLength = bwtInc->bwt->textLength;

		BWTIncFree(mmPool, bwtInc);



        //Building Reversed BWT
		printf("Building Reversed BWT..\n");

		bwtInc = BWTIncConstructFromPacked(mmPool, RevPackedDNAFileName, ShowProgress, 
										   TargetNBit, InitialMaxBuildSize, IncMaxBuildSize);

		printf("Finished constructing Reversed BWT in %u iterations.  ", bwtInc->numberOfIterationDone);
		
		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		printf("Saving BWT..\n");
		BWTSaveBwtCodeAndOcc(bwtInc->bwt, RevBWTCodeFileName, RevBWTOccValueFileName);
		printf("Finished saving BWT.  ");
		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		textLength = bwtInc->bwt->textLength;

		BWTIncFree(mmPool, bwtInc);

	}

	// Load BWT
	if (BuildSaValue || BuildSaIndex) {

		printf("Loading BWT...\n");

		bwt = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, NULL, NULL, NULL, NULL);
        //Use BWT to build the hash table

		printf("Finished loading BWT.  ");

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		textLength = bwt->textLength;

	} else if (BuildHOT) {
		printf("Loading BWT and SA Values...\n");

		bwt = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, SaValueFileName, NULL, NULL, NULL);
        //Use BWT to build the hash table

		printf("Finished loading BWT and SA Values.  ");

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

		textLength = bwt->textLength;

    }
	

	if (BuildSaValue) {

		printf("Building SA value...\n");
		
		if (ShowProgress) {
			BWTGenerateSaValue(mmPool, bwt, SaValueFreq, bwt->textLength / SaValueFreq / 10);
		} else {
			BWTGenerateSaValue(mmPool, bwt, SaValueFreq, 0);
		}
		BWTSaveSaValue(bwt, SaValueFileName);

		printf("Finished building SA value.  ");

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

	}
	
	//Build High Occurrences Pattern Hash Table
	if (BuildHOT) {
		printf("Building High-Occ Hash Table...\n");
		ho_initialize();
        //printf("# of Occ = %u\n",ho_ttlOccurrence);
		unsigned int l = 0;
		unsigned int r = bwt->textLength;
		ExtractionHighOccPattern(bwt,HashPatternLength,l,r);
		//printf("# of SA = %u\n",ho_count());
        //printf("# of Occ = %u\n",ho_ttlOccurrence);
		ho_buildandpopulate(bwt);
		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Finished.\nElapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");
		ho_free();
    }

	if (BuildSaIndex) {

		printf("Building SA index...\n");
		
		BWTGenerateSaRangeTable(bwt, SaIndexNumOfChar, SaIndexFileName);

		printf("Finished building SA index.  ");

		elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
		printf("Elapsed time = ");
		printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
		totalElapsedTime += elapsedTime;
		printf("\n");

	}

	// Free BWT
	if (BuildSaValue || BuildSaIndex) {
		BWTFree(mmPool, bwt);
	}

	// Finished all construction tasks
	printf("Index building is completed.\n");
	totalElapsedTime = getElapsedTime(startTime);
	printf("Total elapsed time = ");
	printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, totalElapsedTime);
	printf("\n");

	//MMMasterPrintReport(stdout, FALSE, FALSE, FALSE);
	if (BuildSaValue) {
		//fprintf(stdout, "Number of char   :  %u\n", textLength);
		//fprintf(stdout, "Bit per char     :  %.2f\n", (float)MMMasterMaxTotalByteDispatched() * BITS_IN_BYTE / textLength);
		//printf("\n");
	}

	MMPoolFree(mmPool);

	iniparser_freedict(programInput);

	return 0;

}

dictionary *ParseInput(int argc, char** argv) {

	dictionary *programInput;
	char t1[3] = "-c";	// specify that this is a boolean type parameter
	char t2[3] = "-U";	// specify that this is a boolean type parameter
	char *d[2];

	d[0] = t1;
	d[1] = t2;
	
	programInput = paraparser_load(argc, argv, 2, d);	// 2 boolean type parameters

	// Get database name
	if (!iniparser_find_entry(programInput, "argument:1")) {
		PrintHelp();
		exit(1);
	}
	iniparser_copystring(programInput, "argument:1", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
	if (strlen(DatabaseName) + 4 > MAX_FILENAME_LEN) {
		PrintHelp();
		exit(1);
	}

	// Get FASTA file name
	iniparser_copystring(programInput, "argument:2", FASTAFileName, DatabaseName, MAX_FILENAME_LEN);
	if (strlen(FASTAFileName) > MAX_FILENAME_LEN) {
		PrintHelp();
		exit(1);
	}


	// Whether confirmation is needed
	Confirmation = iniparser_find_entry(programInput, "parameter:-c");

	MaskLowerCase = iniparser_find_entry(programInput, "parameter:-U");

	return programInput;

}

void ParseIniFile(char *iniFileName) {

	dictionary *ini;

	//printf("Loading %s ..", iniFileName);
	ini = iniparser_load(iniFileName, FALSE);
	if (ini == NULL) {
	//	printf("not found.\n");
		return;
	}
	//printf("done.\n");

	// BuildTasks parameters
	ParseFASTA = iniparser_getboolean(ini, "BuildTasks:ParseFASTA", ParseFASTA);
	BuildBWT = iniparser_getboolean(ini, "BuildTasks:BuildBWT", BuildBWT);
	BuildSaValue = iniparser_getboolean(ini, "BuildTasks:BuildSaValue", BuildSaValue);
	BuildLookUp = iniparser_getboolean(ini, "BuildTasks:BuildLookUp", BuildLookUp);
	BuildSaIndex = iniparser_getboolean(ini, "BuildTasks:BuildSaIndex", BuildSaIndex);
	BuildHOT = iniparser_getboolean(ini, "BuildTasks:BuildHOT", BuildHOT);

	// Display parameters
	ShowProgress = iniparser_getboolean(ini, "Display:ShowProgress", ShowProgress);

	// Parse FASTA parameters
	FASTARandomSeed = iniparser_getint(ini, "ParseFASTA:RandomSeed", FASTARandomSeed);
	if (FASTARandomSeed == 0) {
		FASTARandomSeed = getRandomSeed();
	}

	// Build BWT parameters
	OccValueFreq = iniparser_getint(ini, "BuildBWT:OccValueFreq", OccValueFreq);
	TargetNBit = (float)iniparser_getdouble(ini, "BuildBWT:TargetNBit", TargetNBit);
	InitialMaxBuildSize = iniparser_getint(ini, "BuildBWT:InitialMaxBuildSize", InitialMaxBuildSize);
	IncMaxBuildSize = iniparser_getint(ini, "BuildBWT:IncMaxBuildSize", IncMaxBuildSize);

	// Build SA value parameters
	SaValueFreq = iniparser_getint(ini, "BuildSAValue:SaValueFreq", SaValueFreq);

	// Build SA index parameters
	SaIndexNumOfChar = iniparser_getint(ini, "BuildSAIndex:SaIndexNumOfChar", SaIndexNumOfChar);

    // Build Look Up parameters
	LookUpTableSize = iniparser_getint(ini, "BuildLookUp:LookUpTableSize", LookUpTableSize);
    ReversedLookUpTableSize = iniparser_getint(ini, "BuildLookUp:ReversedLookUpTableSize", ReversedLookUpTableSize);
	

	// Build High Occurrences Pattern Hash Table parameters
	HashPatternLength = iniparser_getint(ini, "BuildHOT:PatternLength", HashPatternLength);
	HashOccThreshold = iniparser_getint(ini, "BuildHOT:OccurrencesThreshold", HashOccThreshold);
	// Database parameters
	iniparser_copystring(ini, "Database:AnnotationFileName", AnnotationFileName, AnnotationFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:AmbiguityFileName", AmbiguityFileName, AmbiguityFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:PackedDNAFileName", PackedDNAFileName, PackedDNAFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:BWTCodeFileName", BWTCodeFileName, BWTCodeFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:BWTOccValueFileName", BWTOccValueFileName, BWTOccValueFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:SaValueFileName", SaValueFileName, SaValueFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:SaIndexFileName", SaIndexFileName, SaIndexFileName, MAX_FILENAME_LEN);
	
	iniparser_copystring(ini, "Database:RevPackedDNAFileName", RevPackedDNAFileName, RevPackedDNAFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:RevBWTCodeFileName", RevBWTCodeFileName, RevBWTCodeFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:RevBWTOccValueFileName", RevBWTOccValueFileName, RevBWTOccValueFileName, MAX_FILENAME_LEN);
	
	iniparser_copystring(ini, "Database:LookupTableFileName", LookupTableFileName, LookupTableFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:RevLookupTableFileName", RevLookupTableFileName, RevLookupTableFileName, MAX_FILENAME_LEN);
	iniparser_copystring(ini, "Database:HighOccHashTableFileName", HighOccHashTableFileName, HighOccHashTableFileName, MAX_FILENAME_LEN);
	iniparser_freedict(ini);

}

void ProcessIni() {

	ProcessFileName(AnnotationFileName, AnnotationFileName, DatabaseName);
	ProcessFileName(AmbiguityFileName, AmbiguityFileName, DatabaseName);
	ProcessFileName(PackedDNAFileName, PackedDNAFileName, DatabaseName);
	ProcessFileName(RevPackedDNAFileName, RevPackedDNAFileName, DatabaseName);
	ProcessFileName(BWTCodeFileName, BWTCodeFileName, DatabaseName);
	ProcessFileName(RevBWTCodeFileName, RevBWTCodeFileName, DatabaseName);
	ProcessFileName(BWTOccValueFileName, BWTOccValueFileName, DatabaseName);
	ProcessFileName(RevBWTOccValueFileName, RevBWTOccValueFileName, DatabaseName);
	ProcessFileName(SaValueFileName, SaValueFileName, DatabaseName);
	ProcessFileName(SaIndexFileName, SaIndexFileName, DatabaseName);
	
	ProcessFileName(LookupTableFileName, LookupTableFileName, DatabaseName);
	ProcessFileName(RevLookupTableFileName, RevLookupTableFileName, DatabaseName);
	
	ProcessFileName(HighOccHashTableFileName, HighOccHashTableFileName, DatabaseName);

}

void ValidateIni() {

	if (!ParseFASTA && !BuildBWT && !BuildSaValue && !BuildSaIndex && !BuildLookUp && !BuildHOT) {
		fprintf(stderr, "No action is specified!\n");
		exit(1);
	}
	if (BuildLookUp) {
		if (PackedDNAFileName[0] == '\0') {
			fprintf(stderr, "Packed DNA file name is not specified!\n");
			exit(1);
		}
		if (RevPackedDNAFileName[0] == '\0') {
			fprintf(stderr, "Reversed Packed DNA file name is not specified!\n");
			exit(1);
		}
    }
	if (ParseFASTA) {
		if (PackedDNAFileName[0] == '\0') {
			fprintf(stderr, "Packed DNA file name is not specified!\n");
			exit(1);
		}
		if (AnnotationFileName[0] == '\0') {
			fprintf(stderr, "Annotation file name is not specified!\n");
			exit(1);
		}
		if (AmbiguityFileName[0] == '\0') {
			fprintf(stderr, "Ambiguity file name is not specified!\n");
			exit(1);
		}
	}
	if (BuildBWT) {
		if (PackedDNAFileName[0] == '\0') {
			fprintf(stderr, "Packed DNA file is not specified!\n");
			exit(1);
		}
		if (BWTCodeFileName[0] == '\0') {
			fprintf(stderr, "BWT code file name is not specified!\n");
			exit(1);
		}
		if (BWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "BWT Occ value file name is not specified!\n");
			exit(1);
		}
		if (TargetNBit < 2.5) {
			fprintf(stderr, "Target NBit should be at least 2.5!\n");
			exit(1);
		}
	}
	if (BuildSaValue) {
		if (BWTCodeFileName[0] == '\0') {
			fprintf(stderr, "BWT code file is not specified!\n");
			exit(1);
		}
		if (BWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "BWT Occ value file is not specified!\n");
			exit(1);
		}
		if (SaValueFileName[0] == '\0') {
			fprintf(stderr, "SA value file name is not specified!\n");
			exit(1);
		}
		if (SaValueFreq <= 0) {
			fprintf(stderr, "SA value frequency must > 0!\n");
			exit(1);
		}
	}

	if (BuildSaIndex) {
		if (BWTCodeFileName[0] == '\0') {
			fprintf(stderr, "BWT code file is not specified!\n");
			exit(1);
		}
		if (BWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "BWT Occ value file is not specified!\n");
			exit(1);
		}
		if (SaIndexFileName[0] == '\0') {
			fprintf(stderr, "SA index file name is not specified!\n");
			exit(1);
		}
		if (SaIndexNumOfChar <= 0) {
			fprintf(stderr, "SA index number of character must > 0!\n");
			exit(1);
		}
		if (SaIndexNumOfChar > 13) {
			fprintf(stderr, "SA index number of character must <= 13!\n");
			exit(1);
		}
	}

	if (BuildHOT) {
		if (BWTCodeFileName[0] == '\0') {
			fprintf(stderr, "BWT code file is not specified!\n");
			exit(1);
		}
		if (BWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "BWT Occ value file is not specified!\n");
			exit(1);
		}
		if (RevBWTCodeFileName[0] == '\0') {
			fprintf(stderr, "Reversed BWT code file is not specified!\n");
			exit(1);
		}
		if (RevBWTOccValueFileName[0] == '\0') {
			fprintf(stderr, "Reversed BWT Occ value file is not specified!\n");
			exit(1);
		}
    }
}


void PrintIni() {

	char boolean[2];

	boolean[0] = 'N';
	boolean[1] = 'Y';

	/*printf("Parse FASTA file    : %c\n", boolean[ParseFASTA]);
	printf("Build BWT           : %c\n", boolean[BuildBWT]);
	printf("Build SA value      : %c\n", boolean[BuildSaValue]);
	printf("Build SA index      : %c\n", boolean[BuildSaIndex]);
	printf("\n");

	printf("Show progress       : %c\n", boolean[ShowProgress]);
	printf("\n");

	if (ParseFASTA) {
		printf("Parse FASTA :\n");
		printf("Mask lower case         : %c\n", boolean[MaskLowerCase]);
		printf("Random seed             : %u\n", FASTARandomSeed);
		printf("\n");
	}

	if (BuildBWT) {
		printf("Build BWT :\n");
		printf("Target N Bits           : %.2f\n", TargetNBit);
		printf("Occ value frequency     : %u\n", OccValueFreq);
		printf("Initial Max Build Size  : %u    Inc Max Build Size : %u\n", 
				InitialMaxBuildSize, IncMaxBuildSize);
		printf("\n");
	}

	if (BuildSaValue) {
		printf("Build SA value :\n");
		printf("SA value frequency      : %u\n", SaValueFreq);
		printf("\n");
	}

	if (BuildSaIndex) {
		printf("Build SA index :\n");
		printf("SA index no. of char    : %u\n", SaIndexNumOfChar);
		printf("\n");
	}

	printf("Annotation file          : %s\n", AnnotationFileName);
	printf("Ambigurity file          : %s\n", AmbiguityFileName);
	printf("Packed DNA file          : %s\n", PackedDNAFileName);
	printf("BWT Code file            : %s\n", BWTCodeFileName);
	printf("BWT Occ value file       : %s\n", BWTOccValueFileName);
	printf("SA value file            : %s\n", SaValueFileName);
	printf("SA index file            : %s\n", SaIndexFileName);
	printf("\n");
	printf("Reversed Packed DNA file          : %s\n", RevPackedDNAFileName);
	printf("Reversed BWT Code file            : %s\n", RevBWTCodeFileName);
	printf("Reversed BWT Occ value file       : %s\n", RevBWTOccValueFileName);
	printf("\n");
	printf("Look-Up Table file                : %s\n", LookupTableFileName);
	printf("Reversed Look-Up Table file       : %s\n", RevLookupTableFileName);
	printf("\n");
	printf("High Occ Hash Table file          : %s\n", HighOccHashTableFileName);
	printf("\n");*/

}

void PrintShortDesc() {

	/*printf("BWTFormatdb v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
	printf("BWTFormatdb comes with ABSOLUTELY NO WARRENTY.\n");
	printf("BWTFormatdb is free software, and you are welcome to\n");
	printf("redistribute it under certain conditions.\n");
	printf("For details type BWTFormatdb.\n");
	printf("\n");*/

}

void PrintHelp() {

	/*printf("BWTFormatdb v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
	printf("\n");

	printf("This program is free software; you can redistribute it and/or\n");
	printf("modify it under the terms of the GNU General Public License\n");
	printf("as published by the Free Software Foundation; either version 2\n");
	printf("of the License, or (at your option) any later version.\n");
	printf("\n");

	printf("This program is distributed in the hope that it will be useful,\n");
	printf("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
	printf("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
	printf("GNU General Public License for more details.\n");
	printf("\n");

	printf("You should have received a copy of the GNU General Public License\n");
	printf("along with this program; if not, write to the Free Software\n");
	printf("Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.\n");
	printf("\n");*/

	printf("Syntax: 2bwt-builder <sequence file>\n");
	printf("Compile Time: " COMPILE_TIME"\n");
}

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName) {

	char tempChar[MAX_FILENAME_LEN];
	unsigned int i;

	if (inputFileName == NULL) {
		if (outputFileName != inputFileName) {
			outputFileName[0] = '\0';
		}
		return;
	}

	if (strlen(databaseName) + strlen(inputFileName) > MAX_FILENAME_LEN) {
		fprintf(stderr, "File length is too long!\n");
		exit(1);
	}

	strncpy(tempChar, inputFileName, MAX_FILENAME_LEN);

	// locate the *
	for (i=0; i<MAX_FILENAME_LEN; i++) {
		if (tempChar[i] == '*') {
			break;
		}
	}
	if (i<MAX_FILENAME_LEN) {
		tempChar[i] = '\0';
		sprintf(outputFileName, "%s%s%s", tempChar, databaseName, tempChar + i + 1);
	} else {
		sprintf(outputFileName, "%s", tempChar);
	}

}

