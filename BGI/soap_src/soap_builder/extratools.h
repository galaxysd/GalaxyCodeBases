#ifndef _EXTRATOOLS_H_
#define _EXTRATOOLS_H_

#include <stdio.h>
#include <stdlib.h>
#include "MiscUtilities.h"
#include "MemManager.h"
#include "TextConverter.h"
#include "Timing.h"
#include "BWT.h"
#include "HSP.h"
#include "Types.h"
#include <fcntl.h>

typedef struct LOOKUPTABLE_TYPE {
       unsigned int tableSize;
       unsigned int * table;
}LOOKUPTABLE;

typedef struct HASHCELL_TYPE {
    unsigned int count;
    unsigned int index;
}HASHCELL;

typedef struct HASHITEM_TYPE {
    unsigned int l;
    unsigned int r;
    unsigned int occIndex;
}HASHITEM;

typedef unsigned int OCC;

typedef struct HASHTABLE_TYPE {
       unsigned int prime;
       unsigned int a;
       unsigned int b;
       unsigned int tableSize;
       HASHCELL * table;
       HASHITEM * itemList;
       OCC * occList;
}HASHTABLE;



BWT * occBwt;
HASHTABLE * occHashtable;
unsigned int * occCollector;
unsigned int occCollected;

FILE * textPositionFile;
void registerTPFile(FILE * filePtr,unsigned int searchMode);


void registerQIndex(unsigned int queryIndex);
void registerQSection();


void LoadLookupTable(LOOKUPTABLE * lookupTable, const char * fileName, const int tableSize);
void FreeLookupTable(LOOKUPTABLE * lookupTable);
unsigned int LookupSafe(LOOKUPTABLE lookupTable, BWT * bwt,unsigned long long lKey, unsigned long long rKey,unsigned int *l, unsigned int *r);
void LoadHashTable(HASHTABLE * hashTable, const char * fileName);
HASHITEM * HashFind(HASHTABLE * hashTable, unsigned int l,unsigned int r);
void FreeHashTable(HASHTABLE * hashTable);
void RegisterDecoder(BWT * bwt,HASHTABLE * hashTable);

//void OCCClean();
//void OCCProcess(unsigned int l,unsigned int r);
void OCCProcess(unsigned int l,unsigned int r, int chain, unsigned int allele1, unsigned int allele2, HitInfo *hits, unsigned int *numOfHits);

//void CleanDecoder();

double getTextPositionTime();
unsigned int getSARetrieved();
unsigned int getHASHRetrieved();

#endif /*_EXTRATOOLS_H_*/

