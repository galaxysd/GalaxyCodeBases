#include "extratools.h"

//This file includes the implementations of all the extra tools adding to all steps
//e.g. Look Up Table
//     Hash Table
//     All things like those

#define gen_occ_arr(l, r, occ, s) do{ \
	if ((r)-(l)+1 >= 4) {	\
		HASHITEM *item = HashFind(occHashtable,(l),(r));	\
		if (item==NULL) {	\
			unsigned int k;	\
			for (k=(l);k<=(r);k++,occ++) {	\
				occ->tp = BWTSaValue(occBwt,k);	\
				occ->seg = occ->seg&0|s;			\
			}\
		} else {\
			unsigned int k;\
			for (k=0;k<item->(r)-item->(l)+1;k++,occ++) {\
				occ->tp= occHashtable->occList[item->occIndex+k];\
				occ->seg = occ->seg&0|s;			\
			}\
		}\
	} else {\
		unsigned int k;\
		for (k=(l);k<=(r);k++,occ++) {\
			occ->tp = BWTSaValue(occBwt,k);\
			occ->seg = occ->seg&0|s;			\
		}\
	}\
}while(0);

void LoadLookupTable(LOOKUPTABLE * lookupTable, const char * fileName, const int tableSize)  {
    (*lookupTable).tableSize = tableSize;
	unsigned long long NR_TOP = 1 << (tableSize * 2);
	(*lookupTable).table = malloc(sizeof(unsigned) * NR_TOP);
	int fin = open(fileName, O_RDONLY);
	unsigned step = 1048576;
	unsigned int i;
	for (i = 0; i < NR_TOP; i += step) {
		read(fin, (*lookupTable).table + i, step * sizeof(*(*lookupTable).table));
	}
	close(fin);
}

unsigned int LookupSafe(LOOKUPTABLE lookupTable, BWT * bwt,
                        unsigned long long lKey, unsigned long long rKey,
                        unsigned int *l, unsigned int *r) {

	*l = lKey ? lookupTable.table[lKey-1]+1 : 1;
	*r = lookupTable.table[rKey];

	if (*l == bwt->inverseSa0) {
        *l++;
    }

	return *r-*l+1;
}

unsigned int retrieveSA=0,retrieveHASH=0;
double textPositionTime, textPositionTimeTotal = 0;
unsigned int writeQIndex;

double getTextPositionTime() {return textPositionTimeTotal;}
unsigned int getSARetrieved() {return retrieveSA;}
unsigned int getHASHRetrieved() {return retrieveHASH;}

void FreeLookupTable(LOOKUPTABLE * lookupTable) {
     free((*lookupTable).table);
}

void LoadHashTable(HASHTABLE * hashTable, const char * fileName) {
unsigned int ttlOccurrence=0;
unsigned int ttlItem=0;

    FILE *inFile;
    if(!(inFile = fopen(fileName, "r"))) return;
    fread((unsigned int *)&((*hashTable).tableSize),sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&((*hashTable).a),sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&((*hashTable).b),sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&((*hashTable).prime),sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&ttlItem,sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&ttlOccurrence,sizeof(unsigned int),1,inFile);

    //printf("Initializing the hash table..(n=%u)\n",(*hashTable).tableSize);
    (*hashTable).table = (HASHCELL*) malloc(sizeof(HASHCELL)*((*hashTable).tableSize));
    (*hashTable).itemList = (HASHITEM*) malloc(sizeof(HASHITEM)*ttlItem);
    (*hashTable).occList = (OCC*) malloc(sizeof(OCC)*ttlOccurrence);
    //printf("Initialized the hash table..\n");
    unsigned int i;
    for (i=0;i<((*hashTable).tableSize);i++) {
        char mk;
        fread((char *) &mk,1,1,inFile);
        if (mk==0) {
            //Empty cell
            (*hashTable).table[i].index=0;
            (*hashTable).table[i].count=0;
        } else {
            fread((unsigned int *)&((*hashTable).table[i].index),sizeof(unsigned int),1,inFile);
            fread((unsigned int *)&((*hashTable).table[i].count),sizeof(unsigned int),1,inFile);
        }
    }
    for (i=0;i<ttlItem;i++) {
        fread((unsigned int *)&((*hashTable).itemList[i].l),sizeof(unsigned int),1,inFile);
        fread((unsigned int *)&((*hashTable).itemList[i].r),sizeof(unsigned int),1,inFile);
        fread((unsigned int *)&((*hashTable).itemList[i].occIndex),sizeof(unsigned int),1,inFile);
    }
    for (i=0;i<ttlOccurrence;i++) {
        fread((unsigned int *)&((*hashTable).occList[i]),sizeof(unsigned int),1,inFile);
    }
    fclose(inFile);

}

void FreeHashTable(HASHTABLE * hashTable) {
     free((*hashTable).table);
     free((*hashTable).itemList);
     free((*hashTable).occList);
}

unsigned int Hash(HASHTABLE * hashTable, unsigned int key) {
    //g(x)=(ax+b) mod p
    unsigned long long multipleA=(key * (*hashTable).a)% (*hashTable).prime;
    unsigned int g=(unsigned int) (( multipleA + (*hashTable).b ) % (*hashTable).prime);

    //f(x)=g(x) mod 2n
    unsigned int f=g % ((*hashTable).tableSize);

    return f;
}

HASHITEM * HashFind(HASHTABLE * hashTable, unsigned int l,unsigned int r) {
    unsigned int hashedIndex=Hash(hashTable,l);
    unsigned int index=(*hashTable).table[hashedIndex].index;
    unsigned int count=(*hashTable).table[hashedIndex].count;

    if (hashedIndex>=0 && hashedIndex<(*hashTable).tableSize) {
        unsigned int i=0;
        while (i<count && ((*hashTable).itemList[index+i].l!=l || (*hashTable).itemList[index+i].r!=r)) {
            i++;
        }
        if (i>=count) return NULL;
        return &((*hashTable).itemList[index+i]);
    }
    return NULL;
}

void RegisterDecoder(BWT * bwt,HASHTABLE * hashTable) {
     occBwt=bwt;
     occHashtable=hashTable;
     occCollected=0;
     retrieveSA=0;
     retrieveHASH=0;
     textPositionTimeTotal = 0;
     writeQIndex = 0;
     //occCollector = malloc(sizeof (unsigned int) * 1024*1024);
}

/*void CleanDecoder() {
     //OCCClean();
     //free(occCollector);
}*/


unsigned int allOne = 0;
unsigned int OCCSection=0;
void OCCProcess(const unsigned int l,const unsigned int r, const unsigned int aux_info, match_aux_t *o, match_table_t *matches) {
	if (r-l+1 >= 4) {
		//Hash
		HASHITEM *item = HashFind(occHashtable,l,r);
		if (item==NULL) {
			//printf("Hash entry is missing %u-%u!\n",occL,occR);
			unsigned int k;
			unsigned int tp;
			for (k=l;k<=r;k++) {
				orient_pos(aux_info, BWTSaValue(occBwt,k), matches, o);
			}
		} else {
			unsigned int k;
			unsigned int tp;
			for (k=0;k<item->r-item->l+1;k++) {
				orient_pos(aux_info, occHashtable->occList[item->occIndex+k], matches, o);
			}
		}
	} else {
		//SA
		unsigned int k;
		unsigned int tp;
		for (k=l;k<=r;k++) {
			orient_pos(aux_info, BWTSaValue(occBwt,k), matches, o);
		}
	}
}

//Update hit informations
/*
   void OCCProcess(unsigned int l,unsigned int r, ){
//	fprintf(stderr, "OCC err;%d\n", *numOfHits);

unsigned int mask = (unsigned int) (-1);
int nhits = *numOfHits;
if (r-l+1 >= 4) {
//Hash
HASHITEM *item = HashFind(occHashtable,l,r);
if (item==NULL) {
unsigned int k;
unsigned int tp;
for (k=l; k<=r && nhits<MAX_HITS_NUM; k++,nhits++) {
tp = BWTSaValue(occBwt,k);
//					printf ("tp: %d\n", tp);
hits[nhits].pos = tp;
hits[nhits].chrID = -1;
hits[nhits].chain = chain;
hits[nhits].type = 0;
					hits[nhits].mm_n = 0;
				}
			} else if (allele2 == mask) {
				for (k=l;k<=r && nhits<MAX_HITS_NUM;k++,nhits++) {
					tp = BWTSaValue(occBwt,k);
					hits[nhits].pos = tp;
					hits[nhits].chrID = -1;
					hits[nhits].chain = chain;
					hits[nhits].l = allele1;
					hits[nhits].type = 1;
					hits[nhits].mm_n = 1;
				}
			} else {
				for (k=l; k<=r && nhits<MAX_HITS_NUM; k++,nhits++) {
					tp = BWTSaValue(occBwt,k);
					hits[nhits].pos = tp;
					hits[nhits].chrID = -1;
					hits[nhits].chain = chain;
					hits[nhits].l = allele1;
					hits[nhits].r = allele2;
					hits[nhits].type = 2;
					hits[nhits].mm_n = 2;
				}
			}
			retrieveSA+=r-l+1;
		} else {
			unsigned int k;
			unsigned int tp;
			if (allele1 == allele2) {
				for (k=0;k<item->r-item->l+1 && nhits<MAX_HITS_NUM; k++,nhits++) {
					tp = occHashtable->occList[item->occIndex+k];
//					printf ("tp: %d\n", tp);
					hits[nhits].pos = tp;
					hits[nhits].chrID = -1;
					hits[nhits].chain = chain;
					hits[nhits].type = 0;
					hits[nhits].mm_n = 0;
				}
			} else if (allele2 == mask) {
				for (k=0;k<item->r-item->l+1&& nhits<MAX_HITS_NUM; k++,nhits++) {
					tp = occHashtable->occList[item->occIndex+k];
					hits[nhits].pos = tp;
					hits[nhits].chrID = -1;
					hits[nhits].chain = chain;
					hits[nhits].l = allele1;
					hits[nhits].type = 1;
					hits[nhits].mm_n = 1;
				}
			} else {
				for (k=0;k<item->r-item->l+1 && nhits<MAX_HITS_NUM; k++,nhits++) {
					tp = occHashtable->occList[item->occIndex+k];
					hits[nhits].pos = tp;
					hits[nhits].chrID = -1;
					hits[nhits].chain = chain;
					hits[nhits].l = allele1;
					hits[nhits].r = allele2;
					hits[nhits].type = 2;
					hits[nhits].mm_n = 2;
				}
			}
			retrieveHASH+=(*item).r-(*item).l+1;
		}
	} else {
		//SA
		unsigned int k;
		unsigned int tp;
		if (allele1 == allele2){
			for (k=l; k<=r && nhits<MAX_HITS_NUM; k++,nhits++) {
//				printf("k: %u\n", k);
				tp = BWTSaValue(occBwt,k);
//				printf ("tp: %d\n", tp);
				hits[nhits].pos = tp;
				hits[nhits].chrID = -1;
				hits[nhits].chain=chain;
				hits[nhits].type = 0;
				hits[nhits].mm_n = 0;
			}
		} else if (allele2 == mask) {
			for (k=l; k<=r && nhits<MAX_HITS_NUM; k++,nhits++) {
				tp = BWTSaValue(occBwt,k);
				hits[nhits].pos = tp;
				hits[nhits].chrID = -1;
				hits[nhits].chain=chain;
				hits[nhits].l = allele1;
				hits[nhits].type = 1;
				hits[nhits].mm_n = 1;
			}
		} else {
			for (k=l; k<=r && nhits<MAX_HITS_NUM; k++,nhits++) {
				tp = BWTSaValue(occBwt,k);
				hits[nhits].pos = tp;
				hits[nhits].chrID = -1;
				hits[nhits].chain=chain;
				hits[nhits].l = allele1;
				hits[nhits].r = allele2;
				hits[nhits].type = 2;
				hits[nhits].mm_n = 2;
			}
		}
		retrieveSA+=r-l+1;
	}
	*numOfHits = nhits;
//	printf ("nhits: %u\n", nhits);
//	printf ("pos: %u\n", hits[nhits-1].pos);
}
//*/

void registerTPFile(FILE * filePtr,unsigned int searchMode) {
    textPositionFile=filePtr;
    fwrite(&searchMode,sizeof(unsigned int),1,textPositionFile);
	allOne=(1U<<31)-1;
	allOne<<=1;
	allOne+=1;
}

void registerQIndex(unsigned int index) {
    writeQIndex=index;
    OCCSection=0;
}
void registerQSection() {
    if (writeQIndex==0) {
        fwrite(&allOne,sizeof(unsigned int),1,textPositionFile);
    } else {
        OCCSection++;
    }
}
