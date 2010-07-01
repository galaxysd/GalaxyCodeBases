// vim:set autoindent shiftwidth=4 tabstop=4 noexpandtab:
#include "lookupBuilder.h"

int LEN;
int TOP;
ULL NR_TOP;
ULL ALL_MASK;

const unsigned IN_BUF_SIZE = 100u Mibi;
const unsigned step = 1 Mibi;

unsigned char * ibuf;

ULL text_pos; // Text position, 0-based
ULL window; // Last length-LEN characters

int fin_src, fout_int;
long long start_time;

unsigned * otop;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void process(int v) {
	++text_pos;
	//if ((text_pos & ((1ull << 27) - 1)) == 0)
		//printf("Time %4llds:     text_pos=%lluMi\n", time(0) - start_time, text_pos / 1048576);
	window = (window << 2) | v;
	if (text_pos >= -((ULL)LEN)) return; // First LEN-1 characters
	unsigned top = (window & ALL_MASK);
	++otop[top];
}

void gen() {
	// init
	text_pos = -(ULL)LEN;
	window = 0;
	// start
	TRY(lseek(fin_src, 0, SEEK_SET));
	while (1) {
        unsigned int i;
		unsigned nbuf = read(fin_src, ibuf, sizeof(*ibuf) * IN_BUF_SIZE);
		if (nbuf < IN_BUF_SIZE) nbuf -= 2;
		for (i = 0; i < nbuf; ++i) {
			unsigned short c = ibuf[i];
			int j;
			for (j = 6; j >= 0; j -= 2) process((c >> j) & 3);
		}
		if (nbuf < IN_BUF_SIZE) {
			int rem = ibuf[nbuf + 1];
			unsigned short c = ibuf[nbuf];
			int j;
			for (j = 6; j >= 8 - (rem * 2); j -= 2) process((c >> j) & 3);
			break;
		}
	}
}



int buildLookupTable(char * PackedDNAFileName, char * LookupTableFileName, int lookupTableSize) {
    
    ibuf = malloc(sizeof(unsigned char)*IN_BUF_SIZE);
    LEN = lookupTableSize;
    TOP = lookupTableSize;
    start_time = time(0);
    NR_TOP = 1 << TOP * 2;
    ALL_MASK = NR_TOP - 1;
    otop = malloc(sizeof(unsigned)*NR_TOP);
    unsigned int i;


	TRY(fin_src = open(PackedDNAFileName, O_RDONLY ));
	TRY(fout_int = open(LookupTableFileName, O_CREAT | O_WRONLY , 0664));
	gen();
	for (i = 1; i < LEN; ++i) process(0);
	for (i = 1; i < NR_TOP; ++i) otop[i] += otop[i-1];
	for (i = 0; i < NR_TOP; i += step) {
		TRYEQ(write(fout_int, otop + i, step * sizeof(*otop)), step * sizeof(*otop));
	}
	TRY(close(fin_src));
	TRY(close(fout_int));


	free(ibuf);
	free(otop);
}

/*

void process(int v) {
	++text_pos;
	if ((text_pos & ((1ull << 27) - 1)) == 0)
		printf("Time %4llds:     text_pos=%lluMi\n", time(0) - start_time, text_pos / 1048576);
	window = (window << 2) | v;
	if (text_pos >= -((ULL)LEN)) return; // First LEN-1 characters
	//ULL curr = rotation ? ((window & ALL_MASK) >> BOTTOM * 2) | (window << (TOP * 2)) : window;
	unsigned top = (window & ALL_MASK);
	++otop[top];
}

void gen() {
	// init
	text_pos = -(ULL)LEN;
	window = 0;
	// start
	TRY(fseek(fin_src, 0, SEEK_SET));
	while (1) {
        unsigned int i;
		unsigned nbuf = fread(ibuf,sizeof(unsigned char) * IN_BUF_SIZE,1,fin_src);
		if (nbuf < IN_BUF_SIZE) nbuf -= 2;
		for (i = 0; i < nbuf; ++i) {
			unsigned short c = ibuf[i];
			int j;
			for (j = 6; j >= 0; j -= 2) process((c >> j) & 3);
		}
		if (nbuf < IN_BUF_SIZE) {
			int rem = ibuf[nbuf + 1];
			unsigned short c = ibuf[nbuf];
			int j;
			for (j = 6; j >= 8 - (rem * 2); j -= 2) process((c >> j) & 3);
			break;
		}
	}
}



int buildLookupTable(char * PackedDNAFileName, char * LookupTableFileName, int lookupTableSize) {
    
    ibuf = malloc(sizeof(unsigned char)*IN_BUF_SIZE);
    LEN = lookupTableSize;
    TOP = lookupTableSize;
    start_time = time(0);
    NR_TOP = 1 << TOP * 2;
    ALL_MASK = NR_TOP - 1;
    otop = malloc(sizeof(unsigned)*NR_TOP);
    unsigned int i;
	TRY(fin_src = fopen(PackedDNAFileName, "rb"));
	TRY(fout_int = fopen(LookupTableFileName, "wb"));
	printf("MK\n");
	gen();
	for (i = 1; i < LEN; ++i) process(0);
	for (i = 1; i < NR_TOP; ++i) otop[i] += otop[i-1];
	const unsigned step = 1 Mibi;
	for (i = 0; i < NR_TOP; i += step) {
		TRYEQ(fwrite(otop + i,step * sizeof(*otop),1,fout_int), step * sizeof(*otop));
	}
	TRY(fclose(fin_src));
	TRY(fclose(fout_int));
	free(ibuf);
	free(otop);
}*/
