
typedef struct
{
	char *name;
	char *name2;
	char *path;
	char *SCFname;
	int  length;
	char *data;
	char *qual;
	int  finished;
} fasta;

fasta *decodeFastq (char *fname, int *nContigs, long *tB, int qT);
fasta *splitFastq ( fasta *iseg, int inSeg, int tB, int *nReads, int length, int step);
void fastaLC (fasta *seg, int nSeg);
void fastaUC (fasta *seg, int nSeg);

