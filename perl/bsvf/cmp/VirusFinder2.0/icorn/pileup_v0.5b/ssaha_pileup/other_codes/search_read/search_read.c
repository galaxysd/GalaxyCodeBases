/***********************************************************************\
 *                                                                     * 
 *                     PROJECT   ssaha_Assembly                        *
 *                                                                     * 
 *---------------------------------------------------------------------*
 *                  Assembly Program using SSAHA                       *
 *                                                                     * 
 *      Sequence Search and Alignment by the Hashing Algorithm         *
 *                                                                     *
 *                                By                                   *
 *                                                                     *
 *                    Zemin Ning & James C. Mullikin                   *
 *                                                                     *
 *          Copyright (C) 1999-2000 by Genome Research Limited         *
 *                       All rights reserved                           *
 *                                                                     *
 *---------------------------------------------------------------------*
 #######################################################################
 #                                                                     #
 #             <------   LICENCE NOTICE   ------>                      #
 #                                                                     #
 # This is a licensed software developed by Genome Research Limited    #
 # (GRL) for genomic sequence assembling. For both commercial and non- # 
 # commercial purposes, a licence must be obtained from GRL before     #
 # use. Under no circumstances, the users should copy, modify and      #
 # redistribut the software as a whole or even part of it without      #
 # permission from GRL. For more information about the software and    #
 # its ducumentation particularly, please contact either one of the    # 
 # authors or GRL.                                                     #
 #######################################################################
 *---------------------------------------------------------------------*/

#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 40000 
#define nfm 800000
#define nfm_sub 500000
#define Max_N_NameBase 60
#define Max_N_Pair 100
static int *readlength;
static char *dna,**S_Name,**rdnames;
static long *h_dna;
static int *cons_head,*cons_list,*rd_head,*rd_index;
static int *read_rept,*read_done,*hitlink,*hitlist,*replink;
static int *ctgend_RT,*ctgend_LT,*bdgread,*ctgindex,*ctgoffset;
static int *id_bdglist,*id_hitlist,*id_hitlink,*read_occur,*id_occur;
static int *km2_index,*km2_indexRC,*id_km2,*id_km2RC,*readorder,*readmatch;
static int *rd_hitoffset,*rd_hitindex,*rd_hitorder,*rd_hitmatch;
static int *ctg2wgs_index,*dup_reads;
static int *snp_rdindex,*snp_roffset,*snp_rdirect,*snp_coindex,*snp_coffset;
static char *data,*snp_base;
static long n_Entry,n_consEntry;
static long MAX_S_LEN,MAX_S_LEN2,MAX_S_LEN3;

/* SSAS default parameters   */
static int K_LEN=12;
static int K_CON=12;
static int IMOD=0;
static int ISUB=0;
static int NCUT=50;
static int NHIT=10;
static int k_skip=12;
static int all_out=0;
static int tmod=0;
static int lowHit=4;
static int ReadDepth=15;
static int n_coverage=10;
static int insertSize=3000;
static int smod=1;
static int dir_flag=0;
static int replace_flag=1;
static int NUM_SNP=500000;
static int N_outread=400000;
static int fastq_flag=0;
static int release_flag=0;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

static fasta **expp;
fasta *expt;

static char rc_char[500000];
static char rc_sub[5000];

int ReverseComplement(int seqdex)
{
        int i,len;
        char *tp,*dp;
        fasta *seqp;

        seqp=expt+seqdex;
        len=seqp->length;
        memset(rc_sub,'\0',5000);
        dp=rc_sub;      
        tp = seqp->data+len;
        for(i=len;--i>=0;)
        {
                int tmp = *--tp;
                if     (tmp == 't') *dp++ = 'a';
                else if(tmp == 'g') *dp++ = 'c';
                else if(tmp == 'c') *dp++ = 'g';
                else if(tmp == 'a') *dp++ = 't';
                else                *dp++ = tmp;
        }
        return(0);
}

int Reverse_Complement(int seqdex)
{
	int i,len;
	char *tp,*dp;

        len=readlength[seqdex];
        dp=rc_char;	
	tp = dna+h_dna[seqdex]+len;
	for(i=len;--i>=0;)
	{
		int tmp = *--tp;
		if     (tmp == 't') *dp++ = 'a';
		else if(tmp == 'g') *dp++ = 'c';
		else if(tmp == 'c') *dp++ = 'g';
		else if(tmp == 'a') *dp++ = 't';
		else                *dp++ = tmp;
	}
	return(0);
}

int Reverse_Complement_Contig(char c_array[],int num_len)
{
        int i,len;
        char *tp,*dp;

        len=num_len;
        dp=rc_char;
        tp = c_array+len;
        for(i=len;--i>=0;)
        {
                int tmp = *--tp;
                if     (tmp == 't') *dp++ = 'a';
                else if(tmp == 'g') *dp++ = 'c';
                else if(tmp == 'c') *dp++ = 'g';
                else if(tmp == 'a') *dp++ = 't';
                else                *dp++ = tmp;
        }
        return(0);
}


int main(int argc, char **argv)
{
    FILE *namef,*fpOutfast;
    long dataSize,totalBases;
    int i,nSeq,args,nRead,qthresh=0,n_len;
    int ac,n_contig,n_reads,nseq,n_input;
    fasta *seq,*seqp;
    void decodeReadpair(int nSeq);
    void HashFasta_Head(int i, int nSeq);
    void HashFasta_Table(int i, int nSeq);
    void Search_SM(fasta *seq,int nSeq);
    void Read_Index(char **argv,int args,int nRead,int nSeq);
    void Memory_Allocate(int arr);
    char *st,line[2000]={0},Nameout[100]={0};
    char **cmatrix(long nrl,long nrh,long ncl,long nch);
    int  num_out,num_file,out_flag;

    seq=NULL;
    fflush(stdout);
    system("ps aux | grep sreads; date");
    if(argc < 2)
    {
      printf("Usage: %s <read name file> <output fastq head> <tagname file> <fastq files>\n",argv[0]);
      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
       {
         sscanf(argv[++i],"%d",&IMOD); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-sub"))
       {
         sscanf(argv[++i],"%d",&ISUB); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-cut"))
       {
         sscanf(argv[++i],"%d",&NCUT); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-hit"))
       {
         sscanf(argv[++i],"%d",&NHIT); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-len"))
       {
         sscanf(argv[++i],"%d",&K_LEN);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-con"))
       {
         sscanf(argv[++i],"%d",&K_CON);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-skip"))
       {
         sscanf(argv[++i],"%d",&k_skip);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-allout"))
       {
         sscanf(argv[++i],"%d",&all_out);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-tmod"))
       {
         sscanf(argv[++i],"%d",&tmod);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-lowhit"))
       {
         sscanf(argv[++i],"%d",&lowHit);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-cover"))
       {
         sscanf(argv[++i],"%d",&n_coverage);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-insert"))
       {
         sscanf(argv[++i],"%d",&insertSize);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-dir"))
       {
         sscanf(argv[++i],"%d",&dir_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-replace"))
       {
         sscanf(argv[++i],"%d",&replace_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-depth"))
       {
         sscanf(argv[++i],"%d",&ReadDepth);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-set"))
       {
         sscanf(argv[++i],"%d",&N_outread);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-smod"))
       {
         sscanf(argv[++i],"%d",&smod);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-fastq"))
       {
         sscanf(argv[++i],"%d",&fastq_flag);
         args=args+2;
       }
    }

    nseq=0;
    if((namef = fopen(argv[args+2],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      nseq++;
    }
    fclose(namef); 

    if((rd_index= (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - rd_index\n");
      exit(1);
    }
    if((rd_head= (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - rd_head\n");
      exit(1);
    }

    expp = (fasta **) malloc(argc*sizeof(fasta *));
    n_input=0;
    nseq=0;
    
    for(ac=(args+3);ac<argc;ac++)
    {
       seq = decodeFastq ( argv[ac], &nSeq, &totalBases, qthresh);
       if(seq == NULL)
       {
               printf("ERROR fmate: no data found\n");
               exit(1);
       }
       fastaUC(seq,nSeq);
       expp[n_input]=seq;
       for(i=0;i<nSeq;i++)
       {
          rd_head[i+nseq]=n_input;
          rd_index[i+nseq]=i;
       }
       nseq=nseq+nSeq;
       n_input++;
       free(data);
    }

    nSeq=nseq;
    nRead=nseq;
    S_Name=cmatrix(0,nRead+1,0,Max_N_NameBase);
    if((dup_reads= (int *)calloc(nRead,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - dup_reads\n");
       exit(1);
    }

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the read.place files         */
    while(!feof(namef))
    {
      int nPair=0;
      char *ptr,line[500];

      fgets(line,500,namef);
      for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
      {
      }
      printf("colums: %d\n",nPair);
      if(nPair==9)
        release_flag=1;
      else if(nPair==1)
        release_flag=2;
      else
        release_flag=0;
      break; 
    }
    fclose(namef);

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    n_reads=0;
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      n_reads++;
    }
    fclose(namef); 

    printf("reads: %d %s\n",n_reads,argv[args+1]);
    if((ctg2wgs_index= (int *)calloc(n_reads,sizeof(int))) == NULL) 
    {
       printf("ERROR Contig_Hist: calloc - ctg2wgs_index\n");
       exit(1);
    }
    memset(ctg2wgs_index,-1,n_reads*sizeof(int));
    Read_Index(argv,args,n_reads,nSeq);

    n_Entry=1<<(K_LEN<<1);

    printf("Number of shotgun reads  %d \n",nRead);

/*  define the sequence matrix SM     */
    dataSize=1000;
    MAX_S_LEN=(dataSize/K_LEN)+100;
    MAX_S_LEN2=MAX_S_LEN<<1;
    MAX_S_LEN3=MAX_S_LEN+MAX_S_LEN2;
    n_consEntry=1<<(K_CON<<1);


/*  process contigs one by one   */
    ac=args+1;
    printf("file name  %s \n",argv[ac]);
    n_contig=0;
    num_out=0;
    num_file=0;
    out_flag=1;
    fpOutfast = NULL;
    for(i=0;i<n_reads;i++)
    {
       if(((num_out%N_outread)==0)&&(out_flag==1))
       {
         memset(Nameout,'\0',100);
         sprintf(Nameout,"%s_%03d%s",argv[ac],num_file,".fastq");
         if((fpOutfast = fopen(Nameout,"w")) == NULL)
         {
             printf("ERROR main:: reads group file \n");
             exit(1);
         }
         out_flag=0;
         num_file++;
       }
       if(ctg2wgs_index[i]>=0) 
       {
          int nline,k,g,rc;
          int idd=ctg2wgs_index[i];
          n_contig++;
          seqp=expp[rd_head[idd]]+rd_index[idd];
//          printf("name: %d %s\n",idd,seqp->name);
          n_len=seqp->length;
          if(fastq_flag==0)
          {
            fprintf(fpOutfast,">%s\n",seqp->name);
            nline=n_len/60;
            st=seqp->data;
            for(k=0;k<nline;k++)
            {
               for(g=0;g<60;g++,st++)
                  fprintf(fpOutfast,"%c",*st);
               fprintf(fpOutfast,"\n");
            }
            for(g=0;g<(n_len-(nline*60));g++,st++)
               fprintf(fpOutfast,"%c",*st);
            if((n_len%60)!=0)
              fprintf(fpOutfast,"\n"); 
          }
          else
          {
            if(seqp->name2)
              fprintf(fpOutfast,"@%s %s\n",seqp->name,seqp->name2);
            else
              fprintf(fpOutfast,"@%s\n",seqp->name);
            for(rc=0;rc<seqp->length;rc++)
               fprintf(fpOutfast,"%c",seqp->data[rc]);
            fprintf(fpOutfast,"\n");
            if(seqp->name2)
              fprintf(fpOutfast,"+%s %s\n",seqp->name,seqp->name2);
            else
              fprintf(fpOutfast,"+%s\n",seqp->name);
            for(rc=0;rc<seqp->length;rc++)
            {
               if(seqp->finished)
               {
                 if((seqp->data[rc]=='N')||(seqp->data[rc]=='X'))
                   putc(2+041,fpOutfast);
                 else
                   putc(40+041,fpOutfast);
               }
               else
                 putc(seqp->qual[rc]+041,fpOutfast);
            }
            fprintf(fpOutfast,"\n");
          }
          out_flag=1;
          num_out++;
       }
       else
       {
            printf("missed: %d %s\n",i,rdnames[i]);
       }
       if(((num_out%N_outread)==0)&&(out_flag==1))
         fclose(fpOutfast);
    }
    fclose(fpOutfast);

    if(seq){
        free(seq->name);
        free(seq);
        seq = NULL;
    }    
    printf("Job finished for %d contigs!\n",n_contig);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read index    */
/* =============================== */
void Read_Index(char **argv,int args,int nRead,int nSeq)
/* =============================== */
{
     int i,j,k,n_reads=nSeq+nRead;
     FILE *namef;
     fasta *seqp;
     void ArraySort_String(int n,char **Pair_Name,int *brr);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     char **DBname,dupname[Max_N_NameBase];
     char tempct[60],tempc1[60],tempc2[60],tempc3[60],tempc4[60],tempc5[60];
     int temp1,temp2,temp3;
     int i_contig,i_reads,num_rd_find,stopflag;
     int *readIndex;
     int mapindex=0,n_dupreads,num_nodup;

     if((readIndex= (int *)calloc(n_reads,sizeof(int))) == NULL)
     {
       printf("Error Contig_Merge: calloc - readIndex\n");
       exit(1);
     } 
     DBname=cmatrix(0,n_reads,0,Max_N_NameBase);
     rdnames=cmatrix(0,nRead+1,0,Max_N_NameBase);
     i_contig=0;
     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR Memory_Allocate:: reads group file \n");
       exit(1);
     }

     i_reads=0;
     printf("before read: %d %d\n",release_flag,nRead);
     if(release_flag==0)
     {
       while(fscanf(namef,"%s %s %s %d %s %s %d %d",tempc1,tempc2,tempc3,&temp1,tempc4,tempc5,&temp2,&temp3)!=EOF)
       {
         strcpy(rdnames[i_reads],tempc1);
         i_reads++;
       }
     }
     else if(release_flag==1)
     {
       while(fscanf(namef,"%s %s %s %s %d %s %s %d %d",tempct,tempc1,tempc2,tempc3,&temp1,tempc4,tempc5,&temp2,&temp3)!=EOF)
       {
         strcpy(rdnames[i_reads],tempc1);
         i_reads++;
       }
     }
     else if(release_flag==2)
     {
       while(fscanf(namef,"%s",rdnames[i_reads])!=EOF)
       {
//     printf("after read: %d %s\n",i_reads,rdnames[i_reads]);
         i_reads++;
       }
     }
     fclose(namef);

     for(j=0;j<nSeq;j++)
     {
        seqp=expp[rd_head[j]]+rd_index[j];
        strcpy(DBname[j],seqp->name);
        readIndex[j]=j;
     }
     n_reads=nSeq;
     ArraySort_String(n_reads,DBname,readIndex);

     n_dupreads=0;
     memset(dup_reads,0,nSeq*sizeof(int));
     for(i=0;i<n_reads-1;i++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        while((j<n_reads)&&(stopflag==0))
        {
          if(strcmp(DBname[j],DBname[i])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>1)
        {
          for(k=(i+1);k<j;k++)
          {
             dup_reads[readIndex[k]] = 1;
             n_dupreads++;
          }
        }
        i=j-1;
     }

     printf("number of dup reads: %d %d\n",nSeq,n_dupreads);
     num_nodup = nSeq;
     for(j=0;j<nSeq;j++)
     {
        if(dup_reads[j]!=1)
        {
          seqp=expp[rd_head[j]]+rd_index[j];
          strcpy(DBname[j],seqp->name);
        }
        else
        {
          sprintf(dupname,"%s%09d","dupname-",j);
          strcpy(DBname[j],dupname);
        }
        readIndex[j]=j;
     }
     for(j=0;j<nRead;j++)
     {
        strcpy(DBname[j+num_nodup],rdnames[j]);
        readIndex[j+num_nodup]=j+num_nodup;
     }
     n_reads=num_nodup+nRead;
     printf("before sort: %d %d\n",num_nodup,nRead);
     ArraySort_String(n_reads,DBname,readIndex);

     num_rd_find=0;
     mapindex=0;
     for(i=0;i<n_reads-1;i++)
     {
        if(readIndex[i]>=num_nodup)
          mapindex = readIndex[i];
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        while((j<n_reads)&&(stopflag==0))
        {
          if(strcmp(DBname[j],DBname[i])==0)
          {
            if(readIndex[j]>=num_nodup)
              mapindex = readIndex[j];
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>=2)
        {
          for(k=i;k<j;k++)
          {
             if(readIndex[k]<num_nodup)
             {
               ctg2wgs_index[mapindex-num_nodup]=readIndex[k];
               num_rd_find++;
               k=j;
             }
          }
        }
        i=j-1;
     }
       printf("reads found: %d %d\n",nRead,num_rd_find);
}

/*   subroutine to allocatre memory for Assemb_SM    */
/* =============================== */
void Memory_Allocate(int nSeq)
/* =============================== */
{
     if((read_rept= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - read_rept\n");
       exit(1);
     }
     if((read_done= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - read_done\n");
       exit(1);
     }
     if((hitlink= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - hitlink\n");
       exit(1);
     }
     if((hitlist= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - hitlist\n");
       exit(1);
     }
     if((replink= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - replink\n");
       exit(1);
     }
     if((bdgread= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - bdgread\n");
       exit(1);
     }
     if((ctgindex= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ctgindex\n");
       exit(1);
     }
     if((ctgoffset= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ctgoffset\n");
       exit(1);
     }
     if((ctgend_RT= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ctgend_RT\n");
       exit(1);
     }
     if((ctgend_LT= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ctgend_LT\n");
       exit(1);
     }
     if((readorder= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - readorder\n");
       exit(1);
     }
     if((readmatch= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - readmatch\n");
       exit(1);
     }
     if((rd_hitorder= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - rd_hitorder\n");
       exit(1);
     }
     if((rd_hitmatch= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - rd_hitmatch\n");
       exit(1);
     }
     if((rd_hitindex= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - rd_hitindex\n");
       exit(1);
     }
     if((rd_hitoffset= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - rd_hitoffset\n");
       exit(1);
     }
     if((id_hitlist= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - id_link\n");
       exit(1);
     }
     if((id_bdglist= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - id_link\n");
       exit(1);
     }
     if((id_hitlink= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - id_hitlink\n");
       exit(1);
     }
     if((km2_index= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - km2_index\n");
       exit(1);
     }
     if((km2_indexRC= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - km2_indexRC\n");
       exit(1);
     }
     if((id_km2= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - id_km2\n");
       exit(1);
     }
     if((id_km2RC= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - id_km2RC\n");
       exit(1);
     }
     if((read_occur= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - read_occur\n");
       exit(1);
     }
     if((id_occur= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - id_occur\n");
       exit(1);
     }
/*   allocate memory for arrays HEAD and LIST   */
     if((cons_head = (int *)calloc(n_consEntry,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - cons_head\n");
       exit(1);
     }
     if((cons_list = (int *)calloc(n_consEntry,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - cons_list\n");
       exit(1);
     }
     if((snp_rdindex = (int *)calloc(NUM_SNP,sizeof(int))) == NULL)
     {
       printf("ERROR ssaha: calloc - snp_rdindex\n");
       exit(1);
     }
     if((snp_roffset = (int *)calloc(NUM_SNP,sizeof(int))) == NULL)
     {
       printf("ERROR ssaha: calloc - snp_roffset\n");
       exit(1);
     }
     if((snp_rdirect = (int *)calloc(NUM_SNP,sizeof(int))) == NULL)
     {
       printf("ERROR ssaha: calloc - snp_rdirect\n");
       exit(1);
     }
     if((snp_coindex = (int *)calloc(NUM_SNP,sizeof(int))) == NULL)
     {
       printf("ERROR ssaha: calloc - snp_coindex\n");
       exit(1);
     }
     if((snp_coffset = (int *)calloc(NUM_SNP,sizeof(int))) == NULL)
     {
       printf("ERROR ssaha: calloc - snp_coffset\n");
       exit(1);
     }
     if((snp_base = (char *)calloc(NUM_SNP,sizeof(char))) == NULL)
     {
       printf("ERROR ssaha: calloc - snp_base\n");
       exit(1);
     }
}

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Int(int n, int *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Mix(int n, long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Mix3(int n, long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             c=crr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
                crr[i+1]=crr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
             crr[i+1]=c;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);
          SWAP(crr[k],crr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
            SWAP(crr[m],crr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
            SWAP(crr[m+1],crr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
            SWAP(crr[m],crr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          c=crr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
             SWAP(crr[i],crr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          crr[m+1]=crr[j];
          crr[j]=c;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

