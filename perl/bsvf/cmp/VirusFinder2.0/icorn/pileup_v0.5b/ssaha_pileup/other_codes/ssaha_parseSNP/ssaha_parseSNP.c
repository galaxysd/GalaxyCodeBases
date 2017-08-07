/***********************************************************************\
 *                                                                     * 
 *                         PROJECT   ssahaSNP                          *
 *                                                                     * 
 *---------------------------------------------------------------------*
 *                             parse_SNP                               *
 *                                                                     *
 *                                By                                   *
 *                                                                     *
 *                            Zemin Ning                               *
 *                                                                     *
 *          Copyright (C) 2005 by Genome Research Limited              *
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

/****************************************************************************/

#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <sys/wait.h>
#include <sys/signal.h>
#include <errno.h>
#include "fasta.h"

#define MAXLINE 4096
#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 40000 
#define Max_N_NameBase 60
#define Max_N_Pair 100
static char **cell_line;
static char **cell_name;
static char **rdname;
static int *hit_score,*hit_length,*hit_refst,*hit_refed,*hit_rcdex,*hit_quest,*hit_queed,*hit_identy;
static int *ctg_index,*cell_index;
static int *contig_hist,*contig_body,*contig_list,*cigar_list;
static long *contig_head,*cigar_head,sBase;
static char *dataline,*cigar_line;

/* SSAS default parameters   */
static int IMOD=0;
static int ISUB=0;
static int mapNumber=1;
static int copyNumber=20;
static long line_len=0;
static int num_reads=0;
static int num_cline=0;
static int low_qual=1;
static int debug_flag=0;
static int set_identy = 0;
static int set_cover = 100;
static int set_match = 80;
static int pos_input = 200;
static int cons_flag = 0;
static int ctg_input = 0;
static int view_mod = 1;
static int n_contigs = 0;
static int min_length = 20;
static int cell_flag  = 0;
static long genome_offset;
static char strain_name[100];

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *sub;

static char rc_char[500000];

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
    FILE *namef;
    int i,nSeq,args;
    char line[2000]={0};
    char **cmatrix(long nrl,long nrh,long ncl,long nch);
    fasta *seq; 
    void ArraySort_Mix(int n, long int *arr, int *brr);
    void SNP_Consensus(char **argv,int args,int nSeq);
    void Read_Pairs(char **argv,int args,int nLib,int nSeq);
    void Align_Process(char **argv,int args,int nRead);

    seq=NULL;
//    fflush(stdout);
//    system("ps aux | grep eulerSNP; date");
    memset(strain_name,'\0',100);
    if(argc < 2)
    {
      printf("Usage: %s <reads_alignment_file> <cigar_line_file> <ref_sequence> <reads_fastq_file\n",argv[0]);
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
       else if(!strcmp(argv[i],"-copy"))
       {
         sscanf(argv[++i],"%d",&copyNumber);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-map"))
       {
         sscanf(argv[++i],"%d",&mapNumber);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-len"))
       {
         sscanf(argv[++i],"%d",&min_length);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-cover"))
       {
         sscanf(argv[++i],"%d",&set_cover);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-strain"))
       {
         sscanf(argv[++i],"%s",strain_name);
         cell_flag  = 1;
         args=args+2;
       }
       else if(!strcmp(argv[i],"-view"))
       {
         sscanf(argv[++i],"%d",&view_mod);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-pos"))
       {
         sscanf(argv[++i],"%d",&pos_input);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-contig"))
       {
         sscanf(argv[++i],"%d",&ctg_input);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-debug"))
       {
         sscanf(argv[++i],"%d",&debug_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-idt"))
       {
         sscanf(argv[++i],"%d",&set_identy);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-qual"))
       {
         sscanf(argv[++i],"%d",&low_qual);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-match"))
       {
         sscanf(argv[++i],"%d",&set_match);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-cons"))
       {
         sscanf(argv[++i],"%d",&cons_flag);
         args=args+2;
       }
    }


/*  input read alignment info line   */
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args]);
      exit(1);
    }

    num_reads = 0;
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      num_reads++;
    }
    fclose(namef); 

    Align_Process(argv,args,num_reads);
    SNP_Consensus(argv,args,num_reads);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   Subroutine to process alignment information */
/* ====================================================  */
void Align_Process(char **argv,int args,int nRead)
/* ====================================================  */
{
     int i,j,k,n_reads = nRead;
     void ArraySort_Mix(int n,long int *arr,int *brr);
     int  **imatrix(long nrl,long nrh,long ncl,long nch);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
//     char **rdname;
     char **DBname,*ptr,line[500];
     FILE *namef;
     int n_find,*readIndex,idd,stopflag,num_align,refhit1,refhit2;
     int hit_len = 0,readhit1 = 0,readhit2 = 0;
 
     sub = decodeFastq(argv[args+2],&n_contigs,&sBase,0);
     fastaUC(sub,n_contigs);

     n_reads = nRead + n_contigs;
     cell_name = cmatrix(0,100,0,Max_N_NameBase);
     cell_line = cmatrix(0,nRead,0,Max_N_NameBase);
     rdname=cmatrix(0,nRead,0,Max_N_NameBase);
     DBname=cmatrix(0,n_reads,0,Max_N_NameBase);

     if((readIndex= (int *)calloc(n_reads,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - readIndex\n");
       exit(1);
     }
     if((hit_score= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - nRead\n");
       exit(1);
     }
     if((hit_quest= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_quest\n");
       exit(1);
     }
     if((hit_queed= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_queed\n");
       exit(1);
     }
     if((hit_refst= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_refst\n");
       exit(1);
     }
     if((hit_rcdex= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_rcdex\n");
       exit(1);
     }
     if((hit_refed= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_refed\n");
       exit(1);
     }
     if((hit_length= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_length\n");
       exit(1);
     }
     if((hit_identy= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_identy\n");
       exit(1);
     }
     if((ctg_index= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - ctg_index\n");
       exit(1);
     }
     if((cell_index= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - cell_index\n");
       exit(1);
     }
     if((cigar_list= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - cigar_list\n");
       exit(1);
     }
     if((cigar_head= (long *)calloc(nRead,sizeof(long))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - cigar_head\n");
       exit(1);
     }

/*   read the cigar line file   */
     if((namef = fopen(argv[args+1],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

/*   get the maximum cigar line length         */
     line_len=0;
     i = 0;
     while(!feof(namef))
     {
       fgets(line,2000,namef);
       if(feof(namef)) break;
       cigar_list[i] = strlen(line) - 1;
       line_len = line_len + cigar_list[i];
       cigar_head[i] = line_len - cigar_list[i];
       i++;
     }
     fclose(namef);
     if((cigar_line= (char *)calloc((line_len+5000),sizeof(char))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - cigar_line\n");
       exit(1);
     }

/*   read the cigar line file   */
     if((namef = fopen(argv[args+1],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }

/*   get the maximum cigar line length         */
     line_len=0;
     i = 0;
     while(!feof(namef))
     {
       int c_len;
       fgets(line,2000,namef);
       if(feof(namef)) break;
       c_len = strlen(line) - 1;
       for(j=0;j<c_len;j++)
       {
          cigar_line[cigar_head[i]+j] = line[j];
       }
       i++;
     }
     fclose(namef);

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

/*   read the SNP output file         */
     num_align=0;
     while(!feof(namef))
     {
       int nPair=0,len;
       char line2[500],base[500];
      
       fgets(line,500,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       if((strncmp(line,"ALIGNM",6))==0)
       { 
         refhit1 = 0;
         refhit2 = 0;     
         readhit1 = 0;
         readhit2 = 0;     
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
         i=0;
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==0)
            {
            }
            else if(i==1)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              hit_score[num_align]=atoi(ptr);
            }
            else if(i==2)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              strcpy(rdname[num_align],ptr);
            }
            else if(i==3)
             {
              memset(base,'\0',500);
              strcat(base,ptr);
              strcpy(DBname[num_align],ptr);
            }
            else if(i==4)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              readhit1 = atoi(ptr);
            }
            else if(i==5)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              readhit2 = atoi(ptr);
            }
            else if(i==6)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              refhit1 = atoi(ptr);
            }
            else if(i==7)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              refhit2 = atoi(ptr);
            }
            else if(i==8)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              if((*ptr)=='F')
              {
                hit_rcdex[num_align]=0;
                hit_refst[num_align]=refhit1;
                hit_refed[num_align]=refhit2;
                hit_quest[num_align]=readhit1;
                hit_queed[num_align]=readhit2;
              } 
              else
              {
                if((readhit1-readhit2)<0)
                {
                  hit_refst[num_align]=refhit2;
                  hit_refed[num_align]=refhit1;
                  hit_quest[num_align]=readhit1;
                  hit_queed[num_align]=readhit2;
                } 
                else
                {
                  hit_refst[num_align]=refhit1;
                  hit_refed[num_align]=refhit2;
                  hit_quest[num_align]=readhit2;
                  hit_queed[num_align]=readhit1;
                } 
                hit_rcdex[num_align]=1;
              }
            }
            else if(i==9)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              hit_len = atoi(ptr); 
            }
            else if(i==10)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              hit_identy[num_align] = atoi(ptr); 
            }
            else if(i==11)
            {
              int hit_len2;
              memset(base,'\0',500);
              strcat(base,ptr);
              hit_len = hit_len*100;
              hit_len2 =  atoi(ptr); 
              hit_length[num_align] = hit_len/hit_len2; 
            }
            else if(i==12)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              len = strlen(ptr);
              strncpy(cell_line[num_align],ptr,len-1);
              readIndex[num_align] = num_align;
              num_align++;
            }
         }
       }
     }
     fclose(namef);

     memset(ctg_index,-1,num_align*4);
/*   sort out reference head  */
     if(cell_flag == 0)
     {
       ArraySort_String(num_align,cell_line,readIndex);
       
       num_cline = 0;
       idd = 0;
       for(i=0;i<num_align;i++)
       {
          stopflag=0;
          j=i+1;
          while((j<num_align)&&(stopflag==0))
          {
            if(strcmp(cell_line[j],cell_line[i])==0)
            {
              j++;
            }
            else
              stopflag=1;
          }
          for(k=i;k<j;k++)
          {
             idd = readIndex[k];
             cell_index[idd] = num_cline;
          }
          if(num_cline<50)
            strcpy(cell_name[num_cline],cell_line[i]);
          num_cline++;
          i=j-1;
       }
     }
     else
     {
       num_cline = 1;
       strcpy(cell_name[0],strain_name);
     }

     printf("number of cell lines: %d\n",num_cline);
     free(cell_line);

/*   sort out contig name match */
     for(j=0;j<nRead;j++)
     {
        readIndex[j]=j;
     }

     for(j=0;j<n_contigs;j++)
     {
        strcpy(DBname[j+nRead],(sub+j)->name);
        readIndex[j+nRead]=j+nRead;
     } 
     n_reads = nRead + n_contigs; 
     ArraySort_String(n_reads,DBname,readIndex);

     n_find = 0;
     idd = -1;
     for(i=0;i<n_reads;i++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
//     printf("www2: %d %d %s %s\n",i,j,DBname[i],DBname[j]);
        while((j<n_reads)&&(stopflag==0))
        {
          if(strcmp(DBname[j],DBname[i])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        idd = -1;
        if((j-i)>1)
        {
          for(k=i;k<j;k++)
          {
             if(readIndex[k]>=nRead)
             {
               idd = readIndex[k]-nRead;
               k = j;
             }
          }
          if(idd>=0)
          {
            for(k=i;k<j;k++)
            {
               if(readIndex[k]<nRead)
               {
//     printf("xxx: %d %d %s %d %d %d\n",i,k,DBname[i],idd,nSeq,nRead);
                 ctg_index[readIndex[k]] = idd;
                 n_find++;
               }
            }
          }
        }
        i=j-1;
     }
     printf("number of reads aligned to the genome: %d\n",n_find);

}


/*   Subroutine to find the contig/supercontig structure */
/* ====================================================  */
void SNP_Consensus(char **argv,int args,int nSeq)
/* ====================================================  */
{
     int i,j,k,m,c,b;
     long *ref_head,idk=0;
     int *ref_list;
     void ArraySort_Mix(int n,long int *arr,int *brr);
     void ArraySort2_Int2(int n, int *arr, int *brr);
     long offset,st,ed,kk,tBase,head_size;
     int  **imatrix(long nrl,long nrh,long ncl,long nch);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
     int set_score = 0,err_flag = 0;
     fasta *seq,*seqp;
     int num_traces,read_qual[5000],read_qsum[4],read_nhit[4],**m_align;
     char **q_align,**s_align,refe_base[5000],read_base[5000],SNP_base,SNP2_base;
     int num_hits,sum_Q_set;
     int read_qindex[4],read_qsort[4],*read_offset;
     int *cell_cover,*base_cover;
     int **cell_nhit; 
     void Output_Process(fasta *seq,int iseq,int nRead,int *read_offset,char *read_base);

     m_align=imatrix(0,500,0,200);
     cell_nhit=imatrix(0,4,0,num_cline);
     seq = decodeFastq(argv[args+3],&num_traces,&tBase,0);
     fastaUC(seq,num_traces);

     if((cell_cover= (int *)calloc(num_cline,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - cell_cover\n");
       exit(1);
     }
     if((base_cover= (int *)calloc(num_cline,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - base_cover\n");
       exit(1);
     }
     if((contig_hist= (int *)calloc(sBase,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - contig_hist\n");
       exit(1);
     }
     if((contig_list= (int *)calloc(sBase+1,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - contig_list\n");
       exit(1);
     }
     if((contig_head= (long *)calloc(sBase+1,sizeof(long))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - contig_head\n");
       exit(1);
     }
     if((ref_list= (int *)calloc(n_contigs,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ref_list\n");
       exit(1);
     }
     if((ref_head= (long *)calloc(n_contigs,sizeof(long))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ref_head\n");
       exit(1);
     }
     if((read_offset= (int *)calloc(5000,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - read_offset\n");
       exit(1);
     }

     seqp = sub;
     ref_list[0] = seqp->length;
     ref_head[0] = 0L;
     for(i=1;i<n_contigs;i++)
     {
        seqp = sub + i;
        ref_list[i] = seqp->length;
        ref_head[i] = ref_head[i-1] + ref_list[i-1];
     }

     for(i=0;i<nSeq;i++)
     {
        int idt;
        int match_flag = 0, q_flag=1;
        if(abs(hit_refed[i]-hit_refst[i])>=min_length)
          match_flag = 1;
        set_score = 0;
	seqp = seq + i;
	if(hit_quest[i]>=5)
	{
	  int qedge1 = 0;
	  for(j=(hit_quest[i]-5);j<(hit_quest[i]-1);j++)
	  {
	     if(seqp->qual[j]>=30)
	       qedge1++;
	  }
	  if(qedge1==4)
	    q_flag=0;
	}
	if((seqp->length-hit_queed[i])>=5)
	{
	  int qedge2 = 0;
	  for(j=hit_queed[i];j<(hit_queed[i]+4);j++)
	  {
	     if(seqp->qual[j]>=30)
	       qedge2++;
	  }
	  if(qedge2==4)
	    q_flag = 0;
	}
//	if(q_flag == 0)
//     printf("here0: %s %d %d\n",rdname[i],hit_quest[i],hit_queed[i]);
        if((q_flag)&&(match_flag)&&(ctg_index[i]>=0))
        {
          offset = ref_head[ctg_index[i]];
          st = hit_refst[i];
          ed = hit_refed[i];
          idt = cell_index[i];
          for(kk=st;kk<=ed;kk++)
          {
             contig_list[offset+kk]++;
          }
        }
     }

     contig_head[0] = 0;
     for(i=1;i<=sBase;i++)
     {
        contig_head[i] = contig_head[i-1] + contig_list[i-1];
        contig_list[i-1] = 0;
     }
     head_size = contig_head[sBase-1]+contig_list[sBase-1]+100;
     contig_list[sBase-1] = 0;

     if((contig_body= (int *)calloc(head_size,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - contig_body\n");
       exit(1);
     }

//     printf("here1: %ld\n",sBase);
     for(i=0;i<nSeq;i++)
     {
        int idt;
        int match_flag = 0,q_flag = 1;
        if(abs(hit_refed[i]-hit_refst[i])>=min_length)
          match_flag = 1;
        set_score = 0;
	seqp = seq + i;
	if(hit_quest[i]>=5)
	{
	  int qedge1 = 0;
	  for(j=(hit_quest[i]-5);j<(hit_quest[i]-1);j++)
	  {
	     if(seqp->qual[j]>=30)
	       qedge1++;
	  }
	  if(qedge1==4)
	    q_flag = 0;
	}
	if((seqp->length-hit_queed[i])>=5)
	{
	  int qedge2 = 0;
	  for(j=hit_queed[i];j<(hit_queed[i]+4);j++)
	  {
	     if(seqp->qual[j]>=30)
	       qedge2++;
	  }
	  if(qedge2==4)
	    q_flag = 0;
	}
        if((q_flag)&&(match_flag)&&(ctg_index[i]>=0))
        {
          offset = ref_head[ctg_index[i]];
          st = hit_refst[i];
          ed = hit_refed[i];
          idt = cell_index[i];
          for(kk=st;kk<=ed;kk++)
          {
             contig_body[contig_head[offset+kk]+contig_list[offset+kk]] = i;
             contig_list[offset+kk]++;
          }
        }
     }

//     pos_input = 13227;

     if((dataline = (char *)calloc(1000,sizeof(char))) == NULL)
     {
       printf("ERROR ssaha: calloc - dataline\n");
       exit(1);
     }

     genome_offset = ref_head[ctg_input] + pos_input; 
     j = contig_list[genome_offset];
     num_hits = j;
     q_align = cmatrix(0,j,0,200);
     s_align = cmatrix(0,j,0,200);
     m_align = imatrix(0,j,0,200);
     for(i=0;i<j;i++)
     {
        memset(q_align[i],'\0',200);
        memset(s_align[i],'\0',200);
        memset(m_align[i],0,200*4);
     }

     err_flag = 0;
     sum_Q_set = 0;
     memset(refe_base,'\0',5000);
     memset(read_qual,0,4*5000);
     for(c=0;c<n_contigs;c++)
     {
        for(b=1;b<(sub+c)->length;b++)
        {
           int max_qsum,max_base,rqual,max_qsum2,out_flag;
           int map_base,max_base2;
           max_base2 = -1;
           genome_offset = ref_head[c] + b;
           pos_input = b;
           ctg_input = c;
//           printf("head: %s %d %d ",(sub+c)->name,b,contig_list[genome_offset]); 
           num_hits = contig_list[genome_offset];
           if((num_hits<=0)||(num_hits>200))
             continue;
           if(num_hits<10)
             sum_Q_set = 30;
           else
             sum_Q_set = 50;
           memset(read_base,'\0',5000);
           memset(cell_cover,0,4*num_cline);
           memset(cell_nhit[0],0,4*num_cline);
           memset(cell_nhit[1],0,4*num_cline);
           memset(cell_nhit[2],0,4*num_cline);
           memset(cell_nhit[3],0,4*num_cline);
           for(j=0;j<contig_list[genome_offset];j++)
           {
              int qst,sst,mpos,m_bases,r_bases,s_bases;
              int nPair = 0,r_mod,r_len,idt;
              char *ptr,line2[2000],line[2000];

              read_offset[j] = 0;
//              memset(q_align[j],'\0',200);
              idk = contig_head[genome_offset]+j;
              i = contig_body[idk];
              idt = cell_index[i];
              seqp = seq + i;
//           printf("name: %s %s %ld\n",seqp->name,cell_name[idt],idk); 
              seqp->qual[0] = 40; 
              memset(line,'\0',2000);
              for(k=0;k<cigar_list[i];k++)
                 line[k] = cigar_line[cigar_head[i]+k];
              strcpy(line2,line);
              qst = hit_quest[i];
              sst = hit_refst[i];
              if(hit_rcdex[i]==0)
                r_bases = qst-1;
              else
                r_bases = hit_queed[i]-1;
              s_bases = sst-1;
              m_bases = 0;
              mpos = 0;
              r_mod = 0;
              for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
              {
              }
              m=0;
              for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),m++)
              {
                 if((m<nPair)&&(m>9))
                 {
                   if((m%2)==0)
                   {
                     if(*ptr=='M')
                       r_mod = 1;
                     else if(*ptr=='I')
                       r_mod = 2;
                     else if(*ptr=='D')
                       r_mod = 3;
                     else
                       r_mod = 0;
                   }
                   else
                   {
                     r_len = atoi(ptr);
                     if(hit_rcdex[i]==0)
                     {
                       if(r_mod==1)
                       {
                         for(k=0;k<r_len;k++)
                         {
                            if((s_bases+k)==(pos_input-1))
                            {
                              read_base[j] = seqp->data[k+r_bases]; 
                              read_qual[j] = seqp->qual[k+r_bases];
                              refe_base[j] = (sub+c)->data[k+s_bases];
			      read_offset[j] = k+r_bases;
                              mpos = m_bases + k;
                              k = r_len;
                              m = nPair; 
                            }
                         }
                         r_bases = r_bases + r_len;
                         s_bases = s_bases + r_len;
                         m_bases = m_bases + r_len;
                       }
                       else if(r_mod==2)
                       {
                         r_bases = r_bases + r_len;
                         m_bases = m_bases + r_len;
                       }
                       else if(r_mod==3)
                       {
                         for(k=0;k<r_len;k++)
                         {
                            if((s_bases+k)==(pos_input-1))
                            {
                              read_base[j] = '-'; 
                              read_qual[j] = 0;
                              refe_base[j] = (sub+ctg_input)->data[k+s_bases];
                              mpos = m_bases + k; 
                              k = r_len;
                              m = nPair; 
                            }
                         }
                         s_bases = s_bases + r_len;
                         m_bases = m_bases + r_len;
                       }
                     }
                     else
                     {
                       if(r_mod==1)
                       {
                         for(k=0;k<r_len;k++)
                         {
                            if((s_bases+k)==(pos_input-1))
                            {
                              if(seqp->data[r_bases-k]=='A')
                                read_base[j]='T';
                              else if(seqp->data[r_bases-k]=='C')
                                read_base[j]='G';
                              else if(seqp->data[r_bases-k]=='G')
                                read_base[j]='C';
                              else if(seqp->data[r_bases-k]=='T')
                                read_base[j]='A';
                              else
                                read_base[j]=seqp->data[r_bases-k];
                              read_qual[j] = seqp->qual[r_bases-k];
                              refe_base[j] = (sub+ctg_input)->data[k+s_bases];
			      read_offset[j] = r_bases-k;
                              mpos = m_bases + k;
                              k = r_len;
                              m = nPair; 
                            }
                         }
                         r_bases = r_bases - r_len;
                         s_bases = s_bases + r_len;
                         m_bases = m_bases + r_len;
                       }
                       else if(r_mod==2)
                       {
                         r_bases = r_bases - r_len;
                         m_bases = m_bases + r_len;
                       }
                       else if(r_mod==3)
                       {
                         for(k=0;k<r_len;k++)
                         {
                            if((s_bases+k)==(pos_input-1))
                            {
                              read_base[j] = '-'; 
                              read_qual[j] = 0;
                              refe_base[j] = (sub+c)->data[k+s_bases];
                              mpos = m_bases + k; 
                              k = r_len;
                              m = nPair; 
                            }
                         }
                         s_bases = s_bases + r_len;
                         m_bases = m_bases + r_len;
                       }
                     }
                   }
                 }
              }
              if(j==0)
              {
//                printf("%c ",refe_base[j]);
                memset(read_qsum,0,4*4);
                memset(read_nhit,0,4*4);
              }
              if(read_base[j]!='-')
                cell_cover[idt]++;
              if(read_qual[j]>=23)
                rqual = 18;
              else
                rqual = 6;
              if(read_base[j]=='A')
              {
                read_qsum[0] = read_qsum[0] + rqual;
                cell_nhit[0][idt]++; 
                read_nhit[0]++;
              }
              else if(read_base[j]=='C')
              {
                read_qsum[1] = read_qsum[1] + rqual;
                cell_nhit[1][idt]++; 
                read_nhit[1]++;
              }
              else if(read_base[j]=='G')
              {
                read_qsum[2] = read_qsum[2] + rqual;
                cell_nhit[2][idt]++; 
                read_nhit[2]++;
              }
              else if(read_base[j]=='T')
              {
                read_qsum[3] = read_qsum[3] + rqual;
                cell_nhit[3][idt]++; 
                read_nhit[3]++;
              }

//              printf("%c",read_base[j]);

              if(read_qual[j]<23) 
                 read_base[j] = tolower(read_base[j]);
           }
//           printf("\n");
           max_qsum = 0;
           max_base = 0; 
           max_base2 = -1; 
           for(j=0;j<4;j++)
           {
              read_qindex[j] = j;
              read_qsort[j]  = read_qsum[j];
              if(read_qsum[j]>max_qsum)
              {
                max_base = j;
                max_qsum = read_qsum[j];
              }
           }

           ArraySort2_Int2(4,read_qsort,read_qindex);
           max_qsum2 = read_qsort[1];
           if(max_qsum2 <= 0)
             max_base2 = -1;
           else
             max_base2 = read_qindex[1];
           SNP_base = '-';
           SNP2_base = '-';
           if(read_qindex[1]==0)
             SNP2_base = 'A';
           else if(read_qindex[1]==1) 
             SNP2_base = 'C';
           else if(read_qindex[1]==2) 
             SNP2_base = 'G';
           else if(read_qindex[1]==3) 
             SNP2_base = 'T';
           
           if(max_base==0)
             SNP_base = 'A';
           else if(max_base==1)
             SNP_base = 'C';
           else if(max_base==2)
             SNP_base = 'G';
           else if(max_base==3)
             SNP_base = 'T';
           map_base = 0;
           if(refe_base[0]=='A')
             map_base = 0;
           else if(refe_base[0]=='C')
             map_base = 1;
           else if(refe_base[0]=='G')
             map_base = 2;
           else if(refe_base[0]=='T')
             map_base = 3;

           out_flag = 0;
//           printf("cons: %s %d %d %c %c %s\n",(sub+c)->name,b,contig_list[genome_offset],refe_base[0],SNP_base,read_base); 
           if(max_base!=map_base)
           {
             if((read_nhit[max_base]>=1)&&(max_qsum>=sum_Q_set))
             {
               if(max_qsum2>0)
                 printf("SNP1: %s %d %d %c %c %s %d %d\n",(sub+c)->name,b,contig_list[genome_offset],refe_base[0],SNP_base,read_base,read_nhit[max_base2],read_nhit[max_base]); 
               else
                 printf("SNP1: %s %d %d %c %c %s %d %d\n",(sub+c)->name,b,contig_list[genome_offset],refe_base[0],SNP_base,read_base,0,read_nhit[max_base]); 
             }
           }
           else if(max_base2>=0)
           {
             if(max_qsum2>sum_Q_set)
             {
               int max_second = 0, max_cline = 0;
               float rate = 0.0;
               rate = max_qsum2;
               rate = rate/max_qsum;
               for(k=0;k<num_cline;k++)
               {
                  if(cell_nhit[max_base2][k]>max_second)
                  {
                    max_second = cell_nhit[max_base2][k];
                    max_cline  = k;
                  }
               }
               if(((max_second==1)&&(cell_cover[max_cline]<=2))||(max_second>1))
               {
                 if((max_base==map_base)&&(read_nhit[max_base]<200))
                 {
                   if(((read_nhit[max_base]>20)&&(rate>0.1))||(read_nhit[max_base]<=20))
                   {
		     rate = read_nhit[max_base2];
		     rate = rate/contig_list[genome_offset];
		     if(rate>0.25)
		     {
		       int n_hitsq = 0;
		       int n_hitrq = 0;
		       float rate2 = 0.0;

		       for(k=0;k<contig_list[genome_offset];k++)
		       {
		          if(read_base[k] == refe_base[0])
			    n_hitsq++;
			  if(read_base[k] == SNP2_base)
			    n_hitrq++;
		       }
		       rate2 = n_hitrq;
		       rate2 = rate2/n_hitsq;
		       if(rate2 > 0.25)
		       {
                         printf("SNP2: %s %d %d %c %c %s %d %d %d %d\n",(sub+c)->name,b,contig_list[genome_offset],refe_base[0],SNP2_base,read_base,read_nhit[max_base],read_nhit[max_base2],n_hitsq,n_hitrq);
                         if(debug_flag)
                           Output_Process(seq,genome_offset,contig_list[genome_offset],read_offset,read_base);
		       }
		     }
                   }
                 }
                 else
                 {
                   if(max_base2==map_base)
                     printf("SNP3: %s %d %d %c %c %s %d %d\n",(sub+c)->name,b,contig_list[genome_offset],refe_base[0],SNP_base,read_base,read_nhit[max_base2],read_nhit[max_base]); 
                   else
                     printf("SNP4: %s %d %d %c %c/%c %s %d %d/%d\n",(sub+c)->name,b,contig_list[genome_offset],refe_base[0],SNP_base,SNP2_base,read_base,read_nhit[map_base],read_nhit[max_base],read_nhit[max_base2]); 
                 }
               }
             }
           }
        }
     }

}


/*   Subroutine to process alignment information */
/* ====================================================  */
void Output_Process(fasta *seq,int i_offset,int nRead,int *read_offset,char *read_base)
/* ====================================================  */
{      
     int i,j,jj,k,idk,len;
     fasta *seqp;
     void ArraySort_Int2(int n, int *arr, int *brr);
     int read_index[nRead],read_locus[nRead];

     for(j=0;j<nRead;j++)
     {
        idk = contig_head[i_offset]+j;
        i = contig_body[idk];
        seqp = seq + i;
	len = seqp->length;
        if(hit_rcdex[i]==0)
          read_locus[j] = len-read_offset[j];
	else
	  read_locus[j] = read_offset[j]+1;
	read_index[j] = j;
     }
     ArraySort_Int2(nRead,read_locus,read_index);
     for(jj=0;jj<nRead;jj++)
     {
        j = read_index[jj];
        idk = contig_head[i_offset]+j;
        i = contig_body[idk];
        seqp = seq + i;
	len = seqp->length;
        printf("%-25s %c %2d %2d %2d %2d ",seqp->name,read_base[j],hit_score[i],read_offset[j],hit_refed[i]-hit_refst[i],hit_rcdex[i]);
	if(hit_rcdex[i]==0)
	{
  	  for(k=0;k<(len-read_offset[j]);k++)
             printf("%c",' ');
  	  for(k=0;k<len;k++)
	  {
	     if(seqp->qual[k]>=23)
               printf("%c",seqp->data[k]);
  	     else
	       printf("%c",tolower(seqp->data[k]));
	  }
          printf("\n");
	}
	else
	{
  	  for(k=0;k<(read_offset[j]+1);k++)
             printf("%c",' ');
  	  for(k=(len-1);k>0;k--)
	  {
	     if(seqp->qual[k]>=23)
	     {
	       if(seqp->data[k]=='A')
	         printf("%c",'T');
	       else if(seqp->data[k]=='C')
	         printf("%c",'G');
	       else if(seqp->data[k]=='G')
	         printf("%c",'C');
	       else if(seqp->data[k]=='T')
	         printf("%c",'A');
	       else
	         printf("%c",seqp->data[k]);
	     }
  	     else
	     {
	       if(seqp->data[k]=='A')
	         printf("%c",'t');
	       else if(seqp->data[k]=='C')
	         printf("%c",'g');
	       else if(seqp->data[k]=='G')
	         printf("%c",'c');
	       else if(seqp->data[k]=='T')
	         printf("%c",'a');
	       else
	         printf("%c",seqp->data[k]);
	     }
	  }
          printf("\n");
	}
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
void ArraySort_Mix(int n, long int *arr, int *brr)
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

