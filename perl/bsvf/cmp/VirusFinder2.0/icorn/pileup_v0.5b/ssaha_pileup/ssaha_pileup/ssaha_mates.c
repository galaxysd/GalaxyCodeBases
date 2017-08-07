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
static char **cell_name;
static char **rdname;
static int *map_score,*hit_refst,*hit_refed,*hit_quest,*hit_queed;
static int *readIndex,*map_unique,*map_rdpair,*ctg_index,*mate_pairs,*mate_index,*reads2wgs;
static int *mate_insert,*mate_std,*hit_rcdex,*hit_frdex;

/* SSAS default parameters   */
static int IMOD=0;
static int num_reads=0;
static int num_mates=0;
static int insert_size = 800;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *sub;

static char rc_char[500];

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
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    fasta *seq; 
    void ArraySort_Mix(int n, B64_long *arr, int *brr);
    void SNP_Consensus(char **argv,int args,int nSeq);
    void Read_Pairs(char **argv,int args,int nLib,int nSeq);
    void Align_Process(char **argv,int args,int nRead);

    seq=NULL;
//    fflush(stdout);
//    system("ps aux | grep eulerSNP; date");
    if(argc < 2)
    {
      printf("Usage: %s -insert 800 <cigar_line_file> <mates_file> <output_cigar_line> \n",argv[0]);
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
       else if(!strcmp(argv[i],"-insert"))
       {
         sscanf(argv[++i],"%d",&insert_size);
         args=args+2;
       }
    }


/*  input cigar info line   */
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

/*  input mates info line   */
    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args]);
      exit(1);
    }

    num_mates = 0;
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      num_mates++;
    }
    fclose(namef); 

    printf("number of reads aligned on to genome: %d\n",num_reads);
    Align_Process(argv,args,num_reads);
//    SNP_Consensus(argv,args,num_reads);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   Subroutine to process alignment information */
/* ====================================================  */
void Align_Process(char **argv,int args,int nRead)
/* ====================================================  */
{
     int i,j,k,n_reads = nRead;
     void ArraySort_Mix(int n,B64_long *arr,int *brr);
     int  **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
//     char **rdname;
     char **ctgname,*ptr,line[2000],RC;
     void ArraySort_Int2(int n, int *arr, int *brr);
     void Mates_Process(char **argv,int args,int nRead);
     FILE *namef,*namef2;
     int *rdpair_index,read_index[2000];
     int rd_forward[2000],rd_reverse[2000];
     B64_long read_offsets[2000];
     int n_find,idd,stopflag,num_align,refhit1,refhit2,num_pairs;
     int readhit1 = 0,readhit2 = 0;
     int insertSize1 = insert_size*(1.0+0.2);
     int insertSize2 = insert_size*(1.0-0.2);

     n_reads = nRead;
     RC = '+';
     cell_name = cmatrix(0,2,0,Max_N_NameBase);
     rdname=cmatrix(0,nRead,0,Max_N_NameBase);
     ctgname=cmatrix(0,nRead,0,Max_N_NameBase);

     if((readIndex= (int *)calloc(n_reads,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - readIndex\n");
       exit(1);
     }
     if((rdpair_index= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - rdpair_index\n");
       exit(1);
     }
     if((map_score= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - map_score\n");
       exit(1);
     }
     if((hit_rcdex= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_rcdex\n");
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
     if((hit_refed= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_refed\n");
       exit(1);
     }
     if((map_unique= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - map_unique\n");
       exit(1);
     }
     if((map_rdpair= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - map_rdpair\n");
       exit(1);
     }
     if((ctg_index= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ctg_index\n");
       exit(1);
     }

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

/*   read the cigar line file         */
     num_align=0;
     while(!feof(namef))
     {
       int nPair=0,len;
       char line2[2000],line3[2000],base[500],score[3];
      
       fgets(line,2000,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       strcpy(line3,line);
       if((strncmp(line,"cigar",5))==0)
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
              memset(score,'\0',3);
              memset(base,'\0',500);
              strcat(base,ptr);
	      len = strlen(base);
	      if(len < 8)
	      {
	        map_score[num_align] = 50;
//	        printf("score: %d %s",len,line3);
	      }
	      else
	      {
	        score[0] = line3[7];
	        score[1] = line3[8];
	        map_score[num_align] = atoi(score);
		if(map_score[num_align]>=30)
	          map_unique[num_align] = 1;
//	        printf("score: %d %s",len,line3);
	      }
            }
            else if(i==1)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              strcpy(rdname[num_align],ptr);
            }
            else if(i==2)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              readhit1 = atoi(ptr);
            }
            else if(i==3)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              readhit2 = atoi(ptr);
            }
            else if(i==4)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
	      RC = *ptr;
	    }
            else if(i==5)
            {
              memset(base,'\0',500);
              strcpy(ctgname[num_align],ptr);
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
              if(RC=='+')
              {
                hit_rcdex[num_align]=0;
                hit_refst[num_align]=refhit1;
                hit_refed[num_align]=refhit2;
                hit_quest[num_align]=readhit1;
                hit_queed[num_align]=readhit2;
              } 
              else
              {
                hit_refst[num_align]=refhit2;
                hit_refed[num_align]=refhit1;
                hit_quest[num_align]=readhit2;
                hit_queed[num_align]=readhit1;
                hit_rcdex[num_align]=1;
              }
              readIndex[num_align] = num_align;
              num_align++;
            }
         }
       }
     }
     fclose(namef);

     Mates_Process(argv,args,num_align);
/*   sort out contig/chromosome idnex */
     n_reads = nRead; 
     ArraySort_String(n_reads,ctgname,readIndex);

     n_find = 0;
     idd = 0;
     for(i=0;i<n_reads;i++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        while((j<n_reads)&&(stopflag==0))
        {
          if(strcmp(ctgname[j],ctgname[i])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>=1)
        {
          for(k=i;k<j;k++)
             ctg_index[readIndex[k]] = idd;
        }
	idd++;
        i=j-1;
     }

/*   sort out read pairs */
     n_reads = nRead;
     num_pairs = 0;
     for(i=0;i<n_reads;i++)
        rdpair_index[i] = i+num_mates;
     for(i=0;i<n_reads;i++)
     {
	int rd = reads2wgs[i];
        readIndex[i] = i;
	rdpair_index[i] = mate_index[rd];
     }

     printf("Started the process: %d\n",nRead);
     ArraySort_Int2(n_reads,rdpair_index,readIndex);

     n_find = 0;
     idd = -1;
     for(i=0;i<n_reads;i++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        while((j<n_reads)&&(stopflag==0))
        {
          if(rdpair_index[j]==rdpair_index[i])
          {
            j++;
          }
          else
            stopflag=1;
        }
        idd = -1;
        if((j-i)>=2)
	{
	  int num_hits = j-i;
	  int stopflag2,m,n_pair;
	  int pair1,pair2;

          n_pair = 0;
          pair1 = 0;
          pair2 = 0;
//     printf("www: %d %d\n",i,j);
	  for(k=0;k<num_hits;k++)
	  {
	     int idt = readIndex[i+k];
	     B64_long nn = ctg_index[idt];
	     read_index[k] = k;
	     read_offsets[k] = (nn<<34) + hit_refst[idt];
	  }
          if(num_hits==2)
          {
            int ctg1 = read_offsets[0]>>34;
            int ctg2 = read_offsets[1]>>34;
            int idt = readIndex[i+read_index[0]];
            idd = readIndex[i+read_index[1]];
            if((ctg1==ctg2)&&(strcmp(rdname[idt],rdname[idd])!=0)&&(hit_rcdex[idt]!=hit_rcdex[idd]))
            {
              int offset1 = hit_refst[idt]; 
              int offset2 = hit_refst[idd];
              int gap_size = abs(offset2-offset1);
              int size_flag = 0;
            
              if(mate_std[reads2wgs[idt]] > 0)
              {
                insertSize1 = mate_insert[reads2wgs[idt]] + mate_std[reads2wgs[idt]];
                insertSize2 = mate_insert[reads2wgs[idt]] - mate_std[reads2wgs[idt]];
              }
              else
              {
                insertSize1 = 1.2*insert_size;
                insertSize2 = 0.8*insert_size;
              }
              if((gap_size <= insertSize1)&&(gap_size>=insertSize2))
                size_flag = 1;
              if(size_flag)
              {
                pair1 = idt;
                pair2 = idd;
                n_pair++;
              }
            }
          }
          else if(num_hits>2)
          {
               int num_forward = 0;
               int num_reverse  = 0;
               int ii,jj;
               for(ii=0;ii<num_hits;ii++)
               {
                  idd = readIndex[i+ii];
                  if(hit_frdex[reads2wgs[idd]] == 1)
                  {
                    rd_forward[num_forward] = idd;
                    num_forward++;
                  }
                  else if(hit_frdex[reads2wgs[idd]] == 2)
                  {
                    rd_reverse[num_reverse] = idd;
                    num_reverse++;
                  }
               }
               for(jj=0;jj<num_reverse;jj++)
               {
                  idd = rd_reverse[jj];
                  for(ii=0;ii<num_forward;ii++)
                  {
                     int idt = rd_forward[ii]; 
                     if((strcmp(rdname[idt],rdname[idd])!=0)&&(hit_rcdex[idt]!=hit_rcdex[idd]))
                     {
                       int offset1 = hit_refst[idt]; 
                       int offset2 = hit_refst[idd];
                       int gap_size = abs(offset2-offset1);
                       int size_flag = 0;

                       if(mate_std[reads2wgs[idt]] > 0)
                       {
                         insertSize1 = mate_insert[reads2wgs[idt]] + mate_std[reads2wgs[idt]];
                         insertSize2 = mate_insert[reads2wgs[idt]] - mate_std[reads2wgs[idt]];
                       }
                       else
                       {
                         insertSize1 = 1.2*insert_size;
                         insertSize2 = 0.8*insert_size;
                       }
                       if((gap_size <= insertSize1)&&(gap_size>=insertSize2))
                         size_flag = 1;
//            printf("www2: %d %d %d %s %s %d %d\n",ii,jj,gap_size,rdname[idd],rdname[idt],hit_refst[idd],hit_refst[idt]);
                       if(size_flag)
                       {
                         pair1 = idt;
                         pair2 = idd;
                         n_pair++;
                       }
                     }
                  }
               }
          }
          if(n_pair==1)
          {
            int newscore;
            int score1 = map_score[pair1];
            int score2 = map_score[pair2];
            if(score1>score2)
              newscore = score1;
            else
              newscore = score2;
            if(newscore<20)
              newscore = 20;
            map_score[pair1] = newscore;
            map_score[pair2] = newscore;
            map_unique[pair1] = 1;
            map_unique[pair2] = 1;
          }
        }
        else if((j-i)==1)
          map_unique[readIndex[i]] = 1;
        i=j-1;
     }
/*   read the cigar line file   */
     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }
/*   read the cigar line file   */
     if((namef2 = fopen(argv[args+2],"w")) == NULL)
     {
       printf("ERROR main:: alignment file 2 \n");
       exit(1);
     }
/*   read the SNP output file         */
     i=0;
     n_find = 0;
     while(!feof(namef))
     {
       fgets(line,2000,namef);
//         printf("%s",line);
       if(feof(namef)) break;
       if((map_unique[i]==1)&&(map_rdpair[i]<=1))
       {
         char score[3] = {0};
	 sprintf(score,"%d,",map_score[i]);
	 line[7] = score[0];
	 line[8] = score[1];
         fprintf(namef2,"%s",line);
	 n_find++;
       }
       i++; 
     }
     printf("number of reads uniquely placed on to genome: %d\n",n_find);
     fclose(namef);
     fclose(namef2);
}


/*   Subroutine to process alignment information */
/* ====================================================  */
void Mates_Process(char **argv,int args,int nRead)
/* ====================================================  */
{
     int i,j,k,n_reads = nRead;
     char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
//     char **rdname;
     char **DBname,**RPname,*ptr,line[2000];
     FILE *namef;
     int *read_index;
     int n_find,idd,stopflag,num_align;
 
     n_reads = nRead + num_mates;
     RPname=cmatrix(0,num_mates,0,Max_N_NameBase);
     DBname=cmatrix(0,n_reads,0,Max_N_NameBase);

     if((read_index= (int *)calloc(n_reads,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - readIndex\n");
       exit(1);
     }
     if((mate_pairs= (int *)calloc(num_mates,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - readIndex\n");
       exit(1);
     }
     if((mate_index= (int *)calloc(num_mates,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - mate_index\n");
       exit(1);
     }
     if((hit_frdex= (int *)calloc(num_mates,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_frdex\n");
       exit(1);
     }
     if((mate_insert= (int *)calloc(num_mates,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - mate_insert\n");
       exit(1);
     }
     if((mate_std= (int *)calloc(num_mates,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - mate_std\n");
       exit(1);
     }
     if((reads2wgs= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ctg_index\n");
       exit(1);
     }

     if((namef = fopen(argv[args+1],"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

/*   read the SNP output file         */
     num_align=0;
     while(!feof(namef))
     {
       int nPair=0;
       char line2[2000],base[500];
      
       fgets(line,2000,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
       {
       }
       i=0;
       for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
       {
          if(i==0)
          {
            memset(base,'\0',500);
            strcat(base,ptr);
            strcpy(DBname[num_align],ptr);
	    mate_pairs[num_align] = -1;
	    read_index[num_align] = num_align;
            num_align++;
          }
          else if(i==1)
          {
            memset(base,'\0',500);
            strcat(base,ptr);
            strcpy(RPname[num_align-1],ptr);
          }
          else if(i==2)
          {
            memset(base,'\0',500);
            strcat(base,ptr);
            mate_insert[num_align-1] = atoi(ptr);
          }
          else if(i==3)
          {
            memset(base,'\0',500);
            strcat(base,ptr);
            mate_std[num_align-1] = atoi(ptr);
          }
       }
     }
     fclose(namef);

     ArraySort_String(num_align,RPname,read_index);

/*   sort out read pairs */
     n_reads = num_align;
     for(i=0;i<n_reads;i++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        while((j<n_reads)&&(stopflag==0))
        {
          if(strcmp(RPname[j],RPname[i])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)==2)
        {
	  mate_pairs[read_index[i+1]] = read_index[i];
	  mate_pairs[read_index[i]] = read_index[i+1];
	  mate_index[read_index[i+1]] = read_index[i];
	  mate_index[read_index[i]] = read_index[i];
          hit_frdex[read_index[i]] = 1;
	  hit_frdex[read_index[i+1]] = 2;
        }
        i=j-1;
     }

/*   sort out read pairs */
     for(i=0;i<num_mates;i++)
     {
        read_index[i] = i;
//	printf("mates: %s %d\n",DBname[i],mate_index[i]);
     }
     
     for(i=0;i<nRead;i++)
     {
        read_index[i+num_mates] = i+num_mates;
	strcpy(DBname[i+num_mates],rdname[i]);
     }
     n_reads = nRead + num_mates;
     ArraySort_String(n_reads,DBname,read_index);

     n_find = 0;
     idd = 0;
     for(i=0;i<n_reads;i++)
     {
        if(read_index[i] < num_mates)
	  idd = read_index[i];
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        while((j<n_reads)&&(stopflag==0))
        {
          if(strcmp(DBname[j],DBname[i])==0)
          {
	    if(read_index[j] < num_mates)
	      idd = read_index[j];
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>=2)
	{
	  for(k=i;k<j;k++)
	  {
	     if(read_index[k] >= num_mates)
	     {
	       reads2wgs[read_index[k]-num_mates]=idd;
	       n_find++;
	     }
	  }
	}
        i=j-1;
     }
     printf("reads found: %d\n",n_find);
}


#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, B64_long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

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
void ArraySort_Mix(int n, B64_long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

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
void ArraySort_Mix3(int n, B64_long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     B64_long a,temp,MIN=7;

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
int     **mmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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

/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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
char    **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch)
{
        B64_long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
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

