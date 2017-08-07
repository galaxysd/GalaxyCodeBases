/***********************************************************************\
*                                                                     * 
*                     PROJECT   ssaha_pileup                          *
*                                                                     * 
*---------------------------------------------------------------------*
*                                                                     * 
*                                By                                   *
*                                                                     *
*                            Zemin Ning                               *
*                                                                     *
*               Copyright (C) 2006 by Genome Research Limited         *
*                         All rights reserved                         *
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
#define Max_N_NameBase 60
#define Max_N_Pair 100
static char **S_Name,**rdnames;
static int *rd_head,*rd_index;
static int *ctg2wgs_index,*dup_reads;

fasta *seq;
/* SSAS default parameters   */
static int fastq_flag =1;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

int main(int argc, char **argv)
{
    FILE *fp,*namef,*fpOutfast;
    long totalBases;
    int i,nSeq,args,nRead,n_len;
    int ac,n_contig,n_reads,nseq,n_input;
    fasta *seqp;
    char *st,line[2000]={0};
    void Read_Index(char **argv,int args,int nRead,int nSeq);
    char **cmatrix(long nrl,long nrh,long ncl,long nch);
    int  num_out,num_file,out_flag;
    fasta *segg;
    long Size_pdata,Size_q_pdata;
    unsigned int *pidata;
    int num_seqque;
    char *pdata;

    seq=NULL;
    fflush(stdout);
    system("ps aux | grep sreads; date");
    if(argc < 2)
    {
      printf("Usage: %s <read_name_file> <input_fastq_file> <output_fastq_file\n",argv[0]);
      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-fastq"))
       {
         sscanf(argv[++i],"%d",&fastq_flag); 
         args=args+2;
       }
    }

    n_input=0;
    nseq=0;
    
    if((fp=fopen(argv[args+1],"rb"))==NULL) printf("Cannot open file\n");
    fseek(fp, 0, SEEK_END);
    Size_q_pdata = ftell(fp) + 1;
    fclose(fp);
    if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
       printf("calloc pdata\n");
    num_seqque = extractFastq(argv[args+1],pdata,Size_q_pdata);
    if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL)
      printf("calloc segg\n");
    if((seq=decodeFastq(argv[args+1],&num_seqque,&totalBases,pdata,Size_q_pdata,segg))==NULL)
      printf("no query data found.\n");
    nSeq = num_seqque;
       
    fastaUC(seq,nSeq);

    nRead=nSeq;
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

/*  process contigs one by one   */
    ac=args+2;
    printf("file name  %s \n",argv[ac]);
    if((fpOutfast = fopen(argv[args+2],"w")) == NULL)
    {
        printf("ERROR main:: reads group file \n");
        exit(1);
    }
    n_contig=0;
    num_out=0;
    num_file=0;
    out_flag=1;
    for(i=0;i<n_reads;i++)
    {
       if(ctg2wgs_index[i]>=0) 
       {
          int nline,k,g,rc;
          int idd=ctg2wgs_index[i];
          n_contig++;
          seqp=seq+idd;
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
              fprintf(fpOutfast,"+\n");
            else
              fprintf(fpOutfast,"+\n");
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
       }
       else
       {
            printf("missed: %d %s\n",i,rdnames[i]);
       }
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
     while(fscanf(namef,"%s",rdnames[i_reads])!=EOF)
     {
//    printf("after read: %d %s\n",i_reads,rdnames[i_reads]);
       i_reads++;
     }
     fclose(namef);

     for(j=0;j<nSeq;j++)
     {
        seqp=seq+j;
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
          seqp=seq+j;
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

