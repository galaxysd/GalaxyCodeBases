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
static int *sm,*head,*qm,*readlength;
static char *dna,**S_Name,**rdnames;
static long *h_dna;
static int *cons_head,*cons_list,*rd_head,*rd_index;
static int *rd_hitoffset,*rd_hitindex,*rd_hitorder,*rd_hitmatch;
static int *ctg_head,*ctg2wgs_index;

/* SSAS default parameters   */
static int IMOD=0;
static int n_SNPs=0;
static int num_subs=0;
static int N_outread=400000;
static int flag_fastq = 1;
static int release_flag=0;
static float rat=0.6;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

static fasta **expp;
fasta *expt;

static char rc_char[500];
static char rc_sub[500];

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
    FILE *fp,*namef;
    long dataSize,totalBases;
    int i,j,nSeq,args,qthresh=0;
    int nseq,kick_flag = 0;
    fasta *seq,*seqp;
    void decodeReadpair(int nSeq);
    void HashFasta_Head(int i, int nSeq);
    void HashFasta_Table(int i, int nSeq);
    void Search_SM(fasta *seq,int nSeq);
    void Assemble_SM(int arr,int brr);
    void Read_Index(char **argv,int args,int nRead,int nSeq);
    void Memory_Allocate(int arr);
    char **cmatrix(long nrl,long nrh,long ncl,long nch);
    fasta *segg;
    long Size_pdata,Size_q_pdata;
    unsigned int *pidata;
    int num_seqque,nline;
    char *pdata,*st,Nameout[100];
    int num_out,num_file,rc,ac,n_input,out_flag;

    seq=NULL;
    if(argc < 2)
    {
      printf("USAGE\n");
      printf("\n");
      printf(" %s <-split 4000> <output fastq head> <input fastq files>\n",argv[0]);
      printf("-split 4000 means cut the file into small patches each with 4000 sequences\n");
      printf("\n");
//      information();
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
       else if(!strcmp(argv[i],"-split"))
       {
         sscanf(argv[++i],"%d",&N_outread);
         args=args+2;
       }
    }

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
    nseq=0;
    nSeq = num_seqque;
    printf("Number of shotgun reads  %d \n",nSeq);

    num_out=0; 
    num_file=0;
    out_flag=1;
    for(i=0;i<nSeq;i++)
    {
       if(((num_out%N_outread)==0)&&(out_flag==1))
       {
         memset(Nameout,'\0',100);
         sprintf(Nameout,"%s_%04d",argv[args],num_file);
         if((namef = fopen(Nameout,"w")) == NULL)
         {
             printf("ERROR main:: reads group file \n");
             exit(1);
         }
         out_flag=0;
         num_file++;
       }
       seqp = seq + i;
       if((seqp->qual != NULL)&&(flag_fastq==1))
       {
         if((seqp->qual != NULL)&&(flag_fastq==1))
         {
           fprintf(namef,"@%s\n",seqp->name);
           for(rc=0;rc<seqp->length;rc++)
              fprintf(namef,"%c",seqp->data[rc]);
           fprintf(namef,"\n");
           fprintf(namef,"+\n");
           for(rc=0;rc<seqp->length;rc++)
           {
              if(seqp->finished)
                putc(40+041,namef);
              else
                putc(seqp->qual[rc]+041,namef);
           }
           fprintf(namef,"\n");
         }
         else
         {
           nline=seqp->length/60;
           st = seqp->data;
           fprintf(namef,">%s\n",seqp->name);
           for(rc=0;rc<nline;rc++)
           {
              for(j=0;j<60;j++,st++)
                 fprintf(namef,"%c",*st);
              fprintf(namef,"\n");
           }
           for(j=0;j<(seqp->length-(nline*60));j++,st++)
              fprintf(namef,"%c",*st);
           if((seqp->length%60)!=0)
             fprintf(namef,"\n");
         } 
         out_flag=1;
         num_out++;
       }
       if(((num_out%N_outread)==0)&&(out_flag==1))
       {
         fclose(namef);
       }
    } 
    if((nSeq%N_outread)!=0)
      fclose(namef);

    if(seq){
        free(seq->name);
        free(seq);
        seq = NULL;
    }    
    printf("Job finished for %d contigs!\n",nSeq);
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
     void ArraySort_String(int n,char Pair_Name[][Max_N_NameBase],int *brr);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     char DBname[n_reads][Max_N_NameBase];
     char tempct[60],tempc1[60],tempc2[60],tempc3[60],tempc4[60],tempc5[60];
     int temp1,temp2,temp3;
     int i_contig,i_reads,num_rd_find,stopflag;
     int *readIndex;
     int mapindex=0;

     if((readIndex= (int *)calloc(n_reads,sizeof(int))) == NULL)
     {
       printf("Error Contig_Merge: calloc - readIndex\n");
       exit(1);
     } 
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
     for(j=0;j<nRead;j++)
     {
        strcpy(DBname[j+nSeq],rdnames[j]);
        readIndex[j+nSeq]=j+nSeq;
     }
     n_reads=nSeq+nRead;
     printf("before sort: %d %d\n",nSeq,nRead);
     ArraySort_String(n_reads,DBname,readIndex);

     num_rd_find=0;
     mapindex=0;
     for(i=0;i<n_reads-1;i++)
     {
        if(readIndex[i]>=nSeq)
          mapindex = readIndex[i];
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        while((j<n_reads)&&(stopflag==0))
        {
          if(strcmp(DBname[j],DBname[i])==0)
          {
            if(readIndex[j]>=nSeq)
              mapindex = readIndex[j];
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>1)
        {
          for(k=i;k<j;k++)
          {
             if(readIndex[k]<nSeq)
             {
               ctg2wgs_index[mapindex-nSeq]=readIndex[k];
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
void s_swap(char Pair_Name[][Max_N_NameBase], int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char Pair_Name[][Max_N_NameBase], int *brr)
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

