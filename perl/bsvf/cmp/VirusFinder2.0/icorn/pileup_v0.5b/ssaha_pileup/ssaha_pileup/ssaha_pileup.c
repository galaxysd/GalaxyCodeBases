/***********************************************************************\
 *                                                                     * 
 *                         PROJECT   ssaha_pileup                      *
 *                                                                     * 
 *---------------------------------------------------------------------*
 *                                                                     *
 *                                By                                   *
 *                                                                     *
 *                            Zemin Ning                               *
 *                                                                     *
 *          Copyright (C) 2008 by Genome Research Limited              *
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
#define Max_N_NameBase 50
#define Max_N_Pair 100
static char **cell_name;
static char **rdname;
static int *map_score,*read_arch,*hit_refst,*hit_refed,*hit_rcdex,*hit_quest,*hit_queed;
static int *readIndex,*ctg_index,*cell_index,*read_filt;
static int *cigar_list,*contig_mapst,*contig_maped,*contig_list,*contig_head;
static B64_long *cigar_head,sBase,qBase;
static char *dataline,*cigar_line;

/* SSAS default parameters   */
static int IMOD=0;
static int ISUB=0;
static int mapNumber=2;
static int copyNumber=20;
static B64_long line_len=0;
static int num_reads=0;
static int num_cline=0;
static int low_qual=23;
static int arch_flag = 1;
static int set_identy = 0;
static int set_Qscore = 30;
static int solexa_flag = 0;
static int file_flag = 1;
static int cover_flag = 1;
static int n_contigs = 0;
static int min_length = 20;
static int phred_flag  = 0;
static int maq_flag  = 0;
static int cell_flag  = 0;
static int trans_flag  = 0;
static int cons_flag  = 0;
static int num_depth = 500;
static B64_long genome_offset;
static char strain_name[100];
static float allele_rate = 0.21;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *sub,*seq;

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
    void ArraySort_Mix(int n, B64_long *arr, int *brr);
    void SNP_Consensus(char **argv,int args,int nSeq);
    void Read_Pairs(char **argv,int args,int nLib,int nSeq);
    void Align_Process(char **argv,int args,int nRead);

    seq=NULL;
//    fflush(stdout);
//    system("ps aux | grep eulerSNP; date");
    memset(strain_name,'\0',100);
    if(argc < 2)
    {
      printf("Usage: %s <cigar_line_file> <ref_sequence> <reads_fastq_file> \n",argv[0]);
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
         sscanf(argv[++i],"%d",&cover_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-arch"))
       {
         sscanf(argv[++i],"%d",&arch_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-phred"))
       {
         sscanf(argv[++i],"%d",&phred_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-maq"))
       {
         file_flag = 0;
	 cons_flag = 1;
         sscanf(argv[++i],"%d",&maq_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-score"))
       {
         sscanf(argv[++i],"%d",&set_Qscore);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-strain"))
       {
         sscanf(argv[++i],"%s",strain_name);
         cell_flag  = 1;
         args=args+2;
       }
       else if(!strcmp(argv[i],"-depth"))
       {
         sscanf(argv[++i],"%d",&num_depth);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-idt"))
       {
         sscanf(argv[++i],"%d",&set_identy);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-trans"))
       {
         sscanf(argv[++i],"%d",&trans_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-qual"))
       {
         sscanf(argv[++i],"%d",&low_qual);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-allele"))
       {
         sscanf(argv[++i],"%f",&allele_rate);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&file_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-cons"))
       {
         sscanf(argv[++i],"%d",&cons_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-solexa"))
       {
         sscanf(argv[++i],"%d",&solexa_flag);
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
     void ArraySort_Mix(int n,B64_long *arr,int *brr);
     int  **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
//     char **rdname;
     char **DBname,*ptr,RC;
     char *line,*st,*ed,*match_text = "M ";
     FILE *fp,*namef;
     B64_long *read_offsets,big_num;
     int n_find,idd,stopflag,num_align,refhit1,refhit2;
     int readhit1 = 0,readhit2 = 0;
     fasta *segg;
     B64_long Size_q_pdata;
     int num_seqque;
     char *pdata;
     void Arch_Align(int nRead);

     if((fp=fopen(argv[args+1],"rb"))==NULL) printf("Cannot open file\n");
       fseek(fp, 0, SEEK_END);
     Size_q_pdata = ftell(fp) + 1;
     fclose(fp);
     if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
       printf("calloc pdata\n");
     num_seqque = extractFastq(argv[args+1],pdata,Size_q_pdata);
     if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL)
       printf("calloc segg\n");
     if((sub=decodeFastq(argv[args+1],&num_seqque,&sBase,pdata,Size_q_pdata,segg))==NULL)
       printf("no query data found.\n");
     n_contigs = num_seqque;
     fastaUC(sub,n_contigs);
     n_reads = nRead + n_contigs;

     if((fp=fopen(argv[args+2],"rb"))==NULL) printf("Cannot open file\n");
       fseek(fp, 0, SEEK_END);
     Size_q_pdata = ftell(fp) + 1;
     fclose(fp);
     if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
       printf("calloc pdata\n");
     num_seqque = extractFastq(argv[args+2],pdata,Size_q_pdata);
     if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL)
       printf("calloc segg\n");
     if((seq=decodeFastq(argv[args+2],&num_seqque,&qBase,pdata,Size_q_pdata,segg))==NULL)
       printf("no query data found.\n");
     fastaUC(seq,num_seqque);


     RC = '+';
     cell_name = cmatrix(0,2,0,Max_N_NameBase);
     rdname=cmatrix(0,nRead,0,Max_N_NameBase);
     DBname=cmatrix(0,n_reads,0,Max_N_NameBase);

     if((line= (char *)calloc(2000,sizeof(char))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - line\n");
       exit(1);
     }
     if((readIndex= (int *)calloc(n_reads,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - readIndex\n");
       exit(1);
     }
     if((map_score= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - map_score\n");
       exit(1);
     }
     if((read_arch= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - read_arch\n");
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
     if((read_filt= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - read_filt\n");
       exit(1);
     }
     if((cigar_list= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - cigar_list\n");
       exit(1);
     }
     if((cigar_head= (B64_long *)calloc(nRead,sizeof(B64_long))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - cigar_head\n");
       exit(1);
     }
     if((read_offsets= (B64_long *)calloc(nRead,sizeof(B64_long))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - read_offsets\n");
       exit(1);
     }
     if((contig_list= (int *)calloc(n_contigs+20,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - contig_list\n");
       exit(1);
     }
     if((contig_head= (int *)calloc(n_contigs+1,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - contig_head\n");
       exit(1);
     }

/*   read the cigar line file   */
     if((namef = fopen(argv[args],"r")) == NULL)
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
     if((namef = fopen(argv[args],"r")) == NULL)
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
//	        printf("score: %d %s",len,line3);
	      }
            }
            else if(i==1)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              strcpy(rdname[num_align],ptr);
//	        printf("name: %s %d\n",rdname[num_align],num_align);
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
              strcat(base,ptr);
              strcpy(DBname[num_align],ptr);
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
              readIndex[num_align] = num_align;
              num_align++;
            }
         }
       }
     }
     fclose(namef);

     memset(ctg_index,-1,num_align*4);
/*   sort out reference head  */
     num_cline = 1;
     strcpy(cell_name[0],"READS");

     printf("number of cell lines: %d\n",num_cline);

/*   sort out contig name match */
     if(nRead!= num_align)
     {
       printf("Number of reads not the same. Job stopped! %d %d\n",nRead,num_align);
       exit(1);
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
                 ctg_index[readIndex[k]] = idd;
                 n_find++;
               }
            }
          }
        }
        i=j-1;
     }
     printf("number of reads aligned to the genome: %d\n",n_find);

/*   force arch alignment   */
     if(arch_flag&&solexa_flag)
       Arch_Align(nRead);

/*   sort out match offsets on contigs/chromosomes */
     for(j=0;j<nRead;j++)
     {
        B64_long nn;

	if(ctg_index[j]>=0)
	  nn = ctg_index[j];
	else
	  nn = n_contigs+10;
        readIndex[j]=j;
	read_offsets[j] = (nn<<40) + hit_refst[j]; 
     }
     ArraySort_Mix(nRead,read_offsets,readIndex);

     big_num = 1L<<30;
     for(i=0;i<nRead;i++)
     {
        stopflag=0;
        j=i+1;
        while((j<nRead)&&(stopflag==0))
        {
          if((read_offsets[j]-read_offsets[i])<big_num)
          {
            j++;
          }
          else
            stopflag=1;
        }
	if((j-i)>0)
        {
	  int idd = (read_offsets[i]>>40);
          contig_list[idd] = j-i;
//	  printf("offset: %d %d %d\n",idd,contig_list[idd],idd);
        }
        i=j-1;
     }
     free(read_offsets);
}


/*   Subroutine to find the contig/supercontig structure */
/* ====================================================  */
void SNP_Consensus(char **argv,int args,int nSeq)
/* ====================================================  */
{
     int i,j,k,m,c,b;
     B64_long *ref_head;
     int *ref_list;
     FILE *fp;
     void ArraySort_Mix(int n,B64_long *arr,int *brr);
     void ArraySort2_Int2(int n, int *arr, int *brr);
     B64_long st,ed,kk,tBase,map_size;
     int  **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
     void Map_Memory(int nSeq);
     int set_score = 0,err_flag = 0;
     fasta *seqp;
     int *read_qual,read_qsum[4],read_nhit[6],read_lowq[6],read_scor[4];
     char *refe_base,*read_base,*maq_base,*score_base,*maq_score,*maq_quals,SNP_base,SNP2_base;
     int num_hits,num_maps,sum_Q_set,ifactor;
     int read_qindex[4],read_qsort[4],*read_offset,*map_index,size_index;
     int *cell_cover,*base_cover;
     int **cell_nhit,num_mapreads; 
     fasta *segg;
     B64_long Size_q_pdata,num_mapbases;
     int num_seqque,pos_input,ctg_input,score_offset;
     char *pdata;
     char *ptr,line2[2000],line[2000],score[3] = {0};

     cell_nhit=imatrix(0,4,0,num_cline);

/*     for(j=0;j<nSeq;j++)
     {
        int idk = readIndex[j];
     }    */
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
     if((ref_list= (int *)calloc(n_contigs,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ref_list\n");
       exit(1);
     }
     if((ref_head= (B64_long *)calloc(n_contigs,sizeof(B64_long))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ref_head\n");
       exit(1);
     }
     map_size = sBase+n_contigs+1;
     if((contig_mapst= (int *)calloc(map_size,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - contig_mapst\n");
       exit(1);
     }
     if((contig_maped= (int *)calloc(map_size,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - contig_maped\n");
       exit(1);
     }

     seqp = sub;
     ref_list[0] = seqp->length;
     ref_head[0] = 0L;
     contig_head[0] = 0;
     if(trans_flag)
       min_length = 12;
     for(i=1;i<n_contigs;i++)
     {
        seqp = sub + i;
        ref_list[i] = seqp->length;
        ref_head[i] = ref_head[i-1] + ref_list[i-1];
        contig_head[i] = contig_head[i-1] + contig_list[i-1];
     }

     for(i=0;i<n_contigs;i++)
     {
        int seq_len; 
        seq_len = ((sub+i)->length)+1;
        for(j=0;j<seq_len;j++)
           contig_mapst[ref_head[i]+j] = nSeq;
     }
     for(c=0;c<n_contigs;c++)
     {
        int n_links = contig_list[c];
	B64_long offset;

        for(j=0;j<n_links;j++)
        {
           int match_flag = 0;
	   int q_flag = 1;
	   int idt = contig_head[c]+j; 
     	   i = readIndex[idt];
	   seqp = seq + i;
           if(abs(hit_refed[i]-hit_refst[i])>=min_length)
             match_flag = 1;
	   else
	     read_filt[i] = 1;
           set_score = 0;
	   if(solexa_flag==0)
	     q_flag = 1;
	   else
	   {
	     seqp = seq + i;
	       if(hit_quest[i]>=5)
	       {
	         int qedge1 = 0;
	         for(k=(hit_quest[i]-5);k<(hit_quest[i]-1);k++)
	         {
	            if(seqp->qual[k]>=low_qual)
	              qedge1++;
	         }
	         if(qedge1==4)
	           q_flag = 0;
	       }
               if((seqp->length-hit_queed[i])>=5)
	       {
	         int qedge2 = 0;
	         for(k=hit_queed[i];k<(hit_queed[i]+4);k++)
	         {
	            if(seqp->qual[k]>=low_qual)
	              qedge2++;
	         }
	         if(qedge2==4)
	           q_flag = 0;
	       }
	   }
	   if(!trans_flag)
	   {
	     if(q_flag == 0)
	       read_filt[i] = 1;
	   }
           if((q_flag)&&(map_score[i]>=set_score)&&(match_flag)&&(ctg_index[i]>=0))
           {
	     offset = ref_head[c];
             st = hit_refst[i];
             ed = hit_refed[i];
             for(kk=st;kk<=ed;kk++)
             {
	        if(idt<contig_mapst[offset+kk])
	          contig_mapst[offset+kk] = idt;
	        if(idt>contig_maped[offset+kk])
	          contig_maped[offset+kk] = idt;
             }
           }
        }
     }

     size_index = 0;
     for(c=0;c<n_contigs;c++)
     {
        for(b=1;b<=(sub+c)->length;b++)
        {
           genome_offset = ref_head[c] + b;
           num_maps = contig_maped[genome_offset]-contig_mapst[genome_offset]+1;
	   if(num_maps>size_index)
	     size_index = num_maps;
	}
     }
     num_depth = size_index+10;
     if((read_offset= (int *)calloc(num_depth,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - read_offset\n");
       exit(1);
     }
     if((read_qual= (int *)calloc(num_depth,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - read_qual\n");
       exit(1);
     }
     if((map_index= (int *)calloc(num_depth,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - map_index\n");
       exit(1);
     }
     if((refe_base= (char *)calloc(num_depth,sizeof(char))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - refe_base\n");
       exit(1);
     }
     if((read_base= (char *)calloc(num_depth,sizeof(char))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - read_base\n");
       exit(1);
     }
     if((maq_base= (char *)calloc(num_depth,sizeof(char))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - maq_base\n");
       exit(1);
     }
     if((maq_quals= (char *)calloc(num_depth,sizeof(char))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - maq_quals\n");
       exit(1);
     }
     if((maq_score= (char *)calloc(num_depth,sizeof(char))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - maq_score\n");
       exit(1);
     }
     if((score_base= (char *)calloc(3*num_depth,sizeof(char))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - score_base\n");
       exit(1);
     }

     printf("Maximum read depth: %d\n",num_depth);
     num_mapreads = 0;
     for(i=0;i<nSeq;i++)
     {
        if(read_filt[i] == 0)
	  num_mapreads++;
     }
     if((dataline = (char *)calloc(1000,sizeof(char))) == NULL)
     {
       printf("ERROR ssaha: calloc - dataline\n");
       exit(1);
     }

     if(cons_flag)
     {
       printf("Start consensus process:\n");
       if(file_flag)
       {
         printf("=================================================================================================\n");
         if(phred_flag)
           printf("  refe_name offset N_reads refe_base N_'A' N_'C' N_'G' N_'T' N_'-' N_'N' Q_'A' Q_'C' Q_'G' Q_'T'\n");
	 else
           printf("  refe_name offset N_reads refe_base N_'A' N_'C' N_'G' N_'T' N_'-' N_'N' N_'a' N_'c' N_'g' N_'t'\n");
         printf("=================================================================================================\n");
       }
       else
       {
         if(maq_flag == 0)
	 {
           printf("============================================================\n");
           printf("  refe_name offset N_reads refe_base read_base mapping_score\n");
           printf("============================================================\n");
	 }
       }
     }
     else
     {
       printf("Start SNP calling process:\n");
       if(file_flag)
       {
         printf("==================================================================================================================\n");
         if(phred_flag)
           printf("  refe_name SNP_score offset N_reads refe_base SNP_base N_'A' N_'C' N_'G' N_'T' N_'-' N_'N' Q_'A' Q_'C' Q_'G' Q_'T'\n");
	 else
           printf("  refe_name SNP_score offset N_reads refe_base SNP_base N_'A' N_'C' N_'G' N_'T' N_'-' N_'N' N_'a' N_'c' N_'g' N_'t'\n");
         printf("==================================================================================================================\n");
       }
       else
       {
         printf("===============================================================================\n");
         printf("  refe_name SNP_score offset N_reads refe_base SNP_base read_base mapping_score\n");
         printf("===============================================================================\n");
       }
     }
     err_flag = 0;
     num_mapbases = 0;
     sum_Q_set = set_Qscore;
//     memset(refe_base,'\0',num_depth);
     memset(read_qual,0,4*num_depth);
     for(c=0;c<n_contigs;c++)
     {
        for(b=1;b<=(sub+c)->length;b++)
        {
           int max_qsum,max_base,rqual,max_qsum2,out_flag;
           int map_base,max_base2,ave_score,sum_score;
	   char con;

           max_base2 = -1;
	   ave_score = 0;
	   sum_score = 0;
           genome_offset = ref_head[c] + b;
           pos_input = b;
           ctg_input = c;
	   con = (sub+c)->data[b-1];
           num_maps = contig_maped[genome_offset]-contig_mapst[genome_offset]+1;
           if(num_maps<=0)
	   {
	     if((cons_flag)&&(cover_flag))
	     {
	       if(maq_flag == 0)
	         printf("cons: %s %d %d %c Zero coverage\n",(sub+c)->name,b,0,con);
	       else
	         printf("%-20s %-8d %c       %-8d @  @  @\n",(sub+c)->name,b,con,0);
	     }
             continue;
	   }
           memset(map_index,0,4*num_maps);
	   num_hits = 0;
           for(j=0;j<num_maps;j++)
           {
	      int idk = contig_mapst[genome_offset]+j;
              i = readIndex[idk];
	        
              if((b>hit_refed[i])||(b<hit_refst[i]))
	        continue;
	      if(read_filt[i] == 0)
	      {
	        map_index[num_hits] = i;
                num_hits++;
	      }
	   }
	   if(num_hits<=0)
	   {
	     if(cons_flag)
	     {
	       if(maq_flag == 0)
	         printf("cons: %s %d %d %c Zero coverage\n",(sub+c)->name,b,0,con);
	       else
	         printf("%-20s %-8d %c       %-8d @  @  @\n",(sub+c)->name,b,con,0);
	     }
             continue;
	   }
	   if(solexa_flag)
	   {
	     ifactor = 10;
	     if(num_hits<10)
	       sum_Q_set = 30;
             else
               sum_Q_set = 50;
	   }
	   else
	     ifactor = 5;
           memset(read_base,'\0',num_depth);
	   if(maq_flag)
	   {
             memset(maq_base,'\0',num_depth);
             memset(maq_score,'\0',num_depth);
             memset(maq_quals,'\0',num_depth);
	   }
	   if(file_flag == 0)
             memset(score_base,'\0',3*num_depth);
           memset(cell_cover,0,4*num_cline);
           memset(cell_nhit[0],0,4*num_cline);
           memset(cell_nhit[1],0,4*num_cline);
           memset(cell_nhit[2],0,4*num_cline);
           memset(cell_nhit[3],0,4*num_cline);
	   score_offset = 0;
           for(j=0;j<num_hits;j++)
           {
              int qst,sst,mpos,m_bases,r_bases,s_bases;
              int nPair = 0,r_mod,r_len,idt;

              read_offset[j] = 0;
              i = map_index[j];
              idt = cell_index[i];
              seqp = seq + i;
	      sum_score = sum_score + map_score[i];
	      if(file_flag==0)
	      {
	        memset(score,'\0',3);
		score[0] = map_score[i]/10 + '0';
		score[1] = map_score[i]%10 + '0';
		score_base[score_offset] = score[0];
		score_base[score_offset+1] = score[1];
		score_base[score_offset+2] = ',';
	        score_offset = score_offset + 3; 
	      }
//           printf("name: %d %s %s\n",j,score,score_base); 
//              seqp->qual[0] = 40;
              if(solexa_flag)
                memset(line,'\0',100);
              else
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
//		         if(read_arch[i] > 0)
//			   printf("ext: %d %d %d %s %d %d %d %d\n",b,r_len,read_arch[i],rdname[i],hit_quest[i],hit_queed[i],hit_refst[i],hit_refed[i]);
		         if(read_arch[i] == 1)
			 {
			   if(m == 11)
			     r_len = r_len+2;
			 }
			 else if(read_arch[i] == 2)
			 {
			   if(m == 11) 
			     r_len = r_len+3;
			 }
			 else if(read_arch[i] == 3)
			 {
//	                      printf("www1: %d %d %s %d %d\n",i,read_arch[i],rdname[i],r_len,m);
			   if(hit_queed[i] == seqp->length)
			     r_len = r_len+2;
			 }
			 else if(read_arch[i] == 4)
			 {
	      //                printf("www2: %d %d %s %d %d\n",i,read_arch[i],rdname[i],r_len,m);
			   if(hit_queed[i] == seqp->length)
			     r_len = r_len+3;
			 }
                         for(k=0;k<r_len;k++)
                         {
                            if((s_bases+k)==(pos_input-1))
                            {
                              read_base[j] = seqp->data[k+r_bases]; 
                              read_qual[j] = seqp->qual[k+r_bases];
//                              refe_base[j] = (sub+c)->data[k+s_bases];
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
//                              refe_base[j] = (sub+ctg_input)->data[k+s_bases];
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
		         if(read_arch[i] == 1)
			 {
			   if(m == 11)
			     r_len = r_len+2;
			 }
			 else if(read_arch[i] == 2)
			 {
			   if(m == 11)
			     r_len = r_len+3;
			 }
			 else if(read_arch[i] == 3)
			 {
			   if(hit_queed[i] == seqp->length)
			     r_len = r_len+2;
			 }
			 else if(read_arch[i] == 4)
			 {
			   if(hit_queed[i] == seqp->length)
			     r_len = r_len+3;
			 }
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
//                              refe_base[j] = (sub+ctg_input)->data[k+s_bases];
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
//                              refe_base[j] = (sub+c)->data[k+s_bases];
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
                memset(read_scor,0,4*4);
                memset(read_nhit,0,6*4);
                memset(read_lowq,0,6*4);
              }
//	      if(b==129)
//	        printf("www: %d %d %s %c\n",i,read_arch[i],rdname[i],read_base[j]);
              if(read_base[j]!='-')
                cell_cover[idt]++;
              if(read_qual[j]>=low_qual)
                rqual = 18;
              else
                rqual = 6;
              if(read_base[j]=='A')
              {
                read_qsum[0] = read_qsum[0] + rqual;
                cell_nhit[0][idt]++;
		if(read_qual[j]>=low_qual)
                  read_scor[0] = read_scor[0] + map_score[i];
		else
		{
                  read_scor[0] = read_scor[0] + map_score[i]/2;
		  if(phred_flag==0)
                    read_lowq[0]++;
		}
		if(phred_flag)
		  read_lowq[0] = read_lowq[0] + read_qual[j];
                read_nhit[0]++;
              }
              else if(read_base[j]=='C')
              {
                read_qsum[1] = read_qsum[1] + rqual;
                cell_nhit[1][idt]++; 
		if(read_qual[j]>=low_qual)
                  read_scor[1] = read_scor[1] + map_score[i];
		else
		{
                  read_scor[1] = read_scor[1] + map_score[i]/2;
		  if(phred_flag==0)
                    read_lowq[1]++;
		}
                read_nhit[1]++;
		if(phred_flag)
		  read_lowq[1] = read_lowq[1] + read_qual[j];
              }
              else if(read_base[j]=='G')
              {
                read_qsum[2] = read_qsum[2] + rqual;
                cell_nhit[2][idt]++; 
		if(read_qual[j]>=low_qual)
                  read_scor[2] = read_scor[2] + map_score[i];
		else
		{
                  read_scor[2] = read_scor[2] + map_score[i]/2;
		  if(phred_flag==0)
                    read_lowq[2]++;
		}
		if(phred_flag)
		  read_lowq[2] = read_lowq[2] + read_qual[j];
                read_nhit[2]++;
              }
              else if(read_base[j]=='T')
              {
                read_qsum[3] = read_qsum[3] + rqual;
                cell_nhit[3][idt]++; 
		if(read_qual[j]>=low_qual)
                  read_scor[3] = read_scor[3] + map_score[i];
		else
		{
                  read_scor[3] = read_scor[3] + map_score[i]/2;
		  if(phred_flag==0)
                    read_lowq[3]++;
		}
		if(phred_flag)
		  read_lowq[3] = read_lowq[3] + read_qual[j];
                read_nhit[3]++;
              }
	      else if(read_base[j]=='-')
	      {
                read_nhit[4]++;
	      }
	      else
	      {
                read_nhit[5]++;
	      }

//              printf("%c",read_base[j]);
              if(maq_flag)
	      {
	        int mp_score;
                if(hit_rcdex[i]==0)
		{
		  if(read_base[j] == con)
		    maq_base[j] = ',';
		  else
		    maq_base[j] = read_base[j];
		}
		else
		{
		  if(read_base[j] == con)
		    maq_base[j] = '.';
		  else
		    maq_base[j] = tolower(read_base[j]);
		}
		if(map_score[i] <= 0)
		  mp_score = 0;
		else
		  mp_score = (2*map_score[i] -1);
		if(mp_score >= 93)
		  mp_score = 93;
		maq_score[j] = mp_score+041;
		maq_quals[j] = read_qual[j]+041;
//		putc(read_qual[j]+041,maq_quals[j]);
	      }
              if(read_qual[j]<low_qual) 
                 read_base[j] = tolower(read_base[j]);
           }
//           printf("\n");
           max_qsum = 0;
           max_base = 0; 
           max_base2 = -1;
	   ave_score = sum_score/num_hits;
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
	   refe_base[0] = (sub+c)->data[b-1];
           if(refe_base[0]=='A')
             map_base = 0;
           else if(refe_base[0]=='C')
             map_base = 1;
           else if(refe_base[0]=='G')
             map_base = 2;
           else if(refe_base[0]=='T')
             map_base = 3;

           out_flag = 0;
	   if(cons_flag)
	   {
	     if(file_flag == 0)
	     {
	       if(maq_flag == 0)
	         printf("cons: %s %d %d %c %s %s\n",(sub+c)->name,b,num_hits,refe_base[0],read_base,score_base);
	       else
	         printf("%-20s %-8d %c       %-8d @%s  @%s  @%s\n",(sub+c)->name,b,refe_base[0],num_hits,maq_base,maq_quals,maq_score);
	     }
	     else
               printf("cons::%05d %s %-9d %-7d %c     %-6d %-6d %-6d %-6d %-6d %-6d %-3d %-3d %-3d %-3d\n",c,(sub+c)->name,b,num_hits,refe_base[0],read_nhit[0],read_nhit[1],read_nhit[2],read_nhit[3],read_nhit[4],read_nhit[5],read_lowq[0],read_lowq[1],read_lowq[2],read_lowq[3]); 
	     num_mapbases++;
	   }
//           printf("cons: %s %d %d %c %c %s %d %d\n",(sub+c)->name,b,num_hits,refe_base[0],SNP_base,read_base,max_qsum,sum_Q_set); 
           else
	   {
             char Base2;
             if(max_base!=map_base)
             {
               if((read_nhit[max_base]>=1)&&(max_qsum>=sum_Q_set))
               {
                 int snp_score = read_scor[max_base]/ifactor;
		 int hit_depth;
		 if(snp_score>99)
		   snp_score = 99;
		 if(max_qsum2>0)
		 {
		   hit_depth = read_nhit[max_base2];
		   if(read_nhit[max_base2]==1)
		   {
		     if((max_qsum2==6)&&(max_base2!=map_base))
		       hit_depth = 0;
		   }
                   if(max_base2 == 0)
                     Base2 = 'A';
                   else if(max_base2 == 1)
                     Base2 = 'C';
                   else if(max_base2 == 2)
                     Base2 = 'G';
                   else if(max_base2 == 3)
                     Base2 = 'T';
                   else
                     Base2 = 'N';
		 }
		 else
		   hit_depth = 0;
		 if(file_flag == 1)
		 {
	  	   if(hit_depth==0)
                     printf("SNP_hom: %s %-2d %-9d %-7d %c %c   %-6d %-6d %-6d %-6d %-6d %-6d %-3d %-3d %-3d %-3d\n",(sub+c)->name,snp_score,b,num_hits,refe_base[0],SNP_base,read_nhit[0],read_nhit[1],read_nhit[2],read_nhit[3],read_nhit[4],read_nhit[5],read_lowq[0],read_lowq[1],read_lowq[2],read_lowq[3]); 
		   else
                     printf("SNP_hez: %s %-2d %-9d %-7d %c %c/%c %-6d %-6d %-6d %-6d %-6d %-6d %-3d %-3d %-3d %-3d\n",(sub+c)->name,snp_score,b,num_hits,refe_base[0],SNP_base,Base2,read_nhit[0],read_nhit[1],read_nhit[2],read_nhit[3],read_nhit[4],read_nhit[5],read_lowq[0],read_lowq[1],read_lowq[2],read_lowq[3]); 
		 }
		 else
		 {
	  	   if(hit_depth==0)
                     printf("SNP_hom: %s %d %d %d %c %c %s %d %d %s\n",(sub+c)->name,snp_score,b,num_hits,refe_base[0],SNP_base,read_base,hit_depth,read_nhit[max_base],score_base); 
		   else
                     printf("SNP_hez: %s %d %d %d %c %c/%c %s %d %d %s\n",(sub+c)->name,snp_score,b,num_hits,refe_base[0],SNP_base,Base2,read_base,hit_depth,read_nhit[max_base],score_base); 
		 }
               }
	       if(mapNumber==1)
   	       {
                 if((read_nhit[max_base]==1)&&(max_qsum>=18))
                   printf("SNP4: %s %d %d %d %c %c %s %d %d %s\n",(sub+c)->name,ave_score,b,num_hits,refe_base[0],SNP_base,read_base,read_nhit[max_base],read_nhit[max_base],score_base); 
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
	         if(((max_second==1)&&(cell_cover[max_cline]<=2))||(max_second>1)||(read_nhit[max_base2]>=2))
                 {
                   if((max_base==map_base)&&(read_nhit[max_base]<num_depth))
                   {
		     if(((read_nhit[max_base]>20)&&(rate>0.1))||(read_nhit[max_base]<=20)||(read_nhit[max_base2]>=2))
                     {
		       int n_hitsq = 0;
		       int n_hitrq = 0;
		       float rate2 = 0.0;

	  	       rate = read_nhit[max_base2];
		       rate = rate/num_hits;

		       for(k=0;k<num_hits;k++)
		       {
		          if(read_base[k] == refe_base[0])
		            n_hitsq++;
    			  if(read_base[k] == SNP2_base)
			    n_hitrq++;
		       }
		       rate2 = n_hitrq;
		       rate2 = rate2/n_hitsq;
//		       if((rate>0.25)||(read_nhit[max_base2]>=2))
//		       if(rate>0.21)
//		       if((rate2 > 0.25)||(read_nhit[max_base2]>=2))
//		         if(rate2 > 0.12)
		       if((rate>0.21)||(rate2 > allele_rate))
		       {
		         if((rate2 > 0.12)||(rate2 > allele_rate))
		         {
                           int snp_score = read_scor[max_base2]/ifactor;
		           if(snp_score>99)
		             snp_score = 99;
                           if(max_base2 == 0)
                             Base2 = 'A';
                           else if(max_base2 == 1)
                             Base2 = 'C';
                           else if(max_base2 == 2)
                             Base2 = 'G';
                           else if(max_base2 == 3)
                             Base2 = 'T';
                           else
                             Base2 = 'N';
			   if(file_flag == 1)
			   {
			     if(SNP2_base==Base2)
                               printf("SNP_hez: %s %-2d %-9d %-7d %c %c/%c %-6d %-6d %-6d %-6d %-6d %-6d %-3d %-3d %-3d %-3d\n",(sub+c)->name,snp_score,b,num_hits,refe_base[0],SNP2_base,refe_base[0],read_nhit[0],read_nhit[1],read_nhit[2],read_nhit[3],read_nhit[4],read_nhit[5],read_lowq[0],read_lowq[1],read_lowq[2],read_lowq[3]);
			     else
                               printf("SNP_hez: %s %-2d %-9d %-7d %c %c/%c %-6d %-6d %-6d %-6d %-3d %-3d %-3d %-3d\n",(sub+c)->name,snp_score,b,num_hits,refe_base[0],SNP2_base,Base2,refe_base[0],read_nhit[0],read_nhit[1],read_nhit[2],read_nhit[3],read_nhit[4],read_nhit[5],read_lowq[0],read_lowq[1],read_lowq[2],read_lowq[3]);
			   }
			   else
			   {
			     if(SNP2_base==Base2)
                               printf("SNP_hez: %s %d %d %d %c %c/%c %s %d %d %s\n",(sub+c)->name,snp_score,b,num_hits,refe_base[0],SNP2_base,refe_base[0],read_base,read_nhit[max_base],read_nhit[max_base2],score_base);
			     else
                               printf("SNP_hez: %s %d %d %d %c %c/%c %s %d %d %s\n",(sub+c)->name,snp_score,b,num_hits,refe_base[0],SNP2_base,Base2,read_base,read_nhit[max_base],read_nhit[max_base2],score_base);
			   }
		         }
	               }
                     }
                   }
                 }
               }
	       if((mapNumber==1)&&(read_nhit[max_base2]==1)&&(max_qsum2>=18))
	       {
                 int snp_score = read_scor[max_base2]/read_nhit[max_base2];
		 if(snp_score>99)
		   snp_score = 99;
                 printf("SNP3: %s %d %d %d %c %c %s %d %d %s\n",(sub+c)->name,snp_score,b,num_hits,refe_base[0],SNP2_base,read_base,read_nhit[max_base],read_nhit[max_base2],score_base);
	       }
             }
	   }
        }
     }
     printf("Mapping Reads: %d Bases: %ld\n",num_mapreads,num_mapbases);

}

/*   Subroutine to find the contig/supercontig structure */
/* ====================================================  */
void Arch_Align(int nRead)
/* ====================================================  */
{
     int i,c;
     fasta *seqp;

        for(i=0;i<nRead;i++)
        {
     	     c = ctg_index[i]; 
	     seqp = seq + i;
	     if(hit_quest[i]==3)
	     {
	       if(hit_rcdex[i]==0)
	       {
	         if(seqp->data[hit_quest[i]-3]==(sub+c)->data[hit_refst[i]-3])
		 {
		   hit_quest[i] = hit_quest[i] - 2;
		   hit_refst[i] = hit_refst[i] - 2;
		   read_arch[i] = 1;
		 }
	       }
	       else
	       {
	         if(seqp->data[hit_quest[i]-3]=='A')
		 {
		   if((sub+c)->data[hit_refed[i]+1]=='T')
		   {
		     hit_quest[i] = hit_quest[i] - 2;
		     hit_refed[i] = hit_refed[i] + 2;
		     read_arch[i] = 1;
		   }
		 }
		 else if(seqp->data[hit_quest[i]-3]=='C')
		 {
		   if((sub+c)->data[hit_refed[i]+1]=='G')
		   {
		     hit_quest[i] = hit_quest[i] - 2;
		     hit_refed[i] = hit_refed[i] + 2;
		     read_arch[i] = 1;
		   }
		 }
		 else if(seqp->data[hit_quest[i]-3]=='G')
		 {
		   if((sub+c)->data[hit_refed[i]+1]=='C')
		   {
		     hit_quest[i] = hit_quest[i] - 2;
		     hit_refed[i] = hit_refed[i] + 2;
		     read_arch[i] = 1;
		   }
		 }
		 else if(seqp->data[hit_quest[i]-3]=='T')
		 {
		   if((sub+c)->data[hit_refed[i]+1]=='A')
		   {
		     hit_quest[i] = hit_quest[i] - 2;
		     hit_refed[i] = hit_refed[i] + 2;
		     read_arch[i] = 1;
		   }
		 }
	       }
	     }
	     else if(hit_quest[i]==4)
	     {
	       if(hit_rcdex[i]==0)
	       {
	         if(seqp->data[hit_quest[i]-4]==(sub+c)->data[hit_refst[i]-4])
		 {
		   hit_quest[i] = hit_quest[i] - 3;
		   hit_refst[i] = hit_refst[i] - 3;
		   read_arch[i] = 2;
		 }
	       }
	       else
	       {
	         if(seqp->data[hit_quest[i]-4]=='A')
		 {
		   if((sub+c)->data[hit_refed[i]+2]=='T')
		   {
		     hit_quest[i] = hit_quest[i] - 3;
		     hit_refed[i] = hit_refed[i] + 3;
		     read_arch[i] = 2;
		   }
		 }
		 else if(seqp->data[hit_quest[i]-4]=='C')
		 {
		   if((sub+c)->data[hit_refed[i]+2]=='G')
		   {
		     hit_quest[i] = hit_quest[i] - 3;
		     hit_refed[i] = hit_refed[i] + 3;
		     read_arch[i] = 2;
		   }
		 }
		 else if(seqp->data[hit_quest[i]-4]=='G')
		 {
		   if((sub+c)->data[hit_refed[i]+2]=='C')
		   {
		     hit_quest[i] = hit_quest[i] - 3;
		     hit_refed[i] = hit_refed[i] + 3;
		     read_arch[i] = 2;
		   }
		 }
		 else if(seqp->data[hit_quest[i]-4]=='T')
		 {
		   if((sub+c)->data[hit_refed[i]+2]=='A')
		   {
		     hit_quest[i] = hit_quest[i] - 3;
		     hit_refed[i] = hit_refed[i] + 3;
		     read_arch[i] = 2;
		   }
		 }
	       }
	     }
	     else if((seqp->length-hit_queed[i])==2)
	     {
	       if(hit_rcdex[i]==0)
	       {
	         if(seqp->data[hit_queed[i]+1]==(sub+c)->data[hit_refed[i]+1])
		 {
		   hit_queed[i] = hit_queed[i] + 2;
		   hit_refed[i] = hit_refed[i] + 2;
		   read_arch[i] = 3;
		 }
	       }
	       else
	       {
	         if(seqp->data[hit_queed[i]+1]=='A')
		 {
		   if((sub+c)->data[hit_refst[i]-3]=='T')
		   {
		     hit_queed[i] = hit_queed[i] + 2;
		     hit_refst[i] = hit_refst[i] - 2;
		     read_arch[i] = 3;
		   }
		 }
		 else if(seqp->data[hit_queed[i]+1]=='C')
		 {
		   if((sub+c)->data[hit_refst[i]-3]=='G')
		   {
		     hit_queed[i] = hit_queed[i] + 2;
		     hit_refst[i] = hit_refst[i] - 2;
		     read_arch[i] = 3;
		   }
		 }
		 else if(seqp->data[hit_queed[i]+1]=='G')
		 {
		   if((sub+c)->data[hit_refst[i]-3]=='C')
		   {
		     hit_queed[i] = hit_queed[i] + 2;
		     hit_refst[i] = hit_refst[i] - 2;
		     read_arch[i] = 3;
		   }
		 }
		 else if(seqp->data[hit_queed[i]+1]=='T')
		 {
		   if((sub+c)->data[hit_refst[i]-3]=='A')
		   {
		     hit_queed[i] = hit_queed[i] + 2;
		     hit_refst[i] = hit_refst[i] - 2;
		     read_arch[i] = 3;
		   }
		 }
	       }
	     }
	     else if((seqp->length-hit_queed[i])==3)
	     {
	       if(hit_rcdex[i]==0)
	       {
	         if(seqp->data[hit_queed[i]+2]==(sub+c)->data[hit_refed[i]+2])
		 {
		   hit_queed[i] = hit_queed[i] + 3;
		   hit_refed[i] = hit_refed[i] + 3;
		   read_arch[i] = 4;
		 }
	       }
	       else
	       {
	         if(seqp->data[hit_queed[i]+2]=='A')
		 {
		   if((sub+c)->data[hit_refst[i]-4]=='T')
		   {
		     hit_queed[i] = hit_queed[i] + 3;
		     hit_refst[i] = hit_refst[i] - 3;
		     read_arch[i] = 4;
		   }
		 }
		 else if(seqp->data[hit_queed[i]+2]=='C')
		 {
		   if((sub+c)->data[hit_refst[i]-4]=='G')
		   {
		     hit_queed[i] = hit_queed[i] + 3;
		     hit_refst[i] = hit_refst[i] - 3;
		     read_arch[i] = 4;
		   }
		 }
		 else if(seqp->data[hit_queed[i]+2]=='G')
		 {
		   if((sub+c)->data[hit_refst[i]-4]=='C')
		   {
		     hit_queed[i] = hit_queed[i] + 3;
		     hit_refst[i] = hit_refst[i] - 3;
		     read_arch[i] = 4;
		   }
		 }
		 else if(seqp->data[hit_queed[i]+2]=='T')
		 {
		   if((sub+c)->data[hit_refst[i]-4]=='A')
		   {
		     hit_queed[i] = hit_queed[i] + 3;
		     hit_refst[i] = hit_refst[i] - 3;
		     read_arch[i] = 4;
		   }
		 }
	       }
	     }
//	    printf("xxx2: %d %s %d %d %d %d\n",i,rdname[i],hit_quest[i],hit_queed[i],hit_refst[i],hit_refed[i]);
        }

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

