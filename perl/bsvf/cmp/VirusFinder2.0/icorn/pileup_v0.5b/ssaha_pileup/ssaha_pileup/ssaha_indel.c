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
static int *map_score,*hit_refst,*hit_refed,*hit_rcdex,*hit_quest,*hit_queed;
static int *readIndex,*ctg_index,*cell_index,*read_mask;
static int *ctg_list,*ctg_cover,*ctg_indel;
static B64_long *ctg_head,sBase;
static char *dataline,*cigar_line;
static int *hit_score,*hit_iddex,*indel_st,*indel_size,*hit_offset;

/* SSAS default parameters   */
static int IMOD=0;
static B64_long line_len=0;
static int num_insreads=0;
static int num_delreads=0;
static int num_hitreads=0;
static int num_reads=0;
static int num_cline=0;
static int set_shift = 10;
static int insert_flag = 0;
static int delete_flag = 0;
static int n_contigs = 0;
static float allele_rate = 0.15;
static char strain_name[100];

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *sub;


int main(int argc, char **argv)
{
    FILE *namef;
    int i,nSeq,args;
    char line[2000]={0};
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    fasta *seq; 
    void ArraySort_Mix(int n, B64_long *arr, int *brr);
    void Read_Pairs(char **argv,int args,int nLib,int nSeq);
    void Align_Process(char **argv,int args,int nRead);
    void Cigar_Filter(char **argv,int args,int nRead);

    seq=NULL;
//    fflush(stdout);
//    system("ps aux | grep eulerSNP; date");
    memset(strain_name,'\0',100);
    if(argc < 2)
    {
      printf("Usage: %s [-allele 0.15] <cigar_line_file> <ref_sequence> <cons_file>\n",argv[0]);
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
       else if(!strcmp(argv[i],"-insertion"))
       {
         sscanf(argv[++i],"%d",&insert_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-deletion"))
       {
         sscanf(argv[++i],"%d",&delete_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-shift"))
       {
         sscanf(argv[++i],"%d",&set_shift);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-allele"))
       {
         sscanf(argv[++i],"%f",&allele_rate);
         args=args+2;
       }
    }

    if((delete_flag==0)&&(insert_flag==0))
    {
      printf("Please give parameter for insertion -insertion 1 or deletion -deletion 1 !\n");
      exit(1);
    }
    if((delete_flag)&&(insert_flag))
    {
      printf("Cannot detection insertion and deletion in one go!\n");
      exit(1);
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
    if((read_mask= (int *)calloc(num_reads,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - read_mask\n");
      exit(1);
    }
    Cigar_Filter(argv,args,num_reads);
    Align_Process(argv,args,num_hitreads);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   Subroutine to process alignment information */
/* ====================================================  */
void Cigar_Filter(char **argv,int args,int nRead)
/* ====================================================  */
{
     int i,j,k,num_align,num_hits;
     char ID,line[2000],*ptr,score[6];
     FILE *namef;

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

/*   read the cigar file         */
     num_align=0;
     while(!feof(namef))
     {
       int nPair=0,len;
       char line2[2000],line3[2000],base[50];
      
       fgets(line,2000,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       if((strncmp(line,"cigar",5))==0)
       { 
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
	 if(nPair>=13)
	 {
           i=0;
           for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
           {
              if(i==12)
              {
                memset(base,'\0',50);
                strcat(base,ptr);
                ID = *ptr;
		if(insert_flag)
		{
	          if(ID == 'I')
		  {
                    read_mask[num_align] = 1;
		    num_insreads++;
		  }
		}
		if(delete_flag)
		{
	          if(ID == 'D')
		  {
                    read_mask[num_align] = 1;
		    num_delreads++;
		  }
		}
              }
           }
	 }
       }
       num_align++;
     }
     fclose(namef);

     if(insert_flag)
     {
       num_hitreads = num_insreads;
       printf("number of insertion reads: %d\n",num_insreads);
     }
     if(delete_flag)
     {
       num_hitreads = num_delreads;
       printf("number of deletion reads: %d\n",num_delreads);
     }

/*   read and process the cons file         */
     if((namef = fopen(argv[args+2],"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }
     num_hits = 0;
     while(!feof(namef))
     {
       fgets(line,2000,namef);
       if(feof(namef)) break;
       num_hits++;
     }
     fclose(namef);
}

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
     char **DBname,*ptr,RC,ID;
     char *line,*st,*ed,*match_text = "M ";
     FILE *fp,*namef;
     B64_long *read_offsets,big_num;
     int n_find,idd,stopflag,num_align,refhit1,refhit2;
     int readhit1 = 0,readhit2 = 0;
     fasta *segg,*seqp;
     B64_long Size_q_pdata;
     int num_seqque,size_range[5];
     char *pdata;
     void Read_Cover(char **argv,int args,int nRead);

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
     RC = '+';
     ID = 'I';
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
     if((hit_score= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - hit_score\n");
       exit(1);
     }
     if((indel_st= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - indel_st\n");
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
     if((hit_iddex= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - hit_iddex\n");
       exit(1);
     }
     if((hit_offset= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - hit_offset\n");
       exit(1);
     }
     if((indel_size= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - indel_size\n");
       exit(1);
     }
     if((read_offsets= (B64_long *)calloc(nRead,sizeof(B64_long))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - read_offsets\n");
       exit(1);
     }

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

/*   read the SNP output file         */
     num_align=0;
     num_reads = 0;
     while(!feof(namef))
     {
       int nPair=0,len;
       char line2[2000],line3[2000],base[500],score[3];
      
       fgets(line,2000,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       strcpy(line3,line);
       if(((strncmp(line,"cigar",5))==0)&&(read_mask[num_reads]==1))
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
              hit_refst[num_align]=refhit1;
              hit_refed[num_align]=refhit2;
              if(RC=='+')
              {
                hit_rcdex[num_align]=0;
                hit_quest[num_align]=readhit1;
                hit_queed[num_align]=readhit2;
              } 
              else
              {
                if((readhit1-readhit2)<0)
                {
                  hit_quest[num_align]=readhit1;
                  hit_queed[num_align]=readhit2;
                } 
                else
                {
                  hit_quest[num_align]=readhit2;
                  hit_queed[num_align]=readhit1;
                } 
                hit_rcdex[num_align]=1;
              }
              readIndex[num_align] = num_align;
              num_align++;
            }
            else if(i==9)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              hit_score[num_align-1] = atoi(ptr);
            }
            else if(i==11)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              hit_offset[num_align-1] = atoi(ptr);
            }
            else if(i==12)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              ID = *ptr;
	      if(ID == 'I')
	        hit_iddex[num_align-1] = 0;
	      else
	        hit_iddex[num_align-1] = 1;
            }
            else if(i==13)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              indel_size[num_align-1] = atoi(ptr);
            }
         }
       }
       num_reads++;
     }
     fclose(namef);

     Read_Cover(argv,args,num_align);
     memset(ctg_index,-1,num_align*4);
/*   sort out reference head  */
     num_cline = 1;
     strcpy(cell_name[0],"READS");

     printf("number of cell lines: %d\n",num_cline);

/*   sort out contig name match */
     if(nRead!= num_align)
     {
       printf("Number of reads not the same. Job stopped!\n");
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
                 ctg_index[readIndex[k]] = idd;
                 n_find++;
               }
            }
          }
        }
        i=j-1;
     }
     printf("number of reads aligned to the genome: %d\n",n_find);

/*   sort out match offsets on contigs/chromosomes */
     for(j=0;j<nRead;j++)
     {
        B64_long nn;

	if(ctg_index[j]>=0)
	  nn = ctg_index[j];
	else
	  nn = 50000L;
        readIndex[j]=j;
	read_offsets[j] = (nn<<40) + hit_refst[j]+hit_offset[j]; 
     }
     ArraySort_Mix(nRead,read_offsets,readIndex);

     size_range[1] = 0;
     size_range[2] = 0;
     size_range[3] = 0;
     big_num = 1L<<30;
     n_find = 0;
     for(i=0;i<nRead;i++)
     {
        stopflag=0;
        j=i+1;
        while((j<nRead)&&(stopflag==0))
        {
//          if(((read_offsets[j]-read_offsets[i])>3)&&((read_offsets[j]-read_offsets[i])<30))
          if((read_offsets[j]-read_offsets[i])<set_shift)
          {
            j++;
          }
          else
            stopflag=1;
        }
	if((j-i)>=2)
        {
	  int idd = (read_offsets[i]>>40);
	  int id_size = indel_size[readIndex[i]];
	  int indel_flag = 0;
          for(k=0;k<(j-i);k++)
	  {
	     int idk = readIndex[i+k];
             int ctg_offset = hit_refst[idk]+hit_offset[idk];
	     int id_ctg = ctg_index[idk]; 
	     B64_long gen_offset = ctg_head[id_ctg]+ctg_offset;
             float rate;

	     rate = ctg_cover[gen_offset];
             if(insert_flag)
	       rate = (j-i)/rate;
	     else
	       rate = ctg_indel[gen_offset]/rate;
	     if(rate>allele_rate)
	       indel_flag = 1;
	  }
	  if(indel_flag)
	  {
	    size_range[id_size]++; 
            for(k=0;k<(j-i);k++)
	    {
	       int idk = readIndex[i+k];
               int ctg_offset = hit_refst[idk]+hit_offset[idk];
               int query_offset;
	       int id_ctg = ctg_index[idk]; 
	       B64_long gen_offset = ctg_head[id_ctg]+ctg_offset;

	       if(hit_rcdex[idk]==0)
	         query_offset = hit_quest[idk] + hit_offset[idk];
	       else
	         query_offset = hit_queed[idk] - hit_offset[idk]+1;
	       seqp = sub + idd;
	       if(insert_flag)
  	         printf("Insertion: %d %s %d %d %d %d %s %d %d %d %d\n",n_find,seqp->name,ctg_offset,indel_size[idk],j-i,map_score[idk],rdname[idk],query_offset,hit_rcdex[idk],ctg_cover[gen_offset],j-i);
	       else 
  	         printf("Deletion: %d %s %d %d %d %d %s %d %d %d %d\n",n_find,seqp->name,ctg_offset,indel_size[idk],j-i,map_score[idk],rdname[idk],query_offset,hit_rcdex[idk],ctg_cover[gen_offset],ctg_indel[gen_offset]);
	    }
	    printf("\n");
	    n_find++;
	  }
        }
        i=j-1;
     }
     printf("Total number of indels:       %d\n",n_find);
     printf("Number of single base indels: %d\n",size_range[1]);
     printf("Number of doublet indels:     %d\n",size_range[2]);
     printf("Number of triplet indels:     %d\n",size_range[3]);
     free(read_offsets);
}

/*   Subroutine to process alignment information */
/* ====================================================  */
void Read_Cover(char **argv,int args,int nRead)
/* ====================================================  */
{
     int i,j,k,num_align,num_hits;
     char ID,line[500],*ptr,score[6];
     FILE *namef;

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

     if((ctg_head= (B64_long *)calloc(n_contigs,sizeof(B64_long))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ctg_head=\n");
       exit(1);
     }
     if((ctg_list= (int *)calloc(n_contigs,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ctg_list=\n");
       exit(1);
     }

     ctg_head[0] = 0;
     for(i=0;i<n_contigs;i++)
     {
        ctg_list[i] = (sub+i)->length;
	if(i>0)
  	  ctg_head[i] = ctg_head[i-1]+ctg_list[i-1];
     }
     
     if((ctg_cover= (int *)calloc(sBase,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ctg_cover\n");
       exit(1);
     }
     if((ctg_indel= (int *)calloc(sBase,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - ctg_indel\n");
       exit(1);
     }
/*   read and process the cons file         */
     if((namef = fopen(argv[args+2],"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

     num_hits=0;
     while(!feof(namef))
     {
       int nPair=0,len,id_ctg,locus;
       char line2[500],line3[500],base[100];
      
       fgets(line,500,namef);
       if(feof(namef)) break;
       strcpy(line2,line);
       strcpy(line3,line);
       if((strncmp(line,"cons:",5))==0)
       { 
         for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
         {
         }
	 if(nPair>=14)
	 {
           i=0;
           for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
           {
	      if(i==0)
	      {
                memset(base,'\0',100);
	        strcat(base,ptr);
		len = strlen(ptr);
		if(len<8)
		{
		  printf("Contig number not found! Please use the latest version of ssaha_pileup\n");
		  exit(1);
		}
		else
	        {
		  memset(score,'\0',6);
                  for(j=6;j<11;j++)
		     score[j-6] = line3[j];
                  id_ctg = atoi(score);
		}
              }
              else if(i==2)
              {
                memset(base,'\0',100);
                strcat(base,ptr);
                locus = atoi(ptr);
              }
              else if(i==3)
              {
                memset(base,'\0',100);
                strcat(base,ptr);
                ctg_cover[ctg_head[id_ctg]+locus] = atoi(ptr);
              }
	      else if(i==9)
	      {
	        memset(base,'\0',100);
	        strcat(base,ptr);
                ctg_indel[ctg_head[id_ctg]+locus] = atoi(ptr);
//		if(ctg_indel[ctg_head[id_ctg]+locus]>3)
//		  printf("offset: %d %d %d\n",id_ctg,locus,ctg_indel[ctg_head[id_ctg]+locus]);
	      }
           }
           num_hits++;
	 }
       }
     }
     fclose(namef);
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

