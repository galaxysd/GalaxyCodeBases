/***********************************************************************\
 *                                                                     * 
 *                         PROJECT   ssaha_pileup                      *
 *                                                                     * 
 *---------------------------------------------------------------------*
 *                             view_pileup                             *
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
#define MAX_N_BRG 50000 
#define MAX_N_ROW 40000 
#define Max_N_NameBase 60
#define Max_N_Pair 100
static char **cell_line;
static char **cell_name;
static char **rdname;
static int *hit_refst,*hit_refed,*hit_rcdex,*hit_quest,*hit_queed;
static int *ctg_index,*cell_index;
static int *contig_body,*contig_list,*cigar_list;
static B64_long *contig_head,*cigar_head,sBase;
static char *dataline,*cigar_line;

/* SSAS default parameters   */
static B64_long line_len=0;
static int num_reads=0;
static int low_qual=1;
static int num_cline=0;
static int LISTENQ=1024;
static int sock_port=777777;
static int pos_input = 200;
static int ctg_input = 0;
static int view_mod = 1;
static int n_contigs = 0;
static int min_length = 20;
static B64_long genome_offset;
static int Size_dataline,Size_sendline;

typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *sub;

static char rc_char[5000];

void error(char *msg)
{
    perror(msg);
    exit(1);
}
       
pid_t Fork(void)
{ 
    pid_t pid;
    if ( (pid = fork()) == -1) printf("fork error\n");
    return(pid);
}       
        
void Close(int fd)
{  
    if(close(fd) == -1)
      printf("close error\n");
}

void sig_chld(int signo)
{
    pid_t pid;
    int   stat;

    while((pid = waitpid(-1,&stat, WNOHANG))>0)
      printf("data writing completed in child %d \n",pid);
    return;
}       
        
ssize_t readn(int filedes, void *buff, size_t nbytes);
ssize_t writen(int filedes, const void *buff, size_t nbytes);
ssize_t readline(int filedes, void *buff, size_t maxlen);
        
ssize_t readn(int fd, void *vptr, size_t n)
/* Read "n" bytes from a descriptor. */
{
    size_t        nleft;
    ssize_t       nread;
    char  *ptr;

    ptr = (char*) vptr;
    nleft = n;
    while (nleft > 0)
    {
      if ( (nread = read(fd, ptr, nleft)) < 0)
      {
        if (errno == EINTR)
          nread = 0;              /* and call read() again */
        else
          return(-1);
      } else if (nread == 0)
        break;                            /* EOF */

      nleft -= nread;
      ptr   += nread;
    }
    return(n - nleft);            /* return >= 0 */
}
/* end readn */

ssize_t writen(int fd, const void *vptr, size_t n)
/* Write "n" bytes to a descriptor. */
{
    size_t        nleft;
    ssize_t       nwritten;
    const char    *ptr;

    ptr = (char *) vptr;
    nleft = n;
    while (nleft > 0)
    {
      if ( (nwritten = write(fd, ptr, nleft)) <= 0) {
        if (errno == EINTR)
          nwritten = 0;           /* and call write() again */
        else
          return(-1);                     /* error */
      }

      nleft -= nwritten;
      ptr   += nwritten;
    }
    return(n);
}
/* end writen */


static ssize_t my_read(int fd, char *ptr)
{
    static int    read_cnt = 0;
    static char   *read_ptr;
    static char   read_buf[MAXLINE];


    if (read_cnt <= 0)
    {
    again:
      if ( (read_cnt = read(fd, read_buf, sizeof(read_buf))) < 0)
      {
        if (errno == EINTR)
          goto again;
        return(-1);
      }
      else if (read_cnt == 0) return(0);
      read_ptr = read_buf;
    }

    read_cnt--;
    *ptr = *read_ptr++;
    return(1);
}

ssize_t readline(int fd, void *vptr, size_t maxlen)
{
    int           n, rc;
    char  c, *ptr;

    ptr = (char *) vptr;
    for (n = 1; n < maxlen; n++)
    {
      if ( (rc = my_read(fd, &c)) == 1)
      {
        *ptr++ = c;
//        if ((c == '\n')||(c=='\0')) // modified to stop at end of string too
        if (c == '\n') // modified to stop at end of string too
          break;  /* newline is stored, like fgets() */
      }
      else if (rc == 0)
      {
        if (n == 1)
          return(0);      /* EOF, no data read */
        else
          break;          /* EOF, some data was read */
      } else
        return(-1);               /* error, errno set by read() */
    }

    *ptr = 0;     /* null terminate like fgets() */
    return(n);
}
/* end readline */

ssize_t Readline(int fd, void *ptr, size_t maxlen)
{
    ssize_t               n;
    if ( (n = readline(fd, ptr, maxlen)) < 0)
      error("readline error");
    return(n);
}

void Writen(int fd, void *ptr, size_t nbytes)
{
  if (writen(fd, ptr, nbytes) != nbytes)
    printf("writen error\n");
}

ssize_t Readn(int fd, void *ptr, size_t nbytes)
{
  ssize_t               n;

  if ( (n = readn(fd, ptr, nbytes)) < 0)
    printf("readn error\n");
  return(n);
}

void str_echo(int sockfd)
{
    int     i,start_flag=0,len,numline;
    ssize_t n;
    char    line[MAXLINE],line2[MAXLINE],*ptr,base[20],zero[20]={0};

    numline=0;
    for(;;)
    {
       if((n = Readline(sockfd,line,MAXLINE)) <=0)
       {
//         printf("read out: %d\n",n);
         return;
       }
       len=strlen(line);
       if(strncmp(line,"COMMAND_INPUT",13)==0)
       {
         i=0;
         strcpy(line2,line);
         for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
         {
            if(i==1)
            {
              strcpy(base,zero);
              strcat(base,ptr);
              low_qual = atoi(base);
            }
            else if(i==2)
            {
              strcpy(base,zero);
              strcat(base,ptr);
              pos_input = atoi(base);
            }
            else if(i==3)
            {
              strcpy(base,zero);
              strcat(base,ptr);
              ctg_input = atoi(base);
            }
            else if(i==4)
            {
              strcpy(base,zero);
              strcat(base,ptr);
              view_mod = atoi(base);
            }
         }
//         printf("line: %s %s\n",line,line2);
         Writen(sockfd,line,n);
         return;
       }
       if(start_flag==0)
       {
//         printf("%s",line);
         if(line[0]=='>')
         {
           strncpy(dataline,line,len);
           Size_dataline+=len;
         }
         else
         {
           strncpy(dataline,">unnamedQuery\n",14);
           Size_dataline=14;
           strncat(dataline,line,len);
           Size_dataline+=len;
         }
         start_flag=1;
       }
       else
       {
//         printf("%s",line);
         strncat(dataline,line,len);
         Size_dataline+=len;
       }
//         printf("%d %s",numline,line);
       numline++;
       Writen(sockfd,line,n);
    }
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
    FILE *namef;
    int i,nSeq,args;
    char line[2000]={0};
    char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
    fasta *seq; 
    void ArraySort_Mix(int n, B64_long *arr, int *brr);
    void SNP_View_Process(char **argv,int args,int nSeq);
    void Read_Pairs(char **argv,int args,int nLib,int nSeq);
    void Align_Process(char **argv,int args,int nRead);

    seq=NULL;
//    fflush(stdout);
//    system("ps aux | grep eulerSNP; date");
    if(argc < 2)
    {
      printf("Usage: %s <cigar_line_file> <ref_sequence> <reads_fastq_file>\n",argv[0]);
      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-len"))
       {
         sscanf(argv[++i],"%d",&min_length);
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
       else if(!strcmp(argv[i],"-port"))
       {
         sscanf(argv[++i],"%d",&sock_port);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-qual"))
       {
         sscanf(argv[++i],"%d",&low_qual);
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
    SNP_View_Process(argv,args,num_reads);
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
     char **DBname,*ptr,line[500];
     FILE *namef,*fp;
     int n_find,*readIndex,idd,stopflag,num_align,refhit1,refhit2;
     int readhit1 = 0,readhit2 = 0;
     fasta *segg;
     B64_long Size_q_pdata;
     int num_seqque;
     char *pdata;

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
     cell_name = cmatrix(0,100,0,Max_N_NameBase);
     cell_line = cmatrix(0,nRead,0,Max_N_NameBase);
     rdname=cmatrix(0,nRead,0,Max_N_NameBase);
     DBname=cmatrix(0,n_reads,0,Max_N_NameBase);

     if((readIndex= (int *)calloc(n_reads,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: Align_Process - readIndex\n");
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

     if((namef = fopen(argv[args],"r")) == NULL)
     {
       printf("ERROR main:: reads group file \n");
       exit(1);
     }

/*   read the SNP output file         */
     num_align=0;
     line_len=0;
     while(!feof(namef))
     {
       int nPair=0,len,c_len;
       char line2[2000],base[500];
      
       fgets(line,2000,namef);
       if(feof(namef)) break;
       c_len = strlen(line) - 1;
       for(j=0;j<c_len;j++)
          cigar_line[cigar_head[num_align]+j] = line[j];
       strcpy(line2,line);
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
              memset(base,'\0',500);
              strcat(base,ptr);
              len = strlen(ptr);
              strncpy(cell_line[num_align],ptr,len-2);
              readIndex[num_align] = num_align;
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
              if((*ptr)=='+')
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
              hit_refst[num_align] = atoi(ptr);
            }
            else if(i==7)
            {
              memset(base,'\0',500);
              strcat(base,ptr);
              hit_refed[num_align] = atoi(ptr);
              num_align++;
            }
         }
     }
     fclose(namef);

     memset(ctg_index,-1,num_align*4);
/*   sort out reference head  */
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
     for(i=0;i<num_cline;i++)
        printf("cell lines: %s\n",cell_name[i]);
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
void SNP_View_Process(char **argv,int args,int nSeq)
/* ====================================================  */
{
     int i,j,k,m;
     B64_long *ref_head,idk = 0L;
     int *ref_list;
     void ArraySort_Mix(int n,B64_long *arr,int *brr);
     B64_long offset,st,ed,kk,tBase,head_size;
     int  **imatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     char **cmatrix(B64_long nrl,B64_long nrh,B64_long ncl,B64_long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
     int NUMLINE=5000;
     FILE *namef,*fp;
     int set_score = 0,err_flag = 0;
     fasta *seq,*seqp;
     int num_traces,read_qual[5000],**m_align;
     char **q_align,**s_align,refe_base[5000],read_base[5000];
     int listenfd, connfd, portno,clilen;
     char r2line[MAXLINE],s2line[MAXLINE];
     struct sockaddr_in serv_addr, cli_addr;
     int num_hits;
     pid_t pid;
     int i_offset,i_cover; 
     fasta *segg;
     B64_long Size_q_pdata;
     int num_seqque;
     char *pdata;

     m_align=imatrix(0,500,0,200);
     if((fp=fopen(argv[args+2],"rb"))==NULL) printf("Cannot open file\n");
       fseek(fp, 0, SEEK_END);
     Size_q_pdata = ftell(fp) + 1;
     fclose(fp);
     if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
       printf("calloc pdata\n");
     num_seqque = extractFastq(argv[args+2],pdata,Size_q_pdata);
     if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL)
       printf("calloc segg\n");
     if((seq=decodeFastq(argv[args+2],&num_seqque,&tBase,pdata,Size_q_pdata,segg))==NULL)
       printf("no query data found.\n");
     num_traces = num_seqque; 
     fastaUC(seq,num_traces);

     if((contig_list= (int *)calloc(sBase,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - contig_list\n");
       exit(1);
     }
     if((contig_head= (B64_long *)calloc(sBase,sizeof(B64_long))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - contig_head\n");
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
        int match_flag = 0;
//        if(abs(hit_refed[i]-hit_refst[i])>min_length)
          match_flag = 1;
        set_score = 0;
        if((match_flag)&&(ctg_index[i]>=0))
        {
          offset = ref_head[ctg_index[i]];
          st = hit_refst[i];
          ed = hit_refed[i];
          idt = cell_index[i];
//     printf("here0: %ld %s %d %d\n",offset,rdname[i],hit_refst[i],ctg_index[i]);
          for(kk=st;kk<=ed;kk++)
          {
             contig_list[offset+kk]++;
          }
        }
     }

     contig_head[0] = 0;
     for(i=1;i<sBase;i++)
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
        int match_flag = 0;
//        if(abs(hit_refed[i]-hit_refst[i])>min_length)
          match_flag = 1;
        set_score = 0;
        if((match_flag)&&(ctg_index[i]>=0))
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
     listenfd = socket(AF_INET, SOCK_STREAM, 0);
     if (listenfd < 0)
        error("ERROR opening socket");
     bzero((char *) &serv_addr, sizeof(serv_addr));
     portno = sock_port;
     serv_addr.sin_family = AF_INET;
     serv_addr.sin_addr.s_addr = INADDR_ANY;
     serv_addr.sin_port = htons(portno);
     if(bind(listenfd, (struct sockaddr *) &serv_addr,sizeof(serv_addr)) < 0)
       error("ERROR on binding");
     listen(listenfd,LISTENQ);
//     Signal(SIGCHLD,sig_chld);
     printf("Start SSAHA2 server\n");
     Size_sendline=NUMLINE*200;
  for(;;)
  {
   clilen = sizeof(cli_addr);
   if((connfd = accept(listenfd, (struct sockaddr *) &cli_addr, &clilen))< 0)
     error("ERROR on accept");
   if( (pid = Fork()) == 0)
   {
     Close(listenfd);
     str_echo(connfd);
     printf("New child process:\n");


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
     if(pos_input>(sub+ctg_input)->length)
       err_flag = 1;
     else if(pos_input<=0)
       err_flag = 1;
     else if(j<=0)
     {
       err_flag = 2;
     }
     else if(ctg_input>=n_contigs)
     {
       err_flag = 3;
     }
     for(j=0;(j<contig_list[genome_offset])&&(err_flag==0);j++)
     {
        int qst,sst,mpos,m_bases,r_bases,s_bases;
        int nPair = 0,r_mod,r_len;
        char *ptr,line2[2000],line[2000];
        idk = contig_head[genome_offset]+j;
        i = contig_body[idk];
        seqp = seq + i;
        seqp->qual[0] = 40; 
//        printf("ssss: %-20s %d %d %ld %d %d %d\n",seqp->name,i,idk,genome_offset,j,pos_input,ctg_input);
        r_len = seqp->length;
        r_len = 2*r_len;
        memset(refe_base,'\0',r_len);
        memset(read_base,'\0',r_len);
        memset(read_qual,0,r_len*4);
        memset(line,'\0',2000);
        for(k=0;k<cigar_list[i];k++)
           line[k] = cigar_line[cigar_head[i]+k];
//        strcpy(line,cigar_line[i]);
//        printf("www: %d %s\n",cigar_list[i],line);
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
//        printf("www: %d %d %d %d || %d %d %d\n",hit_quest[i],hit_queed[i],hit_refst[i],hit_refed[i],r_mod,r_len,m_bases);
        for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
        {
        }
        m=0;
        for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),m++)
        {
           if((m<nPair)&&(m>9))
           {
//        printf("kkkk: %d %s %d\n",m,ptr,nPair);
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
                      read_base[k+m_bases] = seqp->data[k+r_bases]; 
                      read_qual[k+m_bases] = seqp->qual[k+r_bases];
                      refe_base[k+m_bases] = (sub+ctg_input)->data[k+s_bases];
                      if((s_bases+k)==pos_input)
                        mpos = m_bases + k; 
                   }
                   r_bases = r_bases + r_len;
                   s_bases = s_bases + r_len;
                   m_bases = m_bases + r_len;
                 }
                 else if(r_mod==2)
                 {
                   for(k=0;k<r_len;k++)
                   {
                      read_base[k+m_bases] = seqp->data[k+r_bases]; 
                      read_qual[k+m_bases] = seqp->qual[k+r_bases];
                      refe_base[k+m_bases] = '-'; 
                   }
                   r_bases = r_bases + r_len;
                   m_bases = m_bases + r_len;
                 }
                 else if(r_mod==3)
                 {
                   for(k=0;k<r_len;k++)
                   {
                      read_base[k+m_bases] = '-'; 
                      read_qual[k+m_bases] = 0;
                      refe_base[k+m_bases] = (sub+ctg_input)->data[k+s_bases];
                      if((s_bases+k)==pos_input)
                        mpos = m_bases + k; 
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
                      read_base[k+m_bases] = seqp->data[r_bases-k]; 
                      read_qual[k+m_bases] = seqp->qual[r_bases-k];
                      refe_base[k+m_bases] = (sub+ctg_input)->data[k+s_bases];
                      if((s_bases+k)==pos_input)
                        mpos = m_bases + k; 
                   }
                   r_bases = r_bases - r_len;
                   s_bases = s_bases + r_len;
                   m_bases = m_bases + r_len;
                 }
                 else if(r_mod==2)
                 {
                   for(k=0;k<r_len;k++)
                   {
                      read_base[k+m_bases] = seqp->data[r_bases-k]; 
                      read_qual[k+m_bases] = seqp->qual[r_bases-k];
                      refe_base[k+m_bases] = '-'; 
                   }
                   r_bases = r_bases - r_len;
                   m_bases = m_bases + r_len;
                 }
                 else if(r_mod==3)
                 {
                   for(k=0;k<r_len;k++)
                   {
                      read_base[k+m_bases] = '-'; 
                      read_qual[k+m_bases] = 0;
                      refe_base[k+m_bases] = (sub+ctg_input)->data[k+s_bases];
                      if((s_bases+k)==pos_input)
                        mpos = m_bases + k; 
                   }
                   s_bases = s_bases + r_len;
                   m_bases = m_bases + r_len;
                 }
               }
             }
           }
        } 
//        printf("www: %d %d %d\n",mpos,m_bases,i);
        if(mpos == 0)
          mpos = m_bases;
        if(j==0)
        {
          printf("%-20s%10d%10d%10d%31s\n","View position:", ctg_input,pos_input,i,"|");
          printf("%-26s%s\n"," ","==================================================================");
          printf("%-30s ",(sub+ctg_input)->name);
          for(k=(mpos-50);k<=(mpos+50);k++)
          {
             if((k<m_bases)&&(k>=0))
             {
               printf("%c",refe_base[k]);
               s_align[j][k-mpos+50] = refe_base[k];  
             }
             else
             {
               printf("%c",' '); 
               s_align[j][k-mpos+50] = ' ';  
             }
          }
//          strcpy(s_align[j],refe_base);
          printf("\n");
          printf("%-26s%s\n"," ","===================================================================");
        }

        if(low_qual==1)
        {
          for(k=(mpos-50);k<=(mpos+50);k++)
          {
             if((k<m_bases)&&(k>=0))
             {
               if(read_qual[k]<23)
               {
                 read_base[k] = tolower(read_base[k]);
                 q_align[j][k-mpos+50] = read_base[k];  
               }
             }       
          }       
        }
        printf("%-20s%5d ",seqp->name,mpos);
        if(hit_rcdex[i]==0)
        {
          for(k=(mpos-50);k<=(mpos+50);k++)
          {
             if((k<m_bases)&&(k>=0))
             {
               printf("%c",read_base[k]);
               q_align[j][k-mpos+50] = read_base[k];  
             }
//             else if(k<0)
             else
             {
               printf("%c",' '); 
               q_align[j][k-mpos+50] = ' ';  
             }
          }
        }
        else
        {
          for(k=(mpos-50);k<=(mpos+50);k++)
          {
             if((k<m_bases)&&(k>=0))
             {
               if(read_base[k]=='A')
               {
                 printf("%c",'T'); 
                 q_align[j][k-mpos+50] = 'T';  
               }
               else if(read_base[k]=='C') 
               {
                 printf("%c",'G'); 
                 q_align[j][k-mpos+50] = 'G';  
               }
               else if(read_base[k]=='G') 
               {
                 printf("%c",'C'); 
                 q_align[j][k-mpos+50] = 'C';  
               }
               else if(read_base[k]=='T') 
               {
                 printf("%c",'A'); 
                 q_align[j][k-mpos+50] = 'A';  
               }
               else if(read_base[k]=='a')
               {
                 printf("%c",'t'); 
                 q_align[j][k-mpos+50] = 't';  
               }
               else if(read_base[k]=='c') 
               {
                 printf("%c",'g'); 
                 q_align[j][k-mpos+50] = 'g';  
               }
               else if(read_base[k]=='g') 
               {
                 printf("%c",'c'); 
                 q_align[j][k-mpos+50] = 'c';  
               }
               else if(read_base[k]=='t') 
               {
                 printf("%c",'a'); 
                 q_align[j][k-mpos+50] = 'a';  
               }
               else
               {
                 printf("%c",read_base[k]); 
                 q_align[j][k-mpos+50] = read_base[k];  
               }
             }
//             else if(k<0)
             else
             {
               printf("%c",' '); 
               q_align[j][k-mpos+50] = ' ';  
             }
          }
        }
        printf("%15s\n",cell_name[cell_index[i]]);
//        printf("%s\n",read_base);
//        printf("%s\n",refe_base);
     }

     if(err_flag ==0)
     {
       bzero(s2line,MAXLINE);
       sprintf(s2line,"%-20s%10d%10d%10d%31s\n","View position:", ctg_input,pos_input,i,"|");
       Writen(connfd,s2line,strlen(s2line));
       if(Readline(connfd,r2line,MAXLINE)==0)
         error("sever terminated prematurely");

       bzero(s2line,MAXLINE);
       sprintf(s2line,"%-26s%s\n"," ","=====================================================================================================");
       Writen(connfd,s2line,strlen(s2line));
       if(Readline(connfd,r2line,MAXLINE)==0)
         error("sever terminated prematurely");

       bzero(s2line,MAXLINE);
       sprintf(s2line,"%-30s %s\n",(sub+ctg_input)->name,s_align[0]);
       Writen(connfd,s2line,strlen(s2line));
       if(Readline(connfd,r2line,MAXLINE)==0)
         error("sever terminated prematurely");

       bzero(s2line,MAXLINE);
       sprintf(s2line,"%-26s%s\n"," ","=====================================================================================================");
       Writen(connfd,s2line,strlen(s2line));
       if(Readline(connfd,r2line,MAXLINE)==0)
         error("sever terminated prematurely");

       for(j=0;j<num_hits;j++)
       {
          idk = contig_head[genome_offset]+j;
          i = contig_body[idk];
          if(abs(hit_refed[i]-hit_refst[i])>min_length)
          {
            seqp = seq+i;
            bzero(s2line,MAXLINE);
            sprintf(s2line,"%-30s %101s%15s\n",seqp->name,q_align[j],cell_name[cell_index[i]]);
            Writen(connfd,s2line,strlen(s2line));
            if(Readline(connfd,r2line,MAXLINE)==0)
              error("sever terminated prematurely");
          }
       }
       printf("Data write finished! %d\n",num_hits);
     }
     else if(err_flag == 1)
     {
       bzero(s2line,MAXLINE);
       sprintf(s2line,"%s\n","WARNING: POSITION GIVEN IS BEYOND THE CONTIG LENGTH!!!");
       Writen(connfd,s2line,strlen(s2line));
       if(Readline(connfd,r2line,MAXLINE)==0)
         error("sever terminated prematurely");
     }
     else if(err_flag == 2)
     {
       bzero(s2line,MAXLINE);
       sprintf(s2line,"%s\n","WARNING: NO READS WERE PLACED ON THIS REAGION!!!");
       Writen(connfd,s2line,strlen(s2line));
       if(Readline(connfd,r2line,MAXLINE)==0)
         error("sever terminated prematurely");
     }
     else if(err_flag == 3)
     {
       bzero(s2line,MAXLINE);
       sprintf(s2line,"%s\n","WARNING: NO SUCH A CONTIG IN THE REFERENCE!!!");
       Writen(connfd,s2line,strlen(s2line));
       if(Readline(connfd,r2line,MAXLINE)==0)
         error("sever terminated prematurely");
     }
     Close(connfd);
     exit(0); 
   }
   close(connfd);
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

