/*

   contigs.cpp        assemble contigs in BASE.

#    Copyright (C) 2015, The University of Hong Kong.
#
#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 3
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  03110-1301, USA.

    Date   : 28th Jan 2016
    Author : Binghang Liu, bhliu@cs.hku.hk
*/

#include "contigs.h"
#include <sys/resource.h>

struct rusage usage;

FILE* fci;
FILE* fcs;
FILE* fcm;
FILE* fcq;
FILE* logs;

int blank_cycle_time=0;

char* types[]={"LowDepth1", "LowDepth2", "HighDepth", "UsedRead1", "ManyIter", "SameSeed", "UsedRead2", "CB", "R2U", "CurrentUse1", "Cut", "CurrentUse2"};

void get_time_rss(double *utime, double *stime, int* maxrss)
{
    if(getrusage(RUSAGE_SELF,&usage)){ return; }
    *utime = 1e-6 * usage.ru_utime.tv_usec + usage.ru_utime.tv_sec;
    *stime = 1e-6 * usage.ru_stime.tv_usec + usage.ru_stime.tv_sec;
    *maxrss = usage.ru_maxrss;
}

void print_maxrss_and_time(){
    
    if(getrusage(RUSAGE_SELF,&usage)){ return; }
    double utime = 1e-6 * usage.ru_utime.tv_usec + usage.ru_utime.tv_sec;
    double stime = 1e-6 * usage.ru_stime.tv_usec + usage.ru_stime.tv_sec;
    fprintf(stderr, "%s:\t", __FILE__);
    fprintf(stderr, "time\t= %lf (usr)\t+ %lf (sys)\t= %lf (sec)\t",utime , stime, utime + stime);
    fprintf(stderr, "maxrss\t= %d\n", usage.ru_maxrss);
}


unsigned long long get_initial_read(unsigned long long startp, unsigned long long pair_num, int flag, int *c, int size)
{
    unsigned long long i;
    int count= *c, f1;
    for(i=startp; i< pair_num; i++)
    {
        count++;
        f1 = is_used_bwt_sparse(i, size);

        if(f1 == 1 - flag )
            continue;

        if(f1 == flag)
        {
            *c = count;
            return i; 
        }
    }
    *c = count;
    return 0;
}

void print_contig(FILE* fcs, FILE* fcq, Contig *ctg, int id, int ctg_count)
{
    fprintf(fcs, ">UC%d %d %d\n", ctg_count, ctg->len, id);//, ctg->seq);
    int i;
    for(i=0; i< ctg->len; i++)
    {
        fprintf(fcs, "%c", ctg->seq[i]);
        if((i+1) %100 == 0)
        {
            fprintf(fcs, "\n");
        }
    }
    if(i%100 != 0)
        fprintf(fcs, "\n");
}

void complemant_seq(char* seq, int len)
{
    int i;
    for(i=0; i<len ; i++)
    {
        seq[i] = complementMap[ seq[i] ];
    }
}

void reverse_seq(char* seq, int len)
{
    char s[len];
    int i,j;
    for(i=0; i<len ; i++)
    {
        s[i] = complementMap[ seq[len-i-1]];
    }
    for(i=0; i<len ; i++)
    {
        seq[i] = s[i];
    }
}

void reverse_contigs(Contig *ctg, int threadId)
{
    char seq[ctg->len];
    int i,j;
    int num = tmpCount[threadId];

    for(i=0; i<ctg->len ; i++)
    {
        seq[i] = complementMap[ ctg->seq[ctg->len-i-1]];
    }

    for(i=0; i<ctg->len ; i++)
    {
        ctg->seq[i] = seq[i];
    }

    if(num == 0)
       return ;

    j=0;
    for(i=0; i< num; i++)
    {
        if(tmpRead[threadId][i].pos <= Read_length )
            continue;
        tmpRead[threadId][j] = tmpRead[threadId][i];
        tmpRead[threadId][j].pos = int(ctg->len) + Read_length - tmpRead[threadId][i].pos;
        tmpRead[threadId][j].strand = 1 - tmpRead[threadId][i].strand;
        j++;
    }
    num = j;

    TmpRead tmp;
    for(j=0; j<= num/2; j++)
    {
        if(j >= num-1-j)
            break;
        tmp = tmpRead[threadId][j];
        tmpRead[threadId][j] = tmpRead[threadId][num-1-j];
        tmpRead[threadId][num-1-j] = tmp;
    }
    tmpCount[threadId] = num;
}

void update_readPos(TmpRead* tmpRead, int initial, int final)
{
    int i;
    for(i=initial; i< final; i++)
    {
        tmpRead[i].pos -= Read_length;
        tmpRead[i].strand = 1 - tmpRead[i].strand;
    }
}

void print_read(FILE* fp, TmpRead* tmpRead, int startp, int count)
{
    int i;
    fprintf(fp, ">C%d %d\n", startp, count);
    for(i=0; i< count; i++)
    {
        fprintf(fp, "%d %d %d %d %d %d\n", i, tmpRead[i].id, tmpRead[i].strand, tmpRead[i].pos, tmpRead[i].iter, tmpRead[i].cons);
    }
}

void print_g_iters(FILE* fp, Iter* iters, int startp, int count)
{
    int i;
    //fprintf(fp, ">C%d %d\n", startp, count);
    for(i=0; i< count; i++)
    {
        fprintf(fp, "SEED%d\t%s\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\n", i, iters[i].seq, double(Read_length*iters[i].depth)/double(Read_length-iters[i].slen+1), iters[i].slen, iters[i].depth, iters[i].len, iters[i].rnum, iters[i].end, iters[i].pe);
    }
}

void print_single_iter(FILE* fp, Iter* iters, int i)
{
    fprintf(fp, "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", i, iters[i].seq, iters[i].slen, iters[i].depth, iters[i].len, iters[i].rnum, iters[i].end, iters[i].pe);
}

void reverse_seq(char* seq)
{
    int i, len = strlen(seq);
    char rc[len];
    for(i=0; i< len; i++)
    {
        rc[i] = complementMap[ seq[len-i-1]];
    }
    for(i=0; i< len; i++)
        seq[i] = rc[i];

}

int is_cons(TmpRead* tmpRead, int count, int startp)
{
    int i,j;
    for(i=0; i< count; i++)
    {
        if(tmpRead[i].id == startp && tmpRead[i].keep == 1)
            return 1;
    }
    return 0;
}

void print_read_seq(FILE* fp, TmpRead* tmpRead, int startp, int count)
{
    fprintf(fp, ">C%d %d\n", startp, count);
    int i,j;
    char seq[101];
    int qual[101];
    for(i=0; i< count; i++)  
    {
        get_id_seq(seq, qual, tmpRead[i].id, tmpRead[i].size);
        if(tmpRead[i].strand == 0)
	{
            reverse_seq(seq);
        }
        for(j=0; j< tmpRead[i].pos; j++)
        {
            fprintf(fp, " ");
        }
        if(tmpRead[i].pos >= 0)
            fprintf(fp, "%s %d %d %d %d %d\n", seq, tmpRead[i].id, tmpRead[i].strand, tmpRead[i].pos, tmpRead[i].iter, tmpRead[i].unique);
        else
            fprintf(fp, "%s %d %d %d %d %d\tCUT:%d\n", seq, tmpRead[i].id, tmpRead[i].strand, tmpRead[i].pos, tmpRead[i].iter, tmpRead[i].unique, - int(tmpRead[i].pos));
    }
}

void clean_uniq(int threadId)
{
    int i;
    for(i=0; i< ctgs[threadId][ ctg_num[threadId] ].len; i++)
        ctgs[threadId][ ctg_num[threadId] ].uniq[i] = '0';
}

void uniqer_single_seq(int threadId)
{
    int end = ctgs[threadId][ ctg_num[threadId] ].len - 1, start, flag;
    unsigned long long depth;
    int step = depth_6_length > MINOVERLAP ? depth_6_length - MINOVERLAP : MINOVERLAP;
    int break_start=0;
    while(end > start + MINOVERLAP)
    {
        depth = 0;
        flag = extend_backwards(ctgs[threadId][ ctg_num[threadId] ].seq, 0, end, tmpbase[threadId], &(depth), &break_start);
        if(flag != -1 && end-flag < 79)
        {
             ctgs[threadId][ ctg_num[threadId] ].uniq[flag] = end - flag + '0';
             end --;
        }else
             end -= step;
    }
}

void uniqer_single_seq_initial(int threadId)
{
    int end = ctgs[threadId][ ctg_num[threadId] ].len - 1, start, flag, break_start;
    unsigned long long depth;
    while(end > start + MINOVERLAP)
    {
        depth = 0;
        flag = extend_backwards(ctgs[threadId][ ctg_num[threadId] ].seq, 0, end, tmpbase[threadId], &(depth), &break_start);
        if(flag != -1 )
        {
             ctgs[threadId][ ctg_num[threadId] ].uniq[flag] = end - flag + '0';
        }
        end --;
    }
}


int is_rich(char* u, int s, int len, int flag)
{
    int i, zero=0;
    int check_win, cret_win;
    if(flag == 1)
    {
        check_win = MINOVERLAP/2;
        cret_win = MINOVERLAP/4;
    }else{
        check_win = 5;
        cret_win = 2;
    }
    //at least 5 28bp.
    for(int i=s+1; i< len  && i< s + check_win; i++)
    {
        if(flag == 1 && u[i] == '0')
            zero++;
        if(flag == 0 && u[i] != '0')
            zero++;
    }
    if(zero > cret_win)
        return 0;
    return 1;
}

int get_next(char* u, int s, int len)
{
    int i;
    for(i=s+1; i<len-1; i++)
    {
        if(u[i] - '0' == MINOVERLAP-1 && is_rich(u, i, len, 1) == 1)
            return i;
    }
    return -1;
}

int get_next_zero(char* u, int s, int len)
{
    int i;
    for(i=s+1; i<len-1; i++)
    {
        if(u[i] == '0' && is_rich(u, i, len, 0) == 1)
            return i;
    }
    return -1;
}

int get_next_one(char* u, int s, int len)
{
    int i;
    for(i=s+1; i<len-1; i++)
    {
        if(u[i] != '0' && is_rich(u, i, len, 1) == 1)
            return i;
    }
    return -1;
}

void remove_island(char* u, int len)
{
    int i, next_i, j, min_island_size = MINOVERLAP;
    for(i=0; i< len ;i++)
    {
        if(u[i] != '0')
        {
            next_i = get_next_zero(u, i, len);
            if(next_i == -1)
            {
                for(j = i; j< len ; j++)
                {
                    if(u[j] == '0')
                        u[j] = '0' + 1;
                }
                i = len;
                break;
            }

            if(next_i - i < min_island_size)
            {
                for(j=i; j< next_i; j++)
                {
                    if(u[j] == '0')
                        continue;
                    u[j] = '0';
                }
            }else{
                for(j=i+1; j< next_i; j++)
                {
                    if(u[j] == '0')
                        u[j] = '0' + 1;
                }
            }
            i = next_i;
        }
    }
}

void keep_repeat_end(char* u, int len)
{
   int i, j, next_i, keep_len = MINOVERLAP;
   for(i=0; i< len ;i++)
   {
       if(u[i] == '0')
       {
           next_i = get_next_one(u, i, len);
           if(next_i == -1)
           {
               for(j = i; j< i+keep_len && j<len ; j++)
               {
                   if(u[j] == '0')
                       u[j] = '0' + 1;
               }
               i=len;
               break;
           }

           if(next_i - i < Read_length -10)
           {
               for(j=i; j< next_i; j++)
               {
                   if(u[j] == '0')
                       u[j] = '0' + 1;
               }
           }else{
               int tmp_keep_len = (next_i - i - 10)/2;
               if(tmp_keep_len > keep_len)
                   tmp_keep_len = keep_len;

               for(j=i; i>keep_len && j< i + tmp_keep_len; j++)
               {
                   if(u[j] == '0')
                       u[j] = '0' + 1;
               }
               for(j=next_i - 1; j>= next_i - tmp_keep_len; j--)
               {
                   if(u[j] == '0')
                       u[j] = '0' + 1;
               }
               for(; j>= i+tmp_keep_len || (i==0 && j> 0); j--)
               {
                   if(u[j] != '0')
                       u[j] = '0';
               }
           }
           i = next_i;
       }
   }
}

void uniqer2region(char* u, int len)
{
    int i, zero = MINOVERLAP, next_i, now_i, j;
    for(i=0; i< len ;i++)
    {
        if(u[i] - '0' == MINOVERLAP-1)
        {
            zero = 0;
        }else{
            next_i = get_next(u, i, len);
            if(next_i == -1)
            {
                for(j = i; j< len ; j++)
                    u[j] = '0';
                i = len;
                break;
            }

            for(j=i+1; j< next_i-1; j++)
            {
                if(u[j] == '0')
                    continue;
                u[j] = '0';
            }

            now_i = next_i;
            char c;
            for(j=now_i-1; j>=i; j--)
            {
                if(now_i - j > Read_length - 1 - (u[now_i] - '0'))
                    break;
                 c = u[now_i] + (now_i - j);
                 if(u[j] > c || u[j] == '0')
                     u[j] = c;
                 else
                     now_i = j;
            }

            i = next_i; 
        }
    }
}

void print_single_uniq(FILE* fo, char* u, char* id, int len)
{
    fprintf(fo, ">%s %d\n", id, len);
    int i, zero = MINOVERLAP, next_i, now_i, j;

    for(i=0; i< len ;i++)
    {
        if(u[i] - '0' == MINOVERLAP-1)
        {
            zero=0;
            fprintf(fo, "%c", u[i]);
            if((i+1)%100 == 0)
                fprintf(fo, "\n");
        }else
        {
            next_i = get_next(u, i, len);
            if(next_i == -1)
            {
                for(j = i; j< len ; j++)
                {
                    fprintf(fo, "%c", '0');
                    if((j+1)%100 == 0)
                        fprintf(fo, "\n");
                }
                i = len;
                break;
            }

            for(j=i+1; j< next_i-1; j++)
            {
                if(u[j] == '0')
                    continue;
                u[j] = '0';
            }

            now_i = next_i;
            char c;
            for(j=now_i-1; j>=i; j--)
            {
                if(now_i - j > Read_length - 1 - (u[now_i] - '0'))
                {
                    break;
                }
                c = u[now_i] + now_i - j;
                if(u[j] > c || u[j] == '0')
                {
                    u[j] = c;
                }else{
                    now_i = j;
                }
            }

            for(j=i; j<= next_i; j++)
            {
                if(u[j] > '0')
                    fprintf(fo, "%c", u[j]);
                else
                    fprintf(fo, "%c", '0');
                if((j+1)%100 == 0)
                    fprintf(fo, "\n");
            }
            i = next_i;
        }
    }
    if(i%100 != 0)
        fprintf(fo, "\n");
}

void print_ctgs(FILE* fcs, FILE* fcq, FILE* fcm, FILE* flog)
{
    int i, j, k;
    for(i=0; i< threadNum; i++)
    {
        for(j=0; j< ctg_num[i]; j++)
        {
            fprintf(flog, "UC%d %d %d %d %s,%s %d,%d %d,%d %d\n", g_count, ctgs[i][j].id, ctgs[i][j].s_len, ctgs[i][j].s_depth, types[ctgs[i][j].back_end], types[ctgs[i][j].forw_end], 
                                                                  ctgs[i][j].back_iter, ctgs[i][j].forw_iter, ctgs[i][j].back_rc, ctgs[i][j].forw_rc, ctgs[i][j].len);
            fprintf(fcs, ">UC%d %d %lld\n", g_count, ctgs[i][j].len, ctgs[i][j].id);
            fprintf(fcm, ">UC%d\n", g_count);
            g_count++;
            for(k=0; k< ctgs[i][j].len; k++)
            {
                fprintf(fcs, "%c", ctgs[i][j].seq[k]);
                if(ctgs[i][j].uniq[k] != '0')
                    fprintf(fcm, "%c", ctgs[i][j].seq[k]);
                else
                    fprintf(fcm, "N");
                if((k+1) %100 == 0)
                {
                    fprintf(fcs, "\n");
                    fprintf(fcm, "\n");
                }
            }
            if(k%100 != 0)
            {
                fprintf(fcs, "\n");
                fprintf(fcm, "\n");
            }
            ctgs[i][j].len = 0;
            ctgs[i][j].seq[0] = '\0';
        }
        ctg_num[i] = 0;
    }
}

void print_ctg(FILE* fcs, FILE* fcq, FILE* flog)
{
    int i, j, k;
    for(i=0; i< threadNum; i++)
    {
        for(j=0; j< ctg_num[i]; j++)
        {
            fprintf(flog, "UC%d %d %d %d %s,%s %d,%d %d,%d %d\n", g_count, ctgs[i][j].id, ctgs[i][j].s_len, ctgs[i][j].s_depth, types[ctgs[i][j].back_end], types[ctgs[i][j].forw_end], 
                                                                  ctgs[i][j].back_iter, ctgs[i][j].forw_iter, ctgs[i][j].back_rc, ctgs[i][j].forw_rc, ctgs[i][j].len);
            fprintf(fcs, ">UC%d %d %lld\n", g_count, ctgs[i][j].len, ctgs[i][j].id);
            g_count++;
            for(k=0; k< ctgs[i][j].len; k++)
            {
                fprintf(fcs, "%c", ctgs[i][j].seq[k]);
                if((k+1) %100 == 0)
                {
                    fprintf(fcs, "\n");
                }
            }
            if(k%100 != 0)
                fprintf(fcs, "\n");
            ctgs[i][j].len = 0;
            ctgs[i][j].seq[0] = '\0';
        }
        ctg_num[i] = 0;
    }
}

unsigned long long get_next_bwt_id(unsigned long long startp, int size, unsigned long long *c)
{
    unsigned long long i;
    int count= *c, f1, f2, cons=1;

    if(startp % 2== 0)
    {
        count++;
        if(is_used_bwt_sparse(startp+1, size) == 0)
            return startp+1;
    }else{
        startp --;
    }
    
    for(i=startp + 2* gthreadNum; i< pair_num; i+= 2*gthreadNum)
    {
        count++;
        if(is_used_bwt_sparse(i, size) == 0)
        {
            *c = count;
            return i; 
        }
        count++;

        if(is_used_bwt_sparse(i+1, size) == 0)
        {
            *c = count;
            return i+1; 
        }
    }
    *c = count;
    return 0;
}

int get_next_in_conf(Id* conf, int conf_num, unsigned long long* id, int *size, int *count)
{
    int i=*count+1, f1=0, cons=0;
    for(; i< conf_num; i++)
    {
        if(conf[i].c == 0)
            continue;
        if(is_used_bwt_sparse(conf[i].c, *size) == 0 )
        {
            *count = i;
            *id = conf[i].c;
            conf[i].c = 0;
            return 1; 
        }
        conf[i].c = 0;
    }
    *count = i;
    return 0;
}

void slim_conf(Id* conf, int* conf_num)
{
    int count = *conf_num;
    int i, j;
    for(i=0, j=0; i< count; i++)
    {
        if(conf[i].c == 0 || conf[j-1].c == conf[i].c)
            continue;
        conf[j].c = conf[i].c;
        conf[i].c = 0;
        j++;
    }
    *conf_num = j;
}

void add_conflict_buffer(unsigned long long start_id, int size, int threadId, int count)
{
    if(conf_size[threadId] == buffer_size)
    {
        //if(count < 0)
        //{
        //    fprintf(stderr, "buffer Error: thread %d\n", threadId);
        //    exit(1);
        //}
        conf_buffer[threadId][count+1].c = (ULL)start_id;
        conf_buffer[threadId][count+1].size = (ULL)size;
    }else{
        conf_buffer[threadId][ conf_size[threadId] ].c = (ULL)start_id;
        conf_buffer[threadId][ conf_size[threadId] ].size = (ULL)size;
        conf_size[threadId]++;
    }
}

void get_next_id(int threadId, unsigned long long* id, int* size, int *count)
{
    unsigned long long start_id;
    int start_size, prev_count;
    if(*id != 0)
    {
        start_id = *id;
        start_size = *size;
    }else{
        start_id = current_id[threadId].c;
        start_size = current_id[threadId].size;
    }
    //find in conf_buffer.
    if(conf_size[threadId] >= buffer_size - 1)
    {
        prev_count = *count + 1;
        int find_flag = get_next_in_conf(conf_buffer[threadId], conf_size[threadId], &start_id, &start_size, count);
        endInf[threadId].c += *count - prev_count;
        if(*count == conf_size[threadId])
        {
            slim_conf(conf_buffer[threadId], &(conf_size[threadId]));
            *count = -1;
        }
        if(find_flag == 1)
        {
            *id = start_id;
            *size = start_size;
            return ;
        }
        start_id = current_id[threadId].c;
        start_size = current_id[threadId].size;
    }
    
    //find in bwt.
    *id = get_next_bwt_id(start_id, start_size, &(endInf[threadId].c));
    *size = start_size;

    if(*id == 0 && conf_size[threadId] > 0)
    {
        prev_count = *count + 1;
        int find_flag = get_next_in_conf(conf_buffer[threadId], conf_size[threadId], &start_id, &start_size, count);
        endInf[threadId].c += *count - prev_count;

        if(*count == conf_size[threadId])
        {
            slim_conf(conf_buffer[threadId], &(conf_size[threadId]));
            *count = -1;
        }

        if(find_flag == 1)
        {
            *id = start_id;
            *size = start_size;
            return ;
        }
        
    }
    current_id[threadId].c = *id;
    current_id[threadId].size = *size;
    
}

int layout_single_contig_thread(unsigned long long startp, int size, int threadId)
{
    clock_t start_time, now_time;

if(DEBUG)
{  
    int id = 20900; //4302;//64416; //1520; 
    if(startp > id)
    {
        print_ctg(fcs, fcq, logs);
        exit(4);
    }

    if(startp != id)
        return 0;
}
    ctgs[threadId][ctg_num[threadId]].id = startp;
    initial_contig_seed(stdout, startp, size, threadId);

if(DEBUG)
{
    fprintf(stdout, "Starting from %d", startp);
    fprintf(stdout, " %d,%d\n", g_iters[threadId][0].slen, g_iters[threadId][0].depth);
}
    g_iters[threadId][0].rnum =0; g_iters[threadId][0].end = 0; g_iters[threadId][0].pe=0; g_iters[threadId][0].len=0;

    if(g_iters[threadId][0].depth == 0)
    {
        ctgs[threadId][ctg_num[threadId]].len=0; ctgs[threadId][ctg_num[threadId]].seq[0]='\0';
        return 0;
    }
    double utime0=0, utime1=0, stime0=0, stime1=0;
    int maxrss0=0, maxrss1=0;

if(DEBUG)
{
    fprintf(stderr, "For contig starting from %lld\n", startp);
    get_time_rss(&utime0, &stime0, &maxrss0);
}
    tmpCount[threadId] = 0;
    iterCount[threadId] = 0;

    four_iter[threadId] = 0;
    four_iter_checked_read[threadId] = 0;
    four_iter_used_read[threadId] = 0;

    ctgs[threadId][ctg_num[threadId]].s_len = g_iters[threadId][0].slen;
    ctgs[threadId][ctg_num[threadId]].s_depth = g_iters[threadId][0].depth;

    int flag = backward_extending_bwtpe(stdout, 1, threadId);
if(DEBUG)
{
    get_time_rss(&utime1, &stime1, &maxrss1);
    fprintf(stderr, "utime diff %.4f, stime diff %.4f, maxrss diff %d\n", utime1 - utime0, stime1 - stime0, maxrss1 - maxrss0);
}
    if(flag == 1)
    {
        ctgs[threadId][ctg_num[threadId]].len=0; ctgs[threadId][ctg_num[threadId]].seq[0]='\0';
        return 1;
    }
    int plus_iter = iterCount[threadId] - 1;
    ctgs[threadId][ctg_num[threadId]].back_end = g_iters[threadId][ plus_iter ].end;
    ctgs[threadId][ctg_num[threadId]].back_iter = iterCount[threadId];

    reverse_contigs(&(ctgs[threadId][ctg_num[threadId]]), threadId);
    ctgs[threadId][ctg_num[threadId]].back_rc = tmpCount[threadId];

if(DEBUG)    fprintf(stdout, "\nReverse\n");

    reverse_com(g_iters[threadId][0].seq, g_iters[threadId][ iterCount[threadId] ].seq, g_iters[threadId][0].slen); 
    unsigned long long depth =0;

    set_SAranges(g_iters[threadId][ iterCount[threadId] ].seq, 0, g_iters[threadId][0].slen-1, g_iters[threadId][ iterCount[threadId] ].b, &(depth));
    g_iters[threadId][ iterCount[threadId] ].slen = g_iters[threadId][0].slen;
   
    if(depth < MAXDEPTH) 
        g_iters[threadId][ iterCount[threadId] ].depth = depth;
    else
        g_iters[threadId][ iterCount[threadId] ].depth = MAXDEPTH;

    g_iters[threadId][ iterCount[threadId] ].rnum =0; g_iters[threadId][ iterCount[threadId] ].end = 0; g_iters[threadId][ iterCount[threadId] ].pe=0; g_iters[threadId][ iterCount[threadId] ].len=0;

    flag = backward_extending_bwtpe(stdout, 1, threadId);
if(DEBUG)
{
    get_time_rss(&utime1, &stime1, &maxrss1);
    fprintf(stderr, "utime diff %.4f, stime diff %.4f, maxrss diff %d\n", utime1 - utime0, stime1 - stime0, maxrss1 - maxrss0);
}
    if(flag == 1)
    {
        ctgs[threadId][ctg_num[threadId]].len=0; ctgs[threadId][ctg_num[threadId]].seq[0]='\0';
        return 1;
    }
    ctgs[threadId][ctg_num[threadId]].forw_end = g_iters[threadId][ iterCount[threadId] - 1 ].end;
    ctgs[threadId][ctg_num[threadId]].forw_iter = iterCount[threadId] - plus_iter;
    ctgs[threadId][ctg_num[threadId]].forw_rc = tmpCount[threadId] - ctgs[threadId][ctg_num[threadId]].back_rc;
    
    if(tmpCount[threadId] > 0)
        update_readPos(tmpRead[threadId], 0, tmpCount[threadId]);

    /*
    if(tmpCount[threadId] > 0)
    {
         print_read(fci, tmpRead[threadId], startp, tmpCount[threadId]);
    }
    */

if(DEBUG)
{    
    if(iterCount[threadId] > 0)
    {
	    fprintf(stdout, ">C%d\n", startp);
        //print_single_iter(stdout, g_iters, plus_iter); print_single_iter(stdout, g_iters, iter-1);
        print_g_iters(stdout, g_iters[threadId], startp, iterCount[threadId]);
    }
}   
    //fflush(stdout);    
    if(ctgs[threadId][ctg_num[threadId]].len > 0)
        complemant_seq(ctgs[threadId][ctg_num[threadId]].seq, ctgs[threadId][ctg_num[threadId]].len);

    if(ctgs[threadId][ctg_num[threadId]].len >=Read_length && ctgs[threadId][ctg_num[threadId]].len <= 1.5 * Read_length)
    {
        if(is_cons(tmpRead[threadId], tmpCount[threadId], startp) == 0 
		|| (g_iters[threadId][ iterCount[threadId] - 1 ].end == 6 && g_iters[threadId][ plus_iter ].end == 6)
		|| (ctgs[threadId][ctg_num[threadId]].len <= 1.1* Read_length && (g_iters[threadId][ iterCount[threadId] - 1 ].end == 6 || g_iters[threadId][ plus_iter ].end == 6))	
	  )
        {
            ctgs[threadId][ctg_num[threadId]].len = 0; ctgs[threadId][ctg_num[threadId]].seq[0] = '\0';
            clean_trimed_reads(ctgs[threadId][ctg_num[threadId]].len, tmpCount[threadId], threadId);
        }
    }

    if(RepeatMask == 1 && ctgs[threadId][ctg_num[threadId]].len >= Read_length)
    {
        clean_uniq(threadId);
        uniqer_single_seq_initial(threadId);
        remove_island(ctgs[threadId][ctg_num[threadId]].uniq, ctgs[threadId][ctg_num[threadId]].len);
    }

if(DEBUG) fprintf(stdout, "Finished, length=%d\n", ctgs[threadId][ctg_num[threadId]].len);
    free_thread_tmpRead(threadId, 0);

if(DEBUG)
{
    get_time_rss(&utime1, &stime1, &maxrss1);
    fprintf(stderr, "utime diff %.4f, stime diff %.4f, maxrss diff %d\n", utime1 - utime0, stime1 - stime0, maxrss1 - maxrss0);
}
    return 0;
}

void extending_thread(int threadId)
{
    if(threadEnd[threadId] == 1)
        return ;
    unsigned long long start_id=0;
    clock_t start_time, now_time;
    int size=0, count=-1;
    get_next_id(threadId, &start_id, &size, &count);

    while(start_id > 0 )
    {
        ctgs[threadId][ ctg_num[threadId] ].len = 0;

        int has_conflict = layout_single_contig_thread(start_id, size, threadId);
        
        endInf[threadId].uc ++;
        endInf[threadId].len += ctgs[threadId][ ctg_num[threadId] ].len;

        if(ctgs[threadId][ ctg_num[threadId] ].len > 2* Read_length)
            endInf[threadId].c = 0;
        else if(has_conflict == 0){
            endInf[threadId].c ++;

            if(endInf[threadId].c > count_pair_thread )//|| (endInf[threadId].uc > count_pair_thread && endInf[threadId].len/endInf[threadId].uc <= MINOVERLAP))
            {
                break;
            }
        }
        if(gthreadNum == 1 && start_id > count_pair)
            break;
        if(has_conflict == 1)
        {
            add_conflict_buffer(start_id, size, threadId, count);
        }else if(ctgs[threadId][ ctg_num[threadId] ].len >= Read_length)
        {
            endInf[threadId].ulen += ULL(ctgs[threadId][ ctg_num[threadId] ].len);
            ctg_num[threadId] ++;
            if(ctg_num[threadId] >= buffer_size)
                break;
        }
        get_next_id(threadId, &start_id, &size, &count);
    }
    if(start_id == 0)
    {
        threadEnd[threadId] = 1;
        fprintf(stderr, "Thread: %d finished.\n", threadId);
    }else{
        fprintf(stderr, "Thread%d buffer full\n", threadId);
    }
}

void* threadRoutine(void* threadId_p)
{
    int threadId = *((int*) threadId_p);
    //check the status of each thread.
    while (1)
    {
        if (signal[threadId] == 0)
	    {
	        usleep(1);
	        continue;
	    }else if (signal[threadId] == 1)
	    {
            extending_thread(threadId);
            signal[threadId] = 0;
	    }else
  	        break; //unkown signal.
    }
    
}

static void creatThrds ( pthread_t * threads , int* threadId)
{
    unsigned char i;
    int temp;
    for(i=0; i< threadNum; i++)
    {
        threadId[i]=i;
        signal[i] = 0;
    }
    for ( i = 0; i < threadNum; i++ )
    {
        if ( ( temp = pthread_create ( threads + i, NULL, threadRoutine, (void*)(threadId + i) ) ) != 0 )
        {
            fprintf ( stderr, "Create threads failed.\n" );
            exit ( 1 );
        }
    }

    fprintf ( stderr, "%d thread(s) initialized.\n", threadNum );
}

void sendSignal(int SIG)
{
    int i;
    for (i=0; i<threadNum; i++)
    {
        signal[i] = SIG;
    }

    //checking the signal of each thread.
    if(SIG == 1)
    {
    while (1)
    {
        usleep(1);

	for (i=0; i<threadNum; i++)
	{
	    if (signal[i] != 0)
	        break;
        }
        
	if (i == threadNum)
	    break;
    }
    }else if(SIG == 2)
    {
        for(i=0; i< threadNum; i++)
        {
            pthread_join ( threads[i], NULL );
        }
    }
    //after return, all the singal=0.
}

void initial_current(int size)
{
    int i, j;
    for(i=0, j=0; i< threadNum; i++, j+=2)
    {
        current_id[i].c = j;
        current_id[i].size = size;
        threadEnd[i] = 0;
    }
 }

int is_all_end()
{
    int i, all_end = 1;
    for(i=0; i< threadNum; i++)
    {
        if(threadEnd[i] == 0)
        {
            all_end = 0;
            break;
        }
    }
    return all_end;
}
int len_dif_zero;
int is_terminate(int to_single)
{
    int i;
    unsigned long long total_c=0, total_used=0, total_len=0, total_usedc=0, total_usedlen = 0;

    for(i=0; i< threadNum; i++)
    {
        total_c += endInf[i].c;
        total_used += endInf[i].ur;
        total_len += endInf[i].len;
        total_usedc += endInf[i].uc;
        total_usedlen += endInf[i].ulen;
        if(gthreadNum != 1)
        {
            endInf[i].len = 0; endInf[i].uc = 0;
        }
    }
    
    double current_read_ratio = double(total_used)/double(pair_num)*100.0;
    if(total_usedc == 0)
        return 1;
    int avg_len = total_len/total_usedc;
    double len_dif = 1;
    double read_ratio_dif = current_read_ratio - read_ratio;

    if(ass_len > 0)
        len_dif = double(total_usedlen - ass_len)/double(ass_len);

    if(gthreadNum == 1)
    {
        if(current_id[0].c < count_pair)
        {
            fprintf(stderr, "curr id %lld and count pair %lld\n", current_id[0].c, count_pair);
            blank_cycle_time ++;
            if(blank_cycle_time < 10 || read_ratio_dif > 0.0001)
                return 0;
        }
        blank_cycle_time = 0;
        if(avg_len >= old_avg_len && avg_len <= MINOVERLAP)
            is_platform ++;
        else if(avg_len < old_avg_len)
            is_platform = 0;

	if(len_dif < MINLENDIF)
		len_dif_zero++;
	else
		len_dif_zero=0;
        fprintf(stderr, "Stat: %lld %d %f %f %lld %d %lld %f %d %d\n", count_pair, total_c, current_read_ratio, read_ratio_dif, total_len, total_usedc, ass_len, len_dif, avg_len, len_dif_zero);
        if(read_ratio_dif < MINRATEDIF || (len_dif < MINLENDIF && len_dif_zero >= 2) || total_c > increase || is_platform > 2)//|| (is_platform == 1 && avg_len < old_avg_len && avg_len < Read_length/2))
        {
            fprintf(stderr, "\nStop extension: %f < %f, or %d > %d, or platform=%d \n", read_ratio_dif, MINRATEDIF, total_c, increase/2, is_platform);
            return 1;
        }
        count_pair += increase;
        endInf[0].len = 0; endInf[0].uc = 0;
	read_ratio = current_read_ratio;
	ass_len = total_usedlen;
	old_avg_len = avg_len;
    }else{
        if(avg_len >= old_avg_len && avg_len <= MINOVERLAP)
            is_platform ++;
        else if(avg_len < old_avg_len)
            is_platform = 0;

        fprintf(stderr, "Stat: %d %lld %lld %d %d %lld %f %f %f\n", total_c, total_used, total_len, total_usedc, avg_len, total_usedlen, len_dif, current_read_ratio, read_ratio_dif);
        if((read_ratio_dif < MINRATEDIF && (current_read_ratio > 50 || read_ratio_dif == 0)) || is_platform > 2 || (read_ratio_dif < MINRATEDIF*threadNum && is_platform >=2) || len_dif < MINLENDIF)
            return 1;
    
        if(current_read_ratio > 80.5 || (read_ratio_dif < 0.1 && current_read_ratio > 50))// || to_single == 1)
        {
            for(i=1; i< threadNum; i++)
            {
                threadEnd[i] = 1;
                endInf[i].c = 0;
            }
            gthreadNum = 1;
            fprintf(stderr, "Starting single thread for finishing assembly\n");
            count_pair_thread = count_pair/2;

            if( current_id[0].c > count_pair)
                count_pair = (current_id[0].c / count_pair + 2) * count_pair;
            is_platform = 0;
        }
	read_ratio = current_read_ratio;
	old_avg_len = avg_len;
        ass_len = total_usedlen;
    }
    return 0;
}

void layout_contigs_bwt_multi(char* prefix, unsigned long long* contig_num, int method)
{
    int c=0, size = 0, i, j;
    clock_t start_time, now_time;
    char name[255];
    int max_node_num = MAXDEPTH * Read_length;
    len_dif_zero = 0; blank_cycle_time = 0;
    count_pair = pair_num * 1.0 / double(Expect_depth * 5); // G/500/2, with expect N50 >1k. then one window, one genome.

    count_pair_thread = count_pair/threadNum;
    gthreadNum = threadNum;
    increase = count_pair;

    buffer_size = count_pair_thread /2 / 10 /3; //10 buffers with thread num could fill in all the whole genome(N50=1k).
    if(buffer_size < 10)
        buffer_size = 10;

    fprintf(stderr, "Checking window for termination: %d and perthread %d and buffer size %d\n", count_pair, count_pair_thread, buffer_size);

    //find the starting library.
    while(FLAGS[size] == 0 && size < bwtCount)
        size++;
    if(size == bwtCount)
    {
        fprintf(stderr, "You must at least set one bwt for seeding.\n");
	    exit(6);
    }
    
    start_time = clock();
    
    //define and initialization.
    
    ctgs = (Contig**)calloc(threadNum, sizeof(Contig*));
    ctg_num = (int*)calloc(threadNum, sizeof(int));
    tmpRead = (TmpRead**)calloc(threadNum, sizeof(TmpRead*));
    tmpCount = (int*)calloc(threadNum, sizeof(int));
    g_iters = (Iter**)calloc(threadNum, sizeof(Iter*));
    iterCount = (int*)calloc(threadNum, sizeof(int));
    
    g_seq = (char**)calloc(threadNum, sizeof(char*));
    g_nodes = (Node**)calloc(threadNum, sizeof(Node*));
    g_levs  = (Lev**)calloc(threadNum, sizeof(Lev*));
    g_reads = (ReadId**)calloc(threadNum, sizeof(ReadId*));

    tmpbase = (Base**)calloc(threadNum, sizeof(Base*));
    sbase = (Base***)calloc(threadNum, sizeof(Base**));
    
    threads = (pthread_t*) calloc( threadNum, sizeof(pthread_t) );
    int threadId[threadNum];
    signal = (int*)calloc(threadNum, sizeof(int));
        
    conf_buffer = (Id**)calloc(threadNum, sizeof(Id*));
    conf_size = (int*)calloc(threadNum, sizeof(int));
    current_id = (Id*)calloc(threadNum, sizeof(Id));
    threadEnd = (int*)calloc(threadNum, sizeof(int));
    endInf = (EndInf*)calloc(threadNum, sizeof(EndInf));

    four_iter_used_read = (int*)calloc(threadNum, sizeof(int));
    four_iter_checked_read = (int*)calloc(threadNum, sizeof(int));
    four_iter = (int*)calloc(threadNum, sizeof(int));   
 
    for(i=0; i< threadNum; i++)
    {
        ctgs[i] = (Contig*)calloc(buffer_size+1, sizeof(Contig));
        conf_buffer[i] = (Id*)calloc(buffer_size+1, sizeof(Id));
        conf_size[i]=0;
        ctg_num[i] = 0;
        for(j=0; j< buffer_size; j++)
        {
            ctgs[i][j].seq = (char*)calloc(MaxCtgLen, sizeof(char));   
            //ctgs[i][j].qual = (char*)calloc(MaxCtgLen, sizeof(char));
            ctgs[i][j].uniq = (char*)calloc(MaxCtgLen, sizeof(char));
        }
        
        tmpRead[i] = (TmpRead*)calloc(MaxReadCount, sizeof(TmpRead));
        tmpCount[i] = 0;
        
        g_iters[i] = (Iter*)calloc(MaxIterNum, sizeof(Iter));
        iterCount[i] = 0;
        for(j=0; j< MaxIterNum; j++)
        {
            g_iters[i][j].seq = (char*)calloc(Read_length, sizeof(char));
            g_iters[i][j].b = (Base*)calloc(bwtCount, sizeof(Base));
        }
        
        g_seq[i] = (char*)calloc(Read_length, sizeof(char));
        g_nodes[i] = (Node*)calloc(max_node_num, sizeof(Node));
        g_levs[i] = (Lev*)calloc(Read_length+1, sizeof(Lev));
        g_reads[i] = (ReadId*)calloc(MAXDEPTH, sizeof(ReadId));
        for(j=0; j< max_node_num; j++)
        {
            g_nodes[i][j].seq = (char*) calloc(Read_length, sizeof(char));
            g_nodes[i][j].len = 0;
            g_nodes[i][j].conf = 0;
            g_nodes[i][j].b = (Base*) calloc(bwtCount, sizeof(Base));
        }

        for(j=0; j< Read_length+1; j++)
        {
            g_levs[i][j].c =(unsigned short*) calloc(5, sizeof(unsigned short));
        }
        
        tmpbase[i] = (Base*)calloc(bwtCount, sizeof(Base));
        for(j=0; j< bwtCount; j++)
        {
            tmpbase[i][j].saL = 0; tmpbase[i][j].saR = 0;  tmpbase[i][j].saLL = 0; tmpbase[i][j].saRR = 0;
            tmpbase[i][j].plus = 0; tmpbase[i][j].minus = 0; tmpbase[i][j].depth = 0;
        }

        sbase[i] = (Base**)calloc(4, sizeof(Base*));
        for(j=0; j< 4; j++)
            sbase[i][j] = (Base*)calloc(bwtCount, sizeof(Base));
        endInf[i].c = 0; endInf[i].ur = 0; endInf[i].uc = 0; endInf[i].len = 0; endInf[i].ulen = 0;
    }
    fprintf(stderr, "After initial, check mem:\n");
    print_maxrss_and_time();
    //generate files.
    //sprintf(name, "%s.contig", prefix);
    //fci = fopen(name, "w");
    sprintf(name, "%s.contig.fa", prefix);
    fcs = fopen(name, "w");
    //sprintf(name, "%s.contig.fa.maskRep", prefix);
    //fcm = fopen(name, "w");
    //sprintf(name, "%s.contig.qual", prefix);
    //fcq = fopen(name, "w");
    sprintf(name, "%s.contig.log", prefix);
    logs = fopen(name, "w");
    
    now_time = clock();

    double startTime;
    double elapsedTime = 0, totalElapsedTime = 0, prevTime = 0;
    startTime = setStartTime();

    fprintf(stderr, "Initial time: %f seconds\n", (double)(now_time - start_time)/CLOCKS_PER_SEC );
    creatThrds (threads, threadId);

    //define for termination of assembly.
    initial_current(size);
    int cycle=0, to_single = 0;;
    while(is_all_end() == 0)
    {
        //assembly per thread.
        sendSignal(1);

        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        if(to_single == 0 && prevTime > 0 && elapsedTime > prevTime*Upper_bound)
            to_single = 1;

        totalElapsedTime += elapsedTime;
        if(prevTime == 0 || prevTime > elapsedTime)
            prevTime = elapsedTime;
        now_time = clock();

        if(RepeatMask == 1)
            print_ctgs(fcs, fcq, fcm, logs);
        else
            print_ctg(fcs, fcq, logs); 

        //check terminate.
        if(is_terminate( to_single ) == 1)
        {
            fprintf(stderr, "Extension finished, ReadUsedRate %f, TotalLength %lld, CPU time %f, ", read_ratio, ass_len, (double)(now_time - start_time)/CLOCKS_PER_SEC);
            fprintf(stderr, "Realtime "); printElapsedTime(stderr, FALSE, FALSE, TRUE, 2, totalElapsedTime);
            break;
        }
        fprintf(stderr, "%d cycle buffer finished, ReadUsedRate %f, TotalLength %lld, CPUtime %f, ", cycle, read_ratio, ass_len, (double)(now_time - start_time)/CLOCKS_PER_SEC);
        fprintf(stderr, "Realtime "); printElapsedTime(stderr, FALSE, FALSE, TRUE, 2, elapsedTime);

        cycle++;
        //checking more library.
        if(is_all_end() == 1)
        {
            size++;
            while(FLAGS[size] == 0 && size < bwtCount)
                size++;
            if(size >= bwtCount)
                break;
            fprintf(stderr, "Checking from library %d \n", size);
            initial_current(size);
        }
    }
    *contig_num = g_count;
    //fprintf(stderr, "All extension finished.\n"); 
    //free memory, thread and close file.

    sendSignal(2);
    //fprintf(stderr, "All thread freeed.\n");

    for(i=0; i< threadNum; i++)
    {
        for(j=0; j< buffer_size; j++)
        {
            free(ctgs[i][j].seq);
            //free(ctgs[i][j].qual);
            free(ctgs[i][j].uniq);
        }

        for(j=0; j< MaxIterNum; j++)
        {
            free(g_iters[i][j].seq);
            free(g_iters[i][j].b);
        }

        for(j=0; j< max_node_num; j++)
        {
            free(g_nodes[i][j].seq);
            free(g_nodes[i][j].b);
        }

        for(j=0; j< Read_length+1; j++)
        {
            free(g_levs[i][j].c);
        }
        
        for(j=0; j< 4; j++)
            free(sbase[i][j]);
    }

    for(i=0; i< threadNum; i++)
    {
        free(ctgs[i]);
        free(tmpRead[i]);
        free(conf_buffer[i]);
        
        free(g_iters[i]);
        free(g_seq[i]);
        free(g_reads[i]);
        
        free(g_nodes[i]); 
        free(g_levs[i]); 
        free(tmpbase[i]);
        free(sbase[i]);
    }
    free(conf_buffer);
    free(conf_size);
    free(current_id);
    free(threadEnd);

    free(four_iter_used_read);
    free(four_iter_checked_read);
    free(four_iter);
    
    free(ctgs);
    free(ctg_num);
    free(tmpRead);
    free(tmpCount);
    free(g_iters);
    free(iterCount);
    
    free(g_seq);
    free(g_nodes);
    free(g_levs);
    free(g_reads);

    free(tmpbase); 
    free(sbase);
    free(threads);
    free(signal);
    free(endInf);

    //fclose(fci);
    fclose(fcs);
    //fclose(fcq);
    //fclose(fcm);
    fclose(logs);
}

void contig_assembly(char* prefix, int method)
{
    char name[255];
    sprintf(name, "%s.num", prefix);

    pair_num =0;
    int i;
    for(i=0; i< bwtCount; i++)
    {
        if(FLAGS[i] == 0)
            continue;
        pair_num += idx2BWT[i]->numReads;
    }
    unsigned long long ctg_num = 0;
    fprintf(stderr, "There are totally %d read number for contig assembly\n", pair_num);

    layout_contigs_bwt_multi(prefix, &ctg_num, method);
    fprintf(stderr, "layout unique contig number: %d\n", ctg_num);
    fprintf(stderr, "Layout all contigs finished\n");
}

