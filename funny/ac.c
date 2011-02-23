#include <stdio.h>
#include <unistd.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*
0	Reset
1	FG BRIGHT
5	BG BRIGHT under sh/putty
7	REVERSE FG and BG
4	UNDERLINE under gnome-terminal and putty
8	HIDDENT (sh => fg black, gnome-terminal => fg==bg)
9	Strikethrough under gnome-terminal
*/

enum ANSIColors {
    Rest=0,
    FgBright=1,
    BgBright=5,

};

#define ANSIColor_RESET       0
#define ANSIColor_FGBRIGHT    1
#define ANSIColor_BGBRIGHT    5
#define ANSIColor_REVERSE     7

#define ANSIColor_UNDERLINE   4
#define ANSIColor_HIDDENT     8
#define ANSIColor_DELETELINE  9

#define ANSIColor_BLACK       0
#define ANSIColor_RED         1
#define ANSIColor_GREEN       2
#define ANSIColor_YELLOW      3
#define ANSIColor_BLUE        4
#define ANSIColor_MAGENTA     5	// 品红
#define ANSIColor_CYAN        6
#define ANSIColor_WHITE       7

/* Accroding to line 261 of [glibc.git]/stdio-common/vfprintf.c
( http://sourceware.org/git/?p=glibc.git;a=blob;f=stdio-common/vfprintf.c ),
'k' is not used. Thus, I choose '%k' to hacking on. (another reason has sth. to do with the ass)
So, %k(int), %k(int,int), %k(int,int,int).
Also,
16    printf("[%k]\n");
17    printf("[%k]\n",1);
18    printf("[%%k]\n");
t.c:16: warning: unknown conversion type character ‘k’ in format
t.c:17: warning: unknown conversion type character ‘k’ in format
t.c:17: warning: too many arguments for format
*/
int colorprintf(const char *format, ...) {  // ANSI color is only useful to STDOUT, thus no FILE *stream here.
    const char cESC='\033';  // all that not conversion  specifications is ordinary characters, which are simply copied to the output stream.
    size_t inlen = strlen(format);
    char *newformat = (char *) malloc( inlen );
    char *ksubstr   = (char *) malloc( inlen );
    const char *old_pos;
    char *des_pos, *ksub_pos;
    old_pos=format;
    des_pos=newformat;
    register char c;
    while ( (c = *old_pos++) ) {
        if (c != '%') {
            *des_pos++ = c;
        } else {
            if ( (c = *old_pos++) != '\0' ) {
                if (c != 'k') {
                    *des_pos++ = '%';
                    *des_pos++ = c;   // "%%" will be OK since c == '%' here.
                } else {
                    if ( (c = *old_pos++) != '\0' ) {
                        if (c == '(') {
                            ksub_pos=ksubstr;
                            char flag=0;
                            char last_num_count=0;
                            size_t sublen = 0;
                            while ( (c = *old_pos++) && c != ')' ) {
                                ++sublen;
                                if (c>='0' && c<='9') {
                                    *ksub_pos++ = c;
                                    ++last_num_count;
                                } else {
                                    if ( c==',' ) {
                                        if (last_num_count) {
                                            *ksub_pos++ = ';';
                                        }
                                    } else {
                                        flag=1;
                                    }
                                    last_num_count=0;
                                }
                            }
                            if (last_num_count) {
                                *ksub_pos = '\0';
                            } else {
                                *--ksub_pos = '\0';
                            }
                            if ( c == ')' ) {
                                if ( flag==0 && sublen ) {  // is %k\([^)]+\)
                                    *des_pos++ = cESC;
                                    *des_pos++ = '[';
                                    ksub_pos=ksubstr;
                                    while ( (*des_pos++ = *ksub_pos++) )
                                        ;
                                    --des_pos;
                                    *des_pos++ = 'm';
                                }
                            } else {    // c == '\0', EOL
                                break;
                            }
                        } else {
                            *des_pos++ = '%';
                            *des_pos++ = 'k';
                            //*des_pos++ = c;   // what if c == '%' ?
                            --old_pos;  // make another jump is slower than DEC.
                        }
                    } else {
                        *des_pos++ = '%';
                        *des_pos++ = 'k';
                        break;
                    }
                }
            } else {
                *des_pos++ = '%'; // the printf will ignore the tailing '%';
                break;
            }
        }
    }
    *des_pos = '\0';    // when c == '\0', copy not done.
    free(ksubstr);
    //printf("[%s]->[%s]\n",format,newformat);
    puts(format);
    puts("->");
    puts(newformat);

    va_list arg;
    int done;
    va_start (arg, format);
    done = vfprintf (stdout, newformat, arg);
    va_end (arg);

    free(newformat);
    return done;
}

int main(int argc, char *argv[]) {
/*
    printf("isatty = %d\n", isatty(0));                                
    if ( isatty(STDOUT_FILENO) ) {                                     
        char e='\033';                                                 
        puts("[\033[1;32mWith color.\033[0m]");                        
    } else {                                                           
        puts("Without Color.");                                        
        fputs("STDERR Without Color.", stderr);                        
    }                                                                  
    colorprintf("[%k(1,32)With color.[%s]   [%s]%k(0)\n","t1","t2 t2");
*/
    //colorprintf("%k(1");
    colorprintf(argv[1],argv[2]);
    return 0;
}

