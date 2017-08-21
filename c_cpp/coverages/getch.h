// From http://faq.cprogramming.com/cgi-bin/smartfaq.cgi?answer=1042856625&id=1043284385
// Edited by Hu Xuesong @ Thu Apr 28 CST 2011

#ifndef _GA_GETCH_H
#define _GA_GETCH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <termios.h>
#include <unistd.h>

static inline int mygetch ( void ) {
  int ch;
  struct termios oldt, newt;

  tcgetattr( STDIN_FILENO, &oldt );
  newt = oldt;
  newt.c_lflag &= ~( ICANON | ECHO );
  tcsetattr( STDIN_FILENO, TCSANOW, &newt );
  ch = getchar();
  tcsetattr( STDIN_FILENO, TCSANOW, &oldt );

  return ch;
}

static inline int pressAnyKey (void) {
  if ( !isatty(STDIN_FILENO) )
    return -2;	// # define EOF (-1) in <stdio.h>
// other errno in /usr/include/asm-generic/errno-base.h
  fputs("\nPress any key to continue ... ", stderr);
  //return mygetch();
  int ch = mygetch();
  fputs("\n", stderr);
  return ch;
}

#ifdef __cplusplus
}
#endif

#endif /* getch.h */
