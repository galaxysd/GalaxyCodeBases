// From http://faq.cprogramming.com/cgi-bin/smartfaq.cgi?answer=1042856625&id=1043284385
// Edited by Hu Xuesong @ Thu Apr 28 03:29:52 CST 2011

#include <stdio.h>
#include <termios.h>
#include <unistd.h>

int mygetch ( void ) {
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

int pressAnyKey (void) {
  if ( !isatty(STDIN_FILENO) )
    return -2;	// # define EOF (-1) in <stdio.h>
// other errno in /usr/include/asm-generic/errno-base.h
  fputs("\nPress any key to continue ...\n", stderr);
  return mygetch();
}
