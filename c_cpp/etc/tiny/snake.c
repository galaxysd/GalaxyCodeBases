/* http://bbs.chinaunix.net/thread-1170202-1-1.html
许久没来，贡献一个原创小程序，是unix下用终端控制字符实现的贪吃蛇小游戏，供大家一乐。

编译方法：
AIX下：cc -qcpluscmt -o snake snake.c
sco unix下：cc -o snake snake.c

[ 本帖最后由 forest077 于 2008-6-25 13:23 编辑 ]
*/
/***snake.c***/
#include <stdio.h>
#include <malloc.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/select.h>
#include <termio.h>
#include <fcntl.h>
#include <unistd.h> // ttyname
#include <stdlib.h> // rand, srand
#include <string.h> // strlen

#define SNAKE_INITX 5
#define SNAKE_INITY 5
#define SNAKE_SHAPE '*'
#define SNAKE_INITLEN 8

#define WIN_X1 1
#define WIN_X2 80
#define WIN_Y1 1
#define WIN_Y2 24

#define MAX_LEVEL 20
#define MAX_INTER 200000
#define MIN_INTER 0
#define MAX_RICH 10
#define DEVRATE 5
#define OVER "Game Over!!!"

struct stNode
{
		int x;
		int y;
		char shape;
		struct stNode *next;
};

struct stFood
{
		int x;
		int y;
};

struct stNode *gpstHead,*gpstTail;
struct stFood gastFood[MAX_RICH];
int giLevel=1;
int giRich=1;
int giScore=0;
int giLen=0;

void settty(int iFlag)
{
		int fd;
		struct termio stTerm;

		if((fd = open(ttyname(1),O_RDWR))==-1)		return;
		if(iFlag == 1)
		{
				ioctl(fd,TCGETA,&stTerm);
				stTerm.c_lflag &= ~ICANON;
				stTerm.c_lflag &= ~ECHO;
				stTerm.c_cc[4] = 1;
				stTerm.c_cc[5] = 0;
				stTerm.c_iflag &= ~ISTRIP;
				stTerm.c_cflag |= CS8;
				stTerm.c_cflag &= ~PARENB;
				ioctl(fd,TCSETA,&stTerm);
		}
		else
		{
				ioctl(fd,TCGETA,&stTerm);
				stTerm.c_lflag |= ICANON;
				stTerm.c_lflag |= ECHO;
				stTerm.c_cc[4] = 4;
				stTerm.c_cc[5] = 5;
				stTerm.c_iflag &= ~ISTRIP;
				stTerm.c_cflag |= CS8;
				stTerm.c_cflag &= ~PARENB;
				ioctl(fd,TCSETA,&stTerm);
		}
		close(fd);
}

void vDrawOneNode(struct stNode *pstNode,int iFlag)
{
		printf("\033[%dm\033[40;%dm\033[%d;%d;H%c",
				iFlag,iFlag*3+30,pstNode->y,pstNode->x,pstNode->shape);
		fflush(stdout);
}

void vDrawOneFood(int x,int y)
{
		printf("\033[1m\033[40;36m\033[%d;%d;H%c",y,x,'@');
		fflush(stdout);
}

int iGetDir(int iOriDir)
{
		fd_set rset;
		struct		timeval		hTmo;
		int iRet,iFlag=0;
		char cCh;

		FD_ZERO(&rset);
		FD_SET(0,&rset);
		hTmo.tv_sec=0;
		hTmo.tv_usec=MAX_INTER-(MAX_INTER-MIN_INTER)/MAX_LEVEL*giLevel;

		iRet=select(1,&rset,NULL,NULL,&hTmo);
		if(iRet<=0)
		{
				return(iOriDir);
		}
		for(;;)
		{
				cCh=getchar();
				if(cCh != -1)
				{
						switch(cCh)
						{
						case 27  :
						case 91  :
								iFlag++;
								break;
						case 65  ://UP
						case 66 ://DOWN
						case 67  ://RIGHT
						case 68 ://LEFT
						if(iFlag==2)
								return((!((cCh-0x41)^iOriDir^1))^(cCh-0x41));
						default  :
								return(iOriDir);
						}
				}
		}
}
void vInitScreen()
{
		settty(1);
		printf("\033[?25l\033[2J");
}

void vRestoreScreen()
{
		printf("\033[24;1H\033[1m\033[40;34m\033[?25h");
		settty(0);
}

void vDrawScope()
{
		int i,j;

		for(j=WIN_Y1;j<=WIN_Y2;j+=WIN_Y2-WIN_Y1)
		{
				printf("\033[%d;%dH+",j,WIN_X1);
				for(i=WIN_X1+1;i<WIN_X2;i++)
						printf("-");
				printf("+");
		}
		for(i=WIN_Y1+1;i<WIN_Y2;i++)
				printf("\033[%d;%dH|%*c|\n",i,WIN_X1,WIN_X2-WIN_X1-1,' ');
}

void vCreateSnake()
{
		struct stNode *pstNew;
		int i;

		gpstHead=(struct stNode*)malloc(sizeof(struct stNode));
		gpstHead->x=SNAKE_INITX;
		gpstHead->y=SNAKE_INITY;
		gpstHead->shape=SNAKE_SHAPE;
		gpstHead->next=NULL;
		vDrawOneNode(gpstHead,1);
		gpstTail=gpstHead;
		for(i=1;i<SNAKE_INITLEN;i++)
		{
				pstNew=(struct stNode*)malloc(sizeof(struct stNode));
				pstNew->x=gpstHead->x+1;
				pstNew->y=gpstHead->y;
				pstNew->shape=SNAKE_SHAPE;
				pstNew->next=NULL;
				vDrawOneNode(pstNew,1);
				gpstHead->next=pstNew;
				gpstHead=pstNew;
		}
		return;
}

void vKillSnake()
{
		struct stNode *pstNode;

		for(pstNode=gpstTail;gpstTail!=NULL;)
		{
				gpstTail=pstNode->next;
				free(pstNode);
				pstNode=gpstTail;
		}
}

void vGenFood(int iIdx)
{
		struct stNode *pstNode;
		int i,iFound=0;

		for(;!iFound;)
		{
				iFound=1;
				gastFood[iIdx].x=rand()%(WIN_X2-WIN_X1-1)+WIN_X1+1;
				gastFood[iIdx].y=rand()%(WIN_Y2-WIN_Y1-1)+WIN_Y1+1;
				for(i=0;i<giRich;i++)
				{
						if(i!=iIdx && gastFood[iIdx].x==gastFood[i].x &&
								gastFood[iIdx].y==gastFood[i].y)
						{
								iFound=0;
								break;
						}
				}
				if(!iFound) continue;
				for(pstNode=gpstTail;pstNode!=NULL;pstNode=pstNode->next)
				{
						if(gastFood[iIdx].x==pstNode->x &&
								gastFood[iIdx].y==pstNode->y)
						{
								iFound=0;
								break;
						}
				}
				if(!iFound) continue;
		}
		vDrawOneFood(gastFood[iIdx].x,gastFood[iIdx].y);
}

void vInitFood()
{
		int i;

		srand(getpid());
		for(i=0;i<giRich;i++)		vGenFood(i);
}

int iIsValid(int x,int y)
{
		struct stNode *pstNode;

		if(x<=WIN_X1 || x>=WIN_X2 || y<=WIN_Y1 || y>=WIN_Y2)
				return(0);
		pstNode=gpstTail;
		for(;pstNode!=NULL;)
		{
				if(x==pstNode->x && y==pstNode->y)
						return(0);
				pstNode=pstNode->next;
		}
		return(1);
}

int iEat(int x,int y)
{
		int i;

		for(i=0;i<giRich;i++)
		{
				if(x==gastFood[i].x && y==gastFood[i].y)
				{
						vGenFood(i);
						giScore+=giLevel*10;
						giLen++;
						if(giLevel<MAX_LEVEL)
								if(giLen%DEVRATE==0)
										giLevel++;
						return(1);
				}
		}
		return(0);
}
int main()
{
		int iDir=2,iNextX,iNextY;
		struct stNode *pstNew;
		char sPrompt[80];

		vInitScreen();
		vDrawScope();
		vCreateSnake();
		vInitFood();
		for(;;)
		{
				iDir=iGetDir(iDir);
				iNextX=gpstHead->x+(iDir>>1)*(5-(iDir<<1));
				iNextY=gpstHead->y-(!(iDir>>1))*(1-(iDir<<1));
				if(!iIsValid(iNextX,iNextY))
				{
						printf("\033[%d;%zdH\033[1m\033[40;34m%s\033[0m",
								WIN_Y2-1,(WIN_X1+WIN_X2)/2-strlen(OVER)/2,OVER);
						break;
				}
				pstNew=(struct stNode*)malloc(sizeof(struct stNode));
				pstNew->x=iNextX;
				pstNew->y=iNextY;
				pstNew->shape=SNAKE_SHAPE;
				pstNew->next=NULL;
				gpstHead->next=pstNew;
				gpstHead=pstNew;
				vDrawOneNode(gpstHead,1);
				if(!iEat(iNextX,iNextY))
				{
						vDrawOneNode(gpstHead,1);
						vDrawOneNode(gpstTail,0);
						pstNew=gpstTail;
						gpstTail=pstNew->next;
						free(pstNew);
				}
				sprintf(sPrompt,"Score:%7d Level:%2d",giScore,giLevel);
				printf("\033[%d;%zdH\033[1m\033[40;34m%s\033[0m",
						WIN_Y2,(WIN_X1+WIN_X2)/2-strlen(sPrompt)/2,sPrompt);
		}
		vKillSnake();
		vRestoreScreen();
}
