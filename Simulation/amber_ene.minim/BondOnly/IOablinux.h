/*ab model mer*/
#include<stdio.h>
#include<string.h>
#include <math.h>
#define MEMERR {fprintf(stderr,"\nNot sufficient memory\n");exit(1);}
#define MAXLENGTH 200 
#define ATOM MAXLENGTH 
#define NATOM 110
#define NTYPE 200
#define NRES NATOM 
#define MAXPAIR 229
#define MAXHBPAIR 66

#define K 0.6/300


struct config{
double x;
double y;
double z;
};

struct res{
char type[4];
char name[2][4];
int atomindex[2];
int typeindex[2];
int num;
struct config coor[2];
};


struct res native[MAXLENGTH];
/*struct config polymer[MAXLENGTH];*/

struct config polymer[MAXLENGTH];
struct config polymer1[MAXLENGTH];

