/*ab model mer*/
#include<stdio.h>
#include<string.h>
#include <math.h>
#define MEMERR {fprintf(stderr,"\nNot sufficient memory\n");exit(1);}
#define MAXLENGTH 200 
#define ATOM MAXLENGTH 
#define NATOM 110
#define NTYPE 200
/*#define MAXBOX 374 */

#define NRES NATOM 
/*#define MAXVdW 53 
#define MAXHB 55 */
#define MAXPAIR 229 /*total number of native contacts to be analyzed*/
/*#define MAXNONPAIR 840 */            /*total number of nonnative contacts to be analyzed*/
#define MAXHBPAIR 66 


#define K 0.6/300

/*#define THRESLENGTH 1.94*/

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
double side;
struct config coor[2];
};


struct res native[MAXLENGTH];
struct config polymer[MAXLENGTH];

struct config npolymer[MAXLENGTH];
/*struct config polymer[MAXBOX];*/
struct config polymer1[MAXLENGTH];
struct config centercoor[1];


