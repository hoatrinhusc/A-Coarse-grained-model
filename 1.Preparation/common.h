#include <fstream>
#include <iostream>
#include <cstdio>
#include <iomanip>
#include <cstring>
#include <stdio.h>

using namespace std;
#define PI 3.14159265
#define NUMCA 110
#define NUMCB 90
#define NATOM 773
#define NUMAB 200

struct ATOMIC
{
  double x[NUMCA];
  double y[NUMCA];
  double z[NUMCA];
  int num[NUMCA];
};


