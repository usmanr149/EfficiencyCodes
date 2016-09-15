#ifndef HEADER_FILE
#define HEADER_FILE

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h> 

#define ROWS 20
#define COLS 20
#define dim 7  //declare matrix dimension here

long double eMatrix[ROWS][COLS];
//long double b[ROWS][COLS];
long double xVector[ROWS], bVector[ROWS];
long double xi[40], yi[40],ei[40],xi0[40],yi0[40],ei0[40];
char line[132],*file_name,*out_file,*out_file2;
FILE *inp,*out,*fp1;
int N=3;
int n,count,count0, upper_limit=2,upper_limit_best=0,zig=0;
long double a=0,b=0,c=0,g=0,alpha=0, beta=0,d,e,f;
long double energy, energy0=0, min=100000,record_energy=0,abest=0,bbest=0,cbest=0,gbest=0,alphabest=0,betabest=0,dbest=0,ebest=0,fbest=0;
long double sum=0, sumx=0,sumx2=0,sumx3=0,sumx4=0,f=0,fx=0,fx2=0;
  long double sumj00=0, sumxj01=0,sumxj02=0,sumxj03=0,sumxj11=0,fj=0,fxj=0,fxj2=0,fxj3=0,sumxj12=0,sumxj13=0;
  long double sumxj22=0,sumxj23=0,sumxj33=0;


#endif
