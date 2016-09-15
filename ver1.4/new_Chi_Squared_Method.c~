/************************************************************/
/* Usman Rizwan                                             */ 
/*                                                          */
/*                                                          */
/* A program to solve linear equations via gaussian         */
/* elimination. This program will incorporate scaling       */
/* scaling parameter as well. This program is currently     */
/* in working order. I haven't incorporated varting x0.     */
/************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h> 

#define ROWS 20
#define COLS 20

long double eMatrix[ROWS][COLS];
long double bVector[ROWS];
long double xVector[ROWS];
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
int dim=7;  //declare matrix dimension here

int DescendingOrder(){
  int i,j;
  long double en=0,eff=0,err=0;

  printf("Count from DescendingOrder = %d\n", count0);

  for(i=0;i<count0;i++)
    for(j=i+1;j<count0;j++)
      if (xi0[i]>xi0[j]){
	  en=xi0[i];
	  eff=yi0[i];
	  err=ei0[i];
	  xi0[i]=xi0[j];
	  yi0[i]=yi0[j];
	  ei0[i]=ei0[j];
	  xi0[j]=en;
	  yi0[j]=eff;
	  ei0[j]=err;
      }
}

int printmatrix(){
   //print matrix "eMatrix"
   int k,l;
   printf("Print Matrix called\n");
   for(k=0; k<dim; k++){
       for(l=0; l<dim; l++){
	 printf(" %Lf*X%d + ", eMatrix[k][l],l);
       }
       
       printf(" =%Lf\n", bVector[k]);
       //xVector[k]=0;
   }
}

int ReadCo60()
{
   FILE * fp;

  if((fp=fopen("Co60_abs_eff.dat","r"))==NULL)
    {
      printf("Cannot open Co60_abs_eff.dat.\n");
      //exit(-1);
      return 0;
    }

  printf("Parameters read from Co60_abs_eff.dat\n");

    if(fgets(line,132,fp)!=NULL)
      {
	count=0;
	while(fscanf(fp,"%Lf %Lf %Lf\n",&xi[count],&yi[count],&ei[count])!=EOF)
	  {
	    printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
	printf("Count for Co-60 = %d\n\n\n", count);
      }
}

int ReadCo60_nolog()
{
   FILE * fp;

  if((fp=fopen("Co60_abs_eff_nolog.dat","r"))==NULL)
    {
      printf("Cannot open Co60_abs_eff_nolog.dat.\n");
      //exit(-1);
      return 0;
    }

  printf("Parameters read from Co60_abs_eff.dat\n");

    if(fgets(line,132,fp)!=NULL)
      {
	count=0;
	while(fscanf(fp,"%Lf %Lf %Lf\n",&xi[count],&yi[count],&ei[count])!=EOF)
	  {
	    printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
	printf("Count for Co-60 = %d\n\n\n", count);
      }
}

int ReadEu152()
{
   FILE * fp;

  if((fp=fopen("Eu152_rel_eff.dat","r"))==NULL)
    {
      printf("Cannot open Eu152_rel_eff.dat.\n");
      //exit(-1);
      return 0;
    }

  printf("Parameters read from Eu152_rel_eff.dat\n");

    if(fgets(line,132,fp)!=NULL)
      {
	count=0;
	while(fscanf(fp,"%Lf %Lf %Lf\n",&xi[count],&yi[count],&ei[count])!=EOF)
	  {
	    printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
	printf("Count for Eu-152 = %d\n\n\n", count);
      }
}

int ReadEu152_nolog()
{
   FILE * fp;

  if((fp=fopen("Eu152_rel_eff_nolog.dat","r"))==NULL)
    {
      printf("Cannot open Eu152_rel_eff_nolog.dat.\n");
      //exit(-1);
      return 0;
    }

  printf("Parameters read from Eu152_rel_eff_nolog.dat\n");

    if(fgets(line,132,fp)!=NULL)
      {
	count=0;
	while(fscanf(fp,"%Lf %Lf %Lf\n",&xi[count],&yi[count],&ei[count])!=EOF)
	  {
	    printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
	printf("Count for Eu-152 = %d\n\n\n", count);
      }
}

int ReadCo56()
{
   FILE * fp;

  if((fp=fopen("Co56_rel_eff.dat","r"))==NULL)
    {
      printf("Cannot open Co56_rel_eff.dat.\n");
      //exit(-1);
      return 0;
    }

  printf("Parameters read from Co56_rel_eff.dat\n");

    if(fgets(line,132,fp)!=NULL)
      {
	count=0;
	while(fscanf(fp,"%Lf %Lf %Lf\n",&xi[count],&yi[count],&ei[count])!=EOF)
	  {
	    
	    printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
}
}

int ReadCo56_nolog()
{
   FILE * fp;

  if((fp=fopen("Co56_rel_eff_nolog.dat","r"))==NULL)
    {
      printf("Cannot open Co56_rel_eff_nolog.dat.\n");
      //exit(-1);
      return 0;
    }

  printf("Parameters read from Co56_rel_eff_nolog.dat\n");

    if(fgets(line,132,fp)!=NULL)
      {
	count=0;
	while(fscanf(fp,"%Lf %Lf %Lf\n",&xi[count],&yi[count],&ei[count])!=EOF)
	  {
	    
	    printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
}
}

int ReadBa133()
{
   FILE * fp;
  if((fp=fopen("Ba133_rel_eff.dat","r"))==NULL)
    {
      printf("Cannot open Ba133_abs_eff.dat.\n");
      //exit(-1);
      return 0;
    }

  printf("Parameters read from Ba133_rel_eff.dat\n");

    if(fgets(line,132,fp)!=NULL)
      {
	count=0;
	while(fscanf(fp,"%Lf %Lf %Lf\n",&xi[count],&yi[count],&ei[count])!=EOF)
	  {
	    printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
      }
}

int ReadBa133_nolog()
{
   FILE * fp;
  if((fp=fopen("Ba133_rel_eff_nolog.dat","r"))==NULL)
    {
      printf("Cannot open Ba133_abs_eff_nolog.dat.\n");
      //exit(-1);
      return 0;
    }

  printf("Parameters read from Ba133_rel_eff_nolog.dat\n");

    if(fgets(line,132,fp)!=NULL)
      {
	count=0;
	while(fscanf(fp,"%Lf %Lf %Lf\n",&xi[count],&yi[count],&ei[count])!=EOF)
	  {
	    printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
      }
}

double formMatrix1(){  //This part is without the gamma function.

    int i;

    printf("Count formMatrix1 = %d\n", count);

    for(i=(zig+1);i<count;i++)
	      	{
		  //ei[i]=1;
	      	  //printf("ei[%d] = %Lf\n", i,ei[i]);
	      	  sum+=1/(ei[i]*ei[i]);
	      	  sumx+=xi[i]/(ei[i]*ei[i]);
	      	  sumx2+=xi[i]*xi[i]/(ei[i]*ei[i]);
	      	  sumx3+=xi[i]*xi[i]*xi[i]/(ei[i]*ei[i]);
	      	  sumx4+=xi[i]*xi[i]*xi[i]*xi[i]/(ei[i]*ei[i]);
	      	  f+=yi[i]/(ei[i]*ei[i]);
	      	  fx+=yi[i]*xi[i]/(ei[i]*ei[i]);
	      	  fx2+=yi[i]*xi[i]*xi[i]/(ei[i]*ei[i]);
		  //printf("sum = %Lf\n", sum);
		}
    
  }

double formMatrix2(){  //This is where the gamma function is include.

  int i;
  //sumj00=0;
  //energy0=xi0[upper_limit];
  for(i=0;i<(zig+1);i++)
	      	{
	      sumj00+=1/(ei[i]*ei[i]);;
	      sumxj01+=xi[i]/(ei[i]*ei[i]);
	      sumxj11+=xi[i]*xi[i]/(ei[i]*ei[i]);
	      sumxj02+=(2*energy0*xi[i]-energy0*energy0)/(ei[i]*ei[i]);
	      sumxj03+=(xi[i]-energy0)*(xi[i]-energy0)/(ei[i]*ei[i]);
	      sumxj12+=xi[i]*(2*energy0*xi[i]-energy0*energy0)/(ei[i]*ei[i]);
	      sumxj13+=xi[i]*(xi[i]-energy0)*(xi[i]-energy0)/(ei[i]*ei[i]);
	      sumxj22+=(2*energy0*xi[i]-energy0*energy0)*(2*energy0*xi[i]-energy0*energy0)/(ei[i]*ei[i]);
	      sumxj23+=(2*energy0*xi[i]-energy0*energy0)*(xi[i]-energy0)*(xi[i]-energy0)/(ei[i]*ei[i]);
	      sumxj33+=(xi[i]-energy0)*(xi[i]-energy0)*(xi[i]-energy0)*(xi[i]-energy0)/(ei[i]*ei[i]);
	      fj+=yi[i]/(ei[i]*ei[i]);
	      fxj+=yi[i]*xi[i]/(ei[i]*ei[i]);
	      fxj2+=yi[i]*(2*energy0*xi[i]-energy0*energy0)/(ei[i]*ei[i]);
	      fxj3+=yi[i]*(xi[i]-energy0)*(xi[i]-energy0)/(ei[i]*ei[i]);
	      //printf("sumj00 = %Lf\n", sumj00);
	      	}

}

double MatrixCoefficients(){

  formMatrix1();
  formMatrix2();
	      
      eMatrix[0][0]=sum+sumj00;
      eMatrix[0][1]=sumx+sumxj01;
      eMatrix[0][2]=sumx2+sumxj02;
      eMatrix[0][3]=sumxj03;
      eMatrix[0][4];
      eMatrix[0][5];
      eMatrix[0][6];
      eMatrix[1][0]=eMatrix[0][1];
      eMatrix[1][1]=sumx2+sumxj11;
      eMatrix[1][2]=sumx3+sumxj12;
      eMatrix[1][3]=sumxj13;
      eMatrix[1][4];
      eMatrix[1][5];
      eMatrix[1][6];
      eMatrix[2][0]=eMatrix[0][2];
      eMatrix[2][1]=eMatrix[1][2];
      eMatrix[2][2]=sumx4+sumxj22;
      eMatrix[2][3]=sumxj23;
      eMatrix[2][4];
      eMatrix[2][5];
      eMatrix[2][6];
      eMatrix[3][0]=eMatrix[0][3];
      eMatrix[3][1]=eMatrix[1][3];
      eMatrix[3][2]=eMatrix[2][3];
      eMatrix[3][3]=sumxj33;
      eMatrix[3][4];
      eMatrix[3][5];
      eMatrix[3][6];
      eMatrix[4][0]=eMatrix[0][4];
      eMatrix[4][1]=eMatrix[1][4];
      eMatrix[4][2]=eMatrix[2][4];
      eMatrix[4][3]=eMatrix[3][4];
      eMatrix[4][4]=eMatrix[0][4];
      eMatrix[4][5]=0;
      eMatrix[4][6]=0;
      eMatrix[5][0]=eMatrix[0][5];
      eMatrix[5][1]=eMatrix[1][5];
      eMatrix[5][2]=eMatrix[2][5];
      eMatrix[5][3]=eMatrix[3][5];
      eMatrix[5][4]=0;
      eMatrix[5][5]=eMatrix[0][5];
      eMatrix[5][6]=0;
      eMatrix[6][0]=eMatrix[0][6];
      eMatrix[6][1]=eMatrix[1][6];
      eMatrix[6][2]=eMatrix[2][6];
      eMatrix[6][3]=eMatrix[3][6];
      eMatrix[6][4]=0;
      eMatrix[6][5]=0;
      eMatrix[6][6]=eMatrix[0][6];
      
      bVector[0]=f+fj;
      bVector[1]=fx+fxj;
      bVector[2]=fx2+fxj2;
      bVector[3]=fxj3;
      bVector[4];
      bVector[5];
      bVector[6];

      //printmatrix();
}

double SpecialFunctionj(int k, int l)
{
  int i;
  eMatrix[k][l]=0.;
  eMatrix[k][l+1]=0.;
  eMatrix[k][l+2]=0.;
  eMatrix[k][l+3]=0.;
  bVector[l]=0.;

  for (i=0;i<(zig+1);i++){ 
    //ei[i]=1;
    eMatrix[k][l]+=1/(ei[i]*ei[i]);
    eMatrix[k+1][l]+=xi[i]/(ei[i]*ei[i]);
    eMatrix[k+2][l]+=(2*energy0*xi[i]-energy0*energy0)/(ei[i]*ei[i]);
    eMatrix[k+3][l]+=(xi[i]-energy0)*(xi[i]-energy0)/(ei[i]*ei[i]);
    bVector[l]+=yi[i]/(ei[i]*ei[i]);
  }
  
  for (i=zig+1;i<count;i++){ //This is for the side without gamma.
    //ei[i]=1;
  eMatrix[k][l]+=1/(ei[i]*ei[i]);
  eMatrix[k+1][l]+=xi[i]/(ei[i]*ei[i]);
  eMatrix[k+2][l]+=xi[i]*xi[i]/(ei[i]*ei[i]);
  eMatrix[k+3][l]+=0;
  bVector[l]+=yi[i]/(ei[i]*ei[i]);
  }
  
}

int diagonal(){
   int i, j, k;
   float temp=0;
   //printf("eMatrix[2][2] = %Lf \n", eMatrix[2][2]);
   for(i=0; i<N; i++){
       if(eMatrix[i][i]==0){
	   for(j=0; j<N; j++){
	       if(j==i) continue;
	       if(eMatrix[j][i] !=0 && eMatrix[i][j]!=0){
		   for(k=0; k<dim; k++){
		       temp = eMatrix[j][k];
		       eMatrix[j][k] = eMatrix[i][k];
		       eMatrix[i][k] = temp;
		       printf("*");
		   }
		   temp = bVector[j];
		   bVector[j] = bVector[i];
		   bVector[i] = temp;
		   break;
	       }
	   }
       }
   }
}


int SolveMatrix(){
  int i=0,j=0, k=0,l;

   //process rows
   for(k=0; k<dim; k++){
       for(i=k+1; i<dim; i++){
	   if(eMatrix[k][k]==0){
	       printf("\nSolution does not exist.\n");
	       return;
	   }
	   long double M = eMatrix[i][k] / eMatrix[k][k];
	   for(j=k; j<dim; j++){
	       eMatrix[i][j] = eMatrix[i][j] - M * eMatrix[k][j];
	   }
	   bVector[i] = bVector[i] - M * bVector[k];
       }
   }
   
   //this loop is for backward subtraction
   for(i=dim-1; i>=0; i--){
       long double s = 0;
       for(j = i; j<=dim; j++){
	 
	   s = s+eMatrix[i][j]*xVector[j];
       }
       xVector[i] = (bVector[i] - s) / eMatrix[i][i];
   }

   for(i=0;i<dim;i++)   printf("xVector[%d] = %Lf\n", i, xVector[i]);
	 

   a=xVector[0];
   b=xVector[1];
   c=xVector[2];
   g=xVector[3];
   d=xVector[4];
   e=xVector[5];
   f=xVector[6];
}

int formMatrix()
{
  int i,k,jug;
  
  DescendingOrder();

  zig=0;

  k=upper_limit;
 
  long double energy00=xi0[k];
		  energy0=xi0[k];
		  
		  printf("energy0 = %Lf\n", (powl(10,energy0)));
		  jug=1;
		  
		  do{
		  printf("energy0 = %Lf\n", (powl(10,energy0)));
		  Matrix_initialize();

  ReadCo60();
  for(i=0;i<count;i++)
    {
      
      if(xi[0]>energy0) zig=-1;
	 else if(energy0>=xi[i]) zig=i;
    }
  
  MatrixCoefficients();
  
  zig=0;

  ReadEu152();
  for(i=0;i<count;i++)
    {
      if(xi[0]>energy0) zig=-1;
	 else if(energy0>=xi[i]) 
	   zig=i;
    }
  SpecialFunctionj(0,4);
  MatrixCoefficients();
  
  
  zig=0;

  ReadBa133();
  for(i=0;i<count;i++)
    {
      if(xi[0]>energy0) zig=-1;
	 else if(energy0>=xi[i]) zig=i;
      
    }
  SpecialFunctionj(0,5);
  
  MatrixCoefficients();

  zig=0;

  ReadCo56();
  for(i=0;i<count;i++)
    {
      if(xi[0]>energy0) zig=-1;
	 else if(energy0>=xi[i]) zig=i;

    }

  SpecialFunctionj(0,6);
  MatrixCoefficients();
  diagonal();
  SolveMatrix();
  ChiSquared();
  energy0=log10(powl (10,energy00)+1*jug);
  jug++;}while(energy0<log10(2000));

}

int ChiSquared(){
  int i;
  long double chisquared=0,chisquaredl=0,chisquaredr=0;

  alpha = a - c*energy0*energy0 + g*energy0*energy0;
  beta = b + 2*c*energy0 - 2*g*energy0;
  
  ReadCo60();
  {
  for(i=0;i<count;i++)
    {
      
      if(xi[0]>energy0) zig=-1;
	 else if(energy0>=xi[i]) zig=i;
    }
  
  for(i=0;i<=zig;i++)
    {
      chisquaredl+=(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]))*(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]))/(ei[i]*ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residuall = %Lf\n", yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i])); */
      /* printf("chisquaredl = %Lf\n", chisquaredl); */
    }

  for(i=zig+1;i<count;i++)
    {
      chisquaredr+=(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]))*(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]))/(ei[i]*ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residualr = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i])); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
    }
  }

   ReadEu152();
  {
  for(i=0;i<count;i++)
    {
      
      if(xi[0]>energy0) zig=-1;
	 else if(energy0>=xi[i]) zig=i;
    }
  
  for(i=0;i<=zig;i++)
    {
      chisquaredl+=((yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+d))*(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+d))/(ei[i]*ei[i]));
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /*  printf("residuall = %Lf\n", yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+d)); */
      /*  printf("chisquaredl = %Lf\n", chisquaredl); */
    }

  for(i=zig+1;i<count;i++)
    {
      chisquaredr+=((yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+d))*(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+d))/(ei[i]*ei[i]));
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residualr = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+d)); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
    }
  /* printf("chisquared = %Lf\n", chisquared); */
  
  }

  ReadCo56();
  {
  for(i=0;i<count;i++)
    {
      
      if(xi[0]>energy0) zig=-1;
	 else if(energy0>=xi[i]) zig=i;
    }
  
  for(i=0;i<=zig;i++)
    {
      chisquaredl+=((yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+f))*(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+f))/(ei[i]*ei[i]));
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /*  printf("residuall = %Lf\n", yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+f)); */
      /*  printf("chisquaredl = %Lf\n", chisquaredl); */
    }

  for(i=zig+1;i<count;i++)
    {
      chisquaredr+=((yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+f))*(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+f)))/(ei[i]*ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residualr = %Lf\n", (yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+f))); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
    }
  /* printf("chisquared = %Lf\n", chisquared); */
  
  }

  ReadBa133();
  {
  for(i=0;i<count;i++)
    {
      
      if(xi[0]>energy0) zig=-1;
	 else if(energy0>=xi[i]) zig=i;
    }
  
  for(i=0;i<=zig;i++)
    {
      chisquaredl+=((yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+e))*(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+e))/(ei[i]*ei[i]));
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residuall = %Lf\n",(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+e))); */
      /* printf("chisquaredl = %Lf\n", chisquaredl); */
    }

  for(i=zig+1;i<count;i++)
    {
      chisquaredr+=((yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+e))*(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+e))/(ei[i]*ei[i]));
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residualr = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+e)); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
    }
  /* printf("chisquared = %Lf\n", chisquared); */
  
  }

  chisquared=chisquaredl+chisquaredr;

  printf("chisquared = %Lf\n", chisquared);

  chisquared=sqrt(chisquared/(count0-dim));
  printf("chisquared = %Lf\n", chisquared);


  if(chisquared<min)
    {
      min=chisquared;
      abest=a;
      bbest=b;
      cbest=c;
      alphabest=alpha;
      betabest=beta;
      gbest=g;
      dbest=d;
      ebest=e;
      fbest=f;
      upper_limit_best=upper_limit;
      record_energy=energy0;
    }
  
              char* output;
	      output = (char*) malloc(strlen(out_file+1));
	      strcpy(output,out_file);

	      fp1 = fopen(output,"a");
	      //fprintf(fp1,"%10.6Lf %10.6Lf\n", xi[upper_limit], chisquared);
	      fprintf(fp1,"%10.6Lf %10.6Lf\n", powl(10,energy0), chisquared);
	      //printf("%10.6Lf %10.6Lf\n", powl(10,energy0), chisquared);
	      fclose(fp1);
	      free(output);

	      }


int Matrix_initialize()
{
  int i,j;

for (i = 0; i < ROWS; i++) {
    for (j = 0; j < COLS; j++) {
      eMatrix[i][j] = 0.0;
      bVector[i]=0.0;
      xVector[i]=0.0;
    }
 }

 sumj00=0; sumxj01=0;sumxj02=0;sumxj03=0;sumxj11=0;fj=0;fxj=0;fxj2=0;fxj3=0;sumxj12=0;sumxj13=0;
  sumxj22=0;sumxj23=0;sumxj33=0;
  sum=0; sumx=0;sumx2=0;sumx3=0;sumx4=0;f=0;fx=0;fx2=0;
  sum=0; sumx=0; sumx2=0; sumx3=0; sumx4=0; f=0; fx=0; fx2=0;

}

double ScaleEfficiencies()
{
  int i;
  long double en=0,abseff=0,abserr=0;

  char* output;
	      
  output = (char*) malloc(strlen(out_file2+1));
  strcpy(output,out_file2);
  fp1 = fopen(output,"a");

  ReadCo60_nolog();
  for(i=0; i<count; i++)
    {
      abseff=yi[i];
      abserr=ei[i];
      fprintf(fp1,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], abseff, abserr);
    }

  ReadEu152_nolog();
  for(i=0; i<count; i++)
    {
      abseff=yi[i]/(powl(10,dbest));
      abserr=ei[i]/(powl(10,dbest));
      fprintf(fp1,"%10.6Lf %10.6Lf %10.6Lf\n",xi[i], abseff, abserr);
    }

		    ReadCo56_nolog();
		    for(i=0; i<count; i++)
		      {
			abseff=yi[i]/(powl(10,fbest));
			abserr=ei[i]/(powl(10,fbest));
			fprintf(fp1,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], abseff, abserr);
						    }
		    
		ReadBa133_nolog();
		    for(i=0; i<count; i++)
		      {
			abseff=yi[i]/(powl(10,ebest));
			abserr=ei[i]/(powl(10,ebest));
			fprintf(fp1,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], abseff, abserr);
						    }
	      
	      //fprintf(fp1,"%10.6Lf %10.6Lf\n", xi[upper_limit], chisquared);
	      //fprintf(fp1,"%10.6Lf %10.6Lf\n", powl(10,energy0), chisquared);
	      //printf("%10.6Lf %10.6Lf\n", powl(10,energy0), chisquared);
	      fclose(fp1);
	      free(output);



}

int gnuplot_commands()
{
FILE *pipe = popen("gnuplot","w");
//fprintf(pipe, "set data style lines\n");
//fprintf(pipe, "plot '%s' using 1:2\n",file_name);
 fprintf(pipe, "set xrange[0.1:5000]\n");
 fprintf(pipe, "set yrange[1e-4:0.01]\n");
 fprintf(pipe, "unset key\n");
 fprintf(pipe, "set logscale\n");
 fprintf(pipe, "set xlabel \"Energy [keV]\"\n");
 fprintf(pipe, "set ylabel \"Absolute Efficiency [%%]\"\n");
 fprintf(pipe, "plot 'Usmantest_scaledeff.dat' using 1:2:3 with yerrorbars, x > 0.001 ? (x**(%Lf*log10(x)))*(x**(%Lf))*(10**%Lf) : 1/0, x < %Lf ? (x**(%Lf*log10(x)))*(x**(%Lf))*(10**%Lf) : 1/0\n",cbest,bbest,abest,(powl(10,energy0)+100), gbest,betabest,alphabest);
 fprintf(pipe, "set terminal pdf\n");
 fprintf(pipe, "set output \"Best_fit.pdf\"\n");
 fprintf(pipe,"replot\n");
 fprintf(pipe,"set terminal x11\n");
 printf("Press enter to continue...\n");
close(pipe); 
fflush(pipe);
getchar();
}

int main(int argc, char *argv[]){

  int i,j;

  if(argc!=4)
    {
      printf("\n HelloWorld Filename outfilename_for_CHISquared outfilename_for_scaled_eff  \n");
      //exit(-1);
      return 0;
    }

  file_name=argv[1];
  out_file=argv[2];
  out_file2=argv[3];
  

  if((inp=fopen(argv[1],"r"))==NULL)
    {
      printf("Cannot open %s.\n",file_name);
      //exit(-1);
      return 0;
    }

  printf("Parameters read from %s.\n",file_name);

    if(fgets(line,132,inp)!=NULL)
      {
	count0=0;
	count=0;
	while(fscanf(inp,"%Lf %Lf %Lf\n",&xi0[count],&yi0[count],&ei0[count])!=EOF)
	  {
	    printf ("%d %Lf %Lf %Lf\n", count, xi0[count],yi0[count],ei0[count]);
	    count++;
	    count0++;
	  }
      }

    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < COLS; j++) {
            eMatrix[i][j] = 0.0;
	    bVector[i]=0.0;
	    xVector[i]=0.0;
        }
    }
    formMatrix();
    ScaleEfficiencies();
    printf("min = %Lf\n", min);
  printf("min energy = %Lf\n", powl(10,record_energy));
   
  printf("abest = %Lf\n", abest);
  printf("bbest = %Lf\n", bbest);
  printf("cbest = %Lf\n", cbest);
  printf("alphabest = %Lf\n", alphabest);
  printf("betabest = %Lf\n", betabest);
  printf("gbest = %Lf\n", gbest);
  printf("dbest = %Lf\n", dbest);
  printf("ebest = %Lf\n", ebest);
  printf("fbest = %Lf\n", fbest);

    
    gnuplot_commands();
   printf("The end\n");
}
