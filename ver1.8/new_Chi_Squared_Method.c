/************************************************************/
/* Usman Rizwan                                             */ 
/*                                                          */
/*                                                          */
/* A program to solve linear equations via gaussian         */
/* elimination. This program will incorporate scaling       */
/* scaling parameter as well. This program is currently     */
/* in working order.                                        */
/************************************************************/

#include "foo.h"

long double eMatrix[ROWS][COLS];
//long double b[ROWS][COLS];
long double xVector[ROWS], bVector[ROWS];
long double xi[40], yi[40],ei[40],xi0[40],yi0[40],ei0[40];
char line[132],*file_name,*out_file,*out_file2;
FILE *inp,*out,*fp1;
int N=3;
int n,count,count0, upper_limit=1,upper_limit_best=0;
long double a=0,b=0,c=0,g=0,alpha=0, beta=0,d,e,f;
long double aerror,berror,cerror,derror,eerror,gerror,fierror;
long double energy, energy0=0, min=100000,record_energy=0,abest=0,bbest=0,cbest=0,gbest=0,alphabest=0,betabest=0,dbest=0,ebest=0,fbest=0;
long double sum=0, sumx=0,sumx2=0,sumx3=0,sumx4=0,f=0,fx=0,fx2=0;
  long double sumj00=0, sumxj01=0,sumxj02=0,sumxj03=0,sumxj11=0,fj=0,fxj=0,fxj2=0,fxj3=0,sumxj12=0,sumxj13=0;
  long double sumxj22=0,sumxj23=0,sumxj33=0;

double formMatrix1();
double formMatrix2();

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
	 printf(" %Le\t ", eMatrix[k][l]);
       }
       printf("\n");
       //printf(" =%Lf\n", bVector[k]);
       //xVector[k]=0;
   }
}

int ReadCo60()
{
   FILE * fp;

  if((fp=fopen("Co60_abs_eff.dat","r"))==NULL)
    {
      //printf("Cannot open Co60_abs_eff.dat.\n");
      //exit(-1);
      return 0;
    }

  printf("Parameters read from Co60_abs_eff.dat\n");

    if(fgets(line,132,fp)!=NULL)
      {
	count=0;
	while(fscanf(fp,"%Lf %Lf %Lf\n",&xi[count],&yi[count],&ei[count])!=EOF)
	  {
	    //printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
	//printf("Count for Co-60 = %d\n\n\n", count);
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
	    //printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
	//printf("Count for Co-60 = %d\n\n\n", count);
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
	    //printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
	//printf("Count for Eu-152 = %d\n\n\n", count);
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
	    //printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
	//printf("Count for Eu-152 = %d\n\n\n", count);
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
	    
	    //printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
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
	    
	    //printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
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
	    //printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
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
	    //printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }
	fclose(fp);
      }
}

double formMatrix1(){  //This part is without the gamma function.

    int i;

    printf("Count formMatrix1 = %d\n", count);

    for(i=0;i<count;i++)
      if(xi[i]>=energy0)
      {
	sum+=1/(ei[i]*ei[i]);
	      	  sumx+=xi[i]/(ei[i]*ei[i]);
	      	  sumx2+=xi[i]*xi[i]/(ei[i]*ei[i]);
	      	  sumx3+=xi[i]*xi[i]*xi[i]/(ei[i]*ei[i]);
	      	  sumx4+=xi[i]*xi[i]*xi[i]*xi[i]/(ei[i]*ei[i]);
	      	  f+=yi[i]/(ei[i]*ei[i]);
	      	  fx+=yi[i]*xi[i]/(ei[i]*ei[i]);
	      	  fx2+=yi[i]*xi[i]*xi[i]/(ei[i]*ei[i]);
      }

       
}

double formMatrix2(){  //This is where the gamma function is include.

  int i;
  //sumj00=0;
  //energy0=xi0[upper_limit];
  for(i=0;i<count;i++)
      if(xi[i]<=energy0)
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

double MatrixCoefficients()
{ 

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

}


double SpecialFunctionj(int k, int l)
{
  int i;
  eMatrix[k][l]=0.;
  eMatrix[k][l+1]=0.;
  eMatrix[k][l+2]=0.;
  eMatrix[k][l+3]=0.;
  bVector[l]=0.;

  for(i=0;i<count;i++)
    if(xi[i]<=energy0){ 
      //ei[i]=1;
      eMatrix[k][l]+=1/(ei[i]*ei[i]);
      eMatrix[k+1][l]+=xi[i]/(ei[i]*ei[i]);
      eMatrix[k+2][l]+=(2*energy0*xi[i]-energy0*energy0)/(ei[i]*ei[i]);
      eMatrix[k+3][l]+=(xi[i]-energy0)*(xi[i]-energy0)/(ei[i]*ei[i]);
      bVector[l]+=yi[i]/(ei[i]*ei[i]);
    }
  
  for(i=0;i<count;i++)
    if(xi[i]>=energy0){ //This is for the side without gamma.
      //ei[i]=1;
      eMatrix[k][l]+=1/(ei[i]*ei[i]);
      eMatrix[k+1][l]+=xi[i]/(ei[i]*ei[i]);
      eMatrix[k+2][l]+=xi[i]*xi[i]/(ei[i]*ei[i]);
      eMatrix[k+3][l]+=0;
      bVector[l]+=yi[i]/(ei[i]*ei[i]);
    }
  
}

double formMatrix(){ 


  long double bMatrix[ROWS][COLS];
  int i,j;

  printf("enery0 = %Lf\n", powl(2.71828182845,energy0));

  ReadCo60();
  MatrixCoefficients();

  ReadEu152();
  SpecialFunctionj(0,4);
  MatrixCoefficients();

  ReadBa133();
  SpecialFunctionj(0,5);
  MatrixCoefficients();

  ReadCo56();
  SpecialFunctionj(0,6);
  MatrixCoefficients();

 

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

   //for(i=0;i<dim;i++)   printf("xVector[%d] = %Lf\n", i, xVector[i]);
	 

   a=xVector[0];
   b=xVector[1];
   c=xVector[2];
   g=xVector[3];
   d=xVector[4];
   e=xVector[5];
   f=xVector[6];

}

int Matrix()
{
  int k,jug;
  
  DescendingOrder();

  k=upper_limit;
 
  long double energy00=xi0[k];
  energy0=xi0[k];
  
  printf("energy0 = %Lf\n", (powl(2.71828182845,energy0)));
  jug=1;
  
  do{
    printf("energy0 = %Lf\n", (powl(2.71828182845,energy0)));
    
    Matrix_initialize();
    formMatrix();
    diagonal();
    SolveMatrix();
    ChiSquared();
    energy0=log(powl (2.71828182845,energy00)+1*jug);
    jug++;}while(energy0<log(2000)); 
}

int ChiSquared(){
  int i;
  long double chisquared=0,chisquaredl=0,chisquaredr=0;

  alpha = a - c*energy0*energy0 + g*energy0*energy0;
  beta = b + 2*c*energy0 - 2*g*energy0;
  
  ReadCo60();
  {

    printf("count : %d\n", count);
    printf("energy0 : %Lf\n", energy0);
    

  for(i=0;i<count;i++)
      if(xi[i]<=energy0)
    {
      chisquaredl+=((yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]))*(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]))/(ei[i]*ei[i]));
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residuall = %Lf\n", yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i])); */
      /* printf("residuall/error = %Lf\n", (yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]))/ei[i]); */
      /* printf("chisquaredl = %Lf\n", chisquaredl); */
      /* printf("LLLLLLLLLLLLLLLLLLLL\n"); */
    }

  for(i=0;i<count;i++)
      if(xi[i]>=energy0)
    {
      chisquaredr+=(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]))*(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]))/(ei[i]*ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residuall = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i])); */
      /* printf("residuall/error = %Lf\n", (yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]))/ei[i]); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
      /* printf("***************************\n"); */
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residualr = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i])); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
      
    }
}

   ReadEu152();
   {
     for(i=0;i<count;i++)
    if(xi[i]<=energy0)
    {
      chisquaredl+=((yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+d))*(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+d))/(ei[i]*ei[i]));
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residuall = %Lf\n", yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+d)); */
      /* printf("residuall/error = %Lf\n", (yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+d))/ei[i]); */
      /* printf("chisquaredl = %Lf\n", chisquaredl); */
      /* printf("LLLLLLLLLLLLLLLLLLLL\n"); */
    }

     for(i=0;i<count;i++)
      if(xi[i]>=energy0)
    {
      chisquaredr+=((yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+d))*(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+d))/(ei[i]*ei[i]));
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residuall = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+d)); */
      /* printf("residuall/error = %Lf\n", (yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+d))/ei[i]); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
      /* printf("***************************\n"); */
    }
  /* printf("chisquared = %Lf\n", chisquared); */
  
   }

  ReadCo56();
  {
    for(i=0;i<count;i++)
      if(xi[i]<=energy0)
    {
      chisquaredl+=((yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+f))*(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+f))/(ei[i]*ei[i]));
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residuall = %Lf\n", yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+f)); */
      /* printf("residuall/error = %Lf\n", (yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+f))/ei[i]); */
      /* printf("chisquaredl = %Lf\n", chisquaredl); */
      /* printf("LLLLLLLLLLLLLLLLLLLL\n"); */
    }

    for(i=0;i<count;i++)
      if(xi[i]>=energy0)
	{
      chisquaredr+=((yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+f))*(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+f)))/(ei[i]*ei[i]);
      /* printf("a = %Lf\n", a); */
      /* printf("b = %Lf\n", b); */
      /* printf("c = %Lf\n", c); */
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("yi[i] = %Lf\n", yi[i]); */
      /* printf("xi[i] = %Lf\n", xi[i]); */
      /* printf("fit = %Lf\n", (a + b*xi[i] + c*xi[i]*xi[i]+f)); */
      /* printf("residuall = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+f)); */
      /* printf("residuall/error = %Lf\n", (yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+f))/ei[i]); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
      /* printf("***************************\n"); */
      
    }
  /* printf("chisquared = %Lf\n", chisquared); */
  
  }

  ReadBa133();
  {
  
    for(i=0;i<count;i++)
      if(xi[i]<=energy0)
    {
      chisquaredl+=((yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+e))*(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+e))/(ei[i]*ei[i]));
       /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residuall = %Lf\n", yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+e)); */
      /* printf("residuall/error = %Lf\n", (yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+e))/ei[i]); */
      /* printf("chisquaredl = %Lf\n", chisquaredl); */
      /* printf("LLLLLLLLLLLLLLLLLLLL\n"); */
      
    }

    for(i=0;i<count;i++)
      if(xi[i]>=energy0)
    {
      chisquaredr+=((yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+e))*(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+e))/(ei[i]*ei[i]));
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residuall = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+e)); */
      /* printf("residuall/error = %Lf\n", (yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+e))/ei[i]); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
      /* printf("***************************\n"); */
    }
  /* printf("chisquared = %Lf\n", chisquared); */
  
  }

  chisquared=chisquaredl+chisquaredr;

  printf("chisquared = %10.6Lf\n", chisquared);
  printf("count0 = %d\n",count0);
  //printf("count0 = %d\n",dim);
  printf("chi^2/dof = %10.6Lf\n", chisquared);
  
  min=chisquared;

  /* if(chisquared<min) */
  /*   { */
  /*     min=chisquared; */
      abest=a;
      bbest=b;
      cbest=c;
      alphabest=alpha;
      betabest=beta;
      gbest=g;
      dbest=d;
      ebest=e;
      fbest=f;
  /*     //upper_limit_best=upper_limit; */
      record_energy=energy0;
  /*   } */
  
              char* output;
	      output = (char*) malloc(strlen(out_file+1));
	      strcpy(output,out_file);

	      fp1 = fopen(output,"a");
	      //fprintf(fp1,"%10.6Lf %10.6Lf\n", xi[upper_limit], chisquared);
	      fprintf(fp1,"%10.6Lf %10.6Lf\n", powl(2.71828182845,energy0), chisquared);
	      //printf("%10.6Lf %10.6Lf\n", powl(2.71828182845,energy0), chisquared);
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
      //b[i][j]=0.0;
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
  long double en=0,abseff=0,abserr=0,relerror=0;
  char* output;
	      
  output = (char*) malloc(strlen(out_file2+1));
  strcpy(output,out_file2);
  fp1 = fopen(output,"a");

  ReadCo60_nolog();
  for(i=0; i<count; i++)
    {
      abseff=yi[i];
      abserr=ei[i];
      fprintf(fp1,"%10.6Lf %10.6Lf %10.6Lf blue\n", xi[i], abseff, abserr);
    }

  ReadEu152_nolog();
  for(i=0; i<count; i++)
    {
      abseff=yi[i]/(powl(2.71828182845,dbest));
      relerror=derror;
      abserr=abseff*sqrt(  relerror * relerror + ( ei[i]/yi[i] ) *( ei[i]/yi[i] ) );
      fprintf(fp1,"%10.6Lf %10.6Lf %10.6Lf red\n",xi[i], abseff, abserr);
    }

  ReadCo56_nolog();
  for(i=0; i<count; i++)
    {
      abseff=yi[i]/(powl(2.71828182845,fbest));
      relerror=fierror;
      abserr=abseff*sqrt((relerror*relerror)+(ei[i]/yi[i])*(ei[i]/yi[i]));
      fprintf(fp1,"%10.6Lf %10.6Lf %10.6Lf green\n", xi[i], abseff, abserr);
    }
		    
  ReadBa133_nolog();
  for(i=0; i<count; i++)
    {
      abseff=yi[i]/(powl(2.71828182845,ebest));
      relerror=eerror;
      abserr=abseff*sqrt((relerror*relerror)+(ei[i]/yi[i])*(ei[i]/yi[i]));
      fprintf(fp1,"%10.6Lf %10.6Lf %10.6Lf pink\n", xi[i], abseff, abserr);
    }
	      
  //fprintf(fp1,"%10.6Lf %10.6Lf\n", xi[upper_limit], chisquared);
  //fprintf(fp1,"%10.6Lf %10.6Lf\n", powl(2.71828182845,energy0), chisquared);
  //printf("%10.6Lf %10.6Lf\n", powl(2.71828182845,energy0), chisquared);
  fclose(fp1);
  free(output);



}

double residual(){
  
  int i;
  long double residual;
  FILE *fp2;

  fp2 = fopen("residual.dat","a");

  alpha = abest - cbest*energy0*energy0 + gbest*energy0*energy0;
  beta = bbest + 2*cbest*energy0 - 2*gbest*energy0;
  
  ReadCo60_nolog();
  {

    for(i=0;i<count;i++)
      if(xi[i]<=energy0)
    {
      printf("count : %d\n", count);
      
      residual=(yi[i]-(powl(2.71828182845,alpha)*powl(xi[i],beta)*powl(xi[i],gbest*log(xi[i]))));
      fprintf(fp2,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], residual, ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      printf("residuall = %Lf\n", yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i])); 
    }

    for(i=0;i<count;i++)
      if(xi[i]>=energy0)
    {
      residual=(yi[i]-(powl(2.71828182845,abest)*powl(xi[i],bbest)*powl(xi[i],cbest*log(xi[i]))));
      fprintf(fp2,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], residual, ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residualr = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i])); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
    }
  }

   ReadEu152_nolog();
  {
  
  
    for(i=0;i<count;i++)
      if(xi[i]<=energy0)
    {
      residual=(yi[i]-(powl(2.71828182845,alpha+dbest)*powl(xi[i],beta)*powl(xi[i],gbest*log(xi[i]))));
      fprintf(fp2,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], residual,ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /*  printf("residuall = %Lf\n", yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+d)); */
      /*  printf("chisquaredl = %Lf\n", chisquaredl); */
    }

    for(i=0;i<count;i++)
      if(xi[i]>=energy0)
    {
      residual=(yi[i]-(powl(2.71828182845,abest+dbest)*powl(xi[i],bbest)*powl(xi[i],cbest*log(xi[i]))));
      fprintf(fp2,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], residual,ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residualr = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+d)); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
    }
  /* printf("chisquared = %Lf\n", chisquared); */
  
  }

  ReadCo56_nolog();
  {
 
  
    for(i=0;i<count;i++)
      if(xi[i]<=energy0)
    {
      residual=(yi[i]-(powl(2.71828182845,alpha+fbest)*powl(xi[i],beta)*powl(xi[i],gbest*log(xi[i]))));
      fprintf(fp2,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], residual,ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /*  printf("residuall = %Lf\n", yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+f)); */
      /*  printf("chisquaredl = %Lf\n", chisquaredl); */
    }

    for(i=0;i<count;i++)
      if(xi[i]>=energy0)
    {
      residual=(yi[i]-(powl(2.71828182845,abest+fbest)*powl(xi[i],bbest)*powl(xi[i],cbest*log(xi[i]))));
      fprintf(fp2,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], residual,ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residualr = %Lf\n", (yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+f))); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
    }
  /* printf("chisquared = %Lf\n", chisquared); */
  
  }

  ReadBa133_nolog();
  {
    
    for(i=0;i<count;i++)
      if(xi[i]<=energy0)
    {
      residual=(yi[i]-(powl(2.71828182845,alpha+ebest)*powl(xi[i],beta)*powl(xi[i],gbest*log(xi[i]))));
      fprintf(fp2,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], residual,ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residuall = %Lf\n",(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]+e))); */
      /* printf("chisquaredl = %Lf\n", chisquaredl); */
    }

    for(i=0;i<count;i++)
      if(xi[i]>=energy0)
    {
      residual=(yi[i]-(powl(2.71828182845,abest+ebest)*powl(xi[i],bbest)*powl(xi[i],cbest*log(xi[i]))));
      fprintf(fp2,"%10.6Lf %10.6Lf %10.6Lf\n", xi[i], residual,ei[i]);
      /* printf("ei[%d] = %Lf\n", i, ei[i]); */
      /* printf("residualr = %Lf\n", yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]+e)); */
      /* printf("chisquaredr = %Lf\n", chisquaredr); */
    }
  /* printf("chisquared = %Lf\n", chisquared); */
  
  }
  fclose(fp2);
/* FILE *pipe = popen("gnuplot","w"); */
/* //fprintf(pipe, "set data style lines\n"); */
/* //fprintf(pipe, "plot '%s' using 1:2\n",file_name); */
/*  fprintf(pipe, "set xrange[50:5000]\n"); */
/*  //fprintf(pipe, "set yrange[1e-4:0.01]\n"); */
/*  fprintf(pipe, "unset key\n"); */
/*  fprintf(pipe, "set xtics (50,100,1000,5000)\n"); */
/*  fprintf(pipe, "set logscale\n"); */
/*  fprintf(pipe, "set xlabel \"Energy [keV]\"\n"); */
/*  fprintf(pipe, "set ylabel \"Residual\"\n"); */
/*  //fprintf(pipe, "plot 'Usmantest_scaledeff.dat' using 1:2:3 with yerrorbars, x > 0.001 ? (x**(%Lf*log(x)))*(x**(%Lf))*(exp%Lf) : 1/0, x < %Lf ? (x**(%Lf*log(x)))*(x**(%Lf))*(exp%Lf) : 1/0\n", cbest, bbest, abest, (powl(2.71828182845,energy0)+100), gbest, betabest, alphabest); */


}


int gnuplot_commands()
{
FILE *pipe = popen("gnuplot","w");
//fprintf(pipe, "set data style lines\n");
 fprintf(pipe, "set terminal postscript eps color enhanced 'Helvetica' 20\n");
 fprintf(pipe, "set output 'Best_fit.eps\n");
 //fprintf(pipe, "plot '%s' using 1:2\n",file_name);
 fprintf(pipe, "set xrange[40:5000]\n");
 fprintf(pipe, "set yrange[5e-4:0.015]\n");
 fprintf(pipe, "set key top\n");
 fprintf(pipe, "set xtics (50,100,1000,5000)\n");
 fprintf(pipe, "set logscale\n");
 fprintf(pipe, "set format y '%%.0e'\n");
 fprintf(pipe, "set xlabel \"Energy [keV]\"\n");
 fprintf(pipe, "set ylabel \"Absolute Efficiency\n");
 fprintf(pipe, "plot 'Usmantest_scaledeff.dat' every ::0::1 using 1:2:3 title '{/Arial=12 ^{60}Co}' with yerrorbars pointtype 3 ps 1 linecolor rgb 'red', 'Usmantest_scaledeff.dat' every ::2::12 using 1:2:3 title '{/Arial=12 ^{152}Eu}' with yerrorbars pointtype 13 ps 1 linecolor rgb 'blue', 'Usmantest_scaledeff.dat' every ::13::24 using 1:2:3 title '{/Arial=12 ^{56}Co}' with yerrorbars pointtype 9 ps 1 linecolor rgb 'green', 'Usmantest_scaledeff.dat' every ::25::30 using 1:2:3 title '{/Arial=12 ^{133}Ba}' with yerrorbars pointtype 11 ps 0.8 linecolor rgb 'violet', x < %Le ? exp(%10.6Le)*x**(%10.6Le)*x**(%10.6Le*log(x))*exp(-sqrt(((%10.6Le+%10.6Le*log(x)+%10.6Le*(2*log(x)*log(%10.6Le)-(log(%10.6Le))**2)+%10.6Le*(log(x)-log(%10.6Le))**2+%10.6Le*(log(x))**2+%10.6Le*(log(x))*(2*log(x)*log(%10.6Le)-(log(%10.6Le))**2)+%10.6Le*(2*log(x)*log(%10.6Le)-(log(%10.6Le))**2)**2+%10.6Le*((log(x)-log(%10.6Le))**2)*(2*log(x)*log(%10.6Le)-(log(%10.6Le))**2)+%10.6Le*(log(x)-log(%10.6Le))**4+%10.6Le*log(x)*((log(x)-log(%10.6Le))**2))*%10.6Le)*1.710882)) : 1/0 notitle lt 1 lc 7 lw 1, x < %Le ? exp(%10.6Le)*x**(%10.6Le)*x**(%10.6Le*log(x))*exp(sqrt(((%10.6Le+%10.6Le*log(x)+%10.6Le*(2*log(x)*log(%10.6Le)-(log(%10.6Le))**2)+%10.6Le*(log(x)-log(%10.6Le))**2+%10.6Le*(log(x))**2+%10.6Le*(log(x))*(2*log(x)*log(%10.6Le)-(log(%10.6Le))**2)+%10.6Le*(2*log(x)*log(%10.6Le)-(log(%10.6Le))**2)**2+%10.6Le*((log(x)-log(%10.6Le))**2)*(2*log(x)*log(%10.6Le)-(log(%10.6Le))**2)+%10.6Le*(log(x)-log(%10.6Le))**4+%10.6Le*log(x)*((log(x)-log(%10.6Le))**2))*%10.6Le)*1.710882)) : 1/0 notitle lt 1 lc 7 lw 1, x > %Le ? (x**(%Le*log(x)))*(x**(%Le))*(exp(%Le)) : 1/0  title '{/Arial=12 Best-fit}' lt 1 lc 1 lw 1, x < %Le ? (x**(%Le*log(x)))*(x**(%Le))*(exp(%Le)) : 1/0  notitle lt 1 lc 1 lw 1, x > %Le ? exp(%10.6Le+%10.6Le*log(x)+%10.6Le*(log(x))**2)*exp(sqrt(((%10.6Le+%10.6Le*log(x)+%10.6Le*(log(x))**2+%10.6Le*(log(x))**3+%10.6Le*(log(x))**4)*%10.6Le)*1.710882)) : 1/0 notitle lt 1 lc 7 lw 1, x > %Le ? exp(%10.6Le+%10.6Le*log(x)+%10.6Le*(log(x))**2)*exp(-sqrt(((%10.6Le+%10.6Le*log(x)+%10.6Le*(log(x))**2+%10.6Le*(log(x))**3+%10.6Le*(log(x))**4)*%10.6Le)*1.710882)) : 1/0 title '{/Arial=12 90%% confidence interval}' lt 1 lc 7 lw 1\n",energy0+5, alphabest, betabest, gbest, inverse[0][0], inverse[0][1]+inverse[1][0], inverse[0][2]+inverse[2][0],energy0, energy0, inverse[0][3]+inverse[3][0], energy0, inverse[1][1], inverse[1][2]+inverse[2][1],energy0, energy0, inverse[2][2], energy0, energy0, inverse[2][3]+inverse[3][2], energy0,energy0, energy0, inverse[3][3], energy0, inverse[1][3]+inverse[3][1], energy0, min/(count0-dim),energy0+5, alphabest, betabest, gbest, inverse[0][0], inverse[0][1]+inverse[1][0], inverse[0][2]+inverse[2][0],energy0, energy0, inverse[0][3]+inverse[3][0], energy0, inverse[1][1], inverse[1][2]+inverse[2][1],energy0, energy0, inverse[2][2], energy0, energy0, inverse[2][3]+inverse[3][2], energy0,energy0, energy0, inverse[3][3], energy0, inverse[1][3]+inverse[3][1], energy0, min/(count0-dim), energy0-5, cbest, bbest, abest, energy0+5, gbest, betabest, alphabest, energy0-10, abest, bbest, cbest, inverse[0][0], inverse[0][1]+inverse[1][0], inverse[0][2]+inverse[1][1]+inverse[2][0], inverse[2][1]+inverse[1][2], inverse[2][2],min/(count0-dim), energy0-5, abest, bbest, cbest, inverse[0][0], inverse[0][1]+inverse[1][0], inverse[0][2]+inverse[1][1]+inverse[2][0], inverse[2][1]+inverse[1][2], inverse[2][2],min/(count0-dim));

 printf("Press enter to continue...\n");
 close(pipe); 
 fflush(pipe);
 
}

double error(){

  aerror=sqrt(min*inverse[0][0]/(count0-dim));
  berror=sqrt(min*inverse[1][1]/(count0-dim));
  cerror=sqrt(min*inverse[2][2]/(count0-dim));
  gerror=sqrt(min*inverse[3][3]/(count0-dim));
  derror=sqrt(min*inverse[4][4]/(count0-dim));
  eerror=sqrt(min*inverse[5][5]/(count0-dim));
  fierror=sqrt(min*inverse[6][6]/(count0-dim)); 
}


int SolveMatrix2(long double qMatrix[ROWS][COLS], long double qVector[ROWS], int d){
  int i=0,j=0, k=0;
  
   //process rows
   for(k=0; k<d; k++){
       for(i=k+1; i<d; i++){
	   if(qMatrix[k][k]==0){
	       printf("\nSolution does not exist.\n");
	       return;
	   }
	   long double M = qMatrix[i][k] / qMatrix[k][k];
	   for(j=k; j<d; j++){
	       qMatrix[i][j] = qMatrix[i][j] - M * qMatrix[k][j];
	   }
	   qVector[i] = qVector[i] - M * qVector[k];
       }
   }

   //this loop is for backward subtraction
   for(i=(d-1); i>=0; i--){
       long double s = 0;
       for(j = (i+1); j<=d; j++){
	   s = s+qMatrix[i][j]*xVector[j];
       }
       xVector[i] = (qVector[i] - s) / qMatrix[i][i];
   }

    printf("\n\n\n");

   for(i=0;i<d;i++){
    for(j=0;j<d;j++){
      printf("%Le\t", qMatrix[i][j]);
    }
    printf("%Le\n", qVector[i]);
  }

   for(i=0;i<d;i++)   printf("xVector[%d] = %Le\n", i, xVector[i]);
   long double a=xVector[0];
   long double b=xVector[1];
   long double c=xVector[2];

   FILE *pipe = popen("gnuplot","w");
   fprintf(pipe, "set terminal postscript eps color enhanced 'Helvetica' 20\n");
   fprintf(pipe, "set output 'Chisquared_bestfit.eps\n");
   fprintf(pipe, "set xrange[90:300]\n");
   //fprintf(pipe, "set yrange[-10:60]\n");
   fprintf(pipe, "set xtics add ('120' 120); set xtics add ('230' 230) \n");
   fprintf(pipe, "set xlabel \"E_{0} [keV]\"\n");
   fprintf(pipe, "set ylabel \"{/Symbol c}^{2}_{min.}\"\n");
   fprintf(pipe, "set label '{/Italic c}_{0}=3.16x10^{2}' at 150,3000\n");
   fprintf(pipe, "set label '{/Italic c}_{1}=-4.67' at 150,2800\n");
   fprintf(pipe, "set label '{/Italic c}_{2}=-2.23x10^{-2}' at 150,2600\n");
   fprintf(pipe, "set label '{/Italic c}_{4}=-3.26x10^{-5}' at 150,2400\n");
   fprintf(pipe, "set arrow from 120,0 to 120, %10.6Le+%10.6Le*120+%10.6Le*120*120+%10.6Le*120*120*120 nohead lc 7 lw 1\n", xVector[0], xVector[1], xVector[2], xVector[3]);
   fprintf(pipe, "set arrow from 230,0 to 230, %10.6Le+%10.6Le*230+%10.6Le*230*230+%10.6Le*230*230*230 nohead lc 7 lw 1\n", xVector[0], xVector[1], xVector[2], xVector[3]);
   fprintf(pipe, "plot 'Usmantest.dat' using 1:2 title '{/Symbol c}^{2}', %10.6Le+%10.6Le*x+%10.6Le*x*x+%10.6Le*x*x*x title 'Best-fit' lt 3\n", xVector[0], xVector[1], xVector[2], xVector[3]);
   //fprintf(pipe, "plot 'Usmantest.dat' using 1:2, %Lf+%Lf*x+%Lf*x*x\n", qVector[0], qVector[1],qVector[2]);
   /* fprintf(pipe, "set terminal pdf\n"); */
   /* fprintf(pipe, "set output \"Chisquared_bestfit.pdf\"\n"); */
   /* fprintf(pipe,"replot\n"); */
   /* fprintf(pipe,"set terminal x11\n"); */
   printf("Press enter to continue...\n");
   close(pipe); 
   fflush(pipe);
     
}

int main(int argc, char *argv[]){

  int i,j;
  long double bMatrix[ROWS][COLS];
  long double energy, chi;
  long double qMatrix[ROWS][COLS], qVector[ROWS];

  if(argc!=4)
    {
      printf("\n HelloWorld Filename outfilename_for_CHISquared outfilename_for_scaled_eff  \n");
      //exit(-1);
      return 0;
    }

  file_name=argv[1];
  out_file=argv[2];
  out_file2=argv[3];

  remove(out_file);
  remove(out_file2);

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
    
    

    fclose(inp);

    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < COLS; j++) {
            eMatrix[i][j] = 0.0;
	    bVector[i]=0.0;
	    xVector[i]=0.0;
	    qMatrix[i][j]=0;
	    qVector[i]=0.0;
        }
    }
    Matrix();
   
    Matrix_initialize();

    //printf("********record energy = %Lf************\n", powl(2.71828182845,record_energy));


    /*This is to fit a quadratic fit on the Chusquared vs. Energy and figuring out the minimum*/

    if((inp=fopen(argv[2],"r"))==NULL)
    {
      printf("Cannot open %s.\n",out_file);
      //exit(-1);
      return 0;
    }

  printf("Parameters read from %s.\n",out_file);

  long double sum_=0,sumx_=0,sumx2_=0, sumx3_=0, sumx4_=0, f_=0, fx_=0,fx2_=0;
  long double sumx5_=0, sumx6_=0, fx3_=0;

    if(fgets(line,132,inp)!=NULL)
      {
	while(fscanf(inp,"%Lf %Lf\n",&energy,&chi)!=EOF)
	  {
	    if(energy>120 && energy<230){
	      sum_+=1; //(ei[i]*ei[i]);
	      sumx_+=(energy); //(ei[i]*ei[i]);
	      sumx2_+=(energy)*(energy); ///(ei[i]*ei[i]);
	      sumx3_+=(energy)*(energy)*(energy); //(ei[i]*ei[i]);
	      sumx4_+=(energy)*(energy)*(energy)*(energy);//(ei[i]*ei[i]);
	      sumx5_+=(energy)*(energy)*(energy)*(energy)*(energy);
	      sumx6_+=(energy)*(energy)*(energy)*(energy)*(energy)*(energy);
	      f_+=chi; //(ei[i]*ei[i]);
	      fx_+=chi*(energy); //(ei[i]*ei[i]);
	      fx2_+=chi*(energy)*(energy); //(ei[i]*ei[i]);
	      fx3_+=chi*(energy)*(energy)*(energy);
	      //printf("energy =%Lf\n", energy);
	      //printf("sum =%Lf\n", sum);
	      
	  }
	    //printf("sum =%Lf\n", sum); 

	  }
      }

	fclose(inp);

	/* qMatrix[0][0]=sum_; */
	/* qMatrix[0][1]=sumx_; */
	/* qMatrix[0][2]=sumx2_; */
	/* qMatrix[1][0]=qMatrix[0][1]; */
	/* qMatrix[1][1]=qMatrix[0][2]; */
	/* qMatrix[1][2]=sumx3_; */
	/* qMatrix[2][0]=qMatrix[0][2]; */
	/* qMatrix[2][1]=qMatrix[1][2]; */
	/* qMatrix[2][2]=sumx4_; */

	/* qVector[0]=f_; */
	/* qVector[1]=fx_; */
	/* qVector[2]=fx2_; */

	qMatrix[0][0]=sum_;
	qMatrix[0][1]=sumx_;
	qMatrix[0][2]=sumx2_;
	qMatrix[0][3]=sumx3_;
	qMatrix[1][0]=qMatrix[0][1];
	qMatrix[1][1]=qMatrix[0][2];
	qMatrix[1][2]=qMatrix[0][3];
	qMatrix[1][3]=sumx4_;
	qMatrix[2][0]=qMatrix[0][2];
	qMatrix[2][1]=qMatrix[1][2];
	qMatrix[2][2]=qMatrix[1][3];
	qMatrix[2][3]=sumx5_;
	qMatrix[3][0]=qMatrix[2][1];
	qMatrix[3][1]=qMatrix[2][2];
	qMatrix[3][2]=qMatrix[2][3];
	qMatrix[3][3]=sumx6_;

	qVector[0]=f_;
	qVector[1]=fx_;
	qVector[2]=fx2_;
	qVector[3]=fx3_;
       
	for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      printf("%Le\t", qMatrix[i][j]);
    }
    printf("%Le\n", qVector[i]);
  }

	SolveMatrix2(qMatrix, qVector, 4);	

	for(j=0;j<4;j++)
	  {
	    printf("xVector[%d] = %Lf\n",j, xVector[j]);
	  }

	energy0=(-2*xVector[2]+sqrt(2*2*xVector[2]*xVector[2]-4*3*xVector[3]*xVector[1]))/(2*3*xVector[3]);
	//energy0=-xVector[1]/(2*xVector[2]);
	printf("min energy = %Lf\n", energy0);
	

	/*This is to fit a quadratic fit on the Chusquared vs. Energy and figuring out the minimum*/

	/*This part is to figure out the best fit parameters at the energy where chisquared is minimum*/
	energy0=log(energy0);
	Matrix_initialize();
	formMatrix();
	printmatrix();

	/*This part is to figure out the covariance matrix */
	int k= dim;
	for(i=0;i<dim;i++){
	for(j=0;j<dim;j++){
	  bMatrix[i][j]=eMatrix[i][j];
	}
	}
	cofactor(bMatrix, k);
	
	/*This part is to figure out the covariance matrix */
	diagonal();
	SolveMatrix();
	ChiSquared();
	error();
	/*This part is to figure out the best fit parameters at the energy where chisquared is minimum*/
	
	/*In this step, the d,e and f values are used to scale the relative efficiencies.*/
	ScaleEfficiencies();
	residual();
	/*In this step, the d,e and f values are used to scale the relative efficiencies.*/

	for(i=0;i<dim;i++)   printf("xVector[%d] = %Lf\n", i, xVector[i]);

	energy0=powl(2.71828182845,energy0);

	printf("min Chisquared = %Lf\n", min/(count0 - dim));
	printf("min energy = %Lf\n", powl(2.71828182845,record_energy));
   
	printf("abest = %Lf(%Lf)\n", abest,aerror);
	printf("bbest = %Lf(%Lf)\n", bbest,berror);
	printf("cbest = %Lf(%Lf)\n", cbest,cerror);
	printf("gbest = %Lf(%Lf)\n", gbest,gerror);
	printf("dbest = %Lf(%Lf)\n", dbest,derror);
	printf("ebest = %Lf(%Lf)\n", ebest,eerror);
	printf("fbest = %Lf(%Lf)\n", fbest,fierror);

	printf("alphabest = %Lf\n", alphabest);
	printf("betabest = %Lf\n", betabest);
	gnuplot_commands();
	printf("The end\n");
}
