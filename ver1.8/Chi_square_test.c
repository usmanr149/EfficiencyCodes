/************************************************************/
/* Usman Rizwan                                             */ 
/*                                                          */
/*                                                          */
/* A program to solve linear equations via gaussian         */
/* elimination.                                             */
/* Chi square analysis                                      */
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
long double xi[40], yi[40],ei[40];
char line[132],*file_name,*out_file;
FILE *inp,*out,*fp1;
int N=3;
int n,count, upper_limit=3;
long double a=0,b=0,c=0,g=0,alpha=0, beta=0;
long double energy;

int dim=4;  //declare matrix dimension here

int DescendingOrder(){
  int i,j;
  long double en=0,eff=0,err=0;

  printf("Count = %d\n\n\n", count);

  for(i=0;i<count;i++)
    for(j=i+1;j<count;j++)
      if (xi[i]>xi[j]){
	  en=xi[i];
	  eff=yi[i];
	  err=ei[i];
	  xi[i]=xi[j];
	  yi[i]=yi[j];
	  ei[i]=ei[j];
	  xi[j]=en;
	  yi[j]=eff;
	  ei[j]=err;
      }

  for(i=0;i<count;i++){
    printf("%Lf   %Lf   %Lf\n",xi[i],yi[i],ei[i]);
    }
  
}


int formmatrix()
{

  //FILE *inp;
  long double x,y,e;
  long double sum=0, sumx=0,sumx2=0,sumx3=0,sumx4=0,f=0,fx=0,fx2=0;
  long double sumj00=0, sumxj01=0,sumxj02=0,sumxj03=0,sumxj11=0,fj=0,fxj=0,fxj2=0,fxj3=0,sumxj12=0,sumxj13=0;
  long double sumxj22=0,sumxj23=0,sumxj33=0;
  int ag;
  int i,j,jug;
  int M;

  if(eMatrix==NULL || bVector==NULL || xVector==NULL){ 
       printf("\nNot enough memory to allocate for %d equations.\n", N);
       return(0);
   }

  printf("Parameters read from %s.\n",file_name);

    if(fgets(line,132,inp)!=NULL)
      {
	count=0;
	while(fscanf(inp,"%Lf %Lf %Lf\n",&xi[count],&yi[count],&ei[count])!=EOF)
	  {
	    //printf ("%d %Lf %Lf %Lf\n", count, xi[count],yi[count],ei[count]);
	    count++;
	      }

	      DescendingOrder();
	      
	      long double energy0=0;

	      for(upper_limit=3;upper_limit<(count-3);upper_limit++)
		{
		  energy00=xi[upper_limit];
		  jug=0;
		  while(energy0<xi[upper_limit+1])
		    {
		  energy0=energy00+log10(5)*jug;
		    {	     
	      for(i=upper_limit;i<count;i++)
	      	{
	      	  ei[i]=1;
	      	  //printf("xi[%d] = %Lf\n", i, xi[i]);
	      	  sum+=1/(ei[i]*ei[i]);
	      	  sumx+=xi[i]/(ei[i]*ei[i]);
	      	  sumx2+=xi[i]*xi[i]/(ei[i]*ei[i]);
	      	  sumx3+=xi[i]*xi[i]*xi[i]/(ei[i]*ei[i]);
	      	  sumx4+=xi[i]*xi[i]*xi[i]*xi[i]/(ei[i]*ei[i]);
	      	  f+=yi[i]/(ei[i]*ei[i]);
	      	  fx+=yi[i]*xi[i]/(ei[i]*ei[i]);
	      	  fx2+=yi[i]*xi[i]*xi[i]/(ei[i]*ei[i]);
	      	  //printf("%d :::: sumx2 = %Lf :::::%Lf\n", i, sumx2, xi[i]);
	      	}

	      printf("x0 = %Lf\n", xi[upper_limit]);

	      for(i=0;i<upper_limit;i++)
	      	{
	      ei[i]=1;
	      sumj00+=1/(ei[i]*ei[i]);;
	      sumxj01+=xi[i]/(ei[i]*ei[i]);
	      sumxj11+=xi[i]*xi[i]/(ei[i]*ei[i]);
	      sumxj02+=(2*xi[upper_limit]*xi[i]-xi[upper_limit]*xi[upper_limit])/(ei[i]*ei[i]);
	      sumxj03+=(xi[i]-xi[upper_limit])*(xi[i]-xi[upper_limit])/(ei[i]*ei[i]);
	      sumxj12+=xi[i]*(2*xi[upper_limit]*xi[i]-xi[upper_limit]*xi[upper_limit])/(ei[i]*ei[i]);
	      sumxj13+=xi[i]*(xi[i]-xi[upper_limit])*(xi[i]-xi[upper_limit])/(ei[i]*ei[i]);
	      sumxj22+=(2*xi[upper_limit]*xi[i]-xi[upper_limit]*xi[upper_limit])*(2*xi[upper_limit]*xi[i]-xi[upper_limit]*xi[upper_limit])/(ei[i]*ei[i]);
	      sumxj23+=(2*xi[upper_limit]*xi[i]-xi[upper_limit]*xi[upper_limit])*(xi[i]-xi[upper_limit])*(xi[i]-xi[upper_limit])/(ei[i]*ei[i]);
	      sumxj33+=(xi[i]-xi[upper_limit])*(xi[i]-xi[upper_limit])*(xi[i]-xi[upper_limit])*(xi[i]-xi[upper_limit])/(ei[i]*ei[i]);
	      fj+=yi[i]/(ei[i]*ei[i]);
	      fxj+=yi[i]*xi[i]/(ei[i]*ei[i]);
	      fxj2+=yi[i]*(2*xi[upper_limit]*xi[i]-xi[upper_limit]*xi[upper_limit])/(ei[i]*ei[i]);
	      fxj3+=yi[i]*(xi[i]-xi[upper_limit])*(xi[i]-xi[upper_limit])/(ei[i]*ei[i]);
	      //printf("%d :::: sumxj02 = %Lf \n", i, sumxj02);
	      //getchar();
	      	}
	     
	      eMatrix[0][0]=sum+sumj00;
	      eMatrix[0][1]=sumx+sumxj01;
	      eMatrix[0][2]=sumx2+sumxj02;
	      eMatrix[0][3]=sumxj03;
	      eMatrix[1][0]=eMatrix[0][1];
	      eMatrix[1][1]=sumx2+sumxj11;
	      eMatrix[1][2]=sumx3+sumxj12;
	      eMatrix[1][3]=sumxj13;
	      eMatrix[2][0]=eMatrix[0][2];
	      eMatrix[2][1]=eMatrix[1][2];
	      eMatrix[2][2]=sumx4+sumxj22;
	      eMatrix[2][3]=sumxj23;
	      eMatrix[3][0]=eMatrix[0][3];
	      eMatrix[3][1]=eMatrix[1][3];
	      eMatrix[3][2]=eMatrix[2][3];
	      eMatrix[3][3]=sumxj33;

	      bVector[0]=f+fj;
	      bVector[1]=fx+fxj;
	      bVector[2]=fx2+fxj2;
	      bVector[3]=fxj3;

	      energy=pow (10, xi[upper_limit]); 
	      //printf("xi[upper_limit] = %Lf : energy = %Lf\n", xi[upper_limit], energy);
	      //getchar();
	      //diagonal();
		SolveMatrix();
	      ChiSquared();
	      //gnuplot_commands();
		    }
		    jug++; }
	      //fprintf(fp1,"%.3Lf %Lf %Lf %Lf %Lf\n",energy,k0*log10(energy)+b0,fmin,fmax,(fmax-fmin)/2);

	      
	  

      }
    //printf("Hello, world!\n");
//return 0;
}

void printmatrix(){
   //print matrix "eMatrix"
   int k,l;
   printf("\n");
   for(k=0; k<dim; k++){
       for(l=0; l<dim; l++){
	 printf(" %Lf*X%d + ", eMatrix[k][l],l);
       }
       
       printf(" =%Lf\n", bVector[k]);
	xVector[k]=0;
   }
   printf("\n");
}

void diagonal(){
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
   //cls();
   //count_num_lines();
   //allocmatrix();
   //rewind(InFile);
   //read data from file
   //readmatrix();
   //check if there are 0 on main diagonal and exchange rows in that case

   //process rows
   for(k=0; k<dim; k++){
       for(i=k+1; i<dim; i++){
	   if(eMatrix[k][k]==0){
	       printf("\nSolution does not exist.\n");
	       return;
	   }
	   long double M = eMatrix[i][k] / eMatrix[k][k];
	   /* printf("eMatrix[%d][%d] = %Lf : eMatrix[%d][%d] = %Lf \n",i,k,k,k,eMatrix[i][k],eMatrix[k][k]); */
	   /* printf("M = %Lf\n", M); */
	   for(j=k; j<dim; j++){
	     /* printf("eMatrix[i][k] = %Lf \n eMatrix[k][j] = %Lf \n eMatrix[i][j] = %Lf\n",eMatrix[i][k],eMatrix[k][j],eMatrix[i][j]); */
	   /* printf("M = %Lf\n", M); */
	     //printf("eMatrix[i][j] before = %Lf \n",eMatrix[i][j]);
	       eMatrix[i][j] = eMatrix[i][j] - M * eMatrix[k][j];
	       //printf("eMatrix[i][j] after = %Lf \n",eMatrix[i][j]);
	       /* getchar(); */
	       /* printf("\n \n"); */
	   }
	   bVector[i] = bVector[i] - M * bVector[k];
       }
   }


   printmatrix();
   
   //this loop is for backward subtraction
   for(i=dim-1; i>=0; i--){
       long double s = 0;
       for(j = i; j<=dim; j++){
	 /* printf("s =  %Lf\n", s); */
	   s = s+eMatrix[i][j]*xVector[j];
	   /* printf("eMatrix[i][j] = %Lf : %d : %d \n",eMatrix[i][j], i,j);  */
	   /* printf("********************\n"); */
       }
       /* printf("s =  %Lf\n", s); */
       /* printf("bVector[i] = %Lf : %d \n",bVector[i], i);  */
       xVector[i] = (bVector[i] - s) / eMatrix[i][i];
       /* printf("xVector[i] = %Lf : %d \n",xVector[i], i); */
       /* printf("&&&&&&&&&&&&&&&&&&\n"); */
   }

   for(i=0;i<dim;i++){
       printf("xVector[%d] = %10.14Lf\n", i, xVector[i]);
   }

   printf("\n \n");

   a=xVector[0];
   b=xVector[1];
   c=xVector[2];
   g=xVector[3];
   
}

   
int ChiSquared(){
  int i;
  long double chisquared=0,chisquaredi=0,chisquaredj=0;

  alpha = a - c*xi[upper_limit]*xi[upper_limit] + g*xi[upper_limit]*xi[upper_limit];
  beta = b + 2*c*xi[upper_limit] - 2*g*xi[upper_limit];
  
  /* printf("a = %Lf\n", a); */
  /* printf("b = %Lf\n", b); */
  /* printf("c = %Lf\n", c); */
  /* printf("alpha = %Lf\n", alpha); */
  /* printf("beta = %Lf\n", beta); */
  /* printf("g = %Lf\n", g); */

  for(i=0;i<upper_limit;i++)
    {
      printf("x[i] : y[i] %Lf : %Lf\n",xi[i], yi[i]);
      chisquaredi+=(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]))*(yi[i]-(alpha + beta*xi[i] + g*xi[i]*xi[i]));
       /* printf("chisquaredi = %Lf\n", chisquaredi); */
    }

  for(i=upper_limit;i<count;i++)
    {
      chisquaredj+=(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]))*(yi[i]-(a + b*xi[i] + c*xi[i]*xi[i]));
    }

  chisquared=chisquaredi+chisquaredj;
  printf("chisquared = %Lf\n", chisquared);
  
              char* output;
	      output = (char*) malloc(strlen(out_file+1));
	      strcpy(output,out_file);

	      fp1 = fopen(output,"a");
	      //fprintf(fp1,"%10.6Lf %10.6Lf\n", xi[upper_limit], chisquared);
	      fprintf(fp1,"%10.6Lf %10.6Lf\n", energy, chisquared);

	      fclose(fp1);
	      free(output);

	      //getchar();
}

int gnuplot_commands()
{
FILE *pipe = popen("gnuplot","w");
//fprintf(pipe, "set data style lines\n");
//fprintf(pipe, "plot '%s' using 1:2\n",file_name);
 fprintf(pipe, "set logscale\n");
 fprintf(pipe, "plot 'Experiment.dat' using 1:2, (x**(%Lf*log10(x)))*(x**(%Lf))*(10**%Lf),(x**(%Lf*log10(x)))*(x**(%Lf))*(10**%Lf)\n",c,b,a,g,beta,alpha);
 getchar();
close(pipe); 
fflush(pipe);
}

int main(int argc, char *argv[]){

  int i,j;

  if(argc!=3)
    {
      printf("\n HelloWorld Filename outfilename \n");
      //exit(-1);
      return 0;
    }

  file_name=argv[1];
  out_file=argv[2];

  if((inp=fopen(argv[1],"r"))==NULL)
    {
      printf("Cannot open %s.\n",file_name);
      //exit(-1);
      return 0;
    }

    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < COLS; j++) {
            eMatrix[i][j] = 0.0;
	    bVector[i]=0.0;
	    xVector[i]=0.0;
        }
    }

    //system("gnuplot -p -e \"plot '/home/griffin/Data/GRIFFIN15/test.txt'\"");

    /* long double pi=3.141592653589793238462643383279502884197169399375105820974944592307816406286; */
    /* float pif=3141.592653589793238462643383279502884197169399375105820974944592307816406286; */
    /* printf(" pi = %10.20Lf\n", pi); */
    /*  printf(" pi = %10.20ef\n", pif); */
    
    //gnuplot_commands();
  formmatrix();
  /* printmatrix(); */
  /* diagonal(); */
  /* SolveMatrix(); */
  //ChiSquared();
  
}
