#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h> 

#define ROWS 20
#define COLS 20

/*For calculating Determinant of the Matrix */
double determinant(long double bMatrix[ROWS][COLS],int d)
{
  long double s=1,det=0,b[ROWS][COLS];
  int i,j,m,n,c,l,k;

  if (d==1)
    {
     return (bMatrix[0][0]);
    }
  else
    {
     det=0;
     for (c=0;c<d;c++)
       {
        m=0;
        n=0;
        for (i=0;i<d;i++)
          {
            for (j=0;j<d;j++)
              {
//printf("i = %d\t j = %d\t c = %d\n", i,j,c);
                b[i][j]=0;
                if (i != 0 && j != c)
                 {
                   b[m][n]=bMatrix[i][j];
//printf("bMatrix[%d][%d] = %Lf\n", i,j,bMatrix[i][j]);
                   if (n<(d-2))
                    n++;
                   else
                    {
                     n=0;
                     m++;
                     }
                   }
	      }
	     /* printf("Press any button to continue..."); */
	  }
          det=det + s * (bMatrix[0][c] * determinant(b,d-1));
          s=-1 * s;
       }
    }
  //printf("det = %Le\n", det);
  //getchar();
  return (det);
}


/*Finding transpose of matrix*/ 
void transpose(long double num[ROWS][COLS],long double fac[ROWS][COLS],long double r)
{
  int i,j;
  long double b[ROWS][COLS],inverse[ROWS][COLS],d;
 
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
         b[i][j]=fac[j][i];
        }
    }
  d=determinant(num,r);
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
        inverse[i][j]=b[i][j] / d;
        }
    }
   printf("\n\n\nThe inverse of matrix is : \n");
 
   for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
         printf("\t%Lf",inverse[i][j]);
        }
    printf("\n");
     }
}

 
void cofactor(long double num[ROWS][COLS],long double f)
{
 long double b[ROWS][COLS],fac[ROWS][COLS];
int p,q,m,n,i,j,l,k;
int dim=7;
 for (q=0;q<f;q++)
 {
   for (p=0;p<f;p++)
    {
     m=0;
     n=0;
     for (i=0;i<f;i++)
     {
       for (j=0;j<f;j++)
        {
          if (i != q && j != p)
          {
            b[m][n]=num[i][j];
            if (n<(f-2))
             n++;
            else
             {
               n=0;
               m++;
               }
            }
        }
      }
      fac[q][p]=pow(-1,q + p) * determinant(b,f-1);
    }
  }

  transpose(num,fac,f);
}

int main()
{
  long double eMatrix[ROWS][COLS],d;
  //long double a,b,c,g,d,e,f;
  int dim;
  int i,j;
  FILE * fp;

  printf("-------------------------------------------------------------\n");
  printf("----------------made by C code champ ------------------------\n");
  printf("-------------------------------------------------------------\n");
  printf("\n  C Program to find inverse of Matrix\n\n");
  printf("Enter the order of the Matrix : ");http://www.cprogramming.com/tutorial/c/lesson6.html
  scanf("%d",&dim);
  printf("Enter the elements of %dX%d Matrix : \n",dim,dim);

  if((fp=fopen("Matrix.txt","r"))==NULL)
    {
      printf("Cannot open Matrix.txt.\n");
      //exit(-1);
      return 0;
    }

for(i = 0; i < dim; i++)
  {
      for(j = 0; j < dim; j++) 
      {
  //Use lf format specifier, %c is for character
       if (!fscanf(fp, "%Lf", &eMatrix[i][j])) 
           break;
      // mat[i][j] -= '0'; 
       printf("%Lf\n",eMatrix[i][j]); //Use lf format specifier, \n is for new line
}
}

  d=determinant(eMatrix,dim);
  printf("Determinant of the Matrix = %Lf",d);
  if (d==0)
   printf("\nInverse of Entered Matrix is not possible\n");
  else
   cofactor(eMatrix,dim);
   printf("\n\n**** Thanks for using the program!!! ****\n");
}
 

