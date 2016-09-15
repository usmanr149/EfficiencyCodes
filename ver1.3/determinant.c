#include<stdio.h>
#include<math.h>

/*For calculating Determinant of the Matrix */
double determinant(long double eMatrix[25][25],int dim)
{
  long double s=1,det=0,b[25][25];
  int i,j,m,n,c;
  if (dim==1)
    {
     return (eMatrix[0][0]);
    }
  else
    {
     det=0;
     for (c=0;c<dim;c++)
       {
        m=0;
        n=0;
        for (i=0;i<dim;i++)
          {
            for (j=0;j<dim;j++)
              {
		printf("i = %d\t j = %d\t c = %d\n", i,j,c);
		getchar();
                b[i][j]=0;
                if (i != 0 && j != c)
                 {
                   b[m][n]=eMatrix[i][j];
		   printf("eMatrix[%d][%d] = %Lf\n", i,j,eMatrix[i][j]);
                   if (n<(dim-2))
                    n++;
                   else
                    {
                     n=0;
                     m++;
                     }
                   }
               }
             }
          det=det + s * (eMatrix[0][c] * determinant(b,dim-1));
          s=-1 * s;
          }
    }
 
    return (det);
}

/*Finding transpose of matrix*/ 
void transpose(long double num[25][25],long double fac[25][25],long double r)
{
  int i,j;
  long double b[25][25],inverse[25][25],d;
 
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

 
void cofactor(long double num[25][25],long double f)
{
 long double b[25][25],fac[25][25];
 int p,q,m,n,i,j;
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
  long double eMatrix[25][25],d;
  int dim;
  int i,j;
  printf("-------------------------------------------------------------\n");
  printf("----------------made by C code champ ------------------------\n");
  printf("-------------------------------------------------------------\n");
  printf("\n  C Program to find inverse of Matrix\n\n");
  printf("Enter the order of the Matrix : ");http://www.cprogramming.com/tutorial/c/lesson6.html
  scanf("%d",&dim);
  printf("Enter the elements of %dX%d Matrix : \n",dim,dim);
  for (i=0;i<dim;i++)
    {
     for (j=0;j<dim;j++)
       {
        scanf("%Lf",&eMatrix[i][j]);
        }
    }
  d=determinant(eMatrix,dim);
  printf("Determinant of the Matrix = %Lf",d);
  if (d==0)
   printf("\nInverse of Entered Matrix is not possible\n");
  else
   cofactor(eMatrix,dim);
   printf("\n\n**** Thanks for using the program!!! ****\n");
   //getch();
}
 

