#include "foo.h"

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
void transpose(long double num[ROWS][COLS],long double fac[ROWS][COLS],int r)
{
  int i,j;
  long double b[ROWS][COLS],d;
 
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
	  printf("\t%Le",inverse[i][j]);
        }
      printf("\n");
    }
  error(inverse);

}


void cofactor(long double num[ROWS][COLS],int f)
{
  long double b[ROWS][COLS],fac[ROWS][COLS];
  int p,q,m,n,i,j,l,k;

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

