#include "foo.h"
#define ROW_SIZE 40
#define COL_SIZE 40

/*For matrix multiplication*/

long double Matrix_multiplication(long double Matrix_1[ROW_SIZE][COL_SIZE], int dimR_1, int dimC_1, long double Matrix_2[COL_SIZE][ROW_SIZE], int dimR_2, int dimC_2){

  int a1,a2,a;
  int b1;

  /* printf("dimR_1 = %d\n", dimR_1 ); */
  /* printf("dimC_1 = %d\n", dimC_1 ); */
  /* printf("dimR_2 = %d\n", dimR_2 ); */
  /* printf("dimC_2 = %d\n", dimC_2 ); */

  /* for(a1=0;a1<dimR_2;a1++){ */
  /*   for(b1=0;b1<dimC_2;b1++){ */
  /*     //printf("Matrix_2[%d][%d] = %Lf\t", a1, b1, Matrix_2[a1][b1]); */
  /*     printf("%Lf\t", Matrix_2[a1][b1]); */
  /*     //getchar(); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  /* getchar(); */

  /* for(a1=0;a1<dimR_1;a1++){ */
  /*   for(b1=0;b1<dimC_1;b1++){ */
  /*     //printf("Matrix_1[%d][%d] = %Lf\t", a1, b1, Matrix_1[a1][b1]); */
  /*     printf("%Lf\t", Matrix_1[a1][b1]); */
  /*     //getchar(); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  /* printf("\nI am inside the matrix\n\n"); */
  /* getchar(); */

  

  //initialize
  for(a1=0; a1<ROWS; a1++){
    for(b1=0; b1<COLS; b1++){
      Multiplication_answer[a1][b1] = 0;
    }
  }

  if( dimC_1 - dimR_2 == 0 )
    { 
      for(a1=0; a1<dimR_1; a1++){
	for(b1=0; b1<dimC_2; b1++){
	  for(a=0; a<dimC_1; a++){
	    Multiplication_answer[a1][b1] += Matrix_1[a1][a]*Matrix_2[a][b1];
	    /* printf("Matrix_1[%d][%d] = %10.6Le\n", a1, a, Matrix_1[a1][a] ); */
	    /* printf("Matrix_1[%d][%d] = %10.6Le\n", a, b1, Matrix_2[a][b1] ); */
	  }
	}
      }

      //prints the answer
      for(a1=0;a1<dimR_1;a1++){
      	for(b1=0;b1<dimC_2;b1++){
      	  printf("%10.6Le\t", Multiplication_answer[a1][b1]);
      	}
      	printf("\n");
      }

    }

  else {
    printf("Not Valid\n");
    getchar();
    return 0;
  }

  printf("\nI am inside the matrix\n\n");
  error( Multiplication_answer );

  //return Matrix_answer[ROWS][COLS];

}

