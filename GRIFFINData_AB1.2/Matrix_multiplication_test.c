#include <stdio.h>
#define ROWS 20
#define COLS 20
/*For matrix multiplication*/

int main(){

			       long double Matrix_1[ROWS][COLS];
			       long double Matrix_2[ROWS][COLS];
			       long double Matrix_answer[ROWS][COLS];
			       char line[132];
			       FILE * fp;
			       char c;
			       int a;
			       int b;
			       int a1; 
			       int b2;
			       int a2,b1;
			       int dim;
			       int dimR_1=0;
			       int dimC_1=1;
			       int dimR_2=0;
			       int dimC_2=1;
			       int ROWS1=0;
			       int COLS1=0;
			       long double zig;


			       //Read Matrix 1
			       if((fp=fopen("Matrix_1.txt","r"))==NULL)
				 {
				   printf("Cannot open Sudoku.txt.\n");
				   return 0;
				   //return 0;
				 }
			       
			       while( (c=fgetc(fp)) != EOF ){
				 if(!fscanf(fp,"%Lf", &zig));
				     printf("zig = %Lf\n", zig);
				     if(c != '\n' && dimR_1 == 0) dimC_1++;
				     if(c == '\n') dimR_1++;
				     //getchar();
			       }			       
			       
			       printf("dimR_1 = %d\n", dimR_1);
			       printf("dimC_1 = %d\n", dimC_1);
			       getchar();
			       
			       //fclear(fp);
			       rewind(fp);
			       for(a = 0 ; a < dimR_1 ; a++){
				 for(b=0 ; b < dimC_1 ; b++){
				     if(!fscanf(fp,"%Lf", &Matrix_1[a][b]));
				     //printf("zig = %Lf\n", Matrix_1[a][b]);
				 }
			       }
			       for(a = 0 ; a < dimR_1 ; a++){
				   for(b=0 ; b < dimC_1 ; b++){
				     printf("%Lf \t", Matrix_1[a][b]);
					    }
				       printf("\n");
				 }
			       fclose(fp);

			       //Read Matrix 2
			       if((fp=fopen("Matrix_2.txt","r"))==NULL)
				 {
				   printf("Cannot open Sudoku.txt.\n");
				   return 0;
				   //return 0;
				 }
			       
			       while( (c=fgetc(fp)) != EOF ){
				 if(!fscanf(fp,"%Lf", &zig));
				     printf("zig = %Lf\n", zig);
				     if(c != '\n' && dimR_2 == 0) dimC_2++;
				     if(c == '\n') dimR_2++;
				     //getchar();
			       }			       
			       
			       printf("dimR_2 = %d\n", dimR_2);
			       printf("dimC_2 = %d\n", dimC_2);
			       getchar();
			       
			       rewind(fp);
			       for(a = 0 ; a < dimR_2 ; a++){
				 for(b=0 ; b < dimC_2 ; b++){
				     if(!fscanf(fp,"%Lf", &Matrix_2[a][b]));
				     //printf("zig = %Lf\n", Matrix_2[a][b]);
				 }
			       }
			       for(a = 0 ; a < dimR_2 ; a++){
				   for(b=0 ; b < dimC_2 ; b++){
				     printf("%Lf \t", Matrix_2[a][b]);
					    }
				       printf("\n");
				 }
			       fclose(fp);
			      

			       printf("diff. = %d \n", dimC_1 - dimR_2);
			       //check and then
			       if( dimC_1 - dimR_2 == 0 )
				 { 
				   for(a1=0; a1<dimR_1; a1++){
				     for(b1=0; b1<dimC_2; b1++){
				       for(a=0; a<dimC_1; a++){
					 Matrix_answer[a1][b1] += Matrix_1[a1][a]*Matrix_2[a][b1];
					 printf("a1 = %d \n", a1);
					 printf("b1 = %d \n", b1);
					 printf("a = %d \n", a);
					 printf("Matrix[%d][%d] = %Lf \n", a1, b1, Matrix_1[a1][b1]);
					 getchar();
				       }
				     }
				   }
				   for(a1=0;a1<dimR_1;a1++){
				     for(b1=0;b1<dimC_2;b1++){
				       printf("%Lf\t", Matrix_answer[a1][b1]);
				     }
				     printf("\n");
				   }
				 }
			       else printf("Not Valid\n");							       			 
}
