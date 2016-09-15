#include <stdio.h>
#include <string.h>  //new header file

int main(){

  double d =123456.7891234;
  char name[20];
  int ar[3][3];
  char in[3]={'x', 'y', 'z'};
  char new[30];
  int a,b,c;

  for(a=0;a<2;a++)
    for(b=0;b<2;b++)
      ar[a][b]=5;

  for(b=0;b<2;b++){
    for(c=0;c<2;c++) {
	printf("%d \t", ar[b][c] );
}
    printf("\n");
  }

  for(a=0;a<2;a++){
    for(b=0;b<2;b++){
      for(c=0;c<2;c++){ 
	printf("%d*%c", ar[b][c], in[a] );
	if(c!=1) printf("+");
      }
      if(b!=1){
	printf("+");
      }
      printf("\n");
    }
    if (a !=1 )
      printf("+");
  }
  printf("\n");

  /* printf("Enter name: ");  */
  /* scanf("%s",name);  */
  /* printf("Your name is %s\n",name); */
  /* printf("%f*name\n", d); */
  return 0; 
}
