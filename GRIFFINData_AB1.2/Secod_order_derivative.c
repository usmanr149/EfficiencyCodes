#define Extended_ROWS 40
#define Extended_COLS 40

#include "foo.h"

/*First-order derivative matrix*/

long double Second_order_der(long double Matrix_1[int dimR_1][int dimC_1], long double Matrix_2[dimR_2][dimC_2]){
  
  long double F_primeprime[40][10][40];
   
  ReadCo60();
   for(i=0; i<count;i++)
     {
       if( xi[i] <= energy0 ){
	           F_primeprime[0][0][i] = 0;
	 F_primeprime[0][1][i] = 0;
	 F_primeprime[0][2][i] = 0;
	 F_primeprime[0][3][i] = 0;
	 F_primeprime[0][4][i] = 0;
	 F_primeprime[0][5][i]= 0;
	 F_primeprime[0][6][i] = 0;
	 F_primeprime[0][7][i] = 0;
	           F_primeprime[1][0][i] = 0;
	 F_primeprime[1][1][i] = 0;
	 F_primeprime[1][2][i] = 0;
	 F_primeprime[1][3][i] = 0;
	 F_primeprime[1][4][i] = 0;
	 F_primeprime[1][5][i]= 0;
	 F_primeprime[1][6][i] = 0;
	 F_primeprime[1][7][i] = 0;
	           F_primeprime[2][0][i] = 0;
	 F_primeprime[2][1][i] = 0;
	 F_primeprime[2][2][i] = 0;
	 F_primeprime[2][3][i] = 0;
	 F_primeprime[2][4][i] = 0;
	 F_primeprime[2][5][i]= 0;
	 F_primeprime[2][6][i] = 0;
	 F_primeprime[2][7][i] = 2*(xi[i]-energy0);
	           F_primeprime[3][0][i] = 0;
	 F_primeprime[3][1][i] = 0;
	 F_primeprime[3][2][i] = 0;
	 F_primeprime[3][3][i] = 0;
	 F_primeprime[3][4][i] = 0;
	 F_primeprime[3][5][i]= 0;
	 F_primeprime[3][6][i] = 0;
	 F_primeprime[3][7][i] = -2*(xi[i]-energy0); 
		 F_primeprime[4][0][i] = 0;
	 F_primeprime[4][1][i] = 0;
	 F_primeprime[4][2][i] = 0;
	 F_primeprime[4][3][i] = 0;
	 F_primeprime[4][4][i] = 0;
	 F_primeprime[4][5][i]= 0;
	 F_primeprime[4][6][i] = 0;
	 F_primeprime[4][7][i] = 0;
		  F_primeprime[5][0][i] = 0;
	 F_primeprime[5][1][i] = 0;
	 F_primeprime[5][2][i] = 0;
	 F_primeprime[5][3][i] = 0;
	 F_primeprime[5][4][i] = 0;
	 F_primeprime[5][5][i]= 0;
	 F_primeprime[5][6][i] = 0;
	 F_primeprime[5][7][i] = 0;
		 F_primeprime[7][0][i] = 0;
	 F_primeprime[6][1][i] = 0;
	 F_primeprime[6][2][i] = 0;
	 F_primeprime[6][3][i] = 0;
	 F_primeprime[6][4][i] = 0;
	 F_primeprime[6][5][i]= 0;
	 F_primeprime[6][6][i] = 0;
	 F_primeprime[6][7][i] = 0; 
		 F_primeprime[7][0][i] = 0;
	 F_primeprime[7][1][i] = 0;
	 F_primeprime[7][2][i] = 2*(xi[i]-energy0);
	 F_primeprime[7][3][i] = -2*(xi[i]-energy0);
	 F_primeprime[7][4][i] = 0;
	 F_primeprime[7][5][i]= 0;
	 F_primeprime[7][6][i] = 0;
	 F_primeprime[7][7][i] = -2*cbest+2*gamma;
	 
       }
   
       if( xi[i] > energy0 ){
	 
	 for(a=0; a<dim; a++)
	   for(b=0;b<dim; b++)
	     F_primeprimeprimeprime[a][b][i] = 0;
       }
     }

       jasper=count;



   ReadEu152();
   for(i=0; i<count;i++)
     {
       if( xi[i] <= energy0 )
	 {
                   F_primeprime[0][0][jasper+i] = 0;
	 F_primeprime[0][1][jasper+i] = 0;
	 F_primeprime[0][2][jasper+i] = 0;
	 F_primeprime[0][3][jasper+i] = 0;
	 F_primeprime[0][4][jasper+i] = 0;
	 F_primeprime[0][5][jasper+i]= 0;
	 F_primeprime[0][6][jasper+i] = 0;
	 F_primeprime[0][7][jasper+i] = 0;
	           F_primeprime[1][0][jasper+i] = 0;
	 F_primeprime[1][1][jasper+i] = 0;
	 F_primeprime[1][2][jasper+i] = 0;
	 F_primeprime[1][3][jasper+i] = 0;
	 F_primeprime[1][4][jasper+i] = 0;
	 F_primeprime[1][5][jasper+i]= 0;
	 F_primeprime[1][6][jasper+i] = 0;
	 F_primeprime[1][7][jasper+i] = 0;
	           F_primeprime[2][0][jasper+i] = 0;
	 F_primeprime[2][1][jasper+i] = 0;
	 F_primeprime[2][2][jasper+i] = 0;
	 F_prim[2][3][jasper+i] = 0;
	 F_primeprime[2][4][jasper+i] = 0;
	 F_primeprime[2][5][jasper+i]= 0;
	 F_primeprime[2][6][jasper+i] = 0;
	 F_primeprime[2][7][jasper+i] = 2*(xi[i]-energy0);
	           F_primeprime[3][0][i] = 0;
	 F_primeprime[3][1][jasper+i] = 0;
	 F_primeprime[3][2][jasper+i] = 0;
	 F_primeprime[3][3][jasper+i] = 0;
	 F_primeprime[3][4][jasper+i] = 0;
	 F_primeprime[3][5][jasper+i]= 0;
	 F_primeprime[3][6][jasper+i] = 0;
	 F_primeprime[3][7][jasper+i] = -2*(xi[i]-energy0); 
		 F_primeprime[4][0][jasper+i] = 0;
	 F_primeprime[4][1][jasper+i] = 0;
	 F_primeprime[4][2][jasper+i] = 0;
	 F_primeprime[4][3][jasper+i] = 0;
	 F_primeprime[4][4][jasper+i] = 0;
	 F_primeprime[4][5][jasper+i]= 0;
	 F_primeprime[4][6][jasper+i] = 0;
	 F_primeprime[4][7][jasper+i] = 0;
		  F_primeprime[5][0][jasper+i] = 0;
	 F_primeprime[5][1][jasper+i] = 0;
	 F_primeprime[5][2][jasper+i] = 0;
	 F_primeprime[5][3][jasper+i] = 0;
	 F_primeprime[5][4][jasper+i] = 0;
	 F_primeprime[5][5][jasper+i]= 0;
	 F_primeprime[5][6][jasper+i] = 0;
	 F_primeprime[5][7][jasper+i] = 0;
		 F_primeprime[6][0][jasper+i] = 0;
	 F_primeprime[6][1][jasper+i] = 0;
	 F_primeprime[6][2][jasper+i] = 0;
	 F_primeprime[6][3][jasper+i] = 0;
	 F_primeprime[6][4][jasper+i] = 0;
	 F_primeprime[6][5][jasper+i]= 0;
	 F_primeprime[6][6][jasper+i] = 0;
	 F_primeprime[6][7][jasper+i] = 0; 
		 F_primeprime[7][0][jasper+i] = 0;
	 F_primeprime[7][1][jasper+i] = 0;
	 F_primeprime[7][2][jasper+i] = 2*(xi[i]-energy0);
	 F_primeprime[7][3][jasper+i] = -2*(xi[i]-energy0);
	 F_primeprime[7][4][jasper+i] = 0;
	 F_primeprime[7][5][jasper+i]= 0;
	 F_primeprime[7][6][jasper+i] = 0;
	 F_primeprime[7][7][jasper+i] = -2*cbest+2*gamma;
	 }
   
       if( xi[i] > energy0 ){
	 for(a=0; a<dim; a++)
	   for(b=0;b<dim; b++)
	     F_primeprimeprimeprime[a][b][jasper+i] = 0;
       }
     }

   jasper=jasper+count;

   ReadBa133();
   for(i=0; i<count;i++)
     {
       if( xi[i] <= energy0 )
	 {
                   F_primeprime[0][0][jasper+i] = 0;
	 F_primeprime[0][1][jasper+i] = 0;
	 F_primeprime[0][2][jasper+i] = 0;
	 F_primeprime[0][3][jasper+i] = 0;
	 F_primeprime[0][4][jasper+i] = 0;
	 F_primeprime[0][5][jasper+i]= 0;
	 F_primeprime[0][6][jasper+i] = 0;
	 F_primeprime[0][7][jasper+i] = 0;
	           F_primeprime[1][0][jasper+i] = 0;
	 F_primeprime[1][1][jasper+i] = 0;
	 F_primeprime[1][2][jasper+i] = 0;
	 F_primeprime[1][3][jasper+i] = 0;
	 F_primeprime[1][4][jasper+i] = 0;
	 F_primeprime[1][5][jasper+i]= 0;
	 F_primeprime[1][6][jasper+i] = 0;
	 F_primeprime[1][7][jasper+i] = 0;
	           F_primeprime[2][0][jasper+i] = 0;
	 F_primeprime[2][1][jasper+i] = 0;
	 F_primeprime[2][2][jasper+i] = 0;
	 F_primeprime[2][3][jasper+i] = 0;
	 F_primeprime[2][4][jasper+i] = 0;
	 F_primeprime[2][5][jasper+i]= 0;
	 F_primeprime[2][6][jasper+i] = 0;
	 F_primeprime[2][7][jasper+i] = 2*(xi[i]-energy0);
	           F_primeprime[3][0][i] = 0;
	 F_primeprime[3][1][jasper+i] = 0;
	 F_primeprime[3][2][jasper+i] = 0;
	 F_primeprime[3][3][jasper+i] = 0;
	 F_primeprime[3][4][jasper+i] = 0;
	 F_primeprime[3][5][jasper+i]= 0;
	 F_primeprime[3][6][jasper+i] = 0;
	 F_primeprime[3][7][jasper+i] = -2*(xi[i]-energy0); 
		 F_primeprime[4][0][jasper+i] = 0;
	 F_primeprime[4][1][jasper+i] = 0;
	 F_primeprime[4][2][jasper+i] = 0;
	 F_primeprime[4][3][jasper+i] = 0;
	 F_primeprime[4][4][jasper+i] = 0;
	 F_primeprime[4][5][jasper+i]= 0;
	 F_primeprime[4][6][jasper+i] = 0;
	 F_primeprime[4][7][jasper+i] = 0;
		  F_primeprime[5][0][jasper+i] = 0;
	 F_primeprime[5][1][jasper+i] = 0;
	 F_primeprime[5][2][jasper+i] = 0;
	 F_primeprime[5][3][jasper+i] = 0;
	 F_primeprime[5][4][jasper+i] = 0;
	 F_primeprime[5][5][jasper+i]= 0;
	 F_primeprime[5][6][jasper+i] = 0;
	 F_primeprime[5][7][jasper+i] = 0;
		 F_primeprime[6][0][jasper+i] = 0;
	 F_primeprime[6][1][jasper+i] = 0;
	 F_primeprime[6][2][jasper+i] = 0;
	 F_primeprime[6][3][jasper+i] = 0;
	 F_primeprime[6][4][jasper+i] = 0;
	 F_primeprime[6][5][jasper+i]= 0;
	 F_primeprime[6][6][jasper+i] = 0;
	 F_primeprime[6][7][jasper+i] = 0; 
		 F_primeprime[7][0][jasper+i] = 0;
	 F_primeprime[7][1][jasper+i] = 0;
	 F_primeprime[7][2][jasper+i] = 2*(xi[i]-energy0);
	 F_primeprime[7][3][jasper+i] = -2*(xi[i]-energy0);
	 F_primeprime[7][4][jasper+i] = 0;
	 F_primeprime[7][5][jasper+i]= 0;
	 F_primeprime[7][6][jasper+i] = 0;
	 F_primeprime[7][7][jasper+i] = -2*cbest+2*gamma;
	 }
   
       if( xi[i] > energy0 ){
	 for(a=0; a<dim; a++)
	   for(b=0;b<dim; b++)
	     F_primeprimeprimeprime[a][b][jasper+i] = 0;
       }
     }

   jasper=jasper+count;  

      ReadCo56();
   for(i=0; i<count;i++)
     {
      if( xi[i] <= energy0 )
	 {
                   F_primeprime[0][0][jasper+i] = 0;
	 F_primeprime[0][1][jasper+i] = 0;
	 F_primeprime[0][2][jasper+i] = 0;
	 F_primeprime[0][3][jasper+i] = 0;
	 F_primeprime[0][4][jasper+i] = 0;
	 F_primeprime[0][5][jasper+i]= 0;
	 F_primeprime[0][6][jasper+i] = 0;
	 F_primeprime[0][7][jasper+i] = 0;
	           F_primeprime[1][0][jasper+i] = 0;
	 F_primeprime[1][1][jasper+i] = 0;
	 F_primeprime[1][2][jasper+i] = 0;
	 F_primeprime[1][3][jasper+i] = 0;
	 F_primeprime[1][4][jasper+i] = 0;
	 F_primeprime[1][5][jasper+i]= 0;
	 F_primeprime[1][6][jasper+i] = 0;
	 F_primeprime[1][7][jasper+i] = 0;
	           F_primeprime[2][0][jasper+i] = 0;
	 F_primeprime[2][1][jasper+i] = 0;
	 F_primeprime[2][2][jasper+i] = 0;
	 F_primeprime[2][3][jasper+i] = 0;
	 F_primeprime[2][4][jasper+i] = 0;
	 F_primeprime[2][5][jasper+i]= 0;
	 F_primeprime[2][6][jasper+i] = 0;
	 F_primeprime[2][7][jasper+i] = 2*(xi[i]-energy0);
	           F_primeprime[3][0][i] = 0;
	 F_primeprime[3][1][jasper+i] = 0;
	 F_primeprime[3][2][jasper+i] = 0;
	 F_primeprime[3][3][jasper+i] = 0;
	 F_primeprime[3][4][jasper+i] = 0;
	 F_primeprime[3][5][jasper+i]= 0;
	 F_primeprime[3][6][jasper+i] = 0;
	 F_primeprime[3][7][jasper+i] = -2*(xi[i]-energy0); 
		 F_primeprime[4][0][jasper+i] = 0;
	 F_primeprime[4][1][jasper+i] = 0;
	 F_primeprime[4][2][jasper+i] = 0;
	 F_primeprime[4][3][jasper+i] = 0;
	 F_primeprime[4][4][jasper+i] = 0;
	 F_primeprime[4][5][jasper+i]= 0;
	 F_primeprime[4][6][jasper+i] = 0;
	 F_primeprime[4][7][jasper+i] = 0;
		  F_primeprime[5][0][jasper+i] = 0;
	 F_primeprime[5][1][jasper+i] = 0;
	 F_primeprime[5][2][jasper+i] = 0;
	 F_primeprime[5][3][jasper+i] = 0;
	 F_primeprime[5][4][jasper+i] = 0;
	 F_primeprime[5][5][jasper+i]= 0;
	 F_primeprime[5][6][jasper+i] = 0;
	 F_primeprime[5][7][jasper+i] = 0;
		 F_primeprime[6][0][jasper+i] = 0;
	 F_primeprime[6][1][jasper+i] = 0;
	 F_primeprime[6][2][jasper+i] = 0;
	 F_primeprime[6][3][jasper+i] = 0;
	 F_primeprime[6][4][jasper+i] = 0;
	 F_primeprime[6][5][jasper+i]= 0;
	 F_primeprime[6][6][jasper+i] = 0;
	 F_primeprime[6][7][jasper+i] = 0; 
		 F_primeprime[7][0][jasper+i] = 0;
	 F_primeprime[7][1][jasper+i] = 0;
	 F_primeprime[7][2][jasper+i] = 2*(xi[i]-energy0);
	 F_primeprime[7][3][jasper+i] = -2*(xi[i]-energy0);
	 F_primeprime[7][4][jasper+i] = 0;
	 F_primeprime[7][5][jasper+i]= 0;
	 F_primeprime[7][6][jasper+i] = 0;
	 F_primeprime[7][7][jasper+i] = -2*cbest+2*gamma;
	 }
   
       if( xi[i] > energy0 ){
	 for(a=0; a<dim; a++)
	   for(b=0;b<dim; b++)
	     F_primeprimeprimeprime[a][b][jasper+i] = 0;
       }
     }

   jasper=jasper+count;
    
}
