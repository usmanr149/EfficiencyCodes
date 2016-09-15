#define Extended_ROWS 40
#define Extended_COLS 40

#include "foo.h"

/*First-order derivative matrix*/

long double First_order_der(long double Matrix_1[int dimR_1][int dimC_1], long double Matrix_2[dimR_2][dimC_2]){

   ReadCo60();
   for(i=0; i<count;i++)
     {
       if( xi[i] <= energy0 ){
	 F_prime[i][0] += 1;
	 F_prime[i][1] += xi[i];
	 F_prime[i][2] += (2*xi[i]*energy0-energy0*energy0);
	 F_primt[i][3] += (xi[i]-energy0)*(xi[i]-energy0);
	 F_prime[i][4] += 0;
	 F_prime[i][5] += 0;
	 F_prime[i][6] += 0;
	 F_prime[i][7] += cbest*(2*xi[i]-2*energy0)+gammabest*2*(xi[i] - energy0);
       }
   
       if( xi[i] <= energy0 ){
	 F_prime[i][0] = 1;
	 F_prime[i][1] = xi[i];
	 F_prime[i][2] = xi[i]*xi[i];
	 F_primt[i][3] = 0;
	 F_prime[i][4] = 0;
	 F_prime[i][5] = 0;
	 F_prime[i][6] = 0;
	 F_prime[i][7] = 0;
       }
       jasper=count+1;
     }
   ReadEu152();
   for(i=0; i<count;i++)
     {
       if( xi[i] <= energy0 ){
	 F_prime[jasper+i][0] += 1;
	 F_prime[jasper+i][1] += xi[i];
	 F_prime[jasper+i][2] += (2*xi[i]*energy0-energy0*energy0);
	 F_primt[jasper+i][3] += (xi[i]-energy0)*(xi[i]-energy0);
	 F_prime[jasper+i][4] += 1;
	 F_prime[jasper+i][5] += 0;
	 F_prime[jasper+i][6] += 0;
	 F_prime[jasper+i][7] += cbest*(2*xi[i]-2*energy0)+gammabest*2*(xi[i] - energy0);
       }
   
       if( xi[i] <= energy0 ){
	 F_prime[jasper+i][0] = 1;
	 F_prime[jasper+i][1] = xi[i];
	 F_prime[jasper+i][2] = xi[i]*xi[i];
	 F_primt[jasper+i][3] = 0;
	 F_prime[jasper+i][4] = 1;
	 F_prime[jasper+i][5] = 0;
	 F_prime[jasper+i][6] = 0;
	 F_prime[jasper+i][7] = 0;
       }
       jasper=count+1;
     }

   ReadBa133();
   for(i=0; i<count;i++)
     {
       if( xi[i] <= energy0 ){
	 F_prime[jasper+i][0] += 1;
	 F_prime[jasper+i][1] += xi[i];
	 F_prime[jasper+i][2] += (2*xi[i]*energy0-energy0*energy0);
	 F_primt[jasper+i][3] += (xi[i]-energy0)*(xi[i]-energy0);
	 F_prime[jasper+i][4] += 0;
	 F_prime[jasper+i][5] += 1;
	 F_prime[jasper+i][6] += 0;
	 F_prime[jasper+i][7] += cbest*(2*xi[i]-2*energy0)+gammabest*2*(xi[i] - energy0);
       }
   
       if( xi[i] <= energy0 ){
	 F_prime[jasper+i][0] = 1;
	 F_prime[jasper+i][1] = xi[i];
	 F_prime[jasper+i][2] = xi[i]*xi[i];
	 F_primt[jasper+i][3] = 0;
	 F_prime[jasper+i][4] = 0;
	 F_prime[jasper+i][5] = 1;
	 F_prime[jasper+i][6] = 0;
	 F_prime[jasper+i][7] = 0;
       }
       jasper=count+1;
     }
   
      ReadCo56();
   for(i=0; i<count;i++)
     {
       if( xi[i] <= energy0 ){
	 F_prime[jasper+i][0] += 1;
	 F_prime[jasper+i][1] += xi[i];
	 F_prime[jasper+i][2] += (2*xi[i]*energy0-energy0*energy0);
	 F_primt[jasper+i][3] += (xi[i]-energy0)*(xi[i]-energy0);
	 F_prime[jasper+i][4] += 0;
	 F_prime[jasper+i][5] += 0;
	 F_prime[jasper+i][6] += 1;
	 F_prime[jasper+i][7] += cbest*(2*xi[i]-2*energy0)+gammabest*2*(xi[i] - energy0);
       }
   
       if( xi[i] <= energy0 ){
	 F_prime[jasper+i][0] = 1;
	 F_prime[jasper+i][1] = xi[i];
	 F_prime[jasper+i][2] = xi[i]*xi[i];
	 F_primt[jasper+i][3] = 0;
	 F_prime[jasper+i][4] = 0;
	 F_prime[jasper+i][5] = 0;
	 F_prime[jasper+i][6] = 1;
	 F_prime[jasper+i][7] = 0;
       }
       jasper=count+1;
     }
}
