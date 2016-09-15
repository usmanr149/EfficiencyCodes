#include "foo.h"

double Vary_chi(long double abest, long double bbest, long double cbest, long double gbest, long double dbest, long double ebest, , long double fbest, long double min){

  long double chisquaredl = 0, chisquaredr = 0, chisquared = 0;
  long double Chi2;
  long double delta_a=0, delta_b=0, delta_c=0, delta_d=0, delta_g=0, delta_e=0, delta_f=0;
  
    long double E0, zig, change;

   Chi2 = min;
  
   do{
  zig =  powl(2.71828182845,energy0) + (long double)deltaE0;
  printf("zig : %Lf\n", zig );
  E0 = log( zig );

  printf("count : %d\n", count);
  printf("energy0 : %Lf\n", powl(2.71828182845,energy0) );
  printf("E0 : %Lf\n", E0);

  alphabest = (abest+delta_a) - (cbest+delta_c)*E0*E0 + (gbest+delta_g)*E0*E0;
  betabest = (bbest+delta_b) + 2*(cbest+delta_c)*E0 - 2*(gbest+delta_g)*E0;

  ReadCo60();
  {
  for(i=0;i<count;i++)
      if(xi[i]<=energy0)
    {
      chisquaredl+=((yi[i]-(alphabest + betabest*xi[i] + (gbest+delta_g)*xi[i]*xi[i]))*(yi[i]-(alphabest + betabest*xi[i] + (gbest+delta_g)*xi[i]*xi[i]))/(ei[i]*ei[i]));
    }

  for(i=0;i<count;i++)
      if(xi[i]>=energy0)
    {
      chisquaredr+=(yi[i]-((abest+delta_a) + (bbest+delta_b)*xi[i] + (cbest+delta_c)*xi[i]*xi[i]))*(yi[i]-((abest+delta_a) + (bbest+delta_b)*xi[i] + (cbest+delta_c)*xi[i]*xi[i]))/(ei[i]*ei[i]);
    }
}

   ReadEu152();
   {
     for(i=0;i<count;i++)
    if(xi[i]<=energy0)
    {
      chisquaredl+=((yi[i]-(alphabest + betabest*xi[i] + (gbest+delta_g)*xi[i]*xi[i]+(dbest+delta_d)))*(yi[i]-(alphabest + betabest*xi[i] + (gbest+delta_g)*xi[i]*xi[i]+(dbest+delta_d)))/(ei[i]*ei[i]));
    }

     for(i=0;i<count;i++)
      if(xi[i]>=energy0)
    {
      chisquaredr+=((yi[i]-((abest+delta_a) + (bbest+delta_b)*xi[i] + (cbest+delta_c)*xi[i]*xi[i]+(dbest+delta_d)))*(yi[i]-((abest+delta_a) + (bbest+delta_b)*xi[i] + (cbest+delta_c)*xi[i]*xi[i]+(dbest+delta_d)))/(ei[i]*ei[i]));
    }
  
   }

  ReadCo56();
  {
    for(i=0;i<count;i++)
      if(xi[i]<=energy0)
    {
      chisquaredl+=((yi[i]-(alphabest + betabest*xi[i] + (gbest+delta_g)*xi[i]*xi[i]+f))*(yi[i]-(alphabest + betabest*xi[i] + (gbest+delta_g)*xi[i]*xi[i]+(fbest+delta_f)))/(ei[i]*ei[i]));
    }

    for(i=0;i<count;i++)
      if(xi[i]>=energy0)
	{
      chisquaredr+=((yi[i]-((abest+delta_a) + (bbest+delta_b)*xi[i] + (cbest+delta_c)*xi[i]*xi[i]+(fbest+delta_f)))*(yi[i]-((abest+delta_a) + (bbest+delta_b)*xi[i] + (cbest+delta_c)*xi[i]*xi[i]+(fbest+delta_f))))/(ei[i]*ei[i]);
    }
  
  }

  ReadBa133();
  {
  
    for(i=0;i<count;i++)
      if(xi[i]<=energy0)
    {
      chisquaredl+=((yi[i]-(alphabest + betabest*xi[i] + (gbest+delta_g)*xi[i]*xi[i]+(ebest+delta_e)))*(yi[i]-(alphabest + betabest*xi[i] + (gbest+delta_g)*xi[i]*xi[i]+(ebest+delta_e)))/(ei[i]*ei[i]));      
    }

    for(i=0;i<count;i++)
      if(xi[i]>=energy0)
    {
      chisquaredr+=((yi[i]-((abest+delta_a) + (bbest+delta_b)*xi[i] + (cbest+delta_c)*xi[i]*xi[i]+(ebest+delta_e)))*(yi[i]-((abest+delta_a) + (bbest+delta_b)*xi[i] + (cbest+delta_c)*xi[i]*xi[i]+(ebest+delta_e)))/(ei[i]*ei[i]));
    }
  }

  chisquared=chisquaredl+chisquaredr;
  change = min-chisquared;
  if(change < 0)
    {
      change = change * -1;
    }
  delta_a += 0.01
   }(while change < 1)
  
      printf("%10.6Lf \n", delta_a);
      
  //printf("%10.6Lf \t %10.6Lf\n", energy0, chisquared);
  //fprintf(fp2,"%10.6Lf \t %10.6Lf\n", powl(2.71828182845,E0), chisquared);


}
