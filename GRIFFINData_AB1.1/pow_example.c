/* pow example */


#include <stdio.h>
#include <math.h>
 
int main()
{
  double c, d, result;
 
  printf("Enter c and d to calculate c^d\n");
  scanf("%lf%lf", &c, &d);
 
  result = pow(c, d);
 
  printf("%.2lf raised to %.2lf = %.2lf\n", c, d, result);
 
  return 0;
}
