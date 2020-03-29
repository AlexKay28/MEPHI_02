#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>

double f(double x);
double f(double x)
{
  return (4.0 / (1.0 + x*x));
}

double pi_calculate(const int n);
double pi_calculate(const int n)
{
  int i;
  double x,h,sum;
  h = 1.0 / (double) n;
  sum = 0.0;
  for (i=1; i<=n; i++){
    x = h * ((double)i - 0.5);
    sum += f(x);
  }
  return h*sum;
}

int main(int argc, char *argv[])
{	
  double t, mypi;
  long long int t1, t2;

  QueryPerformanceCounter(&t1);
  mypi = pi_calculate(100000000);
  QueryPerformanceCounter(&t2);
  t = 1e-6 * (t2 - t1);

  printf("Time: %lf sec mypi = %14.12lf\n",t,mypi);

  return 0;
}
