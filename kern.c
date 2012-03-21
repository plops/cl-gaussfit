#include <math.h>
#include <stdio.h>
// gcc kern.c -std=c99 -lm -o kern

void makeGaussianKernel(double sigma, double accuracy, int maxRadius)
{
  double PI=atan(-1)*-4;
  int kRadius = (int)ceil(sigma*sqrt(-2*log(accuracy)))+1;
  if (maxRadius < 50) maxRadius = 50;
  // too small maxRadius would result in inaccurate sum.
  if (kRadius > maxRadius) kRadius = maxRadius;
  float kernel[2][kRadius];
  for (int i=0; i<kRadius; i++)               // Gaussian function
    kernel[0][i] = (float)(exp(-0.5*i*i/sigma/sigma));
  if (kRadius < maxRadius && kRadius > 3) {   // edge correction
    double sqrtSlope = 1e34;
    int r = kRadius;
    while (r > kRadius/2) {
      r--;
      double a = sqrt(kernel[0][r])/(kRadius-r);
      if (a < sqrtSlope)
	sqrtSlope = a;
      else
	break;
    }
    for (int r1 = r+2; r1 < kRadius; r1++)
      kernel[0][r1] = (float)((kRadius-r1)*(kRadius-r1)*sqrtSlope*sqrtSlope);
  }
  double sum;      // sum over all kernel elements for normalization
  if (kRadius < maxRadius) {
    sum = kernel[0][0];
    for (int i=1; i<kRadius; i++)
      sum += 2*kernel[0][i];
  } else
    sum = sigma * sqrt(2*PI);
   
  double rsum = 0.5 + 0.5*kernel[0][0]/sum;
  for (int i=0; i<kRadius; i++) {
    double v = (kernel[0][i]/sum);
    kernel[0][i] = (float)v;
    rsum -= v;
    kernel[1][i] = (float)rsum;
    printf("%d sum %f rsum %f\n",i,v,rsum);
  }
  return;
}

int main()
{
  makeGaussianKernel(1.2,1e-4,40);
  return 0;
}
