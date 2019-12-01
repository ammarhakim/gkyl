__global__
void sumArray(int n, double a, double *x, double *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i<n) y[i] = a*x[i]+y[i];
}
