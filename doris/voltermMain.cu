#include <stdio.h>
#include <iostream>

#define REAL double

#define checkCudaErrors(stat)                                                     \
  { cudaErrCheck_((stat), __FILE__, __LINE__); }
void cudaErrCheck_(cudaError_t stat, const char *file, int line) {
  if (stat != cudaSuccess) {
    fprintf(stderr, "CUDA Error: %s %s %d\n", cudaGetErrorString(stat), file,
            line);
    exit(-1);
  }
}

#ifdef INLINE
#include <VlasovSer2x3vP1.cu>
#include <VlasovSer2x3vP1_v1.cu>
#else
// forward declaration in case of separate compilation
extern __host__ __device__ double VlasovVol2x3vSerP1(const double *w, const double *dxv, const double *EM, const double *f,
                                                     double *out);
extern __host__ __device__ double VlasovVol2x3vSerP1_v1(const int nC, const int nC_EM, const double *w,
                                                        const double *dxv, const double *EM,
                                                        const double *f, double *out);
#endif

__global__ void k_volTerm_org(const int nCells, const int nfPerCell, const int nv, const REAL* xcC, const REAL* dx,
                              const REAL* EM, const REAL* fin, REAL* fout)
{
  const int tid = threadIdx.x + blockIdx.x*blockDim.x;
  if (tid >= nCells) return;

  // myfin and myfout is the pointer to the first distribution f[0] of the cell 
  const REAL* myfin = &fin[tid*nfPerCell];
  REAL* myfout = &fout[tid*nfPerCell];

  // need to map the pointer of d_EM[x,y] for each thread
  const REAL* myEM = &EM[(tid / nv) * 24];

  double cflRate = VlasovVol2x3vSerP1(xcC, dx, myEM, myfin, myfout);
  
}
__global__ void k_volTerm(const int nCells, const int nfPerCell, const int nv, const REAL* xcC, const REAL* dx,
                          const REAL* EM, const REAL* fin, REAL* fout)
{
  const unsigned int tid = threadIdx.x + blockIdx.x*blockDim.x;
  if (tid >= nCells) return;
  
  // myfin and myfout is the pointer to the first distribution f[0] of the cell
  const REAL* myfin = &fin[tid];
  REAL* myfout = &fout[tid];
  
  // need to map the pointer of d_EM[x,y] for each thread
  const REAL* myEM = &EM[(tid / nv)];

  double cflRate = VlasovVol2x3vSerP1_v1(nCells, nCells/nv, xcC, dx, myEM, myfin, myfout);
  
}


/* a standalone wrapper for VlasovVol1x2vSerP2(const double *w, const double *dxv, const double *EM, const double *f, double *out)
 w[NDIM]: Cell-center coordinates.
 dxv[NDIM]: Cell spacing.
 EM/f: Input EM-field/distribution function.
 out: Incremented output
 the grid for f is (x, y, vx, vy, vz), vz the fastest
 the grid for EM is (x, y), y the fastest

*/
int main(int argc, char **argv)
{
  const int nx = 8;
  const int ny = 32;
  
  const int nDim = 5;
  const int nvx = 16, nvy = 8, nvz = 32;
  const int nv = nvx*nvy*nvz;
  
  const int nfPerCell = 32;
  const size_t nEM = (size_t)nx*ny*24;
  const size_t nCells = (size_t)nx*ny*nvx*nvy*nvz;
  const size_t nf = nCells * nfPerCell;
  
  REAL *h_EM, *h_fin, *h_fout, *h_fout1;
  REAL *d_EM, *d_fin, *d_fout;

  printf("nx: %d, ny: %d, nDim: %d, nVx: %d, nVy: %d, nVz: %d\n", nx, ny, nDim, nvx, nvy, nvz);
  printf("Total number of cells is %d\n", nx*ny*nvx*nvy*nvz);

  cudaMallocHost((void **) &h_EM, nEM*sizeof(REAL));
  cudaMallocHost((void **) &h_fin, nf *sizeof(REAL));
  cudaMalloc((void **) &d_EM, nEM*sizeof(REAL));
  cudaMalloc((void **) &d_fin, nf*sizeof(REAL));
  cudaMalloc((void **) &d_fout, nf*sizeof(REAL));
  
  // To compare results using two different approach
  cudaMallocHost((void **) &h_fout, nf *sizeof(REAL));
  cudaMallocHost((void **) &h_fout1, nf *sizeof(REAL));
  
  // initialize EM and fin
  //srand(time(NULL));
#pragma omp parallel for 
  for (auto i=0; i<nEM; i++)
  {
    h_EM[i] = (float)(rand() % nEM) / nEM * ((i%2 == 0) ? 1 : -1) ; 
  }
  checkCudaErrors(cudaMemcpy(d_EM, h_EM, nEM*sizeof(REAL), cudaMemcpyHostToDevice));

#pragma omp parallel for 
  for (auto i=0; i<nf; i++)
  {
    h_fin[i] =(float)(rand() % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  checkCudaErrors(cudaMemcpy(d_fin, h_fin, nf*sizeof(REAL), cudaMemcpyHostToDevice));


  int NUM_REPS=2;
  cudaEvent_t start, stop;
  checkCudaErrors( cudaEventCreate(&start) );
  checkCudaErrors( cudaEventCreate(&stop) );
  float ms;

  // assuming xcC and dx are constant for each cell
  REAL h_const[2*nDim];
  for (auto i=0; i<2*nDim; i++) h_const[i] = 0.1 * i;

  REAL *d_const;
  cudaMalloc((void **) &d_const, 2*nDim*sizeof(REAL));
  cudaMemcpy(d_const, h_const, 2*nDim*sizeof(REAL), cudaMemcpyHostToDevice);

  const REAL *xcC = d_const;
  const REAL *dx = &d_const[3];
  
  dim3 const bdim(128);
  dim3 const gdim(nCells/bdim.x + 1);

  //****************************
  // original data structure for distribution is (nCells, nfPerCell) with nfPerCell is the faster dim
  checkCudaErrors(cudaEventRecord(start, 0) );
  for (auto i=0; i<NUM_REPS; i++)
  {
    k_volTerm_org<<< gdim, bdim >>>(nCells, nfPerCell, nv, xcC, dx, d_EM, d_fin, d_fout);
  }

  checkCudaErrors( cudaEventRecord(stop, 0) );
  checkCudaErrors( cudaEventSynchronize(stop) );
  checkCudaErrors( cudaEventElapsedTime(&ms, start, stop) );
  
  auto bandwidth = nf * 2 * sizeof(REAL) * 1e-6 * NUM_REPS / ms;
  printf("Original volTerm kernel %.2f GB/s\n", bandwidth);

  checkCudaErrors(cudaMemcpy(h_fout, d_fout, nf*sizeof(REAL), cudaMemcpyDeviceToHost)); 

  //****************************
  // new data structure for distribution is (nfPerCell, nCells) with nfPerCell is the slower dim  
  // reset d_fout to be 0
  cudaMemset(d_fout, 0, nf *sizeof(REAL));

  REAL *h_fin1;
  REAL *h_EM1;
  cudaMallocHost((void **) &h_fin1, nf *sizeof(REAL));
  cudaMallocHost((void **) &h_EM1, nEM *sizeof(REAL));

  // quick and primitive way of transposing host array
  // could have transposed the device array using cuTensor 
  for (auto i=0; i<nCells; i++) {
    for (auto j=0; j<nfPerCell; j++) { 
      h_fin1[j*nCells + i] = h_fin[i*nfPerCell + j];
    }
  }
  for (auto i=0; i<nx*ny; i++) {
    for (auto j=0; j<24; j++) { 
      h_EM1[j*nx*ny + i] = h_EM[i*24 + j];
    }
  }
           
  checkCudaErrors(cudaMemcpy(d_fin, h_fin1, nf*sizeof(REAL), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_EM, h_EM1, nEM*sizeof(REAL), cudaMemcpyHostToDevice));

  checkCudaErrors(cudaEventRecord(start, 0) );
  for (auto i=0; i<NUM_REPS; i++)
  {
    k_volTerm<<< gdim, bdim >>>(nCells, nfPerCell, nv, xcC, dx, d_EM, d_fin, d_fout);
  }

  checkCudaErrors( cudaEventRecord(stop, 0) );
  checkCudaErrors( cudaEventSynchronize(stop) );
  checkCudaErrors( cudaEventElapsedTime(&ms, start, stop) );

  bandwidth = nf * 2 * sizeof(REAL) * 1e-6 * NUM_REPS / ms;
  printf("New volTerm kernel %.2f GB/s\n", bandwidth);

  checkCudaErrors(cudaMemcpy(h_fin1, d_fout, nf*sizeof(REAL), cudaMemcpyDeviceToHost));
  //checkCudaErrors(cudaMemcpy(h_fout1, d_fout, nf*sizeof(REAL), cudaMemcpyDeviceToHost));

  //****************************
  // Compare results using two data structures/two device functions
  // will first need to transpose h_fin1 to fout1 with the original (nCells, nfPerCell) structure
  for (auto j=0; j<nfPerCell; j++) {
    for (auto i=0; i<nCells; i++)  h_fout1[i*nfPerCell + j] = h_fin1[j*nCells + i];
  }

  
  for (auto i = 0; i < nf; i++)
  {
    if (h_fout1[i] != h_fout[i])
    { std::cout << "wrong results: new h_fout1[" << i<< "] = " << h_fout1[i] << " "
                << " correct h_fout[" << i << "] = " << h_fout[i] << std::endl;
      exit(-1);
    }
  }

  std::cout << "TEST PASSED\n";
  cudaFreeHost(h_fout); cudaFreeHost(h_fout1);
  cudaFreeHost(h_fin); cudaFreeHost(h_fin1);
  cudaFree(d_EM); cudaFree(d_fin); cudaFree(d_fout);
  cudaFree(d_const);
  
  

  return 0;
 
}

