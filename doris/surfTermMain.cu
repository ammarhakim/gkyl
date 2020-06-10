#include <stdio.h>
#include <iostream>
#include <algorithm>

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

// forward declaration in case of separate compilation
__host__ __device__ double VlasovSurfElcMag2x3vSer_VX_P1(const double *wl, const double *wr,
                                                         const double *dxvl, const double *dxvr,
                                                         const double amax, const double *EM,
                                                         const double *fl, const double *fr,
                                                         double *outl, double *outr);

__host__ __device__ double VlasovSurfElcMag2x3vSer_VX_P1_v1(const int nCell, const double *wl, const double *wr,
                                                            const double *dxvl, const double *dxvr,
                                                            const double amax, const double *EM,
                                                            const double *fl, const double *fr,
                                                            double *outl, double *outr);
    
__host__ __device__ double VlasovSurfElcMag2x3vSer_VX_P1_v2(const int nCell, const int nSM,  const double *wl,
                                                            const double *wr,
                                                            const double *dxvl, const double *dxvr,
                                                            const double amax, const double *EM,
                                                            const double *fl, const double *fr,
                                                            double *outl, double *outr);

__global__ void k_surfTerm_org(const int nCells, const int nfPerCell, const int nv, const REAL* xcC, const REAL* dx,
                              const REAL* EM, const REAL* fin, REAL* fout)
{
  const int tid = threadIdx.x + blockIdx.x*blockDim.x;
  if (tid >= nCells) return;

  // finL, finC, finR and pointers to the first distribution on the left, center and right cells 
  const REAL *finC = &fin[tid*nfPerCell];

  // perodic BC
  const REAL* finL = (tid > 0) ? &fin[(tid-1)*nfPerCell] : &fin[(nCells-1)*nfPerCell];
  const REAL* finR = (tid < nCells-1) ? &fin[(tid+1)*nfPerCell] : &fin[0];
  
  REAL* foutC = &fout[tid*nfPerCell];

  // need to map the pointer of d_EM[x,y] for each thread
  const REAL* myEM = &EM[(tid / nv) * 24];

  // should be numComponents, but want to avoid dynamic memory alloc
  double dummy[200];
  /* left (of C) surface update. use dummy in place of fRhsOutL so that only current cell updated.
     eq->surfTerm(dir, dummy, dummy, xcL, xcC, dx, dx, 0., idxL, idxC, fInL, fInC, dummy, fRhsOutC);
     Eq/GkyVlasov.cpp Vlasov::surfTerm(int dir, double *cflL, double *cflR,
                   double *xcL, double *xcR, double *dxL, double *dxR,
                   double maxsOld, int *idxL, int *idxR,
                   double *qInL, double *qInR, double *qRhsOutL, double *qRhsOutR)
  */
  VlasovSurfElcMag2x3vSer_VX_P1(xcC, xcC, dx, dx, 0, myEM, finL, finC, dummy, foutC);

  // right (of C) surface update. use dummy in place of fRhsOutR so that only current cell (C) is updated.
  // eq->surfTerm(dir, dummy, dummy, xcC, xcR, dx, dx, 0., idxC, idxR, fInC, fInR, fRhsOutC, dummy);

  VlasovSurfElcMag2x3vSer_VX_P1(xcC, xcC, dx, dx, 0, myEM, finC, finR, foutC, dummy);
}

__global__ void k_surfTerm_v1(const int nCells, const int nfPerCell, const int nv, const REAL* xcC, const REAL* dx,
                          const REAL* EM, const REAL* fin, REAL* fout)
{
  const unsigned int tid = threadIdx.x + blockIdx.x*blockDim.x;
  if (tid >= nCells) return;
  
  // finL, finC, finR and pointers to the first distribution on the left, center and right cells
  // note nfPerCell is passed as 1 in this case
  const REAL *finC = &fin[tid];

  // perodic BC
  const REAL* finL = (tid > 0) ? &fin[tid-1] : &fin[nCells-1];
  const REAL* finR = (tid < nCells-1) ? &fin[tid+1] : &fin[0];

  REAL* foutC = &fout[tid];

  // need to map the pointer of d_EM[x,y] for each thread
  const REAL* myEM = &EM[(tid / nv) * 24];

  // left (of C) surface update. use dummy in place of fRhsOutL so that only current cell updated.
  //   eq->surfTerm(dir, dummy, dummy, xcL, xcC, dx, dx, 0., idxL, idxC, fInL, fInC, dummy, fRhsOutC);
  VlasovSurfElcMag2x3vSer_VX_P1_v1(nCells, xcC, xcC, dx, dx, 0, myEM, finL, finC, NULL, foutC);

  // right (of C) surface update. use dummy in place of fRhsOutR so that only current cell (C) is updated.
  // eq->surfTerm(dir, dummy, dummy, xcC, xcR, dx, dx, 0., idxC, idxR, fInC, fInR, fRhsOutC, dummy);
  VlasovSurfElcMag2x3vSer_VX_P1_v1(nCells, xcC, xcC, dx, dx, 0, myEM, finC, finR, foutC, NULL);
}

__global__ void
__launch_bounds__(256, 1)
k_surfTerm_v2(const int nCells, const int nfPerCell, const int nv, const REAL* xcC, const REAL* dx,
                              const REAL* EM, const REAL* fin, REAL* fout)
{
  const unsigned int tid = threadIdx.x + blockIdx.x*blockDim.x;
  if (tid >= nCells) return;

  // index of the share memory fast dimension is [0:blockDim+2]
  const unsigned int si = threadIdx.x + 1;
  const unsigned int sstride = blockDim.x + 2;
  
  extern __shared__ REAL s_data[];
  
  // each thread loads 32 f's to shared memory
  for (auto i=0; i<nfPerCell; i++) s_data[i*sstride + si] = fin[i*nCells + tid];

  // fill in left halo in shared memory
  if (threadIdx.x == 0 )
  {
    if (tid > 0) {
      for (auto i=0; i<nfPerCell; i++) s_data[i*sstride + si - 1]  = fin[i*nCells + tid - 1];
    }
    else {  // period BC
      for (auto i=0; i<nfPerCell; i++) s_data[i*sstride + si - 1]  = fin[i*nCells + nCells - 1];
    }
  }
  // fill in right halo in shared memory
  // need to take care of cases when nCells % blockDim.x !=0
  // normally "if (tid >= nCells) return" would be ok unless shared memory is used 
  if (nCells % blockDim.x > 0) {
    if (blockIdx.x == gridDim.x - 1 && threadIdx.x == (nCells % blockDim.x) - 1) {
      for (auto i=0; i<nfPerCell; i++) s_data[i*sstride + si + 1]  = fin[i*nCells];
    }
    if (blockIdx.x < gridDim.x - 1 && threadIdx.x == blockDim.x - 1) {
      for (auto i=0; i<nfPerCell; i++) s_data[i*sstride + si + 1]  = fin[i*nCells + tid + 1];
    }
  }
  if (threadIdx.x == blockDim.x - 1)
  {
    if (tid < nCells-1) {
      for (auto i=0; i<nfPerCell; i++) s_data[i*sstride + si + 1]  = fin[i*nCells + tid + 1];
    }
    else { // period BC
        for (auto i=0; i<nfPerCell; i++) s_data[i*sstride + si + 1]  = fin[i*nCells];
    }
  }
  __syncthreads();
  
  // finL, finC, finR and pointers to the first distribution on the left, center and right cells
  const REAL* finC = &s_data[si];
  const REAL* finL = &s_data[si-1];
  const REAL* finR = &s_data[si+1];

  REAL* foutC = &fout[tid];

  // need to map the pointer of d_EM[x,y] for each thread
  const REAL* myEM = &EM[(tid / nv) * 24];

  // should be numComponents, but want to avoid dynamic memory alloc
  // left (of C) surface update. use dummy in place of fRhsOutL so that only current cell updated.
  //   eq->surfTerm(dir, dummy, dummy, xcL, xcC, dx, dx, 0., idxL, idxC, fInL, fInC, dummy, fRhsOutC);
  VlasovSurfElcMag2x3vSer_VX_P1_v2(nCells, sstride, xcC, xcC, dx, dx, 0, myEM, finL, finC, NULL, foutC);
  
  // right (of C) surface update. use dummy in place of fRhsOutR so that only current cell (C) is updated.
  // eq->surfTerm(dir, dummy, dummy, xcC, xcR, dx, dx, 0., idxC, idxR, fInC, fInR, fRhsOutC, dummy);
  VlasovSurfElcMag2x3vSer_VX_P1_v2(nCells, sstride, xcC, xcC, dx, dx, 0, myEM, finC, finR, foutC, NULL);
}



/* a standalone wrapper for VlasovSurfElcMag2x3vSer_VX_P1
 w[NDIM]: Cell-center coordinates.
 dxv[NDIM]: Cell spacing.
 EM/f: Input EM-field/distribution function.
 out: Incremented output
 the grid for f is (x, y, vx, vy, vz), vz the fastest
 the grid for EM is (x, y), y the fastest

*/
int main(int argc, char **argv)
{
  const int nx = 31;
  const int ny = 31;
  
  const int nDim = 5;
  const int nvx = 6, nvy = 6, nvz = 6;
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


  int NUM_REPS=10;
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
  const REAL *dx = &d_const[nDim];
  
  dim3 const bdim(256);
  dim3 const gdim(nCells/bdim.x + 1);

  // warm up the kernel for performance testing
  k_surfTerm_org<<< gdim, bdim >>>(nCells, nfPerCell, nv, xcC, dx, d_EM, d_fin, d_fout);
  cudaMemset(d_fout, 0, nf *sizeof(REAL));
  
  //****************************
  // original data structure for distribution is (nCells, nfPerCell) with nfPerCell is the faster dim
  checkCudaErrors(cudaEventRecord(start, 0) );
  for (auto i=0; i<NUM_REPS; i++)
  {
    k_surfTerm_org<<< gdim, bdim >>>(nCells, nfPerCell, nv, xcC, dx, d_EM, d_fin, d_fout);
  }

  checkCudaErrors( cudaEventRecord(stop, 0) );
  checkCudaErrors( cudaEventSynchronize(stop) );
  checkCudaErrors( cudaEventElapsedTime(&ms, start, stop) );
  
  auto bandwidth = nf * 2 * sizeof(REAL) * 1e-6 * NUM_REPS / ms;
  printf("Original surfTerm kernel %.2f GB/s\n", bandwidth);

  checkCudaErrors(cudaMemcpy(h_fout, d_fout, nf*sizeof(REAL), cudaMemcpyDeviceToHost)); 

  //****************************
  // V1 new data structure for distribution is (nfPerCell, nCells) with nfPerCell is the slower dim  
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
#if 0
  for (auto i=0; i<nx*ny; i++) {
    for (auto j=0; j<24; j++) {
      h_EM1[j*nx*ny + i] = h_EM[i*24 + j];
    }
  }
#endif
  
  // swap the pointers so i can reuse h_fin1 buffer
  std::swap(h_fin1, h_fin);
  dim3 const bdim_v1(256);
  dim3 const gdim_v1(nCells/bdim_v1.x + 1);
  
  checkCudaErrors(cudaMemcpy(d_fin, h_fin, nf*sizeof(REAL), cudaMemcpyHostToDevice));
  //  checkCudaErrors(cudaMemcpy(d_EM, h_EM1, nEM*sizeof(REAL), cudaMemcpyHostToDevice));

  checkCudaErrors(cudaEventRecord(start, 0) );
  for (auto i=0; i<NUM_REPS; i++)
  {
    k_surfTerm_v1<<< gdim_v1, bdim_v1 >>>(nCells, nfPerCell, nv, xcC, dx, d_EM, d_fin, d_fout);
  }

  checkCudaErrors( cudaEventRecord(stop, 0) );
  checkCudaErrors( cudaEventSynchronize(stop) );
  checkCudaErrors( cudaEventElapsedTime(&ms, start, stop) );

  bandwidth = nf * 2 * sizeof(REAL) * 1e-6 * NUM_REPS / ms;
  printf("New surfTerm kernel v1 %.2f GB/s\n", bandwidth);

  checkCudaErrors(cudaMemcpy(h_fin1, d_fout, nf*sizeof(REAL), cudaMemcpyDeviceToHost));

  //****************************
  // Compare results using two data structures/two device functions
  // will first need to transpose h_fin1 to fout1 with the original (nCells, nfPerCell) structure
  for (auto j=0; j<nfPerCell; j++) {
    for (auto i=0; i<nCells; i++)  h_fout1[i*nfPerCell + j] = h_fin1[j*nCells + i];
  }
  
  for (auto i = 0; i < nf; i++)
  {
    if (h_fout1[i] != h_fout[i])
    { std::cout << "wrong results V1: h_fout1[" << i<< "] = " << h_fout1[i] << " "
                << " correct h_fout[" << i << "] = " << h_fout[i] << std::endl;
      exit(-1);
    }
  }

  //checkCudaErrors(cudaDeviceSetCacheConfig(cudaFuncCachePreferEqual));
  //checkCudaErrors(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
  
  //****************************
  // V2 
  cudaMemset(d_fout, 0, nf *sizeof(REAL));

  dim3 const bdim_v2(128);
  dim3 const gdim_v2(nCells/bdim_v2.x + 1);

  // shared memry for save left and right 
  const int SMsize = (bdim_v2.x + 2) * nfPerCell * sizeof(REAL);
  
  checkCudaErrors(cudaEventRecord(start, 0) );
  for (auto i=0; i<NUM_REPS; i++)
  {
    k_surfTerm_v2<<< gdim_v2, bdim_v2, SMsize, 0 >>>(nCells, nfPerCell, nv, xcC, dx, d_EM, d_fin, d_fout);
  }
  checkCudaErrors(cudaGetLastError());
  

  checkCudaErrors( cudaEventRecord(stop, 0) );
  checkCudaErrors( cudaEventSynchronize(stop) );
  checkCudaErrors( cudaEventElapsedTime(&ms, start, stop) );

  bandwidth = nf * 2 * sizeof(REAL) * 1e-6 * NUM_REPS / ms;
  printf("New surfTerm kernel v2 %.2f GB/s\n", bandwidth);

  checkCudaErrors(cudaMemcpy(h_fin1, d_fout, nf*sizeof(REAL), cudaMemcpyDeviceToHost));

  //****************************
  // Compare results using two data structures/two device functions
  for (auto j=0; j<nfPerCell; j++) {
    for (auto i=0; i<nCells; i++)  h_fout1[i*nfPerCell + j] = h_fin1[j*nCells + i];
  }

  for (auto i = 0; i < nf; i++)
  {
    if (h_fout1[i] != h_fout[i])
    { std::cout << "wrong results V2: h_fout1[" << i<< "] = " << h_fout1[i] << " "
                << " correct h_fout[" << i << "] = " << h_fout[i] << std::endl;
      int ix,iy,ivx,ivy,ivz;
      ix = i/(ny*nvx*nvy*nvz*nfPerCell);
      iy = (i - ix*(ny*nvx*nvy*nvz*nfPerCell)) / (nvx*nvy*nvz*nfPerCell);
      ivx = (i - ix*(ny*nvx*nvy*nvz*nfPerCell) - iy * (nvx*nvy*nvz*nfPerCell)) / (nvy*nvz*nfPerCell);
      ivy = (i - ix*(ny*nvx*nvy*nvz*nfPerCell) - iy * (nvx*nvy*nvz*nfPerCell) - ivx * nvy*nvz*nfPerCell) / (nvz*nfPerCell);
      ivz = (i - ix*(ny*nvx*nvy*nvz*nfPerCell) - iy * (nvx*nvy*nvz*nfPerCell) - ivx * nvy*nvz*nfPerCell - ivy * nvz*nfPerCell) / nfPerCell;
      
      std::cout << ix << " " << iy << " " << ivx << " " << ivy << " " << ivz
                << " " << i%nfPerCell << std::endl;
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

