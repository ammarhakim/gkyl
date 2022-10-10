/* -*- c++ -*- */

// MF 2022/10/09: A temporary .cu file to make waf work (a hack).

#include <cstdio>
#include <GkylDummyHeader.h>

__global__ void cuda_GkylDummyKer(const GkylDummyStruct_t* __restrict__ dummy, int a, int b) {

  // get setup data from GkylDummyStruct_t structure
  const bool zeroF = dummy->zeroF;
  if (zeroF > a)
    printf("condition is true\n");
}

void advanceOnDevice(const int numBlocks, const int numThreads, const GkylDummyStruct_t *dummy) {
  cuda_GkylDummyKer<<<numBlocks, numThreads>>>(dummy, 0, 1);
}

__global__ void dummySetDtKer(GkylDummyStruct_t *dummy, double dt) {
  dummy->dt = dt;
}

void setDt(GkylDummyStruct_t *dummy) {
  dummySetDtKer<<<1, 1>>>(dummy, 1.0);
}
