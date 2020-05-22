#ifndef FPO_KERNELS_H
#define FPO_KERNELS_H

#include <../../Lib/gkyl_ipow.h>

extern "C" {

  void fpoMomsKernelP1(const double *dv, const double *vc, const double *f, double *out);

  void fpoMomsKernelP2(const double *dv, const double *vc, const double *f, double *out);
  
  void fpoDiagKernelP1(const double *dv, const double *vc, const double *f, const double *h, double *out);
  
  void fpoDiagKernelP2(const double *dv, const double *vc, const double *f, const double *h, double *out);

  double fpoDragKernelP1(const double dt, const double *dv, const double *fC, const double *fL, const double *fR, const double *fT, const double *fB, const double *hC, const double *hL, const double *hR, const double *hT, const double *hB, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, double *fOut);

  double fpoDragKernelP2(const double dt, const double *dv, const double *fC, const double *fL, const double *fR, const double *fT, const double *fB, const double *hC, const double *hL, const double *hR, const double *hT, const double *hB, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, double *fOut);

  double fpoDiffKernelP1(const double dt, const double *dv, const double *fTL, const double *fT, const double *fTR, const double *fL, const double *fC, const double *fR, const double *fBL, const double *fB, const double *fBR, const double *gTL, const double *gT, const double *gTR, const double *gL, const double *gC, const double *gR, const double *gBL, const double *gB, const double *gBR, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, double *fOut);

  double fpoDragKernel3xP1(const double dt, const double *dv,
                           const double *f,
                           const double *f1L, const double *f1U,
                           const double *f2L, const double *f2U,
                           const double *f3L, const double *f3U,
                           const double *h,
                           const double *h1L, const double *h1U,
                           const double *h2L, const double *h2U,
                           const double *h3L, const double *h3U,
                           const int isTopEdge, const int isBotEdge,
                           const int isLeftEdge, const int isRightEdge,
                           double *fOut);
}
#endif
