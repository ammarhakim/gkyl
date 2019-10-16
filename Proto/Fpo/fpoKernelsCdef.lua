local ffi  = require "ffi"
ffi.cdef [[

void fpoMomsKernelP1(const double *dv, const double *vc, const double *f, double *out);

  
void fpoDiagKernelP1(const double *dv, const double *vc, const double *f, const double *h, double *out);
  

void fpoDragKernelP1(const double dt, const double *dv, const double *fC, const double *fL, const double *fR, const double *fT, const double *fB, const double *hC, const double *hL, const double *hR, const double *hT, const double *hB, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, double *fOut);


void fpoDiffKernelP1(const double dt, const double *dv, const double *fTL, const double *fT, const double *fTR, const double *fL, const double *fC, const double *fR, const double *fBL, const double *fB, const double *fBR, const double *gTL, const double *gT, const double *gTR, const double *gL, const double *gC, const double *gR, const double *gBL, const double *gB, const double *gBR, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, int cross, double *fOut);


]]
