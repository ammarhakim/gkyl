local ffi  = require "ffi"
ffi.cdef [[

void fpoMomsKernelP1(const double *dv, const double *vc, const double *f, double *out);

void fpoMomsKernelP2(const double *dv, const double *vc, const double *f, double *out);
  
void fpoDiagKernelP1(const double *dv, const double *vc, const double *f, const double *h, double *out);
  
void fpoDiagKernelP2(const double *dv, const double *vc, const double *f, const double *h, double *out);

void fpoDragKernelP1(const double dt, const double *dv, const double *f, const double *fL, const double *fR, const double *fT, const double *fB, const double *h, const double *hL, const double *hR, const double *hT, const double *hB, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, double *fOut);

void fpoDragKernelP2(const double dt, const double *dv, const double *f, const double *fL, const double *fR, const double *fT, const double *fB, const double *h, const double *hL, const double *hR, const double *hT, const double *hB, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, double *fOut);

void fpoDiffKernelP1(const double dt, const double *dv, const double *f, const double *fL, const double *fR, const double *fT, const double *fB, const double *g, const double *gL, const double *gR, const double *gT, const double *gB, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, const int cross, double *fOut);

void fpoDiffKernelP2(const double dt, const double *dv, const double *f, const double *fL, const double *fR, const double *fT, const double *fB, const double *g, const double *gL, const double *gR, const double *gT, const double *gB, const int isTopEdge, const int isBotEdge, const int isLeftEdge, const int isRightEdge, const int cross, double *fOut);

]]
