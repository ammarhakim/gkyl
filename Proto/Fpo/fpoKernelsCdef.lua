local ffi  = require "ffi"
ffi.cdef [[

typedef struct {
  double* C;
  double* xL;
  double* xU;
  double* yL;
  double* yU;
  double* zL;
  double* zU;
} stencil7;

typedef struct {
  double* xLyLzC;
  double* xLyCzL;
  double* xLyCzC;
  double* xLyCzU;
  double* xLyUzC;
  double* xCyLzL;
  double* xCyLzC;
  double* xCyLzU;
  double* xCyCzL;
  double* xCyCzC;
  double* xCyCzU;
  double* xCyUzL;
  double* xCyUzC;
  double* xCyUzU;
  double* xUyLzC;
  double* xUyCzL;
  double* xUyCzC;
  double* xUyCzU;
  double* xUyUzC;
} stencil19;

void fpoMomsKernelP1(const double* dv, const double* vc, const double* f, double* out);

void fpoMomsKernelP2(const double* dv, const double* vc, const double* f, double* out);


void fpoDiagKernelP1(const double* dv, const double* vc, const double* f, const double* h, double* out);

void fpoDiagKernelP2(const double* dv, const double* vc, const double* f, const double* h, double* out);


double fpoDragKernel3xP1(const double dt, const double* dv,
                         const stencil7* fStencil, const stencil7* hStencil,
                         double* fOut);

void fpoDiffSurfXLSer3xP1(const double dt, const double* dv,
                          const double* fLCC,
                          const double* fLLC,
                          const double* fLUC,
                          const double* fLCL,
                          const double* fLCU,
                          const double* fCCC,
                          const double* fCLC,
                          const double* fCUC,
                          const double* fCCL,
                          const double* fCCU,
                          const double* gLCC,
                          const double* gLLC,
                          const double* gLUC,
                          const double* gLCL,
                          const double* gLCU,
                          const double* gCCC,
                          const double* gCLC,
                          const double* gCUC,
                          const double* gCCL,
                          const double* gCCU,
                          double *fOut);

void fpoDiffSurfXUSer3xP1(const double dt, const double* dv,
                          const double* fCCC,
                          const double* fCLC,
                          const double* fCUC,
                          const double* fCCL,
                          const double* fCCU,
                          const double* fUCC,
                          const double* fULC,
                          const double* fUUC,
                          const double* fUCL,
                          const double* fUCU,
                          const double* gCCC,
                          const double* gCLC,
                          const double* gCUC,
                          const double* gCCL,
                          const double* gCCU,
                          const double* gUCC,
                          const double* gULC,
                          const double* gUUC,
                          const double* gUCL,
                          const double* gUCU,
                          double *fOut);

void fpoDiffSurfYLSer3xP1(const double dt, const double* dv,
                          const double* fCLC,
                          const double* fLLC,
                          const double* fULC,
                          const double* fCLL,
                          const double* fCLU,
                          const double* fCCC,
                          const double* fLCC,
                          const double* fUCC,
                          const double* fCCL,
                          const double* fCCU,
                          const double* gCLC,
                          const double* gLLC,
                          const double* gULC,
                          const double* gCLL,
                          const double* gCLU,
                          const double* gCCC,
                          const double* gLCC,
                          const double* gUCC,
                          const double* gCCL,
                          const double* gCCU,
                          double* fOut);

void fpoDiffSurfYUSer3xP1(const double dt, const double* dv,
                          const double* fCCC,
                          const double* fLCC,
                          const double* fUCC,
                          const double* fCCL,
                          const double* fCCU,
                          const double* fCUC,
                          const double* fLUC,
                          const double* fUUC,
                          const double* fCUL,
                          const double* fCUU,
                          const double* gCCC,
                          const double* gLCC,
                          const double* gUCC,
                          const double* gCCL,
                          const double* gCCU,
                          const double* gCUC,
                          const double* gLUC,
                          const double* gUUC,
                          const double* gCUL,
                          const double* gCUU,
                          double* fOut);

void fpoDiffSurfZLSer3xP1(const double dt, const double* dv,
                          const double* fCCL,
                          const double* fLCL,
                          const double* fUCL,
                          const double* fCLL,
                          const double* fCUL,
                          const double* fCCC,
                          const double* fLCC,
                          const double* fUCC,
                          const double* fCLC,
                          const double* fCUC,
                          const double* gCCL,
                          const double* gLCL,
                          const double* gUCL,
                          const double* gCLL,
                          const double* gCUL,
                          const double* gCCC,
                          const double* gLCC,
                          const double* gUCC,
                          const double* gCLC,
                          const double* gCUC,
                          double* fOut);

void fpoDiffSurfZUSer3xP1(const double dt, const double* dv,
                          const double* fCCC,
                          const double* fLCC,
                          const double* fUCC,
                          const double* fCLC,
                          const double* fCUC,
                          const double* fCCU,
                          const double* fLCU,
                          const double* fUCU,
                          const double* fCLU,
                          const double* fCUU,
                          const double* gCCC,
                          const double* gLCC,
                          const double* gUCC,
                          const double* gCLC,
                          const double* gCUC,
                          const double* gCCU,
                          const double* gLCU,
                          const double* gUCU,
                          const double* gCLU,
                          const double* gCUU,
                          double* fOut);

double fpoDiffVolSer3xP1(const double dt, const double* dv,
                         const double* fCCC,
                         const double* gCCC,
                         double* fOut);

]]
