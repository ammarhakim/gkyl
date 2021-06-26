#ifndef FPO_KERNELS_H
#define FPO_KERNELS_H

//#include <../../Lib/gkyl_ipow.h>

extern "C" {

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

  void vlasov_fpo_moms_3x_ser_p1(const double* dv, const double* vc, const double* f, double* out);
  double vlasov_fpo_momM0mod_3x_ser_p1(const double* dv, const double* vc, const double* f,
                                       int isXloEdge, int isXupEdge,
                                       int isYloEdge, int isYupEdge,
                                       int isZloEdge, int isZupEdge);
  double vlasov_fpo_momM2_3x_ser_p1(const double* dv, const double* vc, const double* f);
  
  void vlasov_fpo_diag_3x_ser_p1(const double* dv, const double* vc, const double* f, const double* h, double* out);

  double vlasov_fpo_drag_cell_3x_ser_p1(const double dt, const double *dv,
                                        const stencil7* fStencil, const stencil7* hStencil,
                                        const int isXloEdge, const int isXupEdge,
                                        const int isYloEdge, const int isYupEdge,
                                        const int isZloEdge, const int isZupEdge,
                                        double* fOut);

  void vlasov_fpo_diffI_surfxL_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  void vlasov_fpo_diffI_surfxU_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  void vlasov_fpo_diffI_surfyL_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  void vlasov_fpo_diffI_surfyU_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  void vlasov_fpo_diffI_surfzL_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  void vlasov_fpo_diffI_surfzU_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  double vlasov_fpo_diffI_vol_3x_ser_p1(const double dt, const double* dv,
                                       const double* fIn,
                                       const double* gIn,
                                       double *fOut);
  

  void vlasov_fpo_diffII_surfxL_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  void vlasov_fpo_diffII_surfxU_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  void vlasov_fpo_diffII_surfyL_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  void vlasov_fpo_diffII_surfyU_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  void vlasov_fpo_diffII_surfzL_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  void vlasov_fpo_diffII_surfzU_3x_ser_p1(const double dt, const double* dv,
                                        const double* fLC,
                                        const double* fL1L, const double* fL1U,
                                        const double* fL2L, const double* fL2U,
                                        const double* fUC,
                                        const double* fU1L, const double* fU1U,
                                        const double* fU2L, const double* fU2U,
                                        const double* gLC,
                                        const double* gL1L, const double* gL1U,
                                        const double* gL2L, const double* gL2U,
                                        const double* gUC,
                                        const double* gU1L, const double* gU1U,
                                        const double* gU2L, const double* gU2U,
                                        const int isBoundary,
                                        double* fOut);

  double vlasov_fpo_diffII_vol_3x_ser_p1(const double dt, const double* dv,
                                       const double* fIn,
                                       const double* gIn,
                                       double *fOut);
  
}
#endif
