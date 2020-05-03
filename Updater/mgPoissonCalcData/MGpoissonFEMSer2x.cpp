#include <MGpoissonModDecl.h> 
 
void MGpoissonFEMProlong2xSer_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF5[0] += 0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[1]; 
  fldF5[0] += 0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[1]; 
  fldF5[0] += 0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[1]; 
  fldF5[0] += 0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxDirichletLyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF3[0] += fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxDirichletLyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF3[0] += fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxDirichletLyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF3[0] += fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxNeumannLyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF3[0] += fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxNeumannLyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF3[0] += fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxNeumannLyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF3[0] += fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxRobinLyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF3[0] += fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxRobinLyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF3[0] += fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxRobinLyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldF0[0] += 0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF3[0] += fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxDirichletUyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxDirichletUyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxDirichletUyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxNeumannUyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxNeumannUyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxNeumannUyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxRobinUyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxRobinUyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_LxRobinUyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF2[0] += 0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF1[1] += fldC[1]; 
  fldF3[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxDirichletLyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxDirichletLyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxDirichletLyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxNeumannLyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxNeumannLyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxNeumannLyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxRobinLyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxRobinLyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxRobinLyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldF0[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF5[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF2[0] += 0.25*fldC[0]; 
  fldF6[0] += 0.5*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxDirichletUyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[3]+0.25*fldC[2]+0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[2]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[3]+0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF1[1] += fldC[2]; 
  fldF0[2] += 0.5*fldC[3]+0.5*fldC[2]; 
  fldF0[3] += fldC[3]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[2]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[2]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxDirichletUyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[3]+0.25*fldC[2]+0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[2]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[3]+0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF1[1] += fldC[2]; 
  fldF0[2] += 0.5*fldC[3]+0.5*fldC[2]; 
  fldF0[3] += fldC[3]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[2]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[2]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxDirichletUyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[3]+0.25*fldC[2]+0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[2]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[3]+0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF1[1] += fldC[2]; 
  fldF0[2] += 0.5*fldC[3]+0.5*fldC[2]; 
  fldF0[3] += fldC[3]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[2]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[2]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxNeumannUyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[3]+0.25*fldC[2]+0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[2]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[3]+0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF1[1] += fldC[2]; 
  fldF0[2] += 0.5*fldC[3]+0.5*fldC[2]; 
  fldF0[3] += fldC[3]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[2]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[2]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxNeumannUyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[3]+0.25*fldC[2]+0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[2]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[3]+0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF1[1] += fldC[2]; 
  fldF0[2] += 0.5*fldC[3]+0.5*fldC[2]; 
  fldF0[3] += fldC[3]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[2]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[2]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxNeumannUyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[3]+0.25*fldC[2]+0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[2]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[3]+0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF1[1] += fldC[2]; 
  fldF0[2] += 0.5*fldC[3]+0.5*fldC[2]; 
  fldF0[3] += fldC[3]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[2]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[2]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxRobinUyDirichlet_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[3]+0.25*fldC[2]+0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[2]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[3]+0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF1[1] += fldC[2]; 
  fldF0[2] += 0.5*fldC[3]+0.5*fldC[2]; 
  fldF0[3] += fldC[3]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[2]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[2]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxRobinUyNeumann_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[3]+0.25*fldC[2]+0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[2]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[3]+0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF1[1] += fldC[2]; 
  fldF0[2] += 0.5*fldC[3]+0.5*fldC[2]; 
  fldF0[3] += fldC[3]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[2]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[2]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMProlong2xSer_UxRobinUyRobin_P1(const double *fldC, double **fldF) 
{ 
  // fldC: coarse-grid field.
  // fldF: fine-grid field in cells pointed to by the stencil.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldF0[0] += 0.25*fldC[3]+0.25*fldC[2]+0.25*fldC[1]+0.25*fldC[0]; 
  fldF1[0] += 0.5*fldC[2]+0.5*fldC[0]; 
  fldF4[0] += 0.5*fldC[1]+0.5*fldC[0]; 
  fldF7[0] += fldC[0]; 
  fldF0[1] += 0.5*fldC[3]+0.5*fldC[1]; 
  fldF4[1] += fldC[1]; 
  fldF1[1] += fldC[2]; 
  fldF0[2] += 0.5*fldC[3]+0.5*fldC[2]; 
  fldF0[3] += fldC[3]; 
  fldF10[0] += 0.5*fldC[0]; 
  fldF2[0] += 0.25*fldC[2]+0.25*fldC[0]; 
  fldF2[1] += 0.5*fldC[2]; 
  fldF5[0] += 0.25*fldC[1]+0.25*fldC[0]; 
  fldF8[0] += 0.5*fldC[0]; 
  fldF5[1] += 0.5*fldC[1]; 
  fldF11[0] += 0.25*fldC[0]; 
}

void MGpoissonFEMRestrict2xSer_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF2[1]+fldF1[1]+0.5*fldF0[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF2[1]+fldF1[1]+0.5*fldF0[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF2[1]+fldF1[1]+0.5*fldF0[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxDirichletLyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldC[0] += fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxDirichletLyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldC[0] += fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxDirichletLyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldC[0] += fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxNeumannLyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldC[0] += fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxNeumannLyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldC[0] += fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxNeumannLyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldC[0] += fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxRobinLyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldC[0] += fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxRobinLyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldC[0] += fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxRobinLyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];

  fldC[0] += fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxDirichletUyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF1[1]+0.5*fldF0[1]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxDirichletUyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF1[1]+0.5*fldF0[1]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxDirichletUyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF1[1]+0.5*fldF0[1]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxNeumannUyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF1[1]+0.5*fldF0[1]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxNeumannUyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF1[1]+0.5*fldF0[1]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxNeumannUyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF1[1]+0.5*fldF0[1]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxRobinUyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF1[1]+0.5*fldF0[1]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxRobinUyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF1[1]+0.5*fldF0[1]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_LxRobinUyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF3 = fldF[3];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.25*fldF3[0]+0.5*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF1[1]+0.5*fldF0[1]+0.5*fldF1[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxDirichletLyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF4[1]+0.5*fldF0[1]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxDirichletLyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF4[1]+0.5*fldF0[1]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxDirichletLyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF4[1]+0.5*fldF0[1]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxNeumannLyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF4[1]+0.5*fldF0[1]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxNeumannLyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF4[1]+0.5*fldF0[1]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxNeumannLyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF4[1]+0.5*fldF0[1]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxRobinLyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF4[1]+0.5*fldF0[1]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxRobinLyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF4[1]+0.5*fldF0[1]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxRobinLyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF6 = fldF[6];

  fldC[0] += 0.5*fldF6[0]+fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += fldF4[1]+0.5*fldF0[1]+0.5*fldF4[0]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxDirichletUyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]+0.5*fldF0[1]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxDirichletUyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]+0.5*fldF0[1]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxDirichletUyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]+0.5*fldF0[1]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxNeumannUyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]+0.5*fldF0[1]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxNeumannUyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]+0.5*fldF0[1]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxNeumannUyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]+0.5*fldF0[1]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxRobinUyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]+0.5*fldF0[1]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxRobinUyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]+0.5*fldF0[1]+0.25*fldF0[0]; 
}

void MGpoissonFEMRestrict2xSer_UxRobinUyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF0 = fldF[0];
  double *fldF1 = fldF[1];
  double *fldF2 = fldF[2];
  double *fldF4 = fldF[4];
  double *fldF5 = fldF[5];
  double *fldF7 = fldF[7];
  double *fldF8 = fldF[8];
  double *fldF10 = fldF[10];
  double *fldF11 = fldF[11];

  fldC[0] += 0.5*fldF8[0]+fldF7[0]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF2[0]+0.25*fldF11[0]+0.5*fldF10[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[1] += 0.5*fldF5[1]+fldF4[1]+0.5*fldF0[1]+0.25*fldF5[0]+0.5*fldF4[0]+0.25*fldF0[0]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]+0.5*fldF0[1]+0.25*fldF0[0]; 
}

void MGpoissonFEM_DGtoFEM_2xSer_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 
  double *femFldUy = femOut[2]; 
  double *femFldUxUy = femOut[3]; 

  femFld[0] += 0.375*dgFld[3]-0.2165063509461096*dgFld[2]-0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFldUx[0] += (-0.375*dgFld[3])-0.2165063509461096*dgFld[2]+0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFldUy[0] += (-0.375*dgFld[3])+0.2165063509461096*dgFld[2]-0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFldUxUy[0] += 0.375*dgFld[3]+0.2165063509461096*dgFld[2]+0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_2xSer_Lx_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 
  double *femFldUy = femOut[2]; 
  double *femFldUxUy = femOut[3]; 

  femFld[0] += 0.75*dgFld[3]-0.4330127018922193*dgFld[2]-0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFldUx[0] += (-0.375*dgFld[3])-0.2165063509461096*dgFld[2]+0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFldUy[0] += (-0.75*dgFld[3])+0.4330127018922193*dgFld[2]-0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFldUxUy[0] += 0.375*dgFld[3]+0.2165063509461096*dgFld[2]+0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_2xSer_Ux_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 
  double *femFldUy = femOut[2]; 
  double *femFldUxUy = femOut[3]; 

  femFld[0] += 0.375*dgFld[3]-0.2165063509461096*dgFld[2]-0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFld[1] += (-0.75*dgFld[3])-0.4330127018922193*dgFld[2]+0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFldUy[0] += (-0.375*dgFld[3])+0.2165063509461096*dgFld[2]-0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFldUy[1] += 0.75*dgFld[3]+0.4330127018922193*dgFld[2]+0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_2xSer_Ly_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 
  double *femFldUy = femOut[2]; 
  double *femFldUxUy = femOut[3]; 

  femFld[0] += 0.75*dgFld[3]-0.4330127018922193*dgFld[2]-0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFldUx[0] += (-0.75*dgFld[3])-0.4330127018922193*dgFld[2]+0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFldUy[0] += (-0.375*dgFld[3])+0.2165063509461096*dgFld[2]-0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFldUxUy[0] += 0.375*dgFld[3]+0.2165063509461096*dgFld[2]+0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_2xSer_Uy_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 
  double *femFldUy = femOut[2]; 
  double *femFldUxUy = femOut[3]; 

  femFld[0] += 0.375*dgFld[3]-0.2165063509461096*dgFld[2]-0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFldUx[0] += (-0.375*dgFld[3])-0.2165063509461096*dgFld[2]+0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFld[1] += (-0.75*dgFld[3])+0.4330127018922193*dgFld[2]-0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFldUx[1] += 0.75*dgFld[3]+0.4330127018922193*dgFld[2]+0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_2xSer_LxLy_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 
  double *femFldUy = femOut[2]; 
  double *femFldUxUy = femOut[3]; 

  femFld[0] += 1.5*dgFld[3]-0.8660254037844386*dgFld[2]-0.8660254037844386*dgFld[1]+0.5*dgFld[0]; 
  femFldUx[0] += (-0.75*dgFld[3])-0.4330127018922193*dgFld[2]+0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFldUy[0] += (-0.75*dgFld[3])+0.4330127018922193*dgFld[2]-0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFldUxUy[0] += 0.375*dgFld[3]+0.2165063509461096*dgFld[2]+0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_2xSer_LxUy_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 
  double *femFldUy = femOut[2]; 
  double *femFldUxUy = femOut[3]; 

  femFld[0] += 0.75*dgFld[3]-0.4330127018922193*dgFld[2]-0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFldUx[0] += (-0.375*dgFld[3])-0.2165063509461096*dgFld[2]+0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFld[1] += (-1.5*dgFld[3])+0.8660254037844386*dgFld[2]-0.8660254037844386*dgFld[1]+0.5*dgFld[0]; 
  femFldUx[1] += 0.75*dgFld[3]+0.4330127018922193*dgFld[2]+0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_2xSer_UxLy_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 
  double *femFldUy = femOut[2]; 
  double *femFldUxUy = femOut[3]; 

  femFld[0] += 0.75*dgFld[3]-0.4330127018922193*dgFld[2]-0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFld[1] += (-1.5*dgFld[3])-0.8660254037844386*dgFld[2]+0.8660254037844386*dgFld[1]+0.5*dgFld[0]; 
  femFldUy[0] += (-0.375*dgFld[3])+0.2165063509461096*dgFld[2]-0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFldUy[1] += 0.75*dgFld[3]+0.4330127018922193*dgFld[2]+0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 

}

void MGpoissonFEM_DGtoFEM_2xSer_UxUy_P1(const double *dgFld, double **femOut) 
{ 
  // dgFld:  DG (modal) field coefficients.
  // femOut: FEM (nodal) field coefficients.

  double *femFld = femOut[0]; 
  double *femFldUx = femOut[1]; 
  double *femFldUy = femOut[2]; 
  double *femFldUxUy = femOut[3]; 

  femFld[0] += 0.375*dgFld[3]-0.2165063509461096*dgFld[2]-0.2165063509461096*dgFld[1]+0.125*dgFld[0]; 
  femFld[1] += (-0.75*dgFld[3])-0.4330127018922193*dgFld[2]+0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFld[2] += (-0.75*dgFld[3])+0.4330127018922193*dgFld[2]-0.4330127018922193*dgFld[1]+0.25*dgFld[0]; 
  femFld[3] += 1.5*dgFld[3]+0.8660254037844386*dgFld[2]+0.8660254037844386*dgFld[1]+0.5*dgFld[0]; 

}

void MGpoissonFEMDampedJacobi2xSer_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((16.0*rhoUy[0]+4.0*rhoUxUy[0]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((12.0*phiUy[0]+3.0*phiUxUy[0]+3.0*phiUxLy[0]-6.0*phiUx[0]-24.0*phiPrevC[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+((-6.0*phiUy[0])+3.0*phiUxUy[0]+3.0*phiUxLy[0]+12.0*phiUx[0]-24.0*phiPrevC[0]-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+24.0*phiPrevC[0]*rdx2SqVol[1]+24.0*phiPrevC[0]*rdx2SqVol[0])/(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((-2.0*bcVals[2])+phiUx[0]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = -(1.0*((2.0*bcVals[2]+(phiPrevC[0]-1.0*phiUx[0])*bcVals[1]-2.0*bcVals[0]*phiPrevC[0])*omega-1.0*phiPrevC[0]*bcVals[1]+2.0*bcVals[0]*phiPrevC[0]))/(bcVals[1]-2.0*bcVals[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((15.0*rdx2SqVol[0]-3.0*rdx2SqVol[1])*bcVals[5]+(3.0*phiLy[1]+12.0*phiUy[0]-24.0*phiPrevC[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+((-6.0*phiUy[0])-24.0*phiPrevC[0]-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+24.0*phiPrevC[0]*rdx2SqVol[1]+24.0*phiPrevC[0]*rdx2SqVol[0])/(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[5]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((30.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*bcVals[5]+(3.0*phiLy[1]+15.0*phiUy[0]-30.0*phiPrevC[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+((-3.0*phiUy[0])-12.0*phiPrevC[0]-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+30.0*phiPrevC[0]*rdx2SqVol[1]+12.0*phiPrevC[0]*rdx2SqVol[0])/(30.0*rdx2SqVol[1]+12.0*rdx2SqVol[0]); 
  phiC[1] = ((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((54.0*rdx2SqVol[1]+54.0*rdx2SqVol[0])*bcVals[5]+((-30.0*phiPrevC[1])+3.0*phiLy[1]+15.0*phiUy[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*phiPrevC[1]+3.0*rdx2SqVol[0]*phiLy[1]+((-3.0*phiUy[0])-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+30.0*phiPrevC[1]*rdx2SqVol[1]+12.0*rdx2SqVol[0]*phiPrevC[1])/(30.0*rdx2SqVol[1]+12.0*rdx2SqVol[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);

  phiC[0] = (((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*bcVals[4]+(8.0*rhoUy[1]+8.0*rhoLy[1]+32.0*rhoC[1]+32.0*rhoUy[0]+32.0*rhoLy[0]+8.0*rhoLxUy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*bcVals[3])*omega*volFac+((30.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*bcVals[5]+((3.0*phiLy[1]+15.0*phiUy[0]-30.0*phiPrevC[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+((-3.0*phiUy[0])-12.0*phiPrevC[0]-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*bcVals[4]+((6.0*phiLy[1]+24.0*phiUy[0]-48.0*phiPrevC[0]+24.0*phiLy[0]+6.0*phiLxUy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiLy[1]+((-12.0*phiUy[0])-48.0*phiPrevC[0]-12.0*phiLy[0]+6.0*phiLxUy[0]+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0])*bcVals[3])*omega+(30.0*phiPrevC[0]*rdx2SqVol[1]+12.0*phiPrevC[0]*rdx2SqVol[0])*bcVals[4]+(48.0*phiPrevC[0]*rdx2SqVol[1]+48.0*phiPrevC[0]*rdx2SqVol[0])*bcVals[3])/((30.0*rdx2SqVol[1]+12.0*rdx2SqVol[0])*bcVals[4]+(48.0*rdx2SqVol[1]+48.0*rdx2SqVol[0])*bcVals[3]); 
  phiC[1] = (((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*bcVals4R2+(8.0*rhoUy[1]+8.0*rhoLy[1]+32.0*rhoC[1]+32.0*rhoUy[0]+32.0*rhoLy[0]+8.0*rhoLxUy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*bcVals[3]*bcVals[4])*omega*volFac+(((54.0*rdx2SqVol[1]+54.0*rdx2SqVol[0])*bcVals[4]+(96.0*rdx2SqVol[1]+96.0*rdx2SqVol[0])*bcVals[3])*bcVals[5]+(((-30.0*phiPrevC[1])+3.0*phiLy[1]+15.0*phiUy[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*phiPrevC[1]+3.0*rdx2SqVol[0]*phiLy[1]+((-3.0*phiUy[0])-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*bcVals4R2+(((-108.0*phiPrevC[1])+6.0*phiLy[1]+24.0*phiUy[0]+24.0*phiLy[0]+6.0*phiLxUy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*rdx2SqVol[1]-72.0*rdx2SqVol[0]*phiPrevC[1]+6.0*rdx2SqVol[0]*phiLy[1]+((-12.0*phiUy[0])-12.0*phiLy[0]+6.0*phiLxUy[0]+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+((-96.0*phiPrevC[1]*rdx2SqVol[1])-96.0*rdx2SqVol[0]*phiPrevC[1])*bcVals3R2)*omega+(30.0*phiPrevC[1]*rdx2SqVol[1]+12.0*rdx2SqVol[0]*phiPrevC[1])*bcVals4R2+(108.0*phiPrevC[1]*rdx2SqVol[1]+72.0*rdx2SqVol[0]*phiPrevC[1])*bcVals[3]*bcVals[4]+(96.0*phiPrevC[1]*rdx2SqVol[1]+96.0*rdx2SqVol[0]*phiPrevC[1])*bcVals3R2)/((30.0*rdx2SqVol[1]+12.0*rdx2SqVol[0])*bcVals4R2+(108.0*rdx2SqVol[1]+72.0*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(96.0*rdx2SqVol[1]+96.0*rdx2SqVol[0])*bcVals3R2); 

}

void MGpoissonFEMDampedJacobi2xSer_LyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 

}

void MGpoissonFEMDampedJacobi2xSer_LyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((-2.0*bcVals[8])+phiUy[0]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 

}

void MGpoissonFEMDampedJacobi2xSer_LyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = -(1.0*((2.0*bcVals[8]+(phiPrevC[0]-1.0*phiUy[0])*bcVals[7]-2.0*phiPrevC[0]*bcVals[6])*omega-1.0*phiPrevC[0]*bcVals[7]+2.0*phiPrevC[0]*bcVals[6]))/(bcVals[7]-2.0*bcVals[6]); 

}

void MGpoissonFEMDampedJacobi2xSer_UyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((15.0*rdx2SqVol[1]-3.0*rdx2SqVol[0])*bcVals[11]+(3.0*phiLx[1]+3.0*phiUxLy[0]-6.0*phiUx[0]-24.0*phiPrevC[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]+(3.0*phiUxLy[0]+12.0*phiUx[0]-24.0*phiPrevC[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+24.0*phiPrevC[0]*rdx2SqVol[1]+24.0*phiPrevC[0]*rdx2SqVol[0])/(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[11]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_UyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((30.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*bcVals[11]+(3.0*phiLx[1]+3.0*phiUxLy[0]-3.0*phiUx[0]-12.0*phiPrevC[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]+(3.0*phiUxLy[0]+15.0*phiUx[0]-30.0*phiPrevC[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+12.0*phiPrevC[0]*rdx2SqVol[1]+30.0*phiPrevC[0]*rdx2SqVol[0])/(12.0*rdx2SqVol[1]+30.0*rdx2SqVol[0]); 
  phiC[1] = ((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((54.0*rdx2SqVol[1]+54.0*rdx2SqVol[0])*bcVals[11]+((-12.0*phiPrevC[1])+3.0*phiLx[1]+3.0*phiUxLy[0]-3.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]-30.0*rdx2SqVol[0]*phiPrevC[1]+3.0*rdx2SqVol[0]*phiLx[1]+(3.0*phiUxLy[0]+15.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+12.0*phiPrevC[1]*rdx2SqVol[1]+30.0*rdx2SqVol[0]*phiPrevC[1])/(12.0*rdx2SqVol[1]+30.0*rdx2SqVol[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_UyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = (((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*bcVals[10]+(8.0*rhoUx[1]+8.0*rhoLx[1]+32.0*rhoC[1]+8.0*rhoUxLy[0]+32.0*rhoUx[0]+32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*bcVals[9])*omega*volFac+((30.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*bcVals[11]+((3.0*phiLx[1]+3.0*phiUxLy[0]-3.0*phiUx[0]-12.0*phiPrevC[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]+(3.0*phiUxLy[0]+15.0*phiUx[0]-30.0*phiPrevC[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*bcVals[10]+((6.0*phiLx[1]+6.0*phiUxLy[0]-12.0*phiUx[0]-48.0*phiPrevC[0]+24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiLx[1]+(6.0*phiUxLy[0]+24.0*phiUx[0]-48.0*phiPrevC[0]-12.0*phiLy[0]+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0])*bcVals[9])*omega+(12.0*phiPrevC[0]*rdx2SqVol[1]+30.0*phiPrevC[0]*rdx2SqVol[0])*bcVals[10]+(48.0*phiPrevC[0]*rdx2SqVol[1]+48.0*phiPrevC[0]*rdx2SqVol[0])*bcVals[9])/((12.0*rdx2SqVol[1]+30.0*rdx2SqVol[0])*bcVals[10]+(48.0*rdx2SqVol[1]+48.0*rdx2SqVol[0])*bcVals[9]); 
  phiC[1] = (((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*bcVals10R2+(8.0*rhoUx[1]+8.0*rhoLx[1]+32.0*rhoC[1]+8.0*rhoUxLy[0]+32.0*rhoUx[0]+32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*bcVals[9]*bcVals[10])*omega*volFac+(((54.0*rdx2SqVol[1]+54.0*rdx2SqVol[0])*bcVals[10]+(96.0*rdx2SqVol[1]+96.0*rdx2SqVol[0])*bcVals[9])*bcVals[11]+(((-12.0*phiPrevC[1])+3.0*phiLx[1]+3.0*phiUxLy[0]-3.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]-30.0*rdx2SqVol[0]*phiPrevC[1]+3.0*rdx2SqVol[0]*phiLx[1]+(3.0*phiUxLy[0]+15.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*bcVals10R2+(((-72.0*phiPrevC[1])+6.0*phiLx[1]+6.0*phiUxLy[0]-12.0*phiUx[0]+24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*rdx2SqVol[1]-108.0*rdx2SqVol[0]*phiPrevC[1]+6.0*rdx2SqVol[0]*phiLx[1]+(6.0*phiUxLy[0]+24.0*phiUx[0]-12.0*phiLy[0]+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0])*bcVals[9]*bcVals[10]+((-96.0*phiPrevC[1]*rdx2SqVol[1])-96.0*rdx2SqVol[0]*phiPrevC[1])*bcVals9R2)*omega+(12.0*phiPrevC[1]*rdx2SqVol[1]+30.0*rdx2SqVol[0]*phiPrevC[1])*bcVals10R2+(72.0*phiPrevC[1]*rdx2SqVol[1]+108.0*rdx2SqVol[0]*phiPrevC[1])*bcVals[9]*bcVals[10]+(96.0*phiPrevC[1]*rdx2SqVol[1]+96.0*rdx2SqVol[0]*phiPrevC[1])*bcVals9R2)/((12.0*rdx2SqVol[1]+30.0*rdx2SqVol[0])*bcVals10R2+(72.0*rdx2SqVol[1]+108.0*rdx2SqVol[0])*bcVals[9]*bcVals[10]+(96.0*rdx2SqVol[1]+96.0*rdx2SqVol[0])*bcVals9R2); 

}

void MGpoissonFEMDampedJacobi2xSer_LxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((dxC[1]*bcVals[8]+dxC[0]*bcVals[2]-1.0*phiPrevC[0]*dxC[1]-1.0*dxC[0]*phiPrevC[0])*omega+phiPrevC[0]*dxC[1]+dxC[0]*phiPrevC[0])/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_LxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxDirichletLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((-2.0*bcVals[8])-2.0*bcVals[2]+phiUxUy[0]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxNeumannLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*(((2.0*dxC[1]*bcVals[7]-4.0*dxC[1]*bcVals[6]+2.0*dxC[0])*bcVals[8]+(2.0*dxC[1]*bcVals[2]+(phiPrevC[0]-1.0*phiUxUy[0])*dxC[1])*bcVals7R2+(((2.0*phiUxUy[0]-4.0*phiPrevC[0])*dxC[1]-4.0*dxC[1]*bcVals[2])*bcVals[6]+2.0*dxC[0]*bcVals[2]-1.0*dxC[0]*phiUxUy[0]+dxC[0]*phiPrevC[0])*bcVals[7]+4.0*phiPrevC[0]*dxC[1]*bcVals6R2+((-4.0*dxC[0]*bcVals[2])-2.0*dxC[0]*phiPrevC[0])*bcVals[6])*omega-1.0*phiPrevC[0]*dxC[1]*bcVals7R2+(4.0*phiPrevC[0]*dxC[1]*bcVals[6]-1.0*dxC[0]*phiPrevC[0])*bcVals[7]-4.0*phiPrevC[0]*dxC[1]*bcVals6R2+2.0*dxC[0]*phiPrevC[0]*bcVals[6]))/(dxC[1]*bcVals7R2+(dxC[0]-4.0*dxC[1]*bcVals[6])*bcVals[7]+4.0*dxC[1]*bcVals6R2-2.0*dxC[0]*bcVals[6]); 

}

void MGpoissonFEMDampedJacobi2xSer_LxRobinLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxRobinLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);

  phiC[0] = -(1.0*((((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]+2.0*dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[8]+(2.0*dxC[1]+2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+((phiPrevC[0]-1.0*phiUxUy[0])*bcVals[1]-2.0*bcVals[0]*phiPrevC[0])*dxC[1]+(dxC[0]*phiPrevC[0]-1.0*dxC[0]*phiUxUy[0])*bcVals1R2+(2.0*bcVals[0]*dxC[0]*phiUxUy[0]-4.0*bcVals[0]*dxC[0]*phiPrevC[0])*bcVals[1]+4.0*bcVals0R2*dxC[0]*phiPrevC[0])*omega+(2.0*bcVals[0]*phiPrevC[0]-1.0*phiPrevC[0]*bcVals[1])*dxC[1]-1.0*dxC[0]*phiPrevC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*phiPrevC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0]*phiPrevC[0]))/((bcVals[1]-2.0*bcVals[0])*dxC[1]+dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_LxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);
  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*((((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]*bcVals[7]+(8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*bcVals[6]+2.0*dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[8]+(2.0*dxC[1]*bcVals[2]+((phiPrevC[0]-1.0*phiUxUy[0])*bcVals[1]-2.0*bcVals[0]*phiPrevC[0])*dxC[1])*bcVals7R2+((((2.0*phiUxUy[0]-4.0*phiPrevC[0])*bcVals[1]+8.0*bcVals[0]*phiPrevC[0])*dxC[1]-4.0*dxC[1]*bcVals[2])*bcVals[6]+(2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+(dxC[0]*phiPrevC[0]-1.0*dxC[0]*phiUxUy[0])*bcVals1R2+(2.0*bcVals[0]*dxC[0]*phiUxUy[0]-4.0*bcVals[0]*dxC[0]*phiPrevC[0])*bcVals[1]+4.0*bcVals0R2*dxC[0]*phiPrevC[0])*bcVals[7]+(4.0*phiPrevC[0]*bcVals[1]-8.0*bcVals[0]*phiPrevC[0])*dxC[1]*bcVals6R2+((8.0*bcVals[0]*dxC[0]-4.0*dxC[0]*bcVals[1])*bcVals[2]-2.0*dxC[0]*phiPrevC[0]*bcVals1R2+8.0*bcVals[0]*dxC[0]*phiPrevC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0]*phiPrevC[0])*bcVals[6])*omega+(2.0*bcVals[0]*phiPrevC[0]-1.0*phiPrevC[0]*bcVals[1])*dxC[1]*bcVals7R2+((4.0*phiPrevC[0]*bcVals[1]-8.0*bcVals[0]*phiPrevC[0])*dxC[1]*bcVals[6]-1.0*dxC[0]*phiPrevC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*phiPrevC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0]*phiPrevC[0])*bcVals[7]+(8.0*bcVals[0]*phiPrevC[0]-4.0*phiPrevC[0]*bcVals[1])*dxC[1]*bcVals6R2+(2.0*dxC[0]*phiPrevC[0]*bcVals1R2-8.0*bcVals[0]*dxC[0]*phiPrevC[0]*bcVals[1]+8.0*bcVals0R2*dxC[0]*phiPrevC[0])*bcVals[6]))/((bcVals[1]-2.0*bcVals[0])*dxC[1]*bcVals7R2+((8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*bcVals[6]+dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0])*bcVals[7]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals6R2+((-2.0*dxC[0]*bcVals1R2)+8.0*bcVals[0]*dxC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0])*bcVals[6]); 

}

void MGpoissonFEMDampedJacobi2xSer_LxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 
  phiC[1] = ((dxC[1]*bcVals[11]+dxC[0]*bcVals[2]+((-1.0*dxC[1])-1.0*dxC[0])*phiPrevC[1])*omega+(dxC[1]+dxC[0])*phiPrevC[1])/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_LxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 
  phiC[1] = (bcVals[2]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxDirichletUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 
  phiC[1] = (bcVals[2]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUx[1]+8.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+8.0*rhoLy[0]+32.0*rhoC[0])*omega*volFac+(9.0*rdx2SqVol[1]*bcVals[11]+(3.0*phiUxLy[0]-6.0*phiUx[0]-12.0*phiPrevC[0]+6.0*phiLy[0])*rdx2SqVol[1]+(3.0*phiUxLy[0]+6.0*phiUx[0]-6.0*phiPrevC[0]-3.0*phiLy[0])*rdx2SqVol[0])*omega+12.0*phiPrevC[0]*rdx2SqVol[1]+6.0*phiPrevC[0]*rdx2SqVol[0])/(12.0*rdx2SqVol[1]+6.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[11]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((-2.0*bcVals[2])+phiUx[0]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 
  phiC[1] = (2.0*bcVals[11]-2.0*bcVals[2]-1.0*phiPrevC[1]+phiUx[0])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxNeumannUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = ((-2.0*bcVals[2])+phiUx[0]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 
  phiC[1] = (((2.0*dxC[1]*bcVals[10]+4.0*dxC[1]*bcVals[9]-2.0*dxC[0])*bcVals[11]+((-2.0*dxC[1]*bcVals[2])-1.0*dxC[1]*phiPrevC[1]+phiUx[0]*dxC[1])*bcVals10R2+(((-4.0*dxC[1]*bcVals[2])-4.0*dxC[1]*phiPrevC[1]+2.0*phiUx[0]*dxC[1])*bcVals[9]+2.0*dxC[0]*bcVals[2]+dxC[0]*phiPrevC[1]-1.0*dxC[0]*phiUx[0])*bcVals[10]-4.0*dxC[1]*phiPrevC[1]*bcVals9R2+(4.0*dxC[0]*bcVals[2]+2.0*dxC[0]*phiPrevC[1])*bcVals[9])*omega+dxC[1]*phiPrevC[1]*bcVals10R2+(4.0*dxC[1]*phiPrevC[1]*bcVals[9]-1.0*dxC[0]*phiPrevC[1])*bcVals[10]+4.0*dxC[1]*phiPrevC[1]*bcVals9R2-2.0*dxC[0]*phiPrevC[1]*bcVals[9])/(dxC[1]*bcVals10R2+(4.0*dxC[1]*bcVals[9]-1.0*dxC[0])*bcVals[10]+4.0*dxC[1]*bcVals9R2-2.0*dxC[0]*bcVals[9]); 

}

void MGpoissonFEMDampedJacobi2xSer_LxRobinUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUx[1]+8.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+8.0*rhoLy[0]+32.0*rhoC[0])*omega*volFac+(9.0*rdx2SqVol[1]*bcVals[11]+(3.0*phiUxLy[0]-6.0*phiUx[0]-12.0*phiPrevC[0]+6.0*phiLy[0])*rdx2SqVol[1]+(3.0*phiUxLy[0]+6.0*phiUx[0]-6.0*phiPrevC[0]-3.0*phiLy[0])*rdx2SqVol[0])*omega+12.0*phiPrevC[0]*rdx2SqVol[1]+6.0*phiPrevC[0]*rdx2SqVol[0])/(12.0*rdx2SqVol[1]+6.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[11]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_LxRobinUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);

  phiC[0] = -(1.0*((2.0*bcVals[2]+(phiPrevC[0]-1.0*phiUx[0])*bcVals[1]-2.0*bcVals[0]*phiPrevC[0])*omega-1.0*phiPrevC[0]*bcVals[1]+2.0*bcVals[0]*phiPrevC[0]))/(bcVals[1]-2.0*bcVals[0]); 
  phiC[1] = ((((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]-2.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[11]+((-2.0*dxC[1])+2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+((2.0*bcVals[0]-1.0*bcVals[1])*dxC[1]+dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0])*phiPrevC[1]+phiUx[0]*bcVals[1]*dxC[1]-1.0*dxC[0]*phiUx[0]*bcVals1R2+2.0*bcVals[0]*dxC[0]*phiUx[0]*bcVals[1])*omega+((bcVals[1]-2.0*bcVals[0])*dxC[1]-1.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0])*phiPrevC[1])/((bcVals[1]-2.0*bcVals[0])*dxC[1]-1.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_LxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);
  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = -(1.0*((2.0*bcVals[2]+(phiPrevC[0]-1.0*phiUx[0])*bcVals[1]-2.0*bcVals[0]*phiPrevC[0])*omega-1.0*phiPrevC[0]*bcVals[1]+2.0*bcVals[0]*phiPrevC[0]))/(bcVals[1]-2.0*bcVals[0]); 
  phiC[1] = ((((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]*bcVals[10]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals[9]-2.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[11]+((-2.0*dxC[1]*bcVals[2])+(2.0*bcVals[0]-1.0*bcVals[1])*dxC[1]*phiPrevC[1]+phiUx[0]*bcVals[1]*dxC[1])*bcVals10R2+(((-4.0*dxC[1]*bcVals[2])+(8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*phiPrevC[1]+2.0*phiUx[0]*bcVals[1]*dxC[1])*bcVals[9]+(2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+(dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0])*phiPrevC[1]-1.0*dxC[0]*phiUx[0]*bcVals1R2+2.0*bcVals[0]*dxC[0]*phiUx[0]*bcVals[1])*bcVals[10]+(8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*phiPrevC[1]*bcVals9R2+((4.0*dxC[0]*bcVals[1]-8.0*bcVals[0]*dxC[0])*bcVals[2]+(2.0*dxC[0]*bcVals1R2-8.0*bcVals[0]*dxC[0]*bcVals[1]+8.0*bcVals0R2*dxC[0])*phiPrevC[1])*bcVals[9])*omega+(bcVals[1]-2.0*bcVals[0])*dxC[1]*phiPrevC[1]*bcVals10R2+((4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*phiPrevC[1]*bcVals[9]+((-1.0*dxC[0]*bcVals1R2)+4.0*bcVals[0]*dxC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0])*phiPrevC[1])*bcVals[10]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*phiPrevC[1]*bcVals9R2+((-2.0*dxC[0]*bcVals1R2)+8.0*bcVals[0]*dxC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0])*phiPrevC[1]*bcVals[9])/((bcVals[1]-2.0*bcVals[0])*dxC[1]*bcVals10R2+((4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals[9]-1.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0])*bcVals[10]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals9R2+((-2.0*dxC[0]*bcVals1R2)+8.0*bcVals[0]*dxC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0])*bcVals[9]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 
  phiC[1] = ((dxC[1]*bcVals[8]+dxC[0]*bcVals[5]+((-1.0*dxC[1])-1.0*dxC[0])*phiPrevC[1])*omega+(dxC[1]+dxC[0])*phiPrevC[1])/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUy[1]+8.0*rhoC[1]+16.0*rhoUy[0]+4.0*rhoLxUy[0]+8.0*rhoLx[0]+32.0*rhoC[0])*omega*volFac+(9.0*rdx2SqVol[0]*bcVals[5]+(6.0*phiUy[0]-6.0*phiPrevC[0]+3.0*phiLxUy[0]-3.0*phiLx[0])*rdx2SqVol[1]+((-6.0*phiUy[0])-12.0*phiPrevC[0]+3.0*phiLxUy[0]+6.0*phiLx[0])*rdx2SqVol[0])*omega+6.0*phiPrevC[0]*rdx2SqVol[1]+12.0*phiPrevC[0]*rdx2SqVol[0])/(6.0*rdx2SqVol[1]+12.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[5]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_UxDirichletLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUy[1]+8.0*rhoC[1]+16.0*rhoUy[0]+4.0*rhoLxUy[0]+8.0*rhoLx[0]+32.0*rhoC[0])*omega*volFac+(9.0*rdx2SqVol[0]*bcVals[5]+(6.0*phiUy[0]-6.0*phiPrevC[0]+3.0*phiLxUy[0]-3.0*phiLx[0])*rdx2SqVol[1]+((-6.0*phiUy[0])-12.0*phiPrevC[0]+3.0*phiLxUy[0]+6.0*phiLx[0])*rdx2SqVol[0])*omega+6.0*phiPrevC[0]*rdx2SqVol[1]+12.0*phiPrevC[0]*rdx2SqVol[0])/(6.0*rdx2SqVol[1]+12.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[5]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_UxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 
  phiC[1] = (bcVals[8]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_UxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((-2.0*bcVals[8])+phiUy[0]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 
  phiC[1] = ((-2.0*bcVals[8])+2.0*bcVals[5]-1.0*phiPrevC[1]+phiUy[0])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_UxNeumannLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*((2.0*bcVals[8]+(phiPrevC[0]-1.0*phiUy[0])*bcVals[7]-2.0*phiPrevC[0]*bcVals[6])*omega-1.0*phiPrevC[0]*bcVals[7]+2.0*phiPrevC[0]*bcVals[6]))/(bcVals[7]-2.0*bcVals[6]); 
  phiC[1] = -(1.0*(((2.0*dxC[1]*bcVals[7]-4.0*dxC[1]*bcVals[6]-2.0*dxC[0])*bcVals[8]+((-2.0*dxC[1]*bcVals[5])+dxC[1]*phiPrevC[1]-1.0*phiUy[0]*dxC[1])*bcVals7R2+((4.0*dxC[1]*bcVals[5]-4.0*dxC[1]*phiPrevC[1]+2.0*phiUy[0]*dxC[1])*bcVals[6]+2.0*dxC[0]*bcVals[5]-1.0*dxC[0]*phiPrevC[1]+dxC[0]*phiUy[0])*bcVals[7]+4.0*dxC[1]*phiPrevC[1]*bcVals6R2+(2.0*dxC[0]*phiPrevC[1]-4.0*dxC[0]*bcVals[5])*bcVals[6])*omega-1.0*dxC[1]*phiPrevC[1]*bcVals7R2+(4.0*dxC[1]*phiPrevC[1]*bcVals[6]+dxC[0]*phiPrevC[1])*bcVals[7]-4.0*dxC[1]*phiPrevC[1]*bcVals6R2-2.0*dxC[0]*phiPrevC[1]*bcVals[6]))/(dxC[1]*bcVals7R2+((-4.0*dxC[1]*bcVals[6])-1.0*dxC[0])*bcVals[7]+4.0*dxC[1]*bcVals6R2+2.0*dxC[0]*bcVals[6]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxRobinLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 
  phiC[1] = (bcVals[8]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 

}

void MGpoissonFEMDampedJacobi2xSer_UxRobinLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);

  phiC[0] = ((-2.0*bcVals[8])+phiUy[0]-1.0*phiPrevC[0])*omega+phiPrevC[0]; 
  phiC[1] = -(1.0*(((2.0*dxC[0]*bcVals4R2+(4.0*dxC[0]*bcVals[3]-2.0*dxC[1])*bcVals[4]-4.0*dxC[1]*bcVals[3])*bcVals[8]+((-2.0*dxC[0]*bcVals[4])-4.0*dxC[0]*bcVals[3]+2.0*dxC[1])*bcVals[5]+(dxC[0]*phiPrevC[1]-1.0*dxC[0]*phiUy[0])*bcVals4R2+((4.0*dxC[0]*phiPrevC[1]-2.0*dxC[0]*phiUy[0])*bcVals[3]-1.0*dxC[1]*phiPrevC[1]+phiUy[0]*dxC[1])*bcVals[4]+4.0*dxC[0]*phiPrevC[1]*bcVals3R2-2.0*dxC[1]*phiPrevC[1]*bcVals[3])*omega-1.0*dxC[0]*phiPrevC[1]*bcVals4R2+(dxC[1]*phiPrevC[1]-4.0*dxC[0]*phiPrevC[1]*bcVals[3])*bcVals[4]-4.0*dxC[0]*phiPrevC[1]*bcVals3R2+2.0*dxC[1]*phiPrevC[1]*bcVals[3]))/(dxC[0]*bcVals4R2+(4.0*dxC[0]*bcVals[3]-1.0*dxC[1])*bcVals[4]+4.0*dxC[0]*bcVals3R2-2.0*dxC[1]*bcVals[3]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);
  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*((2.0*bcVals[8]+(phiPrevC[0]-1.0*phiUy[0])*bcVals[7]-2.0*phiPrevC[0]*bcVals[6])*omega-1.0*phiPrevC[0]*bcVals[7]+2.0*phiPrevC[0]*bcVals[6]))/(bcVals[7]-2.0*bcVals[6]); 
  phiC[1] = -(1.0*((((2.0*dxC[1]*bcVals[4]+4.0*dxC[1]*bcVals[3])*bcVals[7]+((-4.0*dxC[1]*bcVals[4])-8.0*dxC[1]*bcVals[3])*bcVals[6]-2.0*dxC[0]*bcVals4R2-4.0*dxC[0]*bcVals[3]*bcVals[4])*bcVals[8]+((-2.0*dxC[1]*bcVals[5])+(dxC[1]*phiPrevC[1]-1.0*phiUy[0]*dxC[1])*bcVals[4]+2.0*dxC[1]*phiPrevC[1]*bcVals[3])*bcVals7R2+((4.0*dxC[1]*bcVals[5]+(2.0*phiUy[0]*dxC[1]-4.0*dxC[1]*phiPrevC[1])*bcVals[4]-8.0*dxC[1]*phiPrevC[1]*bcVals[3])*bcVals[6]+(2.0*dxC[0]*bcVals[4]+4.0*dxC[0]*bcVals[3])*bcVals[5]+(dxC[0]*phiUy[0]-1.0*dxC[0]*phiPrevC[1])*bcVals4R2+(2.0*dxC[0]*phiUy[0]-4.0*dxC[0]*phiPrevC[1])*bcVals[3]*bcVals[4]-4.0*dxC[0]*phiPrevC[1]*bcVals3R2)*bcVals[7]+(4.0*dxC[1]*phiPrevC[1]*bcVals[4]+8.0*dxC[1]*phiPrevC[1]*bcVals[3])*bcVals6R2+(((-4.0*dxC[0]*bcVals[4])-8.0*dxC[0]*bcVals[3])*bcVals[5]+2.0*dxC[0]*phiPrevC[1]*bcVals4R2+8.0*dxC[0]*phiPrevC[1]*bcVals[3]*bcVals[4]+8.0*dxC[0]*phiPrevC[1]*bcVals3R2)*bcVals[6])*omega+((-1.0*dxC[1]*phiPrevC[1]*bcVals[4])-2.0*dxC[1]*phiPrevC[1]*bcVals[3])*bcVals7R2+((4.0*dxC[1]*phiPrevC[1]*bcVals[4]+8.0*dxC[1]*phiPrevC[1]*bcVals[3])*bcVals[6]+dxC[0]*phiPrevC[1]*bcVals4R2+4.0*dxC[0]*phiPrevC[1]*bcVals[3]*bcVals[4]+4.0*dxC[0]*phiPrevC[1]*bcVals3R2)*bcVals[7]+((-4.0*dxC[1]*phiPrevC[1]*bcVals[4])-8.0*dxC[1]*phiPrevC[1]*bcVals[3])*bcVals6R2+((-2.0*dxC[0]*phiPrevC[1]*bcVals4R2)-8.0*dxC[0]*phiPrevC[1]*bcVals[3]*bcVals[4]-8.0*dxC[0]*phiPrevC[1]*bcVals3R2)*bcVals[6]))/((dxC[1]*bcVals[4]+2.0*dxC[1]*bcVals[3])*bcVals7R2+(((-4.0*dxC[1]*bcVals[4])-8.0*dxC[1]*bcVals[3])*bcVals[6]-1.0*dxC[0]*bcVals4R2-4.0*dxC[0]*bcVals[3]*bcVals[4]-4.0*dxC[0]*bcVals3R2)*bcVals[7]+(4.0*dxC[1]*bcVals[4]+8.0*dxC[1]*bcVals[3])*bcVals6R2+(2.0*dxC[0]*bcVals4R2+8.0*dxC[0]*bcVals[3]*bcVals[4]+8.0*dxC[0]*bcVals3R2)*bcVals[6]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (((4.0*dxC[1]+4.0*dxC[0])*rhoC[3]+(16.0*dxC[1]+16.0*dxC[0])*rhoC[2]+(4.0*dxC[1]+4.0*dxC[0])*rhoLy[1]+(4.0*dxC[1]+4.0*dxC[0])*rhoLx[1]+(16.0*dxC[1]+16.0*dxC[0])*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*omega*volFac+(((15.0*dxC[1]+12.0*dxC[0])*rdx2SqVol[1]-3.0*rdx2SqVol[0]*dxC[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+(((-6.0*dxC[1])-3.0*dxC[0])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+((3.0*dxC[1]+3.0*dxC[0])*phiLy[1]+(3.0*dxC[1]+3.0*dxC[0])*phiLx[1]+((-24.0*phiPrevC[0])+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1]-24.0*dxC[0]*phiPrevC[0]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+(3.0*rdx2SqVol[0]*dxC[1]+3.0*dxC[0]*rdx2SqVol[0])*phiLy[1]+(3.0*rdx2SqVol[0]*dxC[1]+3.0*dxC[0]*rdx2SqVol[0])*phiLx[1]+((-24.0*phiPrevC[0])-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1]+((-24.0*dxC[0]*phiPrevC[0])-6.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*omega+(24.0*phiPrevC[0]*dxC[1]+24.0*dxC[0]*phiPrevC[0])*rdx2SqVol[1]+24.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1]+24.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])/((24.0*dxC[1]+24.0*dxC[0])*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1]+24.0*dxC[0]*rdx2SqVol[0]); 
  phiC[1] = (bcVals[5]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 
  phiC[2] = (bcVals[11]-1.0*phiPrevC[2])*omega+phiPrevC[2]; 
  phiC[3] = ((dxC[1]*bcVals[11]+dxC[0]*bcVals[5]+((-1.0*dxC[1])-1.0*dxC[0])*phiPrevC[3])*omega+(dxC[1]+dxC[0])*phiPrevC[3])/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = ((20.0*rdx2SqVol[1]*rhoC[3]+80.0*rdx2SqVol[1]*rhoC[2]+(4.0*rdx2SqVol[1]+8.0*rdx2SqVol[0])*rhoLy[1]+20.0*rdx2SqVol[1]*rhoLx[1]+(24.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+24.0*rhoLx[0]+96.0*rhoC[0])*rdx2SqVol[1]+32.0*rdx2SqVol[0]*rhoLy[0]+8.0*rdx2SqVol[0]*rhoLxLy[0]+28.0*rdx2SqVol[0]*rhoLx[0]+112.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+(((-3.0*rdx2SqVol1R2)+27.0*rdx2SqVol[0]*rdx2SqVol[1]+21.0*rdx2SqVol0R2)*bcVals[5]+(3.0*phiLy[1]-3.0*phiLx[1]-12.0*phiPrevC[0]+12.0*phiLy[0]+3.0*phiLxLy[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiLy[1]+24.0*rdx2SqVol[0]*phiLx[1]+((-90.0*phiPrevC[0])+18.0*phiLy[0]+9.0*phiLxLy[0]+3.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+6.0*rdx2SqVol0R2*phiLy[1]+((-42.0*phiPrevC[0])-12.0*phiLy[0]+6.0*phiLxLy[0]+21.0*phiLx[0])*rdx2SqVol0R2)*omega+12.0*phiPrevC[0]*rdx2SqVol1R2+90.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+42.0*phiPrevC[0]*rdx2SqVol0R2)/(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2); 
  phiC[1] = (bcVals[5]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 
  phiC[2] = (((36.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoC[3]+(144.0*rdx2SqVol[1]+112.0*rdx2SqVol[0])*rhoC[2]+(4.0*rdx2SqVol[1]-4.0*rdx2SqVol[0])*rhoLy[1]+(36.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoLx[1]+32.0*rdx2SqVol[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*rdx2SqVol[1]-16.0*rdx2SqVol[0]*rhoLy[0]-4.0*rdx2SqVol[0]*rhoLxLy[0])*omega*volFac+(((-3.0*rdx2SqVol1R2)+54.0*rdx2SqVol[0]*rdx2SqVol[1]+21.0*rdx2SqVol0R2)*bcVals[5]+((-12.0*rdx2SqVol1R2)-90.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*phiPrevC[2]+(3.0*phiLy[1]-9.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]+6.0*phiLx[0])*rdx2SqVol1R2+(12.0*rdx2SqVol[0]*phiLx[1]+(42.0*phiLx[0]-18.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1]-3.0*rdx2SqVol0R2*phiLy[1]+21.0*rdx2SqVol0R2*phiLx[1]+(6.0*phiLy[0]-3.0*phiLxLy[0])*rdx2SqVol0R2)*omega+(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2)*phiPrevC[2])/(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2); 
  phiC[3] = (bcVals[5]-1.0*phiPrevC[3])*omega+phiPrevC[3]; 

}

void MGpoissonFEMDampedJacobi2xSer_UxDirichletUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = ((20.0*rdx2SqVol[1]*rhoC[3]+80.0*rdx2SqVol[1]*rhoC[2]+(4.0*rdx2SqVol[1]+8.0*rdx2SqVol[0])*rhoLy[1]+20.0*rdx2SqVol[1]*rhoLx[1]+(24.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+24.0*rhoLx[0]+96.0*rhoC[0])*rdx2SqVol[1]+32.0*rdx2SqVol[0]*rhoLy[0]+8.0*rdx2SqVol[0]*rhoLxLy[0]+28.0*rdx2SqVol[0]*rhoLx[0]+112.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+(((-3.0*rdx2SqVol1R2)+27.0*rdx2SqVol[0]*rdx2SqVol[1]+21.0*rdx2SqVol0R2)*bcVals[5]+(3.0*phiLy[1]-3.0*phiLx[1]-12.0*phiPrevC[0]+12.0*phiLy[0]+3.0*phiLxLy[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiLy[1]+24.0*rdx2SqVol[0]*phiLx[1]+((-90.0*phiPrevC[0])+18.0*phiLy[0]+9.0*phiLxLy[0]+3.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+6.0*rdx2SqVol0R2*phiLy[1]+((-42.0*phiPrevC[0])-12.0*phiLy[0]+6.0*phiLxLy[0]+21.0*phiLx[0])*rdx2SqVol0R2)*omega+12.0*phiPrevC[0]*rdx2SqVol1R2+90.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+42.0*phiPrevC[0]*rdx2SqVol0R2)/(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2); 
  phiC[1] = (bcVals[5]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 
  phiC[2] = (((36.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoC[3]+(144.0*rdx2SqVol[1]+112.0*rdx2SqVol[0])*rhoC[2]+(4.0*rdx2SqVol[1]-4.0*rdx2SqVol[0])*rhoLy[1]+(36.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoLx[1]+32.0*rdx2SqVol[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*rdx2SqVol[1]-16.0*rdx2SqVol[0]*rhoLy[0]-4.0*rdx2SqVol[0]*rhoLxLy[0])*omega*volFac+(((-3.0*rdx2SqVol1R2)+54.0*rdx2SqVol[0]*rdx2SqVol[1]+21.0*rdx2SqVol0R2)*bcVals[5]+((-12.0*rdx2SqVol1R2)-90.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*phiPrevC[2]+(3.0*phiLy[1]-9.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]+6.0*phiLx[0])*rdx2SqVol1R2+(12.0*rdx2SqVol[0]*phiLx[1]+(42.0*phiLx[0]-18.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1]-3.0*rdx2SqVol0R2*phiLy[1]+21.0*rdx2SqVol0R2*phiLx[1]+(6.0*phiLy[0]-3.0*phiLxLy[0])*rdx2SqVol0R2)*omega+(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2)*phiPrevC[2])/(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2); 
  phiC[3] = (bcVals[5]-1.0*phiPrevC[3])*omega+phiPrevC[3]; 

}

void MGpoissonFEMDampedJacobi2xSer_UxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = ((20.0*rdx2SqVol[0]*rhoC[3]+(28.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[2]+20.0*rdx2SqVol[0]*rhoLy[1]+(8.0*rdx2SqVol[1]+4.0*rdx2SqVol[0])*rhoLx[1]+80.0*rdx2SqVol[0]*rhoC[1]+(28.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+112.0*rhoC[0])*rdx2SqVol[1]+24.0*rdx2SqVol[0]*rhoLy[0]+4.0*rdx2SqVol[0]*rhoLxLy[0]+16.0*rdx2SqVol[0]*rhoLx[0]+96.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+((21.0*rdx2SqVol1R2+27.0*rdx2SqVol[0]*rdx2SqVol[1]-3.0*rdx2SqVol0R2)*bcVals[11]+(6.0*phiLx[1]-42.0*phiPrevC[0]+21.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*rdx2SqVol1R2+(24.0*rdx2SqVol[0]*phiLy[1]+9.0*rdx2SqVol[0]*phiLx[1]+((-90.0*phiPrevC[0])+3.0*phiLy[0]+9.0*phiLxLy[0]+18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-3.0*rdx2SqVol0R2*phiLy[1]+3.0*rdx2SqVol0R2*phiLx[1]+((-12.0*phiPrevC[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol0R2)*omega+42.0*phiPrevC[0]*rdx2SqVol1R2+90.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+12.0*phiPrevC[0]*rdx2SqVol0R2)/(42.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+12.0*rdx2SqVol0R2); 
  phiC[1] = (((28.0*rdx2SqVol[1]+36.0*rdx2SqVol[0])*rhoC[3]+32.0*rdx2SqVol[0]*rhoC[2]+(28.0*rdx2SqVol[1]+36.0*rdx2SqVol[0])*rhoLy[1]+(4.0*rdx2SqVol[0]-4.0*rdx2SqVol[1])*rhoLx[1]+(112.0*rdx2SqVol[1]+144.0*rdx2SqVol[0])*rhoC[1]+((-4.0*rhoLxLy[0])-16.0*rhoLx[0])*rdx2SqVol[1]+32.0*rdx2SqVol[0]*rhoLy[0]+4.0*rdx2SqVol[0]*rhoLxLy[0]+16.0*rdx2SqVol[0]*rhoLx[0]+128.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+((21.0*rdx2SqVol1R2+54.0*rdx2SqVol[0]*rdx2SqVol[1]-3.0*rdx2SqVol0R2)*bcVals[11]+((-42.0*phiPrevC[1])+21.0*phiLy[1]-3.0*phiLx[1]-3.0*phiLxLy[0]+6.0*phiLx[0])*rdx2SqVol1R2+((-90.0*rdx2SqVol[0]*phiPrevC[1])+12.0*rdx2SqVol[0]*phiLy[1]+(42.0*phiLy[0]-18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol0R2*phiPrevC[1]-9.0*rdx2SqVol0R2*phiLy[1]+3.0*rdx2SqVol0R2*phiLx[1]+(6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol0R2)*omega+42.0*phiPrevC[1]*rdx2SqVol1R2+90.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+12.0*rdx2SqVol0R2*phiPrevC[1])/(42.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+12.0*rdx2SqVol0R2); 
  phiC[2] = (bcVals[11]-1.0*phiPrevC[2])*omega+phiPrevC[2]; 
  phiC[3] = (bcVals[11]-1.0*phiPrevC[3])*omega+phiPrevC[3]; 

}

void MGpoissonFEMDampedJacobi2xSer_UxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((30.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*bcVals[11]+(30.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*bcVals[5]+(3.0*phiLy[1]+3.0*phiLx[1]-15.0*phiPrevC[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+((-15.0*phiPrevC[0])-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+15.0*phiPrevC[0]*rdx2SqVol[1]+15.0*phiPrevC[0]*rdx2SqVol[0])/(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0]); 
  phiC[1] = ((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((30.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*bcVals[11]+(24.0*rdx2SqVol[1]+60.0*rdx2SqVol[0])*bcVals[5]+((-15.0*phiPrevC[1])+3.0*phiLy[1]+3.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]-15.0*rdx2SqVol[0]*phiPrevC[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+15.0*phiPrevC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*phiPrevC[1])/(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0]); 
  phiC[2] = ((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((60.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*bcVals[11]+(30.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*bcVals[5]+((-15.0*rdx2SqVol[1])-15.0*rdx2SqVol[0])*phiPrevC[2]+(3.0*phiLy[1]+3.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0])*phiPrevC[2])/(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0]); 
  phiC[3] = ((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((60.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*bcVals[11]+(24.0*rdx2SqVol[1]+60.0*rdx2SqVol[0])*bcVals[5]+((-15.0*rdx2SqVol[1])-15.0*rdx2SqVol[0])*phiPrevC[3]+(3.0*phiLy[1]+3.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0])*phiPrevC[3])/(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxNeumannUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals10R2+((16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[9]+4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals[10]+(16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals9R2+(8.0*dxC[0]*rhoC[3]+32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[9])*omega*volFac+(((30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[10]+(60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+((30.0*rdx2SqVol[0]*dxC[1]-6.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+((-15.0*phiPrevC[0])+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-15.0*phiPrevC[0])-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals10R2+(((108.0*rdx2SqVol[0]*dxC[1]-36.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+(12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+((-90.0*phiPrevC[0])+48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-54.0*phiPrevC[0])-24.0*phiLy[0]+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[9]+(30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[5]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]-15.0*dxC[0]*phiPrevC[0]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-15.0*dxC[0]*phiPrevC[0])-6.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[10]+((96.0*rdx2SqVol[0]*dxC[1]-48.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+(12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+((-120.0*phiPrevC[0])+48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-48.0*phiPrevC[0])-24.0*phiLy[0]+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[5]+(6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]-60.0*dxC[0]*phiPrevC[0]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiPrevC[0])-12.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[9])*omega+(15.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+15.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals10R2+((90.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+54.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals[9]+15.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+15.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals[10]+(120.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+48.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals9R2+(60.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+24.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals[9])/((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[9]); 
  phiC[1] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals10R2+((16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[9]+4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals[10]+(16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals9R2+(8.0*dxC[0]*rhoC[3]+32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[9])*omega*volFac+(((30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[10]+(60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+((24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((-15.0*dxC[1]*phiPrevC[1])+3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]-15.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals10R2+(((144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((-90.0*dxC[1]*phiPrevC[1])+12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]-54.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[9]+(24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+((-15.0*dxC[0]*phiPrevC[1])+3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]-15.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[10]+((192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((-120.0*dxC[1]*phiPrevC[1])+12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]-48.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals9R2+((108.0*dxC[0]*rdx2SqVol[1]+108.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+((-60.0*dxC[0]*phiPrevC[1])+6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]-24.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[9])*omega+(15.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals10R2+((90.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[9]+15.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals[10]+(120.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals9R2+(60.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals[9])/((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[9]); 
  phiC[2] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals10R2+((8.0*dxC[1]*rhoC[3]+32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[9]+4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals[10])*omega*volFac+(((60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+((30.0*rdx2SqVol[0]*dxC[1]-6.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+((-15.0*dxC[1]*rdx2SqVol[1])-15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals10R2+(((48.0*rdx2SqVol[0]*dxC[1]-24.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+((-90.0*dxC[1]*rdx2SqVol[1])-54.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[9]+(30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[5]+((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[10]+((-120.0*dxC[1]*rdx2SqVol[1])-48.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals9R2+((-60.0*dxC[0]*rdx2SqVol[1])-24.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals[9])*omega+(15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals[9])/((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[9]); 
  phiC[3] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals10R2+((8.0*dxC[1]*rhoC[3]+32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[9]+4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals[10])*omega*volFac+(((60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+((24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((-15.0*dxC[1]*rdx2SqVol[1])-15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals10R2+(((96.0*dxC[1]*rdx2SqVol[1]+96.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((-90.0*dxC[1]*rdx2SqVol[1])-54.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[9]+(24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[3]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[10]+((-120.0*dxC[1]*rdx2SqVol[1])-48.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]*bcVals9R2+((120.0*dxC[0]*rdx2SqVol[1]+48.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+((-60.0*dxC[0]*rdx2SqVol[1])-24.0*dxC[0]*rdx2SqVol[0])*phiPrevC[3])*bcVals[9])*omega+(15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[3])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*phiPrevC[3]*bcVals[9])/((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[9]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxRobinUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = ((20.0*rdx2SqVol[0]*rhoC[3]+(28.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[2]+20.0*rdx2SqVol[0]*rhoLy[1]+(8.0*rdx2SqVol[1]+4.0*rdx2SqVol[0])*rhoLx[1]+80.0*rdx2SqVol[0]*rhoC[1]+(28.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+112.0*rhoC[0])*rdx2SqVol[1]+24.0*rdx2SqVol[0]*rhoLy[0]+4.0*rdx2SqVol[0]*rhoLxLy[0]+16.0*rdx2SqVol[0]*rhoLx[0]+96.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+((21.0*rdx2SqVol1R2+27.0*rdx2SqVol[0]*rdx2SqVol[1]-3.0*rdx2SqVol0R2)*bcVals[11]+(6.0*phiLx[1]-42.0*phiPrevC[0]+21.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*rdx2SqVol1R2+(24.0*rdx2SqVol[0]*phiLy[1]+9.0*rdx2SqVol[0]*phiLx[1]+((-90.0*phiPrevC[0])+3.0*phiLy[0]+9.0*phiLxLy[0]+18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-3.0*rdx2SqVol0R2*phiLy[1]+3.0*rdx2SqVol0R2*phiLx[1]+((-12.0*phiPrevC[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol0R2)*omega+42.0*phiPrevC[0]*rdx2SqVol1R2+90.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+12.0*phiPrevC[0]*rdx2SqVol0R2)/(42.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+12.0*rdx2SqVol0R2); 
  phiC[1] = (((28.0*rdx2SqVol[1]+36.0*rdx2SqVol[0])*rhoC[3]+32.0*rdx2SqVol[0]*rhoC[2]+(28.0*rdx2SqVol[1]+36.0*rdx2SqVol[0])*rhoLy[1]+(4.0*rdx2SqVol[0]-4.0*rdx2SqVol[1])*rhoLx[1]+(112.0*rdx2SqVol[1]+144.0*rdx2SqVol[0])*rhoC[1]+((-4.0*rhoLxLy[0])-16.0*rhoLx[0])*rdx2SqVol[1]+32.0*rdx2SqVol[0]*rhoLy[0]+4.0*rdx2SqVol[0]*rhoLxLy[0]+16.0*rdx2SqVol[0]*rhoLx[0]+128.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+((21.0*rdx2SqVol1R2+54.0*rdx2SqVol[0]*rdx2SqVol[1]-3.0*rdx2SqVol0R2)*bcVals[11]+((-42.0*phiPrevC[1])+21.0*phiLy[1]-3.0*phiLx[1]-3.0*phiLxLy[0]+6.0*phiLx[0])*rdx2SqVol1R2+((-90.0*rdx2SqVol[0]*phiPrevC[1])+12.0*rdx2SqVol[0]*phiLy[1]+(42.0*phiLy[0]-18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol0R2*phiPrevC[1]-9.0*rdx2SqVol0R2*phiLy[1]+3.0*rdx2SqVol0R2*phiLx[1]+(6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol0R2)*omega+42.0*phiPrevC[1]*rdx2SqVol1R2+90.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+12.0*rdx2SqVol0R2*phiPrevC[1])/(42.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+12.0*rdx2SqVol0R2); 
  phiC[2] = (bcVals[11]-1.0*phiPrevC[2])*omega+phiPrevC[2]; 
  phiC[3] = (bcVals[11]-1.0*phiPrevC[3])*omega+phiPrevC[3]; 

}

void MGpoissonFEMDampedJacobi2xSer_UxRobinUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);

  phiC[0] = (((4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+((16.0*dxC[0]*bcVals[3]+4.0*dxC[1])*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]+(16.0*dxC[0]*bcVals3R2+8.0*dxC[1]*bcVals[3])*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals3R2+(32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[3])*omega*volFac+(((30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((108.0*dxC[0]*rdx2SqVol[1]-36.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(96.0*dxC[0]*rdx2SqVol[1]-48.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[11]+((30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[4]+(60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[3]-6.0*dxC[1]*rdx2SqVol[1]+30.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]-15.0*dxC[0]*phiPrevC[0]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-15.0*dxC[0]*phiPrevC[0])-6.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+(((12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]-54.0*dxC[0]*phiPrevC[0]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-90.0*dxC[0]*phiPrevC[0])-24.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+((-15.0*phiPrevC[0])+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-15.0*phiPrevC[0])-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]-48.0*dxC[0]*phiPrevC[0]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-120.0*dxC[0]*phiPrevC[0])-24.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals3R2+((6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+((-24.0*phiPrevC[0])+24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-60.0*phiPrevC[0])-12.0*phiLy[0]+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*omega+(15.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+15.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals4R2+((54.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+90.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals[3]+15.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+15.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+120.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+60.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals[3])/((15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]); 
  phiC[1] = (((4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+((8.0*dxC[0]*bcVals[3]+4.0*dxC[1])*rhoC[3]+(32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4])*omega*volFac+(((30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((48.0*dxC[0]*rdx2SqVol[1]-24.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[4])*bcVals[11]+((24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-15.0*dxC[0]*phiPrevC[1])+3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]-15.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+((((-54.0*dxC[0]*phiPrevC[1])+6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]-90.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]+((-15.0*dxC[1]*phiPrevC[1])+3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]-15.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-48.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals3R2+((-24.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[3])*omega+(15.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals4R2+((54.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals[3]+15.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[4]+(48.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals3R2+(24.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[3])/((15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]); 
  phiC[2] = (((4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+((16.0*dxC[0]*bcVals[3]+4.0*dxC[1])*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]+(16.0*dxC[0]*bcVals3R2+8.0*dxC[1]*bcVals[3])*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals3R2+(32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[3])*omega*volFac+(((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(108.0*dxC[1]*rdx2SqVol[1]+108.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[11]+((30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[4]+(60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[3]-6.0*dxC[1]*rdx2SqVol[1]+30.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+((((-54.0*dxC[0]*rdx2SqVol[1])-90.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]+(12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]+((-15.0*dxC[1]*rdx2SqVol[1])-15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+(((-48.0*dxC[0]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]+(12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals3R2+(((-24.0*dxC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*omega+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals[3]+(15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals[3])/((15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]); 
  phiC[3] = (((4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+((8.0*dxC[0]*bcVals[3]+4.0*dxC[1])*rhoC[3]+(32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4])*omega*volFac+(((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((96.0*dxC[0]*rdx2SqVol[1]+96.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[1]*rdx2SqVol[1]+120.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[11]+((24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[3]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+((((-54.0*dxC[0]*rdx2SqVol[1])-90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]-15.0*dxC[1]*rdx2SqVol[1]-15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]+((6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+(((-48.0*dxC[0]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+((-24.0*dxC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*phiPrevC[3])*omega+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[3]*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]*bcVals[4]+((48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*phiPrevC[3])/((15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]); 

}

void MGpoissonFEMDampedJacobi2xSer_UxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *phiPrevC = phiPrev[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phiPrev[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phiPrev[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phiPrev[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phiPrev[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phiPrev[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);
  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = ((((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]+8.0*dxC[1]*bcVals[3]*rhoC[3]+(32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[3])*bcVals10R2+(((16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[4]+32.0*dxC[1]*bcVals[3]*rhoC[3]+(128.0*dxC[1]*rhoC[2]+32.0*dxC[1]*rhoLy[1]+32.0*dxC[1]*rhoLx[1]+128.0*dxC[1]*rhoC[1]+(128.0*rhoLy[0]+32.0*rhoLxLy[0]+128.0*rhoLx[0]+512.0*rhoC[0])*dxC[1])*bcVals[3])*bcVals[9]+(4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+(16.0*dxC[0]*bcVals[3]*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4]+16.0*dxC[0]*bcVals3R2*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals3R2)*bcVals[10]+((16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[4]+32.0*dxC[1]*bcVals[3]*rhoC[3]+(128.0*dxC[1]*rhoC[2]+32.0*dxC[1]*rhoLy[1]+32.0*dxC[1]*rhoLx[1]+128.0*dxC[1]*rhoC[1]+(128.0*rhoLy[0]+32.0*rhoLxLy[0]+128.0*rhoLx[0]+512.0*rhoC[0])*dxC[1])*bcVals[3])*bcVals9R2+((8.0*dxC[0]*rhoC[3]+32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals4R2+(32.0*dxC[0]*bcVals[3]*rhoC[3]+(128.0*dxC[0]*rhoC[2]+32.0*dxC[0]*rhoLy[1]+32.0*dxC[0]*rhoLx[1]+128.0*dxC[0]*rhoC[1]+128.0*dxC[0]*rhoLy[0]+32.0*dxC[0]*rhoLxLy[0]+128.0*dxC[0]*rhoLx[0]+512.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4]+32.0*dxC[0]*bcVals3R2*rhoC[3]+(128.0*dxC[0]*rhoC[2]+32.0*dxC[0]*rhoLy[1]+32.0*dxC[0]*rhoLx[1]+128.0*dxC[0]*rhoC[1]+128.0*dxC[0]*rhoLy[0]+32.0*dxC[0]*rhoLxLy[0]+128.0*dxC[0]*rhoLx[0]+512.0*dxC[0]*rhoC[0])*bcVals3R2)*bcVals[9])*omega*volFac+((((30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[10]+((60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(120.0*dxC[1]*rdx2SqVol[1]-24.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(108.0*dxC[0]*rdx2SqVol[1]-36.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(96.0*dxC[0]*rdx2SqVol[1]-48.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[11]+((30.0*rdx2SqVol[0]*dxC[1]-6.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+((3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+((-15.0*phiPrevC[0])+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-15.0*phiPrevC[0])-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+((-24.0*phiPrevC[0])+24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-60.0*phiPrevC[0])-12.0*phiLy[0]+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((108.0*rdx2SqVol[0]*dxC[1]-36.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+((12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+((-90.0*phiPrevC[0])+48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-54.0*phiPrevC[0])-24.0*phiLy[0]+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((24.0*dxC[1]*phiLy[1]+24.0*dxC[1]*phiLx[1]+((-144.0*phiPrevC[0])+96.0*phiLy[0]+24.0*phiLxLy[0]-48.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+24.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-216.0*phiPrevC[0])-48.0*phiLy[0]+24.0*phiLxLy[0]+96.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+((30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[4]+(60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[3])*bcVals[5]+((3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]-15.0*dxC[0]*phiPrevC[0]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-15.0*dxC[0]*phiPrevC[0])-6.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+((12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]-54.0*dxC[0]*phiPrevC[0]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-90.0*dxC[0]*phiPrevC[0])-24.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+((12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]-48.0*dxC[0]*phiPrevC[0]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-120.0*dxC[0]*phiPrevC[0])-24.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((96.0*rdx2SqVol[0]*dxC[1]-48.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+((12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+((-120.0*phiPrevC[0])+48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-48.0*phiPrevC[0])-24.0*phiLy[0]+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((24.0*dxC[1]*phiLy[1]+24.0*dxC[1]*phiLx[1]+((-192.0*phiPrevC[0])+96.0*phiLy[0]+24.0*phiLxLy[0]-48.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+24.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-192.0*phiPrevC[0])-48.0*phiLy[0]+24.0*phiLxLy[0]+96.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+(((60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[4]+(120.0*dxC[0]*rdx2SqVol[0]-24.0*dxC[0]*rdx2SqVol[1])*bcVals[3])*bcVals[5]+((6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]-60.0*dxC[0]*phiPrevC[0]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiPrevC[0])-12.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+((24.0*dxC[0]*phiLy[1]+24.0*dxC[0]*phiLx[1]-216.0*dxC[0]*phiPrevC[0]+96.0*dxC[0]*phiLy[0]+24.0*dxC[0]*phiLxLy[0]-48.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+24.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-144.0*dxC[0]*phiPrevC[0])-48.0*dxC[0]*phiLy[0]+24.0*dxC[0]*phiLxLy[0]+96.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+((24.0*dxC[0]*phiLy[1]+24.0*dxC[0]*phiLx[1]-192.0*dxC[0]*phiPrevC[0]+96.0*dxC[0]*phiLy[0]+24.0*dxC[0]*phiLxLy[0]-48.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+24.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-192.0*dxC[0]*phiPrevC[0])-48.0*dxC[0]*phiLy[0]+24.0*dxC[0]*phiLxLy[0]+96.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals3R2)*bcVals[9])*omega+((15.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+15.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals[4]+(24.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+60.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((90.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+54.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals[4]+(144.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+216.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+15.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals4R2+(54.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+90.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+120.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((120.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+48.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*phiPrevC[0]*dxC[1]*rdx2SqVol[1]+192.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+24.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+144.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*phiPrevC[0]*rdx2SqVol[1]+192.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[9])/(((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[9]); 
  phiC[1] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]*bcVals10R2+((16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[4]*bcVals[9]+(4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+(8.0*dxC[0]*bcVals[3]*rhoC[3]+(32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4])*bcVals[10]+(16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[4]*bcVals9R2+((8.0*dxC[0]*rhoC[3]+32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals4R2+(16.0*dxC[0]*bcVals[3]*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4])*bcVals[9])*omega*volFac+(((30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[4]*bcVals[10]+(60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[4]*bcVals[9]+(30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(48.0*dxC[0]*rdx2SqVol[1]-24.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4])*bcVals[11]+((24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-15.0*dxC[1]*phiPrevC[1])+3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]-15.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-24.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[3])*bcVals10R2+(((144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-90.0*dxC[1]*phiPrevC[1])+12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]-54.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-144.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1])-216.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[3])*bcVals[9]+((24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals[3])*bcVals[5]+(((-15.0*dxC[0]*phiPrevC[1])+3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]-15.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+(((-54.0*dxC[0]*phiPrevC[1])+6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]-90.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+((-48.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals3R2)*bcVals[10]+((192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-120.0*dxC[1]*phiPrevC[1])+12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]-48.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-192.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1])-192.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[3])*bcVals9R2+(((108.0*dxC[0]*rdx2SqVol[1]+108.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals[3])*bcVals[5]+(((-60.0*dxC[0]*phiPrevC[1])+6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]-24.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+(((-216.0*dxC[0]*phiPrevC[1])+12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]-144.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+((-192.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1])-192.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals3R2)*bcVals[9])*omega+((15.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[4]+(24.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[4]+(144.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals4R2+(54.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[4]+(192.0*dxC[1]*phiPrevC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1]*phiPrevC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals4R2+(216.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*phiPrevC[1]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0]*phiPrevC[1])*bcVals3R2)*bcVals[9])/(((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[9]); 
  phiC[2] = ((((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]+8.0*dxC[1]*bcVals[3]*rhoC[3]+(32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[3])*bcVals10R2+(((8.0*dxC[1]*rhoC[3]+32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[4]+16.0*dxC[1]*bcVals[3]*rhoC[3]+(64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[3])*bcVals[9]+(4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+(16.0*dxC[0]*bcVals[3]*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4]+16.0*dxC[0]*bcVals3R2*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals3R2)*bcVals[10])*omega*volFac+((((60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(108.0*dxC[1]*rdx2SqVol[1]+108.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[11]+((30.0*rdx2SqVol[0]*dxC[1]-6.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+(((-15.0*dxC[1]*rdx2SqVol[1])-15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+(((-24.0*dxC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((48.0*rdx2SqVol[0]*dxC[1]-24.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+(((-90.0*dxC[1]*rdx2SqVol[1])-54.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+(((-144.0*dxC[1]*rdx2SqVol[1])-216.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]+(12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+((30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[4]+(60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[3])*bcVals[5]+(((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+(((-54.0*dxC[0]*rdx2SqVol[1])-90.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]+(12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(((-48.0*dxC[0]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]+(12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+(((-120.0*dxC[1]*rdx2SqVol[1])-48.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals[4]+((-192.0*dxC[1]*rdx2SqVol[1])-192.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals[3])*bcVals9R2+(((-60.0*dxC[0]*rdx2SqVol[1])-24.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals4R2+((-216.0*dxC[0]*rdx2SqVol[1])-144.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals[3]*bcVals[4]+((-192.0*dxC[0]*rdx2SqVol[1])-192.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals3R2)*bcVals[9])*omega+((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*phiPrevC[2]*bcVals[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*phiPrevC[2]*bcVals3R2)*bcVals[9])/(((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[9]); 
  phiC[3] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]*bcVals10R2+((8.0*dxC[1]*rhoC[3]+32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[4]*bcVals[9]+(4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+(8.0*dxC[0]*bcVals[3]*rhoC[3]+(32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4])*bcVals[10])*omega*volFac+((((60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[1]*rdx2SqVol[1]+120.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(96.0*dxC[0]*rdx2SqVol[1]+96.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4])*bcVals[11]+((24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-15.0*dxC[1]*rdx2SqVol[1])-15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-24.0*dxC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiPrevC[3])*bcVals10R2+(((96.0*dxC[1]*rdx2SqVol[1]+96.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-90.0*dxC[1]*rdx2SqVol[1])-54.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-144.0*dxC[1]*rdx2SqVol[1])-216.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiPrevC[3])*bcVals[9]+((24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals[3])*bcVals[5]+(((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[3]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+(((-54.0*dxC[0]*rdx2SqVol[1])-90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*phiPrevC[3]+((6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3])*bcVals[4]+((-48.0*dxC[0]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2*phiPrevC[3])*bcVals[10]+(((-120.0*dxC[1]*rdx2SqVol[1])-48.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]*bcVals[4]+((-192.0*dxC[1]*rdx2SqVol[1])-192.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiPrevC[3])*bcVals9R2+(((120.0*dxC[0]*rdx2SqVol[1]+48.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals[3])*bcVals[5]+((-60.0*dxC[0]*rdx2SqVol[1])-24.0*dxC[0]*rdx2SqVol[0])*phiPrevC[3]*bcVals4R2+((-216.0*dxC[0]*rdx2SqVol[1])-144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*phiPrevC[3]*bcVals[4]+((-192.0*dxC[0]*rdx2SqVol[1])-192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2*phiPrevC[3])*bcVals[9])*omega+((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiPrevC[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiPrevC[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiPrevC[3]*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*phiPrevC[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2*phiPrevC[3])*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*phiPrevC[3]*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiPrevC[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*phiPrevC[3]*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*phiPrevC[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2*phiPrevC[3])*bcVals[9])/(((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[9]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((16.0*rhoUy[0]+4.0*rhoUxUy[0]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((12.0*phiUy[0]+3.0*phiUxUy[0]+3.0*phiUxLy[0]-6.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+((-6.0*phiUy[0])+3.0*phiUxUy[0]+3.0*phiUxLy[0]+12.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0])*omega+24.0*phiC[0]*rdx2SqVol[1]+24.0*phiC[0]*rdx2SqVol[0])/(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiC[0])*omega+phiC[0]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((-2.0*bcVals[2])+phiUx[0]-1.0*phiC[0])*omega+phiC[0]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = -(1.0*((2.0*bcVals[2]+(phiC[0]-1.0*phiUx[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])*omega-1.0*phiC[0]*bcVals[1]+2.0*bcVals[0]*phiC[0]))/(bcVals[1]-2.0*bcVals[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((15.0*rdx2SqVol[0]-3.0*rdx2SqVol[1])*bcVals[5]+(3.0*phiLy[1]+12.0*phiUy[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+((-6.0*phiUy[0])-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0])*omega+24.0*phiC[0]*rdx2SqVol[1]+24.0*phiC[0]*rdx2SqVol[0])/(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[5]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((30.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*bcVals[5]+(3.0*phiLy[1]+15.0*phiUy[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-30.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+((-3.0*phiUy[0])-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol[0])*omega+30.0*phiC[0]*rdx2SqVol[1]+12.0*phiC[0]*rdx2SqVol[0])/(30.0*rdx2SqVol[1]+12.0*rdx2SqVol[0]); 
  phiC[1] = ((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((54.0*rdx2SqVol[1]+54.0*rdx2SqVol[0])*bcVals[5]+(3.0*phiLy[1]-30.0*phiC[1]+15.0*phiUy[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]-12.0*rdx2SqVol[0]*phiC[1]+((-3.0*phiUy[0])-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+30.0*phiC[1]*rdx2SqVol[1]+12.0*rdx2SqVol[0]*phiC[1])/(30.0*rdx2SqVol[1]+12.0*rdx2SqVol[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);

  phiC[0] = (((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*bcVals[4]+(8.0*rhoUy[1]+8.0*rhoLy[1]+32.0*rhoC[1]+32.0*rhoUy[0]+32.0*rhoLy[0]+8.0*rhoLxUy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*bcVals[3])*omega*volFac+((30.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*bcVals[5]+((3.0*phiLy[1]+15.0*phiUy[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-30.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+((-3.0*phiUy[0])-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol[0])*bcVals[4]+((6.0*phiLy[1]+24.0*phiUy[0]+24.0*phiLy[0]+6.0*phiLxUy[0]+6.0*phiLxLy[0]-12.0*phiLx[0]-48.0*phiC[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiLy[1]+((-12.0*phiUy[0])-12.0*phiLy[0]+6.0*phiLxUy[0]+6.0*phiLxLy[0]+24.0*phiLx[0]-48.0*phiC[0])*rdx2SqVol[0])*bcVals[3])*omega+(30.0*phiC[0]*rdx2SqVol[1]+12.0*phiC[0]*rdx2SqVol[0])*bcVals[4]+(48.0*phiC[0]*rdx2SqVol[1]+48.0*phiC[0]*rdx2SqVol[0])*bcVals[3])/((30.0*rdx2SqVol[1]+12.0*rdx2SqVol[0])*bcVals[4]+(48.0*rdx2SqVol[1]+48.0*rdx2SqVol[0])*bcVals[3]); 
  phiC[1] = (((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*bcVals4R2+(8.0*rhoUy[1]+8.0*rhoLy[1]+32.0*rhoC[1]+32.0*rhoUy[0]+32.0*rhoLy[0]+8.0*rhoLxUy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*bcVals[3]*bcVals[4])*omega*volFac+(((54.0*rdx2SqVol[1]+54.0*rdx2SqVol[0])*bcVals[4]+(96.0*rdx2SqVol[1]+96.0*rdx2SqVol[0])*bcVals[3])*bcVals[5]+((3.0*phiLy[1]-30.0*phiC[1]+15.0*phiUy[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]-12.0*rdx2SqVol[0]*phiC[1]+((-3.0*phiUy[0])-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*bcVals4R2+((6.0*phiLy[1]-108.0*phiC[1]+24.0*phiUy[0]+24.0*phiLy[0]+6.0*phiLxUy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiLy[1]-72.0*rdx2SqVol[0]*phiC[1]+((-12.0*phiUy[0])-12.0*phiLy[0]+6.0*phiLxUy[0]+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+((-96.0*phiC[1]*rdx2SqVol[1])-96.0*rdx2SqVol[0]*phiC[1])*bcVals3R2)*omega+(30.0*phiC[1]*rdx2SqVol[1]+12.0*rdx2SqVol[0]*phiC[1])*bcVals4R2+(108.0*phiC[1]*rdx2SqVol[1]+72.0*rdx2SqVol[0]*phiC[1])*bcVals[3]*bcVals[4]+(96.0*phiC[1]*rdx2SqVol[1]+96.0*rdx2SqVol[0]*phiC[1])*bcVals3R2)/((30.0*rdx2SqVol[1]+12.0*rdx2SqVol[0])*bcVals4R2+(108.0*rdx2SqVol[1]+72.0*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(96.0*rdx2SqVol[1]+96.0*rdx2SqVol[0])*bcVals3R2); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiC[0])*omega+phiC[0]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((-2.0*bcVals[8])+phiUy[0]-1.0*phiC[0])*omega+phiC[0]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = -(1.0*((2.0*bcVals[8]+(phiC[0]-1.0*phiUy[0])*bcVals[7]-2.0*phiC[0]*bcVals[6])*omega-1.0*phiC[0]*bcVals[7]+2.0*phiC[0]*bcVals[6]))/(bcVals[7]-2.0*bcVals[6]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((15.0*rdx2SqVol[1]-3.0*rdx2SqVol[0])*bcVals[11]+(3.0*phiLx[1]+3.0*phiUxLy[0]-6.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]+(3.0*phiUxLy[0]+12.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0])*omega+24.0*phiC[0]*rdx2SqVol[1]+24.0*phiC[0]*rdx2SqVol[0])/(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[11]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((30.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*bcVals[11]+(3.0*phiLx[1]+3.0*phiUxLy[0]-3.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]+(3.0*phiUxLy[0]+15.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-30.0*phiC[0])*rdx2SqVol[0])*omega+12.0*phiC[0]*rdx2SqVol[1]+30.0*phiC[0]*rdx2SqVol[0])/(12.0*rdx2SqVol[1]+30.0*rdx2SqVol[0]); 
  phiC[1] = ((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((54.0*rdx2SqVol[1]+54.0*rdx2SqVol[0])*bcVals[11]+(3.0*phiLx[1]-12.0*phiC[1]+3.0*phiUxLy[0]-3.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]-30.0*rdx2SqVol[0]*phiC[1]+(3.0*phiUxLy[0]+15.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+12.0*phiC[1]*rdx2SqVol[1]+30.0*rdx2SqVol[0]*phiC[1])/(12.0*rdx2SqVol[1]+30.0*rdx2SqVol[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = (((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*bcVals[10]+(8.0*rhoUx[1]+8.0*rhoLx[1]+32.0*rhoC[1]+8.0*rhoUxLy[0]+32.0*rhoUx[0]+32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*bcVals[9])*omega*volFac+((30.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*bcVals[11]+((3.0*phiLx[1]+3.0*phiUxLy[0]-3.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]+(3.0*phiUxLy[0]+15.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-30.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((6.0*phiLx[1]+6.0*phiUxLy[0]-12.0*phiUx[0]+24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0]-48.0*phiC[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiLx[1]+(6.0*phiUxLy[0]+24.0*phiUx[0]-12.0*phiLy[0]+6.0*phiLxLy[0]+24.0*phiLx[0]-48.0*phiC[0])*rdx2SqVol[0])*bcVals[9])*omega+(12.0*phiC[0]*rdx2SqVol[1]+30.0*phiC[0]*rdx2SqVol[0])*bcVals[10]+(48.0*phiC[0]*rdx2SqVol[1]+48.0*phiC[0]*rdx2SqVol[0])*bcVals[9])/((12.0*rdx2SqVol[1]+30.0*rdx2SqVol[0])*bcVals[10]+(48.0*rdx2SqVol[1]+48.0*rdx2SqVol[0])*bcVals[9]); 
  phiC[1] = (((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*bcVals10R2+(8.0*rhoUx[1]+8.0*rhoLx[1]+32.0*rhoC[1]+8.0*rhoUxLy[0]+32.0*rhoUx[0]+32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*bcVals[9]*bcVals[10])*omega*volFac+(((54.0*rdx2SqVol[1]+54.0*rdx2SqVol[0])*bcVals[10]+(96.0*rdx2SqVol[1]+96.0*rdx2SqVol[0])*bcVals[9])*bcVals[11]+((3.0*phiLx[1]-12.0*phiC[1]+3.0*phiUxLy[0]-3.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]-30.0*rdx2SqVol[0]*phiC[1]+(3.0*phiUxLy[0]+15.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*bcVals10R2+((6.0*phiLx[1]-72.0*phiC[1]+6.0*phiUxLy[0]-12.0*phiUx[0]+24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiLx[1]-108.0*rdx2SqVol[0]*phiC[1]+(6.0*phiUxLy[0]+24.0*phiUx[0]-12.0*phiLy[0]+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0])*bcVals[9]*bcVals[10]+((-96.0*phiC[1]*rdx2SqVol[1])-96.0*rdx2SqVol[0]*phiC[1])*bcVals9R2)*omega+(12.0*phiC[1]*rdx2SqVol[1]+30.0*rdx2SqVol[0]*phiC[1])*bcVals10R2+(72.0*phiC[1]*rdx2SqVol[1]+108.0*rdx2SqVol[0]*phiC[1])*bcVals[9]*bcVals[10]+(96.0*phiC[1]*rdx2SqVol[1]+96.0*rdx2SqVol[0]*phiC[1])*bcVals9R2)/((12.0*rdx2SqVol[1]+30.0*rdx2SqVol[0])*bcVals10R2+(72.0*rdx2SqVol[1]+108.0*rdx2SqVol[0])*bcVals[9]*bcVals[10]+(96.0*rdx2SqVol[1]+96.0*rdx2SqVol[0])*bcVals9R2); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((dxC[1]*bcVals[8]+dxC[0]*bcVals[2]-1.0*phiC[0]*dxC[1]-1.0*dxC[0]*phiC[0])*omega+phiC[0]*dxC[1]+dxC[0]*phiC[0])/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiC[0])*omega+phiC[0]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxDirichletLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiC[0])*omega+phiC[0]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiC[0])*omega+phiC[0]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((-2.0*bcVals[8])-2.0*bcVals[2]+phiUxUy[0]-1.0*phiC[0])*omega+phiC[0]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxNeumannLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*(((2.0*dxC[1]*bcVals[7]-4.0*dxC[1]*bcVals[6]+2.0*dxC[0])*bcVals[8]+(2.0*dxC[1]*bcVals[2]+(phiC[0]-1.0*phiUxUy[0])*dxC[1])*bcVals7R2+(((2.0*phiUxUy[0]-4.0*phiC[0])*dxC[1]-4.0*dxC[1]*bcVals[2])*bcVals[6]+2.0*dxC[0]*bcVals[2]-1.0*dxC[0]*phiUxUy[0]+dxC[0]*phiC[0])*bcVals[7]+4.0*phiC[0]*dxC[1]*bcVals6R2+((-4.0*dxC[0]*bcVals[2])-2.0*dxC[0]*phiC[0])*bcVals[6])*omega-1.0*phiC[0]*dxC[1]*bcVals7R2+(4.0*phiC[0]*dxC[1]*bcVals[6]-1.0*dxC[0]*phiC[0])*bcVals[7]-4.0*phiC[0]*dxC[1]*bcVals6R2+2.0*dxC[0]*phiC[0]*bcVals[6]))/(dxC[1]*bcVals7R2+(dxC[0]-4.0*dxC[1]*bcVals[6])*bcVals[7]+4.0*dxC[1]*bcVals6R2-2.0*dxC[0]*bcVals[6]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxRobinLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiC[0])*omega+phiC[0]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxRobinLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);

  phiC[0] = -(1.0*((((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]+2.0*dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[8]+(2.0*dxC[1]+2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+((phiC[0]-1.0*phiUxUy[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])*dxC[1]+(dxC[0]*phiC[0]-1.0*dxC[0]*phiUxUy[0])*bcVals1R2+(2.0*bcVals[0]*dxC[0]*phiUxUy[0]-4.0*bcVals[0]*dxC[0]*phiC[0])*bcVals[1]+4.0*bcVals0R2*dxC[0]*phiC[0])*omega+(2.0*bcVals[0]*phiC[0]-1.0*phiC[0]*bcVals[1])*dxC[1]-1.0*dxC[0]*phiC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*phiC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0]*phiC[0]))/((bcVals[1]-2.0*bcVals[0])*dxC[1]+dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);
  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*((((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]*bcVals[7]+(8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*bcVals[6]+2.0*dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[8]+(2.0*dxC[1]*bcVals[2]+((phiC[0]-1.0*phiUxUy[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])*dxC[1])*bcVals7R2+((((2.0*phiUxUy[0]-4.0*phiC[0])*bcVals[1]+8.0*bcVals[0]*phiC[0])*dxC[1]-4.0*dxC[1]*bcVals[2])*bcVals[6]+(2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+(dxC[0]*phiC[0]-1.0*dxC[0]*phiUxUy[0])*bcVals1R2+(2.0*bcVals[0]*dxC[0]*phiUxUy[0]-4.0*bcVals[0]*dxC[0]*phiC[0])*bcVals[1]+4.0*bcVals0R2*dxC[0]*phiC[0])*bcVals[7]+(4.0*phiC[0]*bcVals[1]-8.0*bcVals[0]*phiC[0])*dxC[1]*bcVals6R2+((8.0*bcVals[0]*dxC[0]-4.0*dxC[0]*bcVals[1])*bcVals[2]-2.0*dxC[0]*phiC[0]*bcVals1R2+8.0*bcVals[0]*dxC[0]*phiC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0]*phiC[0])*bcVals[6])*omega+(2.0*bcVals[0]*phiC[0]-1.0*phiC[0]*bcVals[1])*dxC[1]*bcVals7R2+((4.0*phiC[0]*bcVals[1]-8.0*bcVals[0]*phiC[0])*dxC[1]*bcVals[6]-1.0*dxC[0]*phiC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*phiC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0]*phiC[0])*bcVals[7]+(8.0*bcVals[0]*phiC[0]-4.0*phiC[0]*bcVals[1])*dxC[1]*bcVals6R2+(2.0*dxC[0]*phiC[0]*bcVals1R2-8.0*bcVals[0]*dxC[0]*phiC[0]*bcVals[1]+8.0*bcVals0R2*dxC[0]*phiC[0])*bcVals[6]))/((bcVals[1]-2.0*bcVals[0])*dxC[1]*bcVals7R2+((8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*bcVals[6]+dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0])*bcVals[7]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals6R2+((-2.0*dxC[0]*bcVals1R2)+8.0*bcVals[0]*dxC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0])*bcVals[6]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiC[0])*omega+phiC[0]; 
  phiC[1] = ((dxC[1]*bcVals[11]+dxC[0]*bcVals[2]+((-1.0*dxC[1])-1.0*dxC[0])*phiC[1])*omega+(dxC[1]+dxC[0])*phiC[1])/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiC[0])*omega+phiC[0]; 
  phiC[1] = (bcVals[2]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxDirichletUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[2]-1.0*phiC[0])*omega+phiC[0]; 
  phiC[1] = (bcVals[2]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUx[1]+8.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+8.0*rhoLy[0]+32.0*rhoC[0])*omega*volFac+(9.0*rdx2SqVol[1]*bcVals[11]+(3.0*phiUxLy[0]-6.0*phiUx[0]+6.0*phiLy[0]-12.0*phiC[0])*rdx2SqVol[1]+(3.0*phiUxLy[0]+6.0*phiUx[0]-3.0*phiLy[0]-6.0*phiC[0])*rdx2SqVol[0])*omega+12.0*phiC[0]*rdx2SqVol[1]+6.0*phiC[0]*rdx2SqVol[0])/(12.0*rdx2SqVol[1]+6.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[11]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((-2.0*bcVals[2])+phiUx[0]-1.0*phiC[0])*omega+phiC[0]; 
  phiC[1] = (2.0*bcVals[11]-2.0*bcVals[2]-1.0*phiC[1]+phiUx[0])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxNeumannUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = ((-2.0*bcVals[2])+phiUx[0]-1.0*phiC[0])*omega+phiC[0]; 
  phiC[1] = (((2.0*dxC[1]*bcVals[10]+4.0*dxC[1]*bcVals[9]-2.0*dxC[0])*bcVals[11]+((-2.0*dxC[1]*bcVals[2])-1.0*dxC[1]*phiC[1]+phiUx[0]*dxC[1])*bcVals10R2+(((-4.0*dxC[1]*bcVals[2])-4.0*dxC[1]*phiC[1]+2.0*phiUx[0]*dxC[1])*bcVals[9]+2.0*dxC[0]*bcVals[2]+dxC[0]*phiC[1]-1.0*dxC[0]*phiUx[0])*bcVals[10]-4.0*dxC[1]*phiC[1]*bcVals9R2+(4.0*dxC[0]*bcVals[2]+2.0*dxC[0]*phiC[1])*bcVals[9])*omega+dxC[1]*phiC[1]*bcVals10R2+(4.0*dxC[1]*phiC[1]*bcVals[9]-1.0*dxC[0]*phiC[1])*bcVals[10]+4.0*dxC[1]*phiC[1]*bcVals9R2-2.0*dxC[0]*phiC[1]*bcVals[9])/(dxC[1]*bcVals10R2+(4.0*dxC[1]*bcVals[9]-1.0*dxC[0])*bcVals[10]+4.0*dxC[1]*bcVals9R2-2.0*dxC[0]*bcVals[9]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxRobinUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUx[1]+8.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+8.0*rhoLy[0]+32.0*rhoC[0])*omega*volFac+(9.0*rdx2SqVol[1]*bcVals[11]+(3.0*phiUxLy[0]-6.0*phiUx[0]+6.0*phiLy[0]-12.0*phiC[0])*rdx2SqVol[1]+(3.0*phiUxLy[0]+6.0*phiUx[0]-3.0*phiLy[0]-6.0*phiC[0])*rdx2SqVol[0])*omega+12.0*phiC[0]*rdx2SqVol[1]+6.0*phiC[0]*rdx2SqVol[0])/(12.0*rdx2SqVol[1]+6.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[11]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxRobinUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);

  phiC[0] = -(1.0*((2.0*bcVals[2]+(phiC[0]-1.0*phiUx[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])*omega-1.0*phiC[0]*bcVals[1]+2.0*bcVals[0]*phiC[0]))/(bcVals[1]-2.0*bcVals[0]); 
  phiC[1] = ((((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]-2.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[11]+((-2.0*dxC[1])+2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+((2.0*bcVals[0]-1.0*bcVals[1])*dxC[1]+dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0])*phiC[1]+phiUx[0]*bcVals[1]*dxC[1]-1.0*dxC[0]*phiUx[0]*bcVals1R2+2.0*bcVals[0]*dxC[0]*phiUx[0]*bcVals[1])*omega+((bcVals[1]-2.0*bcVals[0])*dxC[1]-1.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0])*phiC[1])/((bcVals[1]-2.0*bcVals[0])*dxC[1]-1.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_LxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);
  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = -(1.0*((2.0*bcVals[2]+(phiC[0]-1.0*phiUx[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])*omega-1.0*phiC[0]*bcVals[1]+2.0*bcVals[0]*phiC[0]))/(bcVals[1]-2.0*bcVals[0]); 
  phiC[1] = ((((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]*bcVals[10]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals[9]-2.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[11]+((-2.0*dxC[1]*bcVals[2])+(2.0*bcVals[0]-1.0*bcVals[1])*dxC[1]*phiC[1]+phiUx[0]*bcVals[1]*dxC[1])*bcVals10R2+(((-4.0*dxC[1]*bcVals[2])+(8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*phiC[1]+2.0*phiUx[0]*bcVals[1]*dxC[1])*bcVals[9]+(2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+(dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0])*phiC[1]-1.0*dxC[0]*phiUx[0]*bcVals1R2+2.0*bcVals[0]*dxC[0]*phiUx[0]*bcVals[1])*bcVals[10]+(8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*phiC[1]*bcVals9R2+((4.0*dxC[0]*bcVals[1]-8.0*bcVals[0]*dxC[0])*bcVals[2]+(2.0*dxC[0]*bcVals1R2-8.0*bcVals[0]*dxC[0]*bcVals[1]+8.0*bcVals0R2*dxC[0])*phiC[1])*bcVals[9])*omega+(bcVals[1]-2.0*bcVals[0])*dxC[1]*phiC[1]*bcVals10R2+((4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*phiC[1]*bcVals[9]+((-1.0*dxC[0]*bcVals1R2)+4.0*bcVals[0]*dxC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0])*phiC[1])*bcVals[10]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*phiC[1]*bcVals9R2+((-2.0*dxC[0]*bcVals1R2)+8.0*bcVals[0]*dxC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0])*phiC[1]*bcVals[9])/((bcVals[1]-2.0*bcVals[0])*dxC[1]*bcVals10R2+((4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals[9]-1.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0])*bcVals[10]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals9R2+((-2.0*dxC[0]*bcVals1R2)+8.0*bcVals[0]*dxC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0])*bcVals[9]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxDirichletLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiC[0])*omega+phiC[0]; 
  phiC[1] = ((dxC[1]*bcVals[8]+dxC[0]*bcVals[5]+((-1.0*dxC[1])-1.0*dxC[0])*phiC[1])*omega+(dxC[1]+dxC[0])*phiC[1])/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxDirichletLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUy[1]+8.0*rhoC[1]+16.0*rhoUy[0]+4.0*rhoLxUy[0]+8.0*rhoLx[0]+32.0*rhoC[0])*omega*volFac+(9.0*rdx2SqVol[0]*bcVals[5]+(6.0*phiUy[0]+3.0*phiLxUy[0]-3.0*phiLx[0]-6.0*phiC[0])*rdx2SqVol[1]+((-6.0*phiUy[0])+3.0*phiLxUy[0]+6.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol[0])*omega+6.0*phiC[0]*rdx2SqVol[1]+12.0*phiC[0]*rdx2SqVol[0])/(6.0*rdx2SqVol[1]+12.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[5]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxDirichletLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoUy[1]+8.0*rhoC[1]+16.0*rhoUy[0]+4.0*rhoLxUy[0]+8.0*rhoLx[0]+32.0*rhoC[0])*omega*volFac+(9.0*rdx2SqVol[0]*bcVals[5]+(6.0*phiUy[0]+3.0*phiLxUy[0]-3.0*phiLx[0]-6.0*phiC[0])*rdx2SqVol[1]+((-6.0*phiUy[0])+3.0*phiLxUy[0]+6.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol[0])*omega+6.0*phiC[0]*rdx2SqVol[1]+12.0*phiC[0]*rdx2SqVol[0])/(6.0*rdx2SqVol[1]+12.0*rdx2SqVol[0]); 
  phiC[1] = (bcVals[5]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxNeumannLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiC[0])*omega+phiC[0]; 
  phiC[1] = (bcVals[8]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxNeumannLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((-2.0*bcVals[8])+phiUy[0]-1.0*phiC[0])*omega+phiC[0]; 
  phiC[1] = ((-2.0*bcVals[8])+2.0*bcVals[5]-1.0*phiC[1]+phiUy[0])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxNeumannLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*((2.0*bcVals[8]+(phiC[0]-1.0*phiUy[0])*bcVals[7]-2.0*phiC[0]*bcVals[6])*omega-1.0*phiC[0]*bcVals[7]+2.0*phiC[0]*bcVals[6]))/(bcVals[7]-2.0*bcVals[6]); 
  phiC[1] = -(1.0*(((2.0*dxC[1]*bcVals[7]-4.0*dxC[1]*bcVals[6]-2.0*dxC[0])*bcVals[8]+((-2.0*dxC[1]*bcVals[5])+dxC[1]*phiC[1]-1.0*phiUy[0]*dxC[1])*bcVals7R2+((4.0*dxC[1]*bcVals[5]-4.0*dxC[1]*phiC[1]+2.0*phiUy[0]*dxC[1])*bcVals[6]+2.0*dxC[0]*bcVals[5]-1.0*dxC[0]*phiC[1]+dxC[0]*phiUy[0])*bcVals[7]+4.0*dxC[1]*phiC[1]*bcVals6R2+(2.0*dxC[0]*phiC[1]-4.0*dxC[0]*bcVals[5])*bcVals[6])*omega-1.0*dxC[1]*phiC[1]*bcVals7R2+(4.0*dxC[1]*phiC[1]*bcVals[6]+dxC[0]*phiC[1])*bcVals[7]-4.0*dxC[1]*phiC[1]*bcVals6R2-2.0*dxC[0]*phiC[1]*bcVals[6]))/(dxC[1]*bcVals7R2+((-4.0*dxC[1]*bcVals[6])-1.0*dxC[0])*bcVals[7]+4.0*dxC[1]*bcVals6R2+2.0*dxC[0]*bcVals[6]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxRobinLyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (bcVals[8]-1.0*phiC[0])*omega+phiC[0]; 
  phiC[1] = (bcVals[8]-1.0*phiC[1])*omega+phiC[1]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxRobinLyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);

  phiC[0] = ((-2.0*bcVals[8])+phiUy[0]-1.0*phiC[0])*omega+phiC[0]; 
  phiC[1] = -(1.0*(((2.0*dxC[0]*bcVals4R2+(4.0*dxC[0]*bcVals[3]-2.0*dxC[1])*bcVals[4]-4.0*dxC[1]*bcVals[3])*bcVals[8]+((-2.0*dxC[0]*bcVals[4])-4.0*dxC[0]*bcVals[3]+2.0*dxC[1])*bcVals[5]+(dxC[0]*phiC[1]-1.0*dxC[0]*phiUy[0])*bcVals4R2+((4.0*dxC[0]*phiC[1]-2.0*dxC[0]*phiUy[0])*bcVals[3]-1.0*dxC[1]*phiC[1]+phiUy[0]*dxC[1])*bcVals[4]+4.0*dxC[0]*phiC[1]*bcVals3R2-2.0*dxC[1]*phiC[1]*bcVals[3])*omega-1.0*dxC[0]*phiC[1]*bcVals4R2+(dxC[1]*phiC[1]-4.0*dxC[0]*phiC[1]*bcVals[3])*bcVals[4]-4.0*dxC[0]*phiC[1]*bcVals3R2+2.0*dxC[1]*phiC[1]*bcVals[3]))/(dxC[0]*bcVals4R2+(4.0*dxC[0]*bcVals[3]-1.0*dxC[1])*bcVals[4]+4.0*dxC[0]*bcVals3R2-2.0*dxC[1]*bcVals[3]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxRobinLyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);
  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*((2.0*bcVals[8]+(phiC[0]-1.0*phiUy[0])*bcVals[7]-2.0*phiC[0]*bcVals[6])*omega-1.0*phiC[0]*bcVals[7]+2.0*phiC[0]*bcVals[6]))/(bcVals[7]-2.0*bcVals[6]); 
  phiC[1] = -(1.0*((((2.0*dxC[1]*bcVals[4]+4.0*dxC[1]*bcVals[3])*bcVals[7]+((-4.0*dxC[1]*bcVals[4])-8.0*dxC[1]*bcVals[3])*bcVals[6]-2.0*dxC[0]*bcVals4R2-4.0*dxC[0]*bcVals[3]*bcVals[4])*bcVals[8]+((-2.0*dxC[1]*bcVals[5])+(dxC[1]*phiC[1]-1.0*phiUy[0]*dxC[1])*bcVals[4]+2.0*dxC[1]*phiC[1]*bcVals[3])*bcVals7R2+((4.0*dxC[1]*bcVals[5]+(2.0*phiUy[0]*dxC[1]-4.0*dxC[1]*phiC[1])*bcVals[4]-8.0*dxC[1]*phiC[1]*bcVals[3])*bcVals[6]+(2.0*dxC[0]*bcVals[4]+4.0*dxC[0]*bcVals[3])*bcVals[5]+(dxC[0]*phiUy[0]-1.0*dxC[0]*phiC[1])*bcVals4R2+(2.0*dxC[0]*phiUy[0]-4.0*dxC[0]*phiC[1])*bcVals[3]*bcVals[4]-4.0*dxC[0]*phiC[1]*bcVals3R2)*bcVals[7]+(4.0*dxC[1]*phiC[1]*bcVals[4]+8.0*dxC[1]*phiC[1]*bcVals[3])*bcVals6R2+(((-4.0*dxC[0]*bcVals[4])-8.0*dxC[0]*bcVals[3])*bcVals[5]+2.0*dxC[0]*phiC[1]*bcVals4R2+8.0*dxC[0]*phiC[1]*bcVals[3]*bcVals[4]+8.0*dxC[0]*phiC[1]*bcVals3R2)*bcVals[6])*omega+((-1.0*dxC[1]*phiC[1]*bcVals[4])-2.0*dxC[1]*phiC[1]*bcVals[3])*bcVals7R2+((4.0*dxC[1]*phiC[1]*bcVals[4]+8.0*dxC[1]*phiC[1]*bcVals[3])*bcVals[6]+dxC[0]*phiC[1]*bcVals4R2+4.0*dxC[0]*phiC[1]*bcVals[3]*bcVals[4]+4.0*dxC[0]*phiC[1]*bcVals3R2)*bcVals[7]+((-4.0*dxC[1]*phiC[1]*bcVals[4])-8.0*dxC[1]*phiC[1]*bcVals[3])*bcVals6R2+((-2.0*dxC[0]*phiC[1]*bcVals4R2)-8.0*dxC[0]*phiC[1]*bcVals[3]*bcVals[4]-8.0*dxC[0]*phiC[1]*bcVals3R2)*bcVals[6]))/((dxC[1]*bcVals[4]+2.0*dxC[1]*bcVals[3])*bcVals7R2+(((-4.0*dxC[1]*bcVals[4])-8.0*dxC[1]*bcVals[3])*bcVals[6]-1.0*dxC[0]*bcVals4R2-4.0*dxC[0]*bcVals[3]*bcVals[4]-4.0*dxC[0]*bcVals3R2)*bcVals[7]+(4.0*dxC[1]*bcVals[4]+8.0*dxC[1]*bcVals[3])*bcVals6R2+(2.0*dxC[0]*bcVals4R2+8.0*dxC[0]*bcVals[3]*bcVals[4]+8.0*dxC[0]*bcVals3R2)*bcVals[6]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxDirichletUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = (((4.0*dxC[1]+4.0*dxC[0])*rhoC[3]+(16.0*dxC[1]+16.0*dxC[0])*rhoC[2]+(4.0*dxC[1]+4.0*dxC[0])*rhoLy[1]+(4.0*dxC[1]+4.0*dxC[0])*rhoLx[1]+(16.0*dxC[1]+16.0*dxC[0])*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*omega*volFac+(((15.0*dxC[1]+12.0*dxC[0])*rdx2SqVol[1]-3.0*rdx2SqVol[0]*dxC[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+(((-6.0*dxC[1])-3.0*dxC[0])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+((3.0*dxC[1]+3.0*dxC[0])*phiLy[1]+(3.0*dxC[1]+3.0*dxC[0])*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*dxC[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0]-24.0*dxC[0]*phiC[0])*rdx2SqVol[1]+(3.0*rdx2SqVol[0]*dxC[1]+3.0*dxC[0]*rdx2SqVol[0])*phiLy[1]+(3.0*rdx2SqVol[0]*dxC[1]+3.0*dxC[0]*rdx2SqVol[0])*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]*dxC[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0]-24.0*dxC[0]*phiC[0])*rdx2SqVol[0])*omega+(24.0*phiC[0]*dxC[1]+24.0*dxC[0]*phiC[0])*rdx2SqVol[1]+24.0*phiC[0]*rdx2SqVol[0]*dxC[1]+24.0*dxC[0]*phiC[0]*rdx2SqVol[0])/((24.0*dxC[1]+24.0*dxC[0])*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1]+24.0*dxC[0]*rdx2SqVol[0]); 
  phiC[1] = (bcVals[5]-1.0*phiC[1])*omega+phiC[1]; 
  phiC[2] = (bcVals[11]-1.0*phiC[2])*omega+phiC[2]; 
  phiC[3] = ((dxC[1]*bcVals[11]+dxC[0]*bcVals[5]+((-1.0*dxC[1])-1.0*dxC[0])*phiC[3])*omega+(dxC[1]+dxC[0])*phiC[3])/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxDirichletUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = ((20.0*rdx2SqVol[1]*rhoC[3]+80.0*rdx2SqVol[1]*rhoC[2]+(4.0*rdx2SqVol[1]+8.0*rdx2SqVol[0])*rhoLy[1]+20.0*rdx2SqVol[1]*rhoLx[1]+(24.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+24.0*rhoLx[0]+96.0*rhoC[0])*rdx2SqVol[1]+32.0*rdx2SqVol[0]*rhoLy[0]+8.0*rdx2SqVol[0]*rhoLxLy[0]+28.0*rdx2SqVol[0]*rhoLx[0]+112.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+(((-3.0*rdx2SqVol1R2)+27.0*rdx2SqVol[0]*rdx2SqVol[1]+21.0*rdx2SqVol0R2)*bcVals[5]+(3.0*phiLy[1]-3.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-12.0*phiC[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiLy[1]+24.0*rdx2SqVol[0]*phiLx[1]+(18.0*phiLy[0]+9.0*phiLxLy[0]+3.0*phiLx[0]-90.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+6.0*rdx2SqVol0R2*phiLy[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+21.0*phiLx[0]-42.0*phiC[0])*rdx2SqVol0R2)*omega+12.0*phiC[0]*rdx2SqVol1R2+90.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+42.0*phiC[0]*rdx2SqVol0R2)/(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2); 
  phiC[1] = (bcVals[5]-1.0*phiC[1])*omega+phiC[1]; 
  phiC[2] = (((36.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoC[3]+(144.0*rdx2SqVol[1]+112.0*rdx2SqVol[0])*rhoC[2]+(4.0*rdx2SqVol[1]-4.0*rdx2SqVol[0])*rhoLy[1]+(36.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoLx[1]+32.0*rdx2SqVol[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*rdx2SqVol[1]-16.0*rdx2SqVol[0]*rhoLy[0]-4.0*rdx2SqVol[0]*rhoLxLy[0])*omega*volFac+(((-3.0*rdx2SqVol1R2)+54.0*rdx2SqVol[0]*rdx2SqVol[1]+21.0*rdx2SqVol0R2)*bcVals[5]+((-12.0*rdx2SqVol1R2)-90.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*phiC[2]+(3.0*phiLy[1]-9.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]+6.0*phiLx[0])*rdx2SqVol1R2+(12.0*rdx2SqVol[0]*phiLx[1]+(42.0*phiLx[0]-18.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1]-3.0*rdx2SqVol0R2*phiLy[1]+21.0*rdx2SqVol0R2*phiLx[1]+(6.0*phiLy[0]-3.0*phiLxLy[0])*rdx2SqVol0R2)*omega+(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2)*phiC[2])/(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2); 
  phiC[3] = (bcVals[5]-1.0*phiC[3])*omega+phiC[3]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxDirichletUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = ((20.0*rdx2SqVol[1]*rhoC[3]+80.0*rdx2SqVol[1]*rhoC[2]+(4.0*rdx2SqVol[1]+8.0*rdx2SqVol[0])*rhoLy[1]+20.0*rdx2SqVol[1]*rhoLx[1]+(24.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+24.0*rhoLx[0]+96.0*rhoC[0])*rdx2SqVol[1]+32.0*rdx2SqVol[0]*rhoLy[0]+8.0*rdx2SqVol[0]*rhoLxLy[0]+28.0*rdx2SqVol[0]*rhoLx[0]+112.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+(((-3.0*rdx2SqVol1R2)+27.0*rdx2SqVol[0]*rdx2SqVol[1]+21.0*rdx2SqVol0R2)*bcVals[5]+(3.0*phiLy[1]-3.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-12.0*phiC[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiLy[1]+24.0*rdx2SqVol[0]*phiLx[1]+(18.0*phiLy[0]+9.0*phiLxLy[0]+3.0*phiLx[0]-90.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+6.0*rdx2SqVol0R2*phiLy[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+21.0*phiLx[0]-42.0*phiC[0])*rdx2SqVol0R2)*omega+12.0*phiC[0]*rdx2SqVol1R2+90.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+42.0*phiC[0]*rdx2SqVol0R2)/(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2); 
  phiC[1] = (bcVals[5]-1.0*phiC[1])*omega+phiC[1]; 
  phiC[2] = (((36.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoC[3]+(144.0*rdx2SqVol[1]+112.0*rdx2SqVol[0])*rhoC[2]+(4.0*rdx2SqVol[1]-4.0*rdx2SqVol[0])*rhoLy[1]+(36.0*rdx2SqVol[1]+28.0*rdx2SqVol[0])*rhoLx[1]+32.0*rdx2SqVol[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*rdx2SqVol[1]-16.0*rdx2SqVol[0]*rhoLy[0]-4.0*rdx2SqVol[0]*rhoLxLy[0])*omega*volFac+(((-3.0*rdx2SqVol1R2)+54.0*rdx2SqVol[0]*rdx2SqVol[1]+21.0*rdx2SqVol0R2)*bcVals[5]+((-12.0*rdx2SqVol1R2)-90.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*phiC[2]+(3.0*phiLy[1]-9.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]+6.0*phiLx[0])*rdx2SqVol1R2+(12.0*rdx2SqVol[0]*phiLx[1]+(42.0*phiLx[0]-18.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1]-3.0*rdx2SqVol0R2*phiLy[1]+21.0*rdx2SqVol0R2*phiLx[1]+(6.0*phiLy[0]-3.0*phiLxLy[0])*rdx2SqVol0R2)*omega+(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2)*phiC[2])/(12.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+42.0*rdx2SqVol0R2); 
  phiC[3] = (bcVals[5]-1.0*phiC[3])*omega+phiC[3]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxNeumannUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = ((20.0*rdx2SqVol[0]*rhoC[3]+(28.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[2]+20.0*rdx2SqVol[0]*rhoLy[1]+(8.0*rdx2SqVol[1]+4.0*rdx2SqVol[0])*rhoLx[1]+80.0*rdx2SqVol[0]*rhoC[1]+(28.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+112.0*rhoC[0])*rdx2SqVol[1]+24.0*rdx2SqVol[0]*rhoLy[0]+4.0*rdx2SqVol[0]*rhoLxLy[0]+16.0*rdx2SqVol[0]*rhoLx[0]+96.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+((21.0*rdx2SqVol1R2+27.0*rdx2SqVol[0]*rdx2SqVol[1]-3.0*rdx2SqVol0R2)*bcVals[11]+(6.0*phiLx[1]+21.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0]-42.0*phiC[0])*rdx2SqVol1R2+(24.0*rdx2SqVol[0]*phiLy[1]+9.0*rdx2SqVol[0]*phiLx[1]+(3.0*phiLy[0]+9.0*phiLxLy[0]+18.0*phiLx[0]-90.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-3.0*rdx2SqVol0R2*phiLy[1]+3.0*rdx2SqVol0R2*phiLx[1]+(3.0*phiLxLy[0]+12.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol0R2)*omega+42.0*phiC[0]*rdx2SqVol1R2+90.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+12.0*phiC[0]*rdx2SqVol0R2)/(42.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+12.0*rdx2SqVol0R2); 
  phiC[1] = (((28.0*rdx2SqVol[1]+36.0*rdx2SqVol[0])*rhoC[3]+32.0*rdx2SqVol[0]*rhoC[2]+(28.0*rdx2SqVol[1]+36.0*rdx2SqVol[0])*rhoLy[1]+(4.0*rdx2SqVol[0]-4.0*rdx2SqVol[1])*rhoLx[1]+(112.0*rdx2SqVol[1]+144.0*rdx2SqVol[0])*rhoC[1]+((-4.0*rhoLxLy[0])-16.0*rhoLx[0])*rdx2SqVol[1]+32.0*rdx2SqVol[0]*rhoLy[0]+4.0*rdx2SqVol[0]*rhoLxLy[0]+16.0*rdx2SqVol[0]*rhoLx[0]+128.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+((21.0*rdx2SqVol1R2+54.0*rdx2SqVol[0]*rdx2SqVol[1]-3.0*rdx2SqVol0R2)*bcVals[11]+(21.0*phiLy[1]-3.0*phiLx[1]-42.0*phiC[1]-3.0*phiLxLy[0]+6.0*phiLx[0])*rdx2SqVol1R2+(12.0*rdx2SqVol[0]*phiLy[1]-90.0*rdx2SqVol[0]*phiC[1]+(42.0*phiLy[0]-18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-9.0*rdx2SqVol0R2*phiLy[1]+3.0*rdx2SqVol0R2*phiLx[1]-12.0*rdx2SqVol0R2*phiC[1]+(6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol0R2)*omega+42.0*phiC[1]*rdx2SqVol1R2+90.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+12.0*rdx2SqVol0R2*phiC[1])/(42.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+12.0*rdx2SqVol0R2); 
  phiC[2] = (bcVals[11]-1.0*phiC[2])*omega+phiC[2]; 
  phiC[3] = (bcVals[11]-1.0*phiC[3])*omega+phiC[3]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxNeumannUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  phiC[0] = ((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((30.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*bcVals[11]+(30.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*bcVals[5]+(3.0*phiLy[1]+3.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-15.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-15.0*phiC[0])*rdx2SqVol[0])*omega+15.0*phiC[0]*rdx2SqVol[1]+15.0*phiC[0]*rdx2SqVol[0])/(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0]); 
  phiC[1] = ((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((30.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*bcVals[11]+(24.0*rdx2SqVol[1]+60.0*rdx2SqVol[0])*bcVals[5]+(3.0*phiLy[1]+3.0*phiLx[1]-15.0*phiC[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]-15.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+15.0*phiC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*phiC[1])/(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0]); 
  phiC[2] = ((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((60.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*bcVals[11]+(30.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*bcVals[5]+((-15.0*rdx2SqVol[1])-15.0*rdx2SqVol[0])*phiC[2]+(3.0*phiLy[1]+3.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0])*phiC[2])/(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0]); 
  phiC[3] = ((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*omega*volFac+((60.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*bcVals[11]+(24.0*rdx2SqVol[1]+60.0*rdx2SqVol[0])*bcVals[5]+((-15.0*rdx2SqVol[1])-15.0*rdx2SqVol[0])*phiC[3]+(3.0*phiLy[1]+3.0*phiLx[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0])*omega+(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0])*phiC[3])/(15.0*rdx2SqVol[1]+15.0*rdx2SqVol[0]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxNeumannUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals10R2+((16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[9]+4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals[10]+(16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals9R2+(8.0*dxC[0]*rhoC[3]+32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[9])*omega*volFac+(((30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[10]+(60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+((30.0*rdx2SqVol[0]*dxC[1]-6.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-15.0*phiC[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-15.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals10R2+(((108.0*rdx2SqVol[0]*dxC[1]-36.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+(12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0]-90.0*phiC[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0]-54.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals[9]+(30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[5]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0]-15.0*dxC[0]*phiC[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0]-15.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals[10]+((96.0*rdx2SqVol[0]*dxC[1]-48.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+(12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0]-120.0*phiC[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0]-48.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[5]+(6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0]-60.0*dxC[0]*phiC[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0]-24.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals[9])*omega+(15.0*phiC[0]*dxC[1]*rdx2SqVol[1]+15.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals10R2+((90.0*phiC[0]*dxC[1]*rdx2SqVol[1]+54.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals[9]+15.0*dxC[0]*phiC[0]*rdx2SqVol[1]+15.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals[10]+(120.0*phiC[0]*dxC[1]*rdx2SqVol[1]+48.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals9R2+(60.0*dxC[0]*phiC[0]*rdx2SqVol[1]+24.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals[9])/((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[9]); 
  phiC[1] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals10R2+((16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[9]+4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals[10]+(16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals9R2+(8.0*dxC[0]*rhoC[3]+32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[9])*omega*volFac+(((30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[10]+(60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+((24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]-15.0*dxC[1]*phiC[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]-15.0*rdx2SqVol[0]*dxC[1]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals10R2+(((144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]-90.0*dxC[1]*phiC[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]-54.0*rdx2SqVol[0]*dxC[1]*phiC[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[9]+(24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]-15.0*dxC[0]*phiC[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]-15.0*dxC[0]*rdx2SqVol[0]*phiC[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[10]+((192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]-120.0*dxC[1]*phiC[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]-48.0*rdx2SqVol[0]*dxC[1]*phiC[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals9R2+((108.0*dxC[0]*rdx2SqVol[1]+108.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+(6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]-60.0*dxC[0]*phiC[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]-24.0*dxC[0]*rdx2SqVol[0]*phiC[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[9])*omega+(15.0*dxC[1]*phiC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals10R2+((90.0*dxC[1]*phiC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[9]+15.0*dxC[0]*phiC[1]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals[10]+(120.0*dxC[1]*phiC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals9R2+(60.0*dxC[0]*phiC[1]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals[9])/((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[9]); 
  phiC[2] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals10R2+((8.0*dxC[1]*rhoC[3]+32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[9]+4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals[10])*omega*volFac+(((60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+((30.0*rdx2SqVol[0]*dxC[1]-6.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+((-15.0*dxC[1]*rdx2SqVol[1])-15.0*rdx2SqVol[0]*dxC[1])*phiC[2]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals10R2+(((48.0*rdx2SqVol[0]*dxC[1]-24.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+((-90.0*dxC[1]*rdx2SqVol[1])-54.0*rdx2SqVol[0]*dxC[1])*phiC[2]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[9]+(30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[5]+((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiC[2]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[10]+((-120.0*dxC[1]*rdx2SqVol[1])-48.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals9R2+((-60.0*dxC[0]*rdx2SqVol[1])-24.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals[9])*omega+(15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiC[2])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals[9])/((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[9]); 
  phiC[3] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals10R2+((8.0*dxC[1]*rhoC[3]+32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[9]+4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals[10])*omega*volFac+(((60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+((24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((-15.0*dxC[1]*rdx2SqVol[1])-15.0*rdx2SqVol[0]*dxC[1])*phiC[3]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals10R2+(((96.0*dxC[1]*rdx2SqVol[1]+96.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((-90.0*dxC[1]*rdx2SqVol[1])-54.0*rdx2SqVol[0]*dxC[1])*phiC[3]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[9]+(24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiC[3]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[10]+((-120.0*dxC[1]*rdx2SqVol[1])-48.0*rdx2SqVol[0]*dxC[1])*phiC[3]*bcVals9R2+((120.0*dxC[0]*rdx2SqVol[1]+48.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+((-60.0*dxC[0]*rdx2SqVol[1])-24.0*dxC[0]*rdx2SqVol[0])*phiC[3])*bcVals[9])*omega+(15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiC[3]*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*phiC[3]*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiC[3])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*phiC[3]*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*phiC[3]*bcVals[9])/((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals10R2+((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[9]+15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals[10]+(120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals9R2+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals[9]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxRobinUyDirichlet_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = ((20.0*rdx2SqVol[0]*rhoC[3]+(28.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[2]+20.0*rdx2SqVol[0]*rhoLy[1]+(8.0*rdx2SqVol[1]+4.0*rdx2SqVol[0])*rhoLx[1]+80.0*rdx2SqVol[0]*rhoC[1]+(28.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+112.0*rhoC[0])*rdx2SqVol[1]+24.0*rdx2SqVol[0]*rhoLy[0]+4.0*rdx2SqVol[0]*rhoLxLy[0]+16.0*rdx2SqVol[0]*rhoLx[0]+96.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+((21.0*rdx2SqVol1R2+27.0*rdx2SqVol[0]*rdx2SqVol[1]-3.0*rdx2SqVol0R2)*bcVals[11]+(6.0*phiLx[1]+21.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0]-42.0*phiC[0])*rdx2SqVol1R2+(24.0*rdx2SqVol[0]*phiLy[1]+9.0*rdx2SqVol[0]*phiLx[1]+(3.0*phiLy[0]+9.0*phiLxLy[0]+18.0*phiLx[0]-90.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-3.0*rdx2SqVol0R2*phiLy[1]+3.0*rdx2SqVol0R2*phiLx[1]+(3.0*phiLxLy[0]+12.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol0R2)*omega+42.0*phiC[0]*rdx2SqVol1R2+90.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+12.0*phiC[0]*rdx2SqVol0R2)/(42.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+12.0*rdx2SqVol0R2); 
  phiC[1] = (((28.0*rdx2SqVol[1]+36.0*rdx2SqVol[0])*rhoC[3]+32.0*rdx2SqVol[0]*rhoC[2]+(28.0*rdx2SqVol[1]+36.0*rdx2SqVol[0])*rhoLy[1]+(4.0*rdx2SqVol[0]-4.0*rdx2SqVol[1])*rhoLx[1]+(112.0*rdx2SqVol[1]+144.0*rdx2SqVol[0])*rhoC[1]+((-4.0*rhoLxLy[0])-16.0*rhoLx[0])*rdx2SqVol[1]+32.0*rdx2SqVol[0]*rhoLy[0]+4.0*rdx2SqVol[0]*rhoLxLy[0]+16.0*rdx2SqVol[0]*rhoLx[0]+128.0*rdx2SqVol[0]*rhoC[0])*omega*volFac+((21.0*rdx2SqVol1R2+54.0*rdx2SqVol[0]*rdx2SqVol[1]-3.0*rdx2SqVol0R2)*bcVals[11]+(21.0*phiLy[1]-3.0*phiLx[1]-42.0*phiC[1]-3.0*phiLxLy[0]+6.0*phiLx[0])*rdx2SqVol1R2+(12.0*rdx2SqVol[0]*phiLy[1]-90.0*rdx2SqVol[0]*phiC[1]+(42.0*phiLy[0]-18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-9.0*rdx2SqVol0R2*phiLy[1]+3.0*rdx2SqVol0R2*phiLx[1]-12.0*rdx2SqVol0R2*phiC[1]+(6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol0R2)*omega+42.0*phiC[1]*rdx2SqVol1R2+90.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+12.0*rdx2SqVol0R2*phiC[1])/(42.0*rdx2SqVol1R2+90.0*rdx2SqVol[0]*rdx2SqVol[1]+12.0*rdx2SqVol0R2); 
  phiC[2] = (bcVals[11]-1.0*phiC[2])*omega+phiC[2]; 
  phiC[3] = (bcVals[11]-1.0*phiC[3])*omega+phiC[3]; 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxRobinUyNeumann_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);

  phiC[0] = (((4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+((16.0*dxC[0]*bcVals[3]+4.0*dxC[1])*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]+(16.0*dxC[0]*bcVals3R2+8.0*dxC[1]*bcVals[3])*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals3R2+(32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[3])*omega*volFac+(((30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((108.0*dxC[0]*rdx2SqVol[1]-36.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(96.0*dxC[0]*rdx2SqVol[1]-48.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[11]+((30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[4]+(60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[3]-6.0*dxC[1]*rdx2SqVol[1]+30.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0]-15.0*dxC[0]*phiC[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0]-15.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals4R2+(((12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0]-54.0*dxC[0]*phiC[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0]-90.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals[3]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-15.0*phiC[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-15.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0]-48.0*dxC[0]*phiC[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0]-120.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals3R2+((6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0]-24.0*phiC[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0]-60.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*omega+(15.0*dxC[0]*phiC[0]*rdx2SqVol[1]+15.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals4R2+((54.0*dxC[0]*phiC[0]*rdx2SqVol[1]+90.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals[3]+15.0*phiC[0]*dxC[1]*rdx2SqVol[1]+15.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[0]*phiC[0]*rdx2SqVol[1]+120.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*phiC[0]*dxC[1]*rdx2SqVol[1]+60.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals[3])/((15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]); 
  phiC[1] = (((4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+((8.0*dxC[0]*bcVals[3]+4.0*dxC[1])*rhoC[3]+(32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4])*omega*volFac+(((30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((48.0*dxC[0]*rdx2SqVol[1]-24.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[4])*bcVals[11]+((24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]-15.0*dxC[0]*phiC[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]-15.0*dxC[0]*rdx2SqVol[0]*phiC[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+(((6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]-54.0*dxC[0]*phiC[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]-90.0*dxC[0]*rdx2SqVol[0]*phiC[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]-15.0*dxC[1]*phiC[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]-15.0*rdx2SqVol[0]*dxC[1]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-48.0*dxC[0]*phiC[1]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals3R2+((-24.0*dxC[1]*phiC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[3])*omega+(15.0*dxC[0]*phiC[1]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals4R2+((54.0*dxC[0]*phiC[1]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals[3]+15.0*dxC[1]*phiC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[4]+(48.0*dxC[0]*phiC[1]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals3R2+(24.0*dxC[1]*phiC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[3])/((15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]); 
  phiC[2] = (((4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+((16.0*dxC[0]*bcVals[3]+4.0*dxC[1])*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]+(16.0*dxC[0]*bcVals3R2+8.0*dxC[1]*bcVals[3])*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals3R2+(32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[3])*omega*volFac+(((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(108.0*dxC[1]*rdx2SqVol[1]+108.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[11]+((30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[4]+(60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[3]-6.0*dxC[1]*rdx2SqVol[1]+30.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiC[2]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+((((-54.0*dxC[0]*rdx2SqVol[1])-90.0*dxC[0]*rdx2SqVol[0])*phiC[2]+(12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]+((-15.0*dxC[1]*rdx2SqVol[1])-15.0*rdx2SqVol[0]*dxC[1])*phiC[2]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+(((-48.0*dxC[0]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0])*phiC[2]+(12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals3R2+(((-24.0*dxC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1])*phiC[2]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*omega+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals[3]+(15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiC[2])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals[3])/((15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]); 
  phiC[3] = (((4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+((8.0*dxC[0]*bcVals[3]+4.0*dxC[1])*rhoC[3]+(32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4])*omega*volFac+(((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((96.0*dxC[0]*rdx2SqVol[1]+96.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[1]*rdx2SqVol[1]+120.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[11]+((24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiC[3]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+((((-54.0*dxC[0]*rdx2SqVol[1])-90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]-15.0*dxC[1]*rdx2SqVol[1]-15.0*rdx2SqVol[0]*dxC[1])*phiC[3]+((6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+(((-48.0*dxC[0]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+((-24.0*dxC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*phiC[3])*omega+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiC[3]*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiC[3]*bcVals[4]+((48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*phiC[3])/((15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+((54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]+15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]); 

}

void MGpoissonFEMDampedGaussSeidel2xSer_UxRobinUyRobin_P1(const double omega, double **dx, const double *bcVals, double **rho, double **phiPrev, double **phi) 
{ 
  // omega:   relaxation parameter.
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phiPrev: (Jacobi-only) iterate cells pointed to by the stencil (only use neighbor cells).
  // phi:     iterate cells pointed to by the stencil (Gauss-Seidel), or cell we are currently updating (Jacobi).

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);
  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = ((((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]+8.0*dxC[1]*bcVals[3]*rhoC[3]+(32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[3])*bcVals10R2+(((16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[4]+32.0*dxC[1]*bcVals[3]*rhoC[3]+(128.0*dxC[1]*rhoC[2]+32.0*dxC[1]*rhoLy[1]+32.0*dxC[1]*rhoLx[1]+128.0*dxC[1]*rhoC[1]+(128.0*rhoLy[0]+32.0*rhoLxLy[0]+128.0*rhoLx[0]+512.0*rhoC[0])*dxC[1])*bcVals[3])*bcVals[9]+(4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+(16.0*dxC[0]*bcVals[3]*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4]+16.0*dxC[0]*bcVals3R2*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals3R2)*bcVals[10]+((16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[4]+32.0*dxC[1]*bcVals[3]*rhoC[3]+(128.0*dxC[1]*rhoC[2]+32.0*dxC[1]*rhoLy[1]+32.0*dxC[1]*rhoLx[1]+128.0*dxC[1]*rhoC[1]+(128.0*rhoLy[0]+32.0*rhoLxLy[0]+128.0*rhoLx[0]+512.0*rhoC[0])*dxC[1])*bcVals[3])*bcVals9R2+((8.0*dxC[0]*rhoC[3]+32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals4R2+(32.0*dxC[0]*bcVals[3]*rhoC[3]+(128.0*dxC[0]*rhoC[2]+32.0*dxC[0]*rhoLy[1]+32.0*dxC[0]*rhoLx[1]+128.0*dxC[0]*rhoC[1]+128.0*dxC[0]*rhoLy[0]+32.0*dxC[0]*rhoLxLy[0]+128.0*dxC[0]*rhoLx[0]+512.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4]+32.0*dxC[0]*bcVals3R2*rhoC[3]+(128.0*dxC[0]*rhoC[2]+32.0*dxC[0]*rhoLy[1]+32.0*dxC[0]*rhoLx[1]+128.0*dxC[0]*rhoC[1]+128.0*dxC[0]*rhoLy[0]+32.0*dxC[0]*rhoLxLy[0]+128.0*dxC[0]*rhoLx[0]+512.0*dxC[0]*rhoC[0])*bcVals3R2)*bcVals[9])*omega*volFac+((((30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[10]+((60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(120.0*dxC[1]*rdx2SqVol[1]-24.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(108.0*dxC[0]*rdx2SqVol[1]-36.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(96.0*dxC[0]*rdx2SqVol[1]-48.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[11]+((30.0*rdx2SqVol[0]*dxC[1]-6.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+((3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-15.0*phiC[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-15.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0]-24.0*phiC[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0]-60.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((108.0*rdx2SqVol[0]*dxC[1]-36.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+((12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0]-90.0*phiC[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0]-54.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((24.0*dxC[1]*phiLy[1]+24.0*dxC[1]*phiLx[1]+(96.0*phiLy[0]+24.0*phiLxLy[0]-48.0*phiLx[0]-144.0*phiC[0])*dxC[1])*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+24.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-48.0*phiLy[0])+24.0*phiLxLy[0]+96.0*phiLx[0]-216.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+((30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[4]+(60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[3])*bcVals[5]+((3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0]-15.0*dxC[0]*phiC[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0]-15.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals4R2+((12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0]-54.0*dxC[0]*phiC[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0]-90.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+((12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0]-48.0*dxC[0]*phiC[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0]-120.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((96.0*rdx2SqVol[0]*dxC[1]-48.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+((12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0]-120.0*phiC[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0]-48.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((24.0*dxC[1]*phiLy[1]+24.0*dxC[1]*phiLx[1]+(96.0*phiLy[0]+24.0*phiLxLy[0]-48.0*phiLx[0]-192.0*phiC[0])*dxC[1])*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+24.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-48.0*phiLy[0])+24.0*phiLxLy[0]+96.0*phiLx[0]-192.0*phiC[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+(((60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[4]+(120.0*dxC[0]*rdx2SqVol[0]-24.0*dxC[0]*rdx2SqVol[1])*bcVals[3])*bcVals[5]+((6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0]-60.0*dxC[0]*phiC[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0]-24.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals4R2+((24.0*dxC[0]*phiLy[1]+24.0*dxC[0]*phiLx[1]+96.0*dxC[0]*phiLy[0]+24.0*dxC[0]*phiLxLy[0]-48.0*dxC[0]*phiLx[0]-216.0*dxC[0]*phiC[0])*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+24.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-48.0*dxC[0]*phiLy[0])+24.0*dxC[0]*phiLxLy[0]+96.0*dxC[0]*phiLx[0]-144.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+((24.0*dxC[0]*phiLy[1]+24.0*dxC[0]*phiLx[1]+96.0*dxC[0]*phiLy[0]+24.0*dxC[0]*phiLxLy[0]-48.0*dxC[0]*phiLx[0]-192.0*dxC[0]*phiC[0])*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+24.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-48.0*dxC[0]*phiLy[0])+24.0*dxC[0]*phiLxLy[0]+96.0*dxC[0]*phiLx[0]-192.0*dxC[0]*phiC[0])*rdx2SqVol[0])*bcVals3R2)*bcVals[9])*omega+((15.0*phiC[0]*dxC[1]*rdx2SqVol[1]+15.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals[4]+(24.0*phiC[0]*dxC[1]*rdx2SqVol[1]+60.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((90.0*phiC[0]*dxC[1]*rdx2SqVol[1]+54.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals[4]+(144.0*phiC[0]*dxC[1]*rdx2SqVol[1]+216.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*phiC[0]*rdx2SqVol[1]+15.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals4R2+(54.0*dxC[0]*phiC[0]*rdx2SqVol[1]+90.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*phiC[0]*rdx2SqVol[1]+120.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((120.0*phiC[0]*dxC[1]*rdx2SqVol[1]+48.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*phiC[0]*dxC[1]*rdx2SqVol[1]+192.0*phiC[0]*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*phiC[0]*rdx2SqVol[1]+24.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*phiC[0]*rdx2SqVol[1]+144.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*phiC[0]*rdx2SqVol[1]+192.0*dxC[0]*phiC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[9])/(((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[9]); 
  phiC[1] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]*bcVals10R2+((16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[4]*bcVals[9]+(4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+(8.0*dxC[0]*bcVals[3]*rhoC[3]+(32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4])*bcVals[10]+(16.0*dxC[1]*rhoC[3]+64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[4]*bcVals9R2+((8.0*dxC[0]*rhoC[3]+32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals4R2+(16.0*dxC[0]*bcVals[3]*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4])*bcVals[9])*omega*volFac+(((30.0*dxC[1]*rdx2SqVol[1]-6.0*rdx2SqVol[0]*dxC[1])*bcVals[4]*bcVals[10]+(60.0*dxC[1]*rdx2SqVol[1]-12.0*rdx2SqVol[0]*dxC[1])*bcVals[4]*bcVals[9]+(30.0*dxC[0]*rdx2SqVol[1]-6.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(48.0*dxC[0]*rdx2SqVol[1]-24.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4])*bcVals[11]+((24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]-15.0*dxC[1]*phiC[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]-15.0*rdx2SqVol[0]*dxC[1]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-24.0*dxC[1]*phiC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[3])*bcVals10R2+(((144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]-90.0*dxC[1]*phiC[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]-54.0*rdx2SqVol[0]*dxC[1]*phiC[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-144.0*dxC[1]*phiC[1]*rdx2SqVol[1])-216.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[3])*bcVals[9]+((24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals[3])*bcVals[5]+((3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]-15.0*dxC[0]*phiC[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]-15.0*dxC[0]*rdx2SqVol[0]*phiC[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+((6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]-54.0*dxC[0]*phiC[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]-90.0*dxC[0]*rdx2SqVol[0]*phiC[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+((-48.0*dxC[0]*phiC[1]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals3R2)*bcVals[10]+((192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+((12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]-120.0*dxC[1]*phiC[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]-48.0*rdx2SqVol[0]*dxC[1]*phiC[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-192.0*dxC[1]*phiC[1]*rdx2SqVol[1])-192.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[3])*bcVals9R2+(((108.0*dxC[0]*rdx2SqVol[1]+108.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals[3])*bcVals[5]+((6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]-60.0*dxC[0]*phiC[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]-24.0*dxC[0]*rdx2SqVol[0]*phiC[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+((12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]-216.0*dxC[0]*phiC[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]-144.0*dxC[0]*rdx2SqVol[0]*phiC[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+((-192.0*dxC[0]*phiC[1]*rdx2SqVol[1])-192.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals3R2)*bcVals[9])*omega+((15.0*dxC[1]*phiC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[4]+(24.0*dxC[1]*phiC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*phiC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[4]+(144.0*dxC[1]*phiC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*phiC[1]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals4R2+(54.0*dxC[0]*phiC[1]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*phiC[1]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*phiC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[4]+(192.0*dxC[1]*phiC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1]*phiC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*phiC[1]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals4R2+(216.0*dxC[0]*phiC[1]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*phiC[1]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0]*phiC[1])*bcVals3R2)*bcVals[9])/(((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[9]); 
  phiC[2] = ((((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]+8.0*dxC[1]*bcVals[3]*rhoC[3]+(32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[3])*bcVals10R2+(((8.0*dxC[1]*rhoC[3]+32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[4]+16.0*dxC[1]*bcVals[3]*rhoC[3]+(64.0*dxC[1]*rhoC[2]+16.0*dxC[1]*rhoLy[1]+16.0*dxC[1]*rhoLx[1]+64.0*dxC[1]*rhoC[1]+(64.0*rhoLy[0]+16.0*rhoLxLy[0]+64.0*rhoLx[0]+256.0*rhoC[0])*dxC[1])*bcVals[3])*bcVals[9]+(4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+(16.0*dxC[0]*bcVals[3]*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4]+16.0*dxC[0]*bcVals3R2*rhoC[3]+(64.0*dxC[0]*rhoC[2]+16.0*dxC[0]*rhoLy[1]+16.0*dxC[0]*rhoLx[1]+64.0*dxC[0]*rhoC[1]+64.0*dxC[0]*rhoLy[0]+16.0*dxC[0]*rhoLxLy[0]+64.0*dxC[0]*rhoLx[0]+256.0*dxC[0]*rhoC[0])*bcVals3R2)*bcVals[10])*omega*volFac+((((60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(108.0*dxC[1]*rdx2SqVol[1]+108.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[11]+((30.0*rdx2SqVol[0]*dxC[1]-6.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+(((-15.0*dxC[1]*rdx2SqVol[1])-15.0*rdx2SqVol[0]*dxC[1])*phiC[2]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+(((-24.0*dxC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1])*phiC[2]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((48.0*rdx2SqVol[0]*dxC[1]-24.0*dxC[1]*rdx2SqVol[1])*bcVals[5]+(((-90.0*dxC[1]*rdx2SqVol[1])-54.0*rdx2SqVol[0]*dxC[1])*phiC[2]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+(((-144.0*dxC[1]*rdx2SqVol[1])-216.0*rdx2SqVol[0]*dxC[1])*phiC[2]+(12.0*dxC[1]*phiLy[1]+12.0*dxC[1]*phiLx[1]+(48.0*phiLy[0]+12.0*phiLxLy[0]-24.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+12.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-24.0*phiLy[0])+12.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+((30.0*dxC[0]*rdx2SqVol[0]-6.0*dxC[0]*rdx2SqVol[1])*bcVals[4]+(60.0*dxC[0]*rdx2SqVol[0]-12.0*dxC[0]*rdx2SqVol[1])*bcVals[3])*bcVals[5]+(((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiC[2]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+(((-54.0*dxC[0]*rdx2SqVol[1])-90.0*dxC[0]*rdx2SqVol[0])*phiC[2]+(12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(((-48.0*dxC[0]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0])*phiC[2]+(12.0*dxC[0]*phiLy[1]+12.0*dxC[0]*phiLx[1]+48.0*dxC[0]*phiLy[0]+12.0*dxC[0]*phiLxLy[0]-24.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+12.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-24.0*dxC[0]*phiLy[0])+12.0*dxC[0]*phiLxLy[0]+48.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+(((-120.0*dxC[1]*rdx2SqVol[1])-48.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals[4]+((-192.0*dxC[1]*rdx2SqVol[1])-192.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals[3])*bcVals9R2+(((-60.0*dxC[0]*rdx2SqVol[1])-24.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals4R2+((-216.0*dxC[0]*rdx2SqVol[1])-144.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals[3]*bcVals[4]+((-192.0*dxC[0]*rdx2SqVol[1])-192.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals3R2)*bcVals[9])*omega+((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*phiC[2]*bcVals[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*phiC[2]*bcVals3R2)*bcVals[9])/(((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[9]); 
  phiC[3] = (((4.0*dxC[1]*rhoC[3]+16.0*dxC[1]*rhoC[2]+4.0*dxC[1]*rhoLy[1]+4.0*dxC[1]*rhoLx[1]+16.0*dxC[1]*rhoC[1]+(16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*dxC[1])*bcVals[4]*bcVals10R2+((8.0*dxC[1]*rhoC[3]+32.0*dxC[1]*rhoC[2]+8.0*dxC[1]*rhoLy[1]+8.0*dxC[1]*rhoLx[1]+32.0*dxC[1]*rhoC[1]+(32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*dxC[1])*bcVals[4]*bcVals[9]+(4.0*dxC[0]*rhoC[3]+16.0*dxC[0]*rhoC[2]+4.0*dxC[0]*rhoLy[1]+4.0*dxC[0]*rhoLx[1]+16.0*dxC[0]*rhoC[1]+16.0*dxC[0]*rhoLy[0]+4.0*dxC[0]*rhoLxLy[0]+16.0*dxC[0]*rhoLx[0]+64.0*dxC[0]*rhoC[0])*bcVals4R2+(8.0*dxC[0]*bcVals[3]*rhoC[3]+(32.0*dxC[0]*rhoC[2]+8.0*dxC[0]*rhoLy[1]+8.0*dxC[0]*rhoLx[1]+32.0*dxC[0]*rhoC[1]+32.0*dxC[0]*rhoLy[0]+8.0*dxC[0]*rhoLxLy[0]+32.0*dxC[0]*rhoLx[0]+128.0*dxC[0]*rhoC[0])*bcVals[3])*bcVals[4])*bcVals[10])*omega*volFac+((((60.0*dxC[1]*rdx2SqVol[1]+24.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(48.0*dxC[1]*rdx2SqVol[1]+120.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(96.0*dxC[0]*rdx2SqVol[1]+96.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4])*bcVals[11]+((24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-15.0*dxC[1]*rdx2SqVol[1])-15.0*rdx2SqVol[0]*dxC[1])*phiC[3]+(3.0*dxC[1]*phiLy[1]+3.0*dxC[1]*phiLx[1]+(12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+3.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-24.0*dxC[1]*rdx2SqVol[1])-60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiC[3])*bcVals10R2+(((96.0*dxC[1]*rdx2SqVol[1]+96.0*rdx2SqVol[0]*dxC[1])*bcVals[5]+(((-90.0*dxC[1]*rdx2SqVol[1])-54.0*rdx2SqVol[0]*dxC[1])*phiC[3]+(6.0*dxC[1]*phiLy[1]+6.0*dxC[1]*phiLx[1]+(24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0])*dxC[1])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLy[1]+6.0*rdx2SqVol[0]*dxC[1]*phiLx[1]+((-12.0*phiLy[0])+6.0*phiLxLy[0]+24.0*phiLx[0])*rdx2SqVol[0]*dxC[1])*bcVals[4]+((-144.0*dxC[1]*rdx2SqVol[1])-216.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiC[3])*bcVals[9]+((24.0*dxC[0]*rdx2SqVol[1]+60.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals[3])*bcVals[5]+(((-15.0*dxC[0]*rdx2SqVol[1])-15.0*dxC[0]*rdx2SqVol[0])*phiC[3]+(3.0*dxC[0]*phiLy[1]+3.0*dxC[0]*phiLx[1]+12.0*dxC[0]*phiLy[0]+3.0*dxC[0]*phiLxLy[0]-6.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+3.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-6.0*dxC[0]*phiLy[0])+3.0*dxC[0]*phiLxLy[0]+12.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals4R2+(((-54.0*dxC[0]*rdx2SqVol[1])-90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*phiC[3]+((6.0*dxC[0]*phiLy[1]+6.0*dxC[0]*phiLx[1]+24.0*dxC[0]*phiLy[0]+6.0*dxC[0]*phiLxLy[0]-12.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLy[1]+6.0*dxC[0]*rdx2SqVol[0]*phiLx[1]+((-12.0*dxC[0]*phiLy[0])+6.0*dxC[0]*phiLxLy[0]+24.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*bcVals[3])*bcVals[4]+((-48.0*dxC[0]*rdx2SqVol[1])-120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2*phiC[3])*bcVals[10]+(((-120.0*dxC[1]*rdx2SqVol[1])-48.0*rdx2SqVol[0]*dxC[1])*phiC[3]*bcVals[4]+((-192.0*dxC[1]*rdx2SqVol[1])-192.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiC[3])*bcVals9R2+(((120.0*dxC[0]*rdx2SqVol[1]+48.0*dxC[0]*rdx2SqVol[0])*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals[3])*bcVals[5]+((-60.0*dxC[0]*rdx2SqVol[1])-24.0*dxC[0]*rdx2SqVol[0])*phiC[3]*bcVals4R2+((-216.0*dxC[0]*rdx2SqVol[1])-144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*phiC[3]*bcVals[4]+((-192.0*dxC[0]*rdx2SqVol[1])-192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2*phiC[3])*bcVals[9])*omega+((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*phiC[3]*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiC[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*phiC[3]*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiC[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*phiC[3]*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*phiC[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2*phiC[3])*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*phiC[3]*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3]*phiC[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*phiC[3]*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*phiC[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2*phiC[3])*bcVals[9])/(((15.0*dxC[1]*rdx2SqVol[1]+15.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(24.0*dxC[1]*rdx2SqVol[1]+60.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals10R2+(((90.0*dxC[1]*rdx2SqVol[1]+54.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(144.0*dxC[1]*rdx2SqVol[1]+216.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals[9]+(15.0*dxC[0]*rdx2SqVol[1]+15.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(54.0*dxC[0]*rdx2SqVol[1]+90.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(48.0*dxC[0]*rdx2SqVol[1]+120.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[10]+((120.0*dxC[1]*rdx2SqVol[1]+48.0*rdx2SqVol[0]*dxC[1])*bcVals[4]+(192.0*dxC[1]*rdx2SqVol[1]+192.0*rdx2SqVol[0]*dxC[1])*bcVals[3])*bcVals9R2+((60.0*dxC[0]*rdx2SqVol[1]+24.0*dxC[0]*rdx2SqVol[0])*bcVals4R2+(216.0*dxC[0]*rdx2SqVol[1]+144.0*dxC[0]*rdx2SqVol[0])*bcVals[3]*bcVals[4]+(192.0*dxC[0]*rdx2SqVol[1]+192.0*dxC[0]*rdx2SqVol[0])*bcVals3R2)*bcVals[9]); 

}

void MGpoissonFEMResidue2xSer_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((16.0*rhoUy[0]+4.0*rhoUxUy[0]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(12.0*phiUy[0]+3.0*phiUxUy[0]+3.0*phiUxLy[0]-6.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+((-6.0*phiUy[0])+3.0*phiUxUy[0]+3.0*phiUxLy[0]+12.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 

}

void MGpoissonFEMResidue2xSer_LxDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[2]; 

}

void MGpoissonFEMResidue2xSer_LxNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 2.0*bcVals[2]-1.0*phiUx[0]+phiC[0]; 

}

void MGpoissonFEMResidue2xSer_LxRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = (2.0*bcVals[2]+(phiC[0]-1.0*phiUx[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])/(bcVals[1]-2.0*bcVals[0]); 

}

void MGpoissonFEMResidue2xSer_UxDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*bcVals[5]+(3.0*phiLy[1]-6.0*phiC[1]+12.0*phiUy[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+12.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiUy[0])-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[5]; 

}

void MGpoissonFEMResidue2xSer_UxNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(6.0*rdx2SqVol[1]+6.0*rdx2SqVol[0])*bcVals[5]+(3.0*phiLy[1]-6.0*phiC[1]+15.0*phiUy[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+12.0*rdx2SqVol[0]*phiC[1]+((-3.0*phiUy[0])-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (-2.0*bcVals[5])+phiC[1]-1.0*phiC[0]; 

}

void MGpoissonFEMResidue2xSer_UxRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = (((4.0*rhoUy[1]+4.0*rhoLy[1]+16.0*rhoC[1]+16.0*rhoUy[0]+16.0*rhoLy[0]+4.0*rhoLxUy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*bcVals[4]+(8.0*rhoUy[1]+8.0*rhoLy[1]+32.0*rhoC[1]+32.0*rhoUy[0]+32.0*rhoLy[0]+8.0*rhoLxUy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*bcVals[3])*volFac+(6.0*rdx2SqVol[1]+6.0*rdx2SqVol[0])*bcVals[5]+((3.0*phiLy[1]-6.0*phiC[1]+15.0*phiUy[0]+12.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+12.0*rdx2SqVol[0]*phiC[1]+((-3.0*phiUy[0])-6.0*phiLy[0]+3.0*phiLxUy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0])*bcVals[4]+((6.0*phiLy[1]-12.0*phiC[1]+24.0*phiUy[0]+24.0*phiLy[0]+6.0*phiLxUy[0]+6.0*phiLxLy[0]-12.0*phiLx[0]-48.0*phiC[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiLy[1]+24.0*rdx2SqVol[0]*phiC[1]+((-12.0*phiUy[0])-12.0*phiLy[0]+6.0*phiLxUy[0]+6.0*phiLxLy[0]+24.0*phiLx[0]-48.0*phiC[0])*rdx2SqVol[0])*bcVals[3])/(18.0*bcVals[4]+36.0*bcVals[3]); 
  resOut[1] = -(1.0*(2.0*bcVals[5]+(phiC[0]-1.0*phiC[1])*bcVals[4]-2.0*phiC[1]*bcVals[3]))/(bcVals[4]+2.0*bcVals[3]); 

}

void MGpoissonFEMResidue2xSer_LyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[8]; 

}

void MGpoissonFEMResidue2xSer_LyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 2.0*bcVals[8]-1.0*phiUy[0]+phiC[0]; 

}

void MGpoissonFEMResidue2xSer_LyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = (2.0*bcVals[8]+(phiC[0]-1.0*phiUy[0])*bcVals[7]-2.0*phiC[0]*bcVals[6])/(bcVals[7]-2.0*bcVals[6]); 

}

void MGpoissonFEMResidue2xSer_UyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*bcVals[11]+(3.0*phiLx[1]+12.0*phiC[1]+3.0*phiUxLy[0]-6.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]-6.0*rdx2SqVol[0]*phiC[1]+(3.0*phiUxLy[0]+12.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[11]; 

}

void MGpoissonFEMResidue2xSer_UyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(6.0*rdx2SqVol[1]+6.0*rdx2SqVol[0])*bcVals[11]+(3.0*phiLx[1]+12.0*phiC[1]+3.0*phiUxLy[0]-3.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]-6.0*rdx2SqVol[0]*phiC[1]+(3.0*phiUxLy[0]+15.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (-2.0*bcVals[11])+phiC[1]-1.0*phiC[0]; 

}

void MGpoissonFEMResidue2xSer_UyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = (((4.0*rhoUx[1]+4.0*rhoLx[1]+16.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*bcVals[10]+(8.0*rhoUx[1]+8.0*rhoLx[1]+32.0*rhoC[1]+8.0*rhoUxLy[0]+32.0*rhoUx[0]+32.0*rhoLy[0]+8.0*rhoLxLy[0]+32.0*rhoLx[0]+128.0*rhoC[0])*bcVals[9])*volFac+(6.0*rdx2SqVol[1]+6.0*rdx2SqVol[0])*bcVals[11]+((3.0*phiLx[1]+12.0*phiC[1]+3.0*phiUxLy[0]-3.0*phiUx[0]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLx[1]-6.0*rdx2SqVol[0]*phiC[1]+(3.0*phiUxLy[0]+15.0*phiUx[0]-6.0*phiLy[0]+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((6.0*phiLx[1]+24.0*phiC[1]+6.0*phiUxLy[0]-12.0*phiUx[0]+24.0*phiLy[0]+6.0*phiLxLy[0]-12.0*phiLx[0]-48.0*phiC[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiLx[1]-12.0*rdx2SqVol[0]*phiC[1]+(6.0*phiUxLy[0]+24.0*phiUx[0]-12.0*phiLy[0]+6.0*phiLxLy[0]+24.0*phiLx[0]-48.0*phiC[0])*rdx2SqVol[0])*bcVals[9])/(18.0*bcVals[10]+36.0*bcVals[9]); 
  resOut[1] = -(1.0*(2.0*bcVals[11]+(phiC[0]-1.0*phiC[1])*bcVals[10]-2.0*phiC[1]*bcVals[9]))/(bcVals[10]+2.0*bcVals[9]); 

}

void MGpoissonFEMResidue2xSer_LxDirichletLyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = -(1.0*(dxC[1]*bcVals[8]+dxC[0]*bcVals[2]-1.0*phiC[0]*dxC[1]-1.0*dxC[0]*phiC[0]))/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMResidue2xSer_LxDirichletLyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[2]; 

}

void MGpoissonFEMResidue2xSer_LxDirichletLyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[2]; 

}

void MGpoissonFEMResidue2xSer_LxNeumannLyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[8]; 

}

void MGpoissonFEMResidue2xSer_LxNeumannLyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 2.0*bcVals[8]+2.0*bcVals[2]-1.0*phiUxUy[0]+phiC[0]; 

}

void MGpoissonFEMResidue2xSer_LxNeumannLyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  resOut[0] = ((2.0*dxC[1]*bcVals[7]-4.0*dxC[1]*bcVals[6]+2.0*dxC[0])*bcVals[8]+(2.0*dxC[1]*bcVals[2]+(phiC[0]-1.0*phiUxUy[0])*dxC[1])*bcVals7R2+(((2.0*phiUxUy[0]-4.0*phiC[0])*dxC[1]-4.0*dxC[1]*bcVals[2])*bcVals[6]+2.0*dxC[0]*bcVals[2]-1.0*dxC[0]*phiUxUy[0]+dxC[0]*phiC[0])*bcVals[7]+4.0*phiC[0]*dxC[1]*bcVals6R2+((-4.0*dxC[0]*bcVals[2])-2.0*dxC[0]*phiC[0])*bcVals[6])/(dxC[1]*bcVals7R2+(dxC[0]-4.0*dxC[1]*bcVals[6])*bcVals[7]+4.0*dxC[1]*bcVals6R2-2.0*dxC[0]*bcVals[6]); 

}

void MGpoissonFEMResidue2xSer_LxRobinLyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[8]; 

}

void MGpoissonFEMResidue2xSer_LxRobinLyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);

  resOut[0] = (((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]+2.0*dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[8]+(2.0*dxC[1]+2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+((phiC[0]-1.0*phiUxUy[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])*dxC[1]+(dxC[0]*phiC[0]-1.0*dxC[0]*phiUxUy[0])*bcVals1R2+(2.0*bcVals[0]*dxC[0]*phiUxUy[0]-4.0*bcVals[0]*dxC[0]*phiC[0])*bcVals[1]+4.0*bcVals0R2*dxC[0]*phiC[0])/((bcVals[1]-2.0*bcVals[0])*dxC[1]+dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0]); 

}

void MGpoissonFEMResidue2xSer_LxRobinLyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);
  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  resOut[0] = (((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]*bcVals[7]+(8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*bcVals[6]+2.0*dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[8]+(2.0*dxC[1]*bcVals[2]+((phiC[0]-1.0*phiUxUy[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])*dxC[1])*bcVals7R2+((((2.0*phiUxUy[0]-4.0*phiC[0])*bcVals[1]+8.0*bcVals[0]*phiC[0])*dxC[1]-4.0*dxC[1]*bcVals[2])*bcVals[6]+(2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+(dxC[0]*phiC[0]-1.0*dxC[0]*phiUxUy[0])*bcVals1R2+(2.0*bcVals[0]*dxC[0]*phiUxUy[0]-4.0*bcVals[0]*dxC[0]*phiC[0])*bcVals[1]+4.0*bcVals0R2*dxC[0]*phiC[0])*bcVals[7]+(4.0*phiC[0]*bcVals[1]-8.0*bcVals[0]*phiC[0])*dxC[1]*bcVals6R2+((8.0*bcVals[0]*dxC[0]-4.0*dxC[0]*bcVals[1])*bcVals[2]-2.0*dxC[0]*phiC[0]*bcVals1R2+8.0*bcVals[0]*dxC[0]*phiC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0]*phiC[0])*bcVals[6])/((bcVals[1]-2.0*bcVals[0])*dxC[1]*bcVals7R2+((8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*bcVals[6]+dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0])*bcVals[7]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals6R2+((-2.0*dxC[0]*bcVals1R2)+8.0*bcVals[0]*dxC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0])*bcVals[6]); 

}

void MGpoissonFEMResidue2xSer_LxDirichletUyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[2]; 
  resOut[1] = -(1.0*(dxC[1]*bcVals[11]+dxC[0]*bcVals[2]+((-1.0*dxC[1])-1.0*dxC[0])*phiC[1]))/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMResidue2xSer_LxDirichletUyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[2]; 
  resOut[1] = phiC[1]-1.0*bcVals[2]; 

}

void MGpoissonFEMResidue2xSer_LxDirichletUyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[2]; 
  resOut[1] = phiC[1]-1.0*bcVals[2]; 

}

void MGpoissonFEMResidue2xSer_LxNeumannUyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoUx[1]+8.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+8.0*rhoLy[0]+32.0*rhoC[0])*volFac+3.0*rdx2SqVol[1]*bcVals[11]+(6.0*phiC[1]+3.0*phiUxLy[0]-6.0*phiUx[0]+6.0*phiLy[0]-12.0*phiC[0])*rdx2SqVol[1]+(3.0*phiUxLy[0]+6.0*phiUx[0]-3.0*phiLy[0]-6.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[11]; 

}

void MGpoissonFEMResidue2xSer_LxNeumannUyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 2.0*bcVals[2]-1.0*phiUx[0]+phiC[0]; 
  resOut[1] = (-2.0*bcVals[11])+2.0*bcVals[2]+phiC[1]-1.0*phiUx[0]; 

}

void MGpoissonFEMResidue2xSer_LxNeumannUyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  resOut[0] = 2.0*bcVals[2]-1.0*phiUx[0]+phiC[0]; 
  resOut[1] = -(1.0*((2.0*dxC[1]*bcVals[10]+4.0*dxC[1]*bcVals[9]-2.0*dxC[0])*bcVals[11]+((-2.0*dxC[1]*bcVals[2])-1.0*dxC[1]*phiC[1]+phiUx[0]*dxC[1])*bcVals10R2+(((-4.0*dxC[1]*bcVals[2])-4.0*dxC[1]*phiC[1]+2.0*phiUx[0]*dxC[1])*bcVals[9]+2.0*dxC[0]*bcVals[2]+dxC[0]*phiC[1]-1.0*dxC[0]*phiUx[0])*bcVals[10]-4.0*dxC[1]*phiC[1]*bcVals9R2+(4.0*dxC[0]*bcVals[2]+2.0*dxC[0]*phiC[1])*bcVals[9]))/(dxC[1]*bcVals10R2+(4.0*dxC[1]*bcVals[9]-1.0*dxC[0])*bcVals[10]+4.0*dxC[1]*bcVals9R2-2.0*dxC[0]*bcVals[9]); 

}

void MGpoissonFEMResidue2xSer_LxRobinUyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoUx[1]+8.0*rhoC[1]+4.0*rhoUxLy[0]+16.0*rhoUx[0]+8.0*rhoLy[0]+32.0*rhoC[0])*volFac+3.0*rdx2SqVol[1]*bcVals[11]+(6.0*phiC[1]+3.0*phiUxLy[0]-6.0*phiUx[0]+6.0*phiLy[0]-12.0*phiC[0])*rdx2SqVol[1]+(3.0*phiUxLy[0]+6.0*phiUx[0]-3.0*phiLy[0]-6.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[11]; 

}

void MGpoissonFEMResidue2xSer_LxRobinUyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);

  resOut[0] = (2.0*bcVals[2]+(phiC[0]-1.0*phiUx[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])/(bcVals[1]-2.0*bcVals[0]); 
  resOut[1] = -(1.0*(((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]-2.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[11]+((-2.0*dxC[1])+2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+((2.0*bcVals[0]-1.0*bcVals[1])*dxC[1]+dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0])*phiC[1]+phiUx[0]*bcVals[1]*dxC[1]-1.0*dxC[0]*phiUx[0]*bcVals1R2+2.0*bcVals[0]*dxC[0]*phiUx[0]*bcVals[1]))/((bcVals[1]-2.0*bcVals[0])*dxC[1]-1.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0]); 

}

void MGpoissonFEMResidue2xSer_LxRobinUyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);
  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  resOut[0] = (2.0*bcVals[2]+(phiC[0]-1.0*phiUx[0])*bcVals[1]-2.0*bcVals[0]*phiC[0])/(bcVals[1]-2.0*bcVals[0]); 
  resOut[1] = -(1.0*(((2.0*bcVals[1]-4.0*bcVals[0])*dxC[1]*bcVals[10]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals[9]-2.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1])*bcVals[11]+((-2.0*dxC[1]*bcVals[2])+(2.0*bcVals[0]-1.0*bcVals[1])*dxC[1]*phiC[1]+phiUx[0]*bcVals[1]*dxC[1])*bcVals10R2+(((-4.0*dxC[1]*bcVals[2])+(8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*phiC[1]+2.0*phiUx[0]*bcVals[1]*dxC[1])*bcVals[9]+(2.0*dxC[0]*bcVals[1]-4.0*bcVals[0]*dxC[0])*bcVals[2]+(dxC[0]*bcVals1R2-4.0*bcVals[0]*dxC[0]*bcVals[1]+4.0*bcVals0R2*dxC[0])*phiC[1]-1.0*dxC[0]*phiUx[0]*bcVals1R2+2.0*bcVals[0]*dxC[0]*phiUx[0]*bcVals[1])*bcVals[10]+(8.0*bcVals[0]-4.0*bcVals[1])*dxC[1]*phiC[1]*bcVals9R2+((4.0*dxC[0]*bcVals[1]-8.0*bcVals[0]*dxC[0])*bcVals[2]+(2.0*dxC[0]*bcVals1R2-8.0*bcVals[0]*dxC[0]*bcVals[1]+8.0*bcVals0R2*dxC[0])*phiC[1])*bcVals[9]))/((bcVals[1]-2.0*bcVals[0])*dxC[1]*bcVals10R2+((4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals[9]-1.0*dxC[0]*bcVals1R2+4.0*bcVals[0]*dxC[0]*bcVals[1]-4.0*bcVals0R2*dxC[0])*bcVals[10]+(4.0*bcVals[1]-8.0*bcVals[0])*dxC[1]*bcVals9R2+((-2.0*dxC[0]*bcVals1R2)+8.0*bcVals[0]*dxC[0]*bcVals[1]-8.0*bcVals0R2*dxC[0])*bcVals[9]); 

}

void MGpoissonFEMResidue2xSer_UxDirichletLyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[8]; 
  resOut[1] = -(1.0*(dxC[1]*bcVals[8]+dxC[0]*bcVals[5]+((-1.0*dxC[1])-1.0*dxC[0])*phiC[1]))/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMResidue2xSer_UxDirichletLyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoUy[1]+8.0*rhoC[1]+16.0*rhoUy[0]+4.0*rhoLxUy[0]+8.0*rhoLx[0]+32.0*rhoC[0])*volFac+3.0*rdx2SqVol[0]*bcVals[5]+(6.0*phiUy[0]+3.0*phiLxUy[0]-3.0*phiLx[0]-6.0*phiC[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiUy[0])+3.0*phiLxUy[0]+6.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[5]; 

}

void MGpoissonFEMResidue2xSer_UxDirichletLyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoUy[1]+8.0*rhoC[1]+16.0*rhoUy[0]+4.0*rhoLxUy[0]+8.0*rhoLx[0]+32.0*rhoC[0])*volFac+3.0*rdx2SqVol[0]*bcVals[5]+(6.0*phiUy[0]+3.0*phiLxUy[0]-3.0*phiLx[0]-6.0*phiC[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiUy[0])+3.0*phiLxUy[0]+6.0*phiLx[0]-12.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[5]; 

}

void MGpoissonFEMResidue2xSer_UxNeumannLyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[8]; 
  resOut[1] = phiC[1]-1.0*bcVals[8]; 

}

void MGpoissonFEMResidue2xSer_UxNeumannLyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 2.0*bcVals[8]-1.0*phiUy[0]+phiC[0]; 
  resOut[1] = 2.0*bcVals[8]-2.0*bcVals[5]+phiC[1]-1.0*phiUy[0]; 

}

void MGpoissonFEMResidue2xSer_UxNeumannLyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  resOut[0] = (2.0*bcVals[8]+(phiC[0]-1.0*phiUy[0])*bcVals[7]-2.0*phiC[0]*bcVals[6])/(bcVals[7]-2.0*bcVals[6]); 
  resOut[1] = ((2.0*dxC[1]*bcVals[7]-4.0*dxC[1]*bcVals[6]-2.0*dxC[0])*bcVals[8]+((-2.0*dxC[1]*bcVals[5])+dxC[1]*phiC[1]-1.0*phiUy[0]*dxC[1])*bcVals7R2+((4.0*dxC[1]*bcVals[5]-4.0*dxC[1]*phiC[1]+2.0*phiUy[0]*dxC[1])*bcVals[6]+2.0*dxC[0]*bcVals[5]-1.0*dxC[0]*phiC[1]+dxC[0]*phiUy[0])*bcVals[7]+4.0*dxC[1]*phiC[1]*bcVals6R2+(2.0*dxC[0]*phiC[1]-4.0*dxC[0]*bcVals[5])*bcVals[6])/(dxC[1]*bcVals7R2+((-4.0*dxC[1]*bcVals[6])-1.0*dxC[0])*bcVals[7]+4.0*dxC[1]*bcVals6R2+2.0*dxC[0]*bcVals[6]); 

}

void MGpoissonFEMResidue2xSer_UxRobinLyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = phiC[0]-1.0*bcVals[8]; 
  resOut[1] = phiC[1]-1.0*bcVals[8]; 

}

void MGpoissonFEMResidue2xSer_UxRobinLyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);

  resOut[0] = 2.0*bcVals[8]-1.0*phiUy[0]+phiC[0]; 
  resOut[1] = ((2.0*dxC[0]*bcVals4R2+(4.0*dxC[0]*bcVals[3]-2.0*dxC[1])*bcVals[4]-4.0*dxC[1]*bcVals[3])*bcVals[8]+((-2.0*dxC[0]*bcVals[4])-4.0*dxC[0]*bcVals[3]+2.0*dxC[1])*bcVals[5]+(dxC[0]*phiC[1]-1.0*dxC[0]*phiUy[0])*bcVals4R2+((4.0*dxC[0]*phiC[1]-2.0*dxC[0]*phiUy[0])*bcVals[3]-1.0*dxC[1]*phiC[1]+phiUy[0]*dxC[1])*bcVals[4]+4.0*dxC[0]*phiC[1]*bcVals3R2-2.0*dxC[1]*phiC[1]*bcVals[3])/(dxC[0]*bcVals4R2+(4.0*dxC[0]*bcVals[3]-1.0*dxC[1])*bcVals[4]+4.0*dxC[0]*bcVals3R2-2.0*dxC[1]*bcVals[3]); 

}

void MGpoissonFEMResidue2xSer_UxRobinLyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);
  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  resOut[0] = (2.0*bcVals[8]+(phiC[0]-1.0*phiUy[0])*bcVals[7]-2.0*phiC[0]*bcVals[6])/(bcVals[7]-2.0*bcVals[6]); 
  resOut[1] = (((2.0*dxC[1]*bcVals[4]+4.0*dxC[1]*bcVals[3])*bcVals[7]+((-4.0*dxC[1]*bcVals[4])-8.0*dxC[1]*bcVals[3])*bcVals[6]-2.0*dxC[0]*bcVals4R2-4.0*dxC[0]*bcVals[3]*bcVals[4])*bcVals[8]+((-2.0*dxC[1]*bcVals[5])+(dxC[1]*phiC[1]-1.0*phiUy[0]*dxC[1])*bcVals[4]+2.0*dxC[1]*phiC[1]*bcVals[3])*bcVals7R2+((4.0*dxC[1]*bcVals[5]+(2.0*phiUy[0]*dxC[1]-4.0*dxC[1]*phiC[1])*bcVals[4]-8.0*dxC[1]*phiC[1]*bcVals[3])*bcVals[6]+(2.0*dxC[0]*bcVals[4]+4.0*dxC[0]*bcVals[3])*bcVals[5]+(dxC[0]*phiUy[0]-1.0*dxC[0]*phiC[1])*bcVals4R2+(2.0*dxC[0]*phiUy[0]-4.0*dxC[0]*phiC[1])*bcVals[3]*bcVals[4]-4.0*dxC[0]*phiC[1]*bcVals3R2)*bcVals[7]+(4.0*dxC[1]*phiC[1]*bcVals[4]+8.0*dxC[1]*phiC[1]*bcVals[3])*bcVals6R2+(((-4.0*dxC[0]*bcVals[4])-8.0*dxC[0]*bcVals[3])*bcVals[5]+2.0*dxC[0]*phiC[1]*bcVals4R2+8.0*dxC[0]*phiC[1]*bcVals[3]*bcVals[4]+8.0*dxC[0]*phiC[1]*bcVals3R2)*bcVals[6])/((dxC[1]*bcVals[4]+2.0*dxC[1]*bcVals[3])*bcVals7R2+(((-4.0*dxC[1]*bcVals[4])-8.0*dxC[1]*bcVals[3])*bcVals[6]-1.0*dxC[0]*bcVals4R2-4.0*dxC[0]*bcVals[3]*bcVals[4]-4.0*dxC[0]*bcVals3R2)*bcVals[7]+(4.0*dxC[1]*bcVals[4]+8.0*dxC[1]*bcVals[3])*bcVals6R2+(2.0*dxC[0]*bcVals4R2+8.0*dxC[0]*bcVals[3]*bcVals[4]+8.0*dxC[0]*bcVals3R2)*bcVals[6]); 

}

void MGpoissonFEMResidue2xSer_UxDirichletUyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*phiC[3]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*phiC[2]+(3.0*phiLy[1]+3.0*phiLx[1]-6.0*phiC[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+12.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[5]; 
  resOut[2] = phiC[2]-1.0*bcVals[11]; 
  resOut[3] = -(1.0*(dxC[1]*bcVals[11]+dxC[0]*bcVals[5]+((-1.0*dxC[1])-1.0*dxC[0])*phiC[3]))/(dxC[1]+dxC[0]); 

}

void MGpoissonFEMResidue2xSer_UxDirichletUyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*phiC[3]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*phiC[2]+(3.0*phiLy[1]+3.0*phiLx[1]-6.0*phiC[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+12.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[5]; 
  resOut[2] = 0.05555555555555555*((8.0*rhoC[3]+32.0*rhoC[2]+8.0*rhoLx[1]+4.0*rhoC[1]+4.0*rhoLx[0]+16.0*rhoC[0])*volFac+6.0*rdx2SqVol[0]*phiC[3]+((-6.0*rdx2SqVol[1])-12.0*rdx2SqVol[0])*phiC[2]+((-3.0*phiLx[1])+3.0*phiLx[0]+6.0*phiC[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiLx[1]+3.0*rdx2SqVol[0]*phiC[1]+(3.0*phiLx[0]-6.0*phiC[0])*rdx2SqVol[0]); 
  resOut[3] = phiC[3]-1.0*bcVals[5]; 

}

void MGpoissonFEMResidue2xSer_UxDirichletUyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*phiC[3]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*phiC[2]+(3.0*phiLy[1]+3.0*phiLx[1]-6.0*phiC[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+12.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[5]; 
  resOut[2] = 0.05555555555555555*((8.0*rhoC[3]+32.0*rhoC[2]+8.0*rhoLx[1]+4.0*rhoC[1]+4.0*rhoLx[0]+16.0*rhoC[0])*volFac+6.0*rdx2SqVol[0]*phiC[3]+((-6.0*rdx2SqVol[1])-12.0*rdx2SqVol[0])*phiC[2]+((-3.0*phiLx[1])+3.0*phiLx[0]+6.0*phiC[0])*rdx2SqVol[1]+6.0*rdx2SqVol[0]*phiLx[1]+3.0*rdx2SqVol[0]*phiC[1]+(3.0*phiLx[0]-6.0*phiC[0])*rdx2SqVol[0]); 
  resOut[3] = phiC[3]-1.0*bcVals[5]; 

}

void MGpoissonFEMResidue2xSer_UxNeumannUyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*phiC[3]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*phiC[2]+(3.0*phiLy[1]+3.0*phiLx[1]-6.0*phiC[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+12.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.05555555555555555*((8.0*rhoC[3]+4.0*rhoC[2]+8.0*rhoLy[1]+32.0*rhoC[1]+4.0*rhoLy[0]+16.0*rhoC[0])*volFac+6.0*rdx2SqVol[1]*phiC[3]+3.0*rdx2SqVol[1]*phiC[2]+(6.0*phiLy[1]-12.0*phiC[1]+3.0*phiLy[0]-6.0*phiC[0])*rdx2SqVol[1]-3.0*rdx2SqVol[0]*phiLy[1]-6.0*rdx2SqVol[0]*phiC[1]+(3.0*phiLy[0]+6.0*phiC[0])*rdx2SqVol[0]); 
  resOut[2] = phiC[2]-1.0*bcVals[11]; 
  resOut[3] = phiC[3]-1.0*bcVals[11]; 

}

void MGpoissonFEMResidue2xSer_UxNeumannUyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*phiC[3]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*phiC[2]+(3.0*phiLy[1]+3.0*phiLx[1]-6.0*phiC[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+12.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (-2.0*bcVals[5])+phiC[1]-1.0*phiC[0]; 
  resOut[2] = (-2.0*bcVals[11])+phiC[2]-1.0*phiC[0]; 
  resOut[3] = (-2.0*bcVals[11])-2.0*bcVals[5]+phiC[3]-1.0*phiC[0]; 

}

void MGpoissonFEMResidue2xSer_UxNeumannUyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  resOut[0] = 0.05555555555555555*((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*phiC[3]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*phiC[2]+(3.0*phiLy[1]+3.0*phiLx[1]-6.0*phiC[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+12.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (-2.0*bcVals[5])+phiC[1]-1.0*phiC[0]; 
  resOut[2] = -(1.0*(2.0*bcVals[11]+(phiC[0]-1.0*phiC[2])*bcVals[10]-2.0*phiC[2]*bcVals[9]))/(bcVals[10]+2.0*bcVals[9]); 
  resOut[3] = -(1.0*((2.0*dxC[1]*bcVals[10]+4.0*dxC[1]*bcVals[9]+2.0*dxC[0])*bcVals[11]+(2.0*dxC[1]*bcVals[5]-1.0*dxC[1]*phiC[3]+phiC[0]*dxC[1])*bcVals10R2+((4.0*dxC[1]*bcVals[5]-4.0*dxC[1]*phiC[3]+2.0*phiC[0]*dxC[1])*bcVals[9]+2.0*dxC[0]*bcVals[5]-1.0*dxC[0]*phiC[3]+dxC[0]*phiC[0])*bcVals[10]-4.0*dxC[1]*phiC[3]*bcVals9R2+(4.0*dxC[0]*bcVals[5]-2.0*dxC[0]*phiC[3])*bcVals[9]))/(dxC[1]*bcVals10R2+(4.0*dxC[1]*bcVals[9]+dxC[0])*bcVals[10]+4.0*dxC[1]*bcVals9R2+2.0*dxC[0]*bcVals[9]); 

}

void MGpoissonFEMResidue2xSer_UxRobinUyDirichlet_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 


  resOut[0] = 0.05555555555555555*((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*phiC[3]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*phiC[2]+(3.0*phiLy[1]+3.0*phiLx[1]-6.0*phiC[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+12.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.05555555555555555*((8.0*rhoC[3]+4.0*rhoC[2]+8.0*rhoLy[1]+32.0*rhoC[1]+4.0*rhoLy[0]+16.0*rhoC[0])*volFac+6.0*rdx2SqVol[1]*phiC[3]+3.0*rdx2SqVol[1]*phiC[2]+(6.0*phiLy[1]-12.0*phiC[1]+3.0*phiLy[0]-6.0*phiC[0])*rdx2SqVol[1]-3.0*rdx2SqVol[0]*phiLy[1]-6.0*rdx2SqVol[0]*phiC[1]+(3.0*phiLy[0]+6.0*phiC[0])*rdx2SqVol[0]); 
  resOut[2] = phiC[2]-1.0*bcVals[11]; 
  resOut[3] = phiC[3]-1.0*bcVals[11]; 

}

void MGpoissonFEMResidue2xSer_UxRobinUyNeumann_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);

  resOut[0] = 0.05555555555555555*((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*phiC[3]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*phiC[2]+(3.0*phiLy[1]+3.0*phiLx[1]-6.0*phiC[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+12.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = -(1.0*(2.0*bcVals[5]+(phiC[0]-1.0*phiC[1])*bcVals[4]-2.0*phiC[1]*bcVals[3]))/(bcVals[4]+2.0*bcVals[3]); 
  resOut[2] = (-2.0*bcVals[11])+phiC[2]-1.0*phiC[0]; 
  resOut[3] = -(1.0*((2.0*dxC[0]*bcVals4R2+(4.0*dxC[0]*bcVals[3]+2.0*dxC[1])*bcVals[4]+4.0*dxC[1]*bcVals[3])*bcVals[11]+(2.0*dxC[0]*bcVals[4]+4.0*dxC[0]*bcVals[3]+2.0*dxC[1])*bcVals[5]+(dxC[0]*phiC[0]-1.0*dxC[0]*phiC[3])*bcVals4R2+(((-4.0*dxC[0]*bcVals[3])-1.0*dxC[1])*phiC[3]+2.0*dxC[0]*phiC[0]*bcVals[3]+phiC[0]*dxC[1])*bcVals[4]+((-4.0*dxC[0]*bcVals3R2)-2.0*dxC[1]*bcVals[3])*phiC[3]))/(dxC[0]*bcVals4R2+(4.0*dxC[0]*bcVals[3]+dxC[1])*bcVals[4]+4.0*dxC[0]*bcVals3R2+2.0*dxC[1]*bcVals[3]); 

}

void MGpoissonFEMResidue2xSer_UxRobinUyRobin_P1(double **dx, const double *bcVals, double **rho, double **phi, double *resOut) 
{ 
  // dx:      cell lengths of cells pointed to by the stencil.
  // bcVals:  values to impose as BCs.
  // rho:     right-side source in the current cell.
  // phi:     iterate cells pointed to by the stencil.
  // resOut:  residue in nodes stored in current cell.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double rdx2SqVol[2]; 
  rdx2SqVol[0] = volFac*4.0/(dxC[0]*dxC[0]); 
  rdx2SqVol[1] = volFac*4.0/(dxC[1]*dxC[1]); 

  double *phiC = phi[0]; 
  double *rhoC = rho[0]; 
  double *phiLx = phi[1]; 
  double *rhoLx = rho[1]; 
  double *phiUx = phi[2]; 
  double *rhoUx = rho[2]; 
  double *phiLy = phi[3]; 
  double *rhoLy = rho[3]; 
  double *phiUy = phi[4]; 
  double *rhoUy = rho[4]; 
  double *phiLxLy = phi[5]; 
  double *rhoLxLy = rho[5]; 
  double *phiLxUy = phi[6]; 
  double *rhoLxUy = rho[6]; 
  double *phiUxLy = phi[7]; 
  double *rhoUxLy = rho[7]; 
  double *phiUxUy = phi[8]; 
  double *rhoUxUy = rho[8]; 

  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);
  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  resOut[0] = 0.05555555555555555*((4.0*rhoC[3]+16.0*rhoC[2]+4.0*rhoLy[1]+4.0*rhoLx[1]+16.0*rhoC[1]+16.0*rhoLy[0]+4.0*rhoLxLy[0]+16.0*rhoLx[0]+64.0*rhoC[0])*volFac+(3.0*rdx2SqVol[1]+3.0*rdx2SqVol[0])*phiC[3]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*phiC[2]+(3.0*phiLy[1]+3.0*phiLx[1]-6.0*phiC[1]+12.0*phiLy[0]+3.0*phiLxLy[0]-6.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[1]+3.0*rdx2SqVol[0]*phiLy[1]+3.0*rdx2SqVol[0]*phiLx[1]+12.0*rdx2SqVol[0]*phiC[1]+((-6.0*phiLy[0])+3.0*phiLxLy[0]+12.0*phiLx[0]-24.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = -(1.0*(2.0*bcVals[5]+(phiC[0]-1.0*phiC[1])*bcVals[4]-2.0*phiC[1]*bcVals[3]))/(bcVals[4]+2.0*bcVals[3]); 
  resOut[2] = -(1.0*(2.0*bcVals[11]+(phiC[0]-1.0*phiC[2])*bcVals[10]-2.0*phiC[2]*bcVals[9]))/(bcVals[10]+2.0*bcVals[9]); 
  resOut[3] = -(1.0*(((2.0*dxC[1]*bcVals[4]+4.0*dxC[1]*bcVals[3])*bcVals[10]+(4.0*dxC[1]*bcVals[4]+8.0*dxC[1]*bcVals[3])*bcVals[9]+2.0*dxC[0]*bcVals4R2+4.0*dxC[0]*bcVals[3]*bcVals[4])*bcVals[11]+(2.0*dxC[1]*bcVals[5]+(phiC[0]*dxC[1]-1.0*dxC[1]*phiC[3])*bcVals[4]-2.0*dxC[1]*bcVals[3]*phiC[3])*bcVals10R2+((4.0*dxC[1]*bcVals[5]+(2.0*phiC[0]*dxC[1]-4.0*dxC[1]*phiC[3])*bcVals[4]-8.0*dxC[1]*bcVals[3]*phiC[3])*bcVals[9]+(2.0*dxC[0]*bcVals[4]+4.0*dxC[0]*bcVals[3])*bcVals[5]+(dxC[0]*phiC[0]-1.0*dxC[0]*phiC[3])*bcVals4R2+(2.0*dxC[0]*phiC[0]*bcVals[3]-4.0*dxC[0]*bcVals[3]*phiC[3])*bcVals[4]-4.0*dxC[0]*bcVals3R2*phiC[3])*bcVals[10]+((-4.0*dxC[1]*phiC[3]*bcVals[4])-8.0*dxC[1]*bcVals[3]*phiC[3])*bcVals9R2+((4.0*dxC[0]*bcVals[4]+8.0*dxC[0]*bcVals[3])*bcVals[5]-2.0*dxC[0]*phiC[3]*bcVals4R2-8.0*dxC[0]*bcVals[3]*phiC[3]*bcVals[4]-8.0*dxC[0]*bcVals3R2*phiC[3])*bcVals[9]))/((dxC[1]*bcVals[4]+2.0*dxC[1]*bcVals[3])*bcVals10R2+((4.0*dxC[1]*bcVals[4]+8.0*dxC[1]*bcVals[3])*bcVals[9]+dxC[0]*bcVals4R2+4.0*dxC[0]*bcVals[3]*bcVals[4]+4.0*dxC[0]*bcVals3R2)*bcVals[10]+(4.0*dxC[1]*bcVals[4]+8.0*dxC[1]*bcVals[3])*bcVals9R2+(2.0*dxC[0]*bcVals4R2+8.0*dxC[0]*bcVals[3]*bcVals[4]+8.0*dxC[0]*bcVals3R2)*bcVals[9]); 

}

