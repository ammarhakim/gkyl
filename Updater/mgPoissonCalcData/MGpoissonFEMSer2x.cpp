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

  double *fldF1 = fldF[1];

  fldC[0] += 0.5*fldF1[0]; 
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
  fldC[1] += 0.5*fldF0[1]; 
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

  double *fldF2 = fldF[2];

  fldC[0] += 0.25*fldF2[0]; 
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
  fldC[1] += 0.5*fldF0[1]; 
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


}

void MGpoissonFEMRestrict2xSer_LxDirichletLyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF1 = fldF[1];

  fldC[0] += 0.5*fldF1[0]; 
}

void MGpoissonFEMRestrict2xSer_LxDirichletLyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF1 = fldF[1];

  fldC[0] += 0.5*fldF1[0]; 
}

void MGpoissonFEMRestrict2xSer_LxNeumannLyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF2 = fldF[2];

  fldC[0] += 0.5*fldF2[0]; 
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

  double *fldF2 = fldF[2];

  fldC[0] += 0.5*fldF2[0]; 
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

  double *fldF1 = fldF[1];

  fldC[0] += 0.5*fldF1[0]; 
}

void MGpoissonFEMRestrict2xSer_LxDirichletUyNeumann_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF1 = fldF[1];

  fldC[0] += 0.5*fldF1[0]; 
  fldC[1] += fldF1[1]+0.5*fldF1[0]; 
}

void MGpoissonFEMRestrict2xSer_LxDirichletUyRobin_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF1 = fldF[1];

  fldC[0] += 0.5*fldF1[0]; 
  fldC[1] += fldF1[1]+0.5*fldF1[0]; 
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
  fldC[1] += 0.5*fldF0[1]; 
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
  fldC[1] += 0.5*fldF0[1]; 
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

  double *fldF2 = fldF[2];

  fldC[0] += 0.25*fldF2[0]; 
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
  fldC[1] += 0.5*fldF0[1]; 
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
  fldC[1] += 0.5*fldF0[1]; 
}

void MGpoissonFEMRestrict2xSer_UxNeumannLyDirichlet_P1(double **fldF, double *fldC) 
{ 
  // fldF: fine-grid field in cells pointed to by the stencil.
  // fldC: coarse-grid field.

  double *fldF2 = fldF[2];

  fldC[0] += 0.25*fldF2[0]; 
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

  double *fldF2 = fldF[2];

  fldC[0] += 0.25*fldF2[0]; 
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
  fldC[1] += 0.5*fldF0[1]; 
  fldC[2] += 0.5*fldF0[2]; 
  fldC[3] += fldF0[3]; 
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
  fldC[1] += 0.5*fldF0[1]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[1]; 
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
  fldC[1] += 0.5*fldF0[1]; 
  fldC[2] += 0.5*fldF0[2]+0.5*fldF2[1]+fldF1[1]+0.25*fldF2[0]+0.5*fldF1[0]+0.25*fldF0[0]; 
  fldC[3] += fldF0[3]+0.5*fldF0[1]; 
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
  fldC[2] += 0.5*fldF0[2]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]; 
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
  fldC[2] += 0.5*fldF0[2]; 
  fldC[3] += fldF0[3]+0.5*fldF0[2]; 
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

void MGpoissonFEM_DGtoFEM_2xSer_LxNonPeriodic_P1(const double *dgFld, double **femOut) 
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

void MGpoissonFEM_DGtoFEM_2xSer_UxNonPeriodic_P1(const double *dgFld, double **femOut) 
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

void MGpoissonFEM_DGtoFEM_2xSer_LyNonPeriodic_P1(const double *dgFld, double **femOut) 
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

void MGpoissonFEM_DGtoFEM_2xSer_UyNonPeriodic_P1(const double *dgFld, double **femOut) 
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

void MGpoissonFEM_DGtoFEM_2xSer_LxNonPeriodicLyNonPeriodic_P1(const double *dgFld, double **femOut) 
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

void MGpoissonFEM_DGtoFEM_2xSer_LxNonPeriodicUyNonPeriodic_P1(const double *dgFld, double **femOut) 
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

void MGpoissonFEM_DGtoFEM_2xSer_UxNonPeriodicLyNonPeriodic_P1(const double *dgFld, double **femOut) 
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

void MGpoissonFEM_DGtoFEM_2xSer_UxNonPeriodicUyNonPeriodic_P1(const double *dgFld, double **femOut) 
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

void MGpoissonFEMproject2xSer_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 
  double *femFldLy = femFld[3]; 
  double *femFldUy = femFld[4]; 
  double *femFldLxLy = femFld[5]; 
  double *femFldLxUy = femFld[6]; 
  double *femFldUxLy = femFld[7]; 
  double *femFldUxUy = femFld[8]; 

  out[0] = 0.1111111111111111*(8.0*femFldUy[0]+2.0*(femFldUxUy[0]+femFldUxLy[0])+8.0*(femFldUx[0]+femFldLy[0])+2.0*(femFldLxUy[0]+femFldLxLy[0])+8.0*femFldLx[0]+32.0*femFldC[0])*volFac; 

}

void MGpoissonFEMproject2xSer_LxNonPeriodic_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 
  double *femFldLy = femFld[3]; 
  double *femFldUy = femFld[4]; 
  double *femFldLxLy = femFld[5]; 
  double *femFldLxUy = femFld[6]; 
  double *femFldUxLy = femFld[7]; 
  double *femFldUxUy = femFld[8]; 

  out[0] = 0.1111111111111111*(4.0*femFldUy[0]+2.0*(femFldUxUy[0]+femFldUxLy[0])+8.0*femFldUx[0]+4.0*femFldLy[0]+16.0*femFldC[0])*volFac; 

}

void MGpoissonFEMproject2xSer_UxNonPeriodic_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 
  double *femFldLy = femFld[3]; 
  double *femFldUy = femFld[4]; 
  double *femFldLxLy = femFld[5]; 
  double *femFldLxUy = femFld[6]; 
  double *femFldUxLy = femFld[7]; 
  double *femFldUxUy = femFld[8]; 

  out[0] = 0.1111111111111111*(2.0*(femFldUy[1]+femFldLy[1])+8.0*(femFldC[1]+femFldUy[0]+femFldLy[0])+2.0*(femFldLxUy[0]+femFldLxLy[0])+8.0*femFldLx[0]+32.0*femFldC[0])*volFac; 
  out[1] = 0.1111111111111111*(4.0*(femFldUy[1]+femFldLy[1])+16.0*femFldC[1]+2.0*(femFldUy[0]+femFldLy[0])+8.0*femFldC[0])*volFac; 

}

void MGpoissonFEMproject2xSer_LyNonPeriodic_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 
  double *femFldLy = femFld[3]; 
  double *femFldUy = femFld[4]; 
  double *femFldLxLy = femFld[5]; 
  double *femFldLxUy = femFld[6]; 
  double *femFldUxLy = femFld[7]; 
  double *femFldUxUy = femFld[8]; 

  out[0] = 0.1111111111111111*(8.0*femFldUy[0]+2.0*femFldUxUy[0]+4.0*femFldUx[0]+2.0*femFldLxUy[0]+4.0*femFldLx[0]+16.0*femFldC[0])*volFac; 

}

void MGpoissonFEMproject2xSer_UyNonPeriodic_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 
  double *femFldLy = femFld[3]; 
  double *femFldUy = femFld[4]; 
  double *femFldLxLy = femFld[5]; 
  double *femFldLxUy = femFld[6]; 
  double *femFldUxLy = femFld[7]; 
  double *femFldUxUy = femFld[8]; 

  out[0] = 0.1111111111111111*(2.0*(femFldUx[1]+femFldLx[1])+8.0*femFldC[1]+2.0*femFldUxLy[0]+8.0*(femFldUx[0]+femFldLy[0])+2.0*femFldLxLy[0]+8.0*femFldLx[0]+32.0*femFldC[0])*volFac; 
  out[1] = 0.1111111111111111*(4.0*(femFldUx[1]+femFldLx[1])+16.0*femFldC[1]+2.0*(femFldUx[0]+femFldLx[0])+8.0*femFldC[0])*volFac; 

}

void MGpoissonFEMproject2xSer_LxNonPeriodicLyNonPeriodic_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 
  double *femFldLy = femFld[3]; 
  double *femFldUy = femFld[4]; 
  double *femFldLxLy = femFld[5]; 
  double *femFldLxUy = femFld[6]; 
  double *femFldUxLy = femFld[7]; 
  double *femFldUxUy = femFld[8]; 

  out[0] = 0.1111111111111111*(4.0*femFldUy[0]+2.0*femFldUxUy[0]+4.0*femFldUx[0]+8.0*femFldC[0])*volFac; 

}

void MGpoissonFEMproject2xSer_LxNonPeriodicUyNonPeriodic_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 
  double *femFldLy = femFld[3]; 
  double *femFldUy = femFld[4]; 
  double *femFldLxLy = femFld[5]; 
  double *femFldLxUy = femFld[6]; 
  double *femFldUxLy = femFld[7]; 
  double *femFldUxUy = femFld[8]; 

  out[0] = 0.1111111111111111*(2.0*femFldUx[1]+4.0*femFldC[1]+2.0*femFldUxLy[0]+8.0*femFldUx[0]+4.0*femFldLy[0]+16.0*femFldC[0])*volFac; 
  out[1] = 0.1111111111111111*(4.0*femFldUx[1]+8.0*femFldC[1]+2.0*femFldUx[0]+4.0*femFldC[0])*volFac; 

}

void MGpoissonFEMproject2xSer_UxNonPeriodicLyNonPeriodic_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 
  double *femFldLy = femFld[3]; 
  double *femFldUy = femFld[4]; 
  double *femFldLxLy = femFld[5]; 
  double *femFldLxUy = femFld[6]; 
  double *femFldUxLy = femFld[7]; 
  double *femFldUxUy = femFld[8]; 

  out[0] = 0.1111111111111111*(2.0*femFldUy[1]+4.0*femFldC[1]+8.0*femFldUy[0]+2.0*femFldLxUy[0]+4.0*femFldLx[0]+16.0*femFldC[0])*volFac; 
  out[1] = 0.1111111111111111*(4.0*femFldUy[1]+8.0*femFldC[1]+2.0*femFldUy[0]+4.0*femFldC[0])*volFac; 

}

void MGpoissonFEMproject2xSer_UxNonPeriodicUyNonPeriodic_P1(double **dx, double **femFld, double *out) 
{ 
  // dx:      cell lengths of cells pointed to by the projection stencil.
  // femFld:  FEM field in cells pointed to by the projection stencil.
  // out:     projection of the FEM field.

  double *dxC  = dx[0]; 

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldLx = femFld[1]; 
  double *femFldUx = femFld[2]; 
  double *femFldLy = femFld[3]; 
  double *femFldUy = femFld[4]; 
  double *femFldLxLy = femFld[5]; 
  double *femFldLxUy = femFld[6]; 
  double *femFldUxLy = femFld[7]; 
  double *femFldUxUy = femFld[8]; 

  out[0] = 0.1111111111111111*(2.0*femFldC[3]+8.0*femFldC[2]+2.0*(femFldLy[1]+femFldLx[1])+8.0*(femFldC[1]+femFldLy[0])+2.0*femFldLxLy[0]+8.0*femFldLx[0]+32.0*femFldC[0])*volFac; 
  out[1] = 0.1111111111111111*(4.0*femFldC[3]+2.0*femFldC[2]+4.0*femFldLy[1]+16.0*femFldC[1]+2.0*femFldLy[0]+8.0*femFldC[0])*volFac; 
  out[2] = 0.1111111111111111*(4.0*femFldC[3]+16.0*femFldC[2]+4.0*femFldLx[1]+2.0*(femFldC[1]+femFldLx[0])+8.0*femFldC[0])*volFac; 
  out[3] = 0.1111111111111111*(8.0*femFldC[3]+4.0*(femFldC[2]+femFldC[1])+2.0*femFldC[0])*volFac; 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = (((4.0*phiUy[0]+phiUxUy[0]+phiUxLy[0]-2.0*phiUx[0]-8.0*phiPrevC[0]+4.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-2.0*phiLx[0])*rdx2SqVol[1]+6.0*rhoC[0]+((-2.0*phiUy[0])+phiUxUy[0]+phiUxLy[0]+4.0*phiUx[0]-8.0*phiPrevC[0]-2.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]+4.0*phiLx[0])*rdx2SqVol[0])*omega+8.0*phiPrevC[0]*rdx2SqVol[1]+8.0*phiPrevC[0]*rdx2SqVol[0])/(8.0*rdx2SqVol[1]+8.0*rdx2SqVol[0]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]-1.0*phiUxLy[0]+2.0*phiUx[0]+4.0*phiPrevC[0]-2.0*phiLy[0])*rdx2SqVol[1]-6.0*rhoC[0]+(phiUy[0]-1.0*phiUxUy[0]-1.0*phiUxLy[0]-4.0*phiUx[0]+4.0*phiPrevC[0]+phiLy[0])*rdx2SqVol[0])*omega-4.0*phiPrevC[0]*rdx2SqVol[1]-4.0*phiPrevC[0]*rdx2SqVol[0]))/(4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]-1.0*phiUxLy[0]+2.0*phiUx[0]+4.0*phiPrevC[0]-2.0*phiLy[0])*bcVals[1]*rdx2SqVol[1]+((phiUy[0]-1.0*phiUxUy[0]-1.0*phiUxLy[0]-4.0*phiUx[0]+4.0*phiPrevC[0]+phiLy[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiPrevC[0])*rdx2SqVol[0])*omega-4.0*phiPrevC[0]*bcVals[1]*rdx2SqVol[1]-4.0*phiPrevC[0]*rdx2SqVol[0]*bcVals[1]+4.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0]))/(4.0*bcVals[1]*rdx2SqVol[1]+4.0*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*rdx2SqVol[0]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*(((rdx2SqVol[1]-5.0*rdx2SqVol[0])*bcVals[5]+((-1.0*phiLy[1])-4.0*phiUy[0]+8.0*phiPrevC[0]-4.0*phiLy[0]-1.0*phiLxUy[0]-1.0*phiLxLy[0]+2.0*phiLx[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-6.0*rhoC[0]+(2.0*phiUy[0]+8.0*phiPrevC[0]+2.0*phiLy[0]-1.0*phiLxUy[0]-1.0*phiLxLy[0]-4.0*phiLx[0])*rdx2SqVol[0])*omega-8.0*phiPrevC[0]*rdx2SqVol[1]-8.0*phiPrevC[0]*rdx2SqVol[0]))/(8.0*rdx2SqVol[1]+8.0*rdx2SqVol[0]); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = -(1.0*(((6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+(6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])+14.0*phiPrevC[0]-7.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]+4.0*phiLx[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-9.0*rdx2SqVol[0]*phiLy[1]-12.0*rhoC[0]+((-5.0*phiUy[0])+40.0*phiPrevC[0]-5.0*phiLy[0]-4.0*phiLxUy[0]-4.0*phiLxLy[0]-4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]+8.0*phiPrevC[0]+2.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol0R2)*omega-14.0*phiPrevC[0]*rdx2SqVol1R2-40.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]-8.0*phiPrevC[0]*rdx2SqVol0R2))/(14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2); 
  phiC[1] = (((24.0*rdx2SqVol[0]*rdx2SqVol[1]+24.0*rdx2SqVol0R2)*bcVals[5]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+(7.0*phiUy[1]-14.0*phiPrevC[1]+7.0*phiLy[1]-1.0*phiLxUy[0]-1.0*phiLxLy[0]+2.0*phiLx[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUy[1]-40.0*rdx2SqVol[0]*phiPrevC[1]+5.0*rdx2SqVol[0]*phiLy[1]-6.0*rhoC[0]+(18.0*phiUy[0]+18.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-2.0*rdx2SqVol0R2*phiUy[1]-8.0*rdx2SqVol0R2*phiPrevC[1]-2.0*rdx2SqVol0R2*phiLy[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLxUy[0]+2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*omega+14.0*phiPrevC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiPrevC[1])/(14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = -(1.0*(((6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+((6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])+14.0*phiPrevC[0]-7.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]+4.0*phiLx[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-9.0*rdx2SqVol[0]*phiLy[1]-12.0*rhoC[0]+((-5.0*phiUy[0])+40.0*phiPrevC[0]-5.0*phiLy[0]-4.0*phiLxUy[0]-4.0*phiLxLy[0]-4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]+8.0*phiPrevC[0]+2.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+(((-4.0*rdx2SqVol[0]*phiUy[1])-2.0*rdx2SqVol[0]*phiLy[1]+((-8.0*phiUy[0])+16.0*phiPrevC[0]-8.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiUy[1]-2.0*rdx2SqVol0R2*phiLy[1]-12.0*rdx2SqVol[0]*rhoC[0]+(4.0*phiUy[0]+16.0*phiPrevC[0]+4.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol0R2)*bcVals[3])*omega+((-14.0*phiPrevC[0]*rdx2SqVol1R2)-40.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]-8.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[4]+((-16.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1])-16.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[3]))/((14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 
  phiC[1] = (((24.0*rdx2SqVol[0]*rdx2SqVol[1]+24.0*rdx2SqVol0R2)*bcVals[5]+((24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+(7.0*phiUy[1]-14.0*phiPrevC[1]+7.0*phiLy[1]-1.0*phiLxUy[0]-1.0*phiLxLy[0]+2.0*phiLx[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUy[1]-40.0*rdx2SqVol[0]*phiPrevC[1]+5.0*rdx2SqVol[0]*phiLy[1]-6.0*rhoC[0]+(18.0*phiUy[0]+18.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-2.0*rdx2SqVol0R2*phiUy[1]-8.0*rdx2SqVol0R2*phiPrevC[1]-2.0*rdx2SqVol0R2*phiLy[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLxUy[0]+2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+(((-8.0*rdx2SqVol[0]*phiUy[1])-16.0*rdx2SqVol[0]*phiPrevC[1])*rdx2SqVol[1]-8.0*rdx2SqVol0R2*phiUy[1]-16.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[3])*omega+(14.0*phiPrevC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[4]+(16.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+16.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[3])/((14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]+((-4.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+4.0*phiPrevC[0]-1.0*phiLxUy[0]+phiLx[0])*rdx2SqVol[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+4.0*phiPrevC[0]-1.0*phiLxUy[0]-2.0*phiLx[0])*rdx2SqVol[0])*omega-4.0*phiPrevC[0]*rdx2SqVol[1]-4.0*phiPrevC[0]*rdx2SqVol[0]))/(4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]+(((-4.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+4.0*phiPrevC[0]-1.0*phiLxUy[0]+phiLx[0])*rdx2SqVol[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+4.0*phiPrevC[0]-1.0*phiLxUy[0]-2.0*phiLx[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiPrevC[0])*rdx2SqVol[1]*bcVals[6])*omega+((-4.0*phiPrevC[0]*rdx2SqVol[1])-4.0*phiPrevC[0]*rdx2SqVol[0])*bcVals[7]+4.0*phiPrevC[0]*rdx2SqVol[1]*bcVals[6]))/((4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0])*bcVals[7]-4.0*rdx2SqVol[1]*bcVals[6]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = (((5.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*bcVals[11]+(phiLx[1]+phiUxLy[0]-2.0*phiUx[0]-8.0*phiPrevC[0]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLx[1]+6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-8.0*phiPrevC[0]-2.0*phiLy[0]+phiLxLy[0]+4.0*phiLx[0])*rdx2SqVol[0])*omega+8.0*phiPrevC[0]*rdx2SqVol[1]+8.0*phiPrevC[0]*rdx2SqVol[0])/(8.0*rdx2SqVol[1]+8.0*rdx2SqVol[0]); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]-8.0*phiPrevC[0]+8.0*phiLy[0]+2.0*phiLxLy[0]-2.0*phiLx[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiUx[1]+9.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]-40.0*phiPrevC[0]+4.0*phiLy[0]+4.0*phiLxLy[0]+5.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-14.0*phiPrevC[0]-4.0*phiLy[0]+2.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol0R2)*omega+8.0*phiPrevC[0]*rdx2SqVol1R2+40.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+14.0*phiPrevC[0]*rdx2SqVol0R2)/(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2); 
  phiC[1] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+((-2.0*phiUx[1])-8.0*phiPrevC[1]-2.0*phiLx[1]+2.0*phiUxLy[0]+8.0*phiLy[0]+2.0*phiLxLy[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUx[1]-40.0*rdx2SqVol[0]*phiPrevC[1]+5.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-8.0*phiLy[0]+phiLxLy[0]+18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiUx[1]-14.0*rdx2SqVol0R2*phiPrevC[1]+7.0*rdx2SqVol0R2*phiLx[1]-6.0*rdx2SqVol[0]*rhoC[0]+((-1.0*phiUxLy[0])+2.0*phiLy[0]-1.0*phiLxLy[0])*rdx2SqVol0R2)*omega+8.0*phiPrevC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+14.0*rdx2SqVol0R2*phiPrevC[1])/(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]-8.0*phiPrevC[0]+8.0*phiLy[0]+2.0*phiLxLy[0]-2.0*phiLx[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiUx[1]+9.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]-40.0*phiPrevC[0]+4.0*phiLy[0]+4.0*phiLxLy[0]+5.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-14.0*phiPrevC[0]-4.0*phiLy[0]+2.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol0R2)*bcVals[10]+(((-2.0*phiUx[1])+2.0*phiLx[1]+2.0*phiUxLy[0]-4.0*phiUx[0]-16.0*phiPrevC[0]+8.0*phiLy[0]+2.0*phiLxLy[0]-4.0*phiLx[0])*rdx2SqVol1R2+(4.0*rdx2SqVol[0]*phiUx[1]+2.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(2.0*phiUxLy[0]+8.0*phiUx[0]-16.0*phiPrevC[0]-4.0*phiLy[0]+2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[9])*omega+(8.0*phiPrevC[0]*rdx2SqVol1R2+40.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+14.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[10]+(16.0*phiPrevC[0]*rdx2SqVol1R2+16.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9])/((8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 
  phiC[1] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+((-2.0*phiUx[1])-8.0*phiPrevC[1]-2.0*phiLx[1]+2.0*phiUxLy[0]+8.0*phiLy[0]+2.0*phiLxLy[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUx[1]-40.0*rdx2SqVol[0]*phiPrevC[1]+5.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-8.0*phiLy[0]+phiLxLy[0]+18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiUx[1]-14.0*rdx2SqVol0R2*phiPrevC[1]+7.0*rdx2SqVol0R2*phiLx[1]-6.0*rdx2SqVol[0]*rhoC[0]+((-1.0*phiUxLy[0])+2.0*phiLy[0]-1.0*phiLxLy[0])*rdx2SqVol0R2)*bcVals[10]+(((-8.0*phiUx[1])-16.0*phiPrevC[1])*rdx2SqVol1R2+((-8.0*rdx2SqVol[0]*phiUx[1])-16.0*rdx2SqVol[0]*phiPrevC[1])*rdx2SqVol[1])*bcVals[9])*omega+(8.0*phiPrevC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+14.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[10]+(16.0*phiPrevC[1]*rdx2SqVol1R2+16.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1])*bcVals[9])/((8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]+6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiPrevC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiPrevC[0])*rdx2SqVol[0])*omega-2.0*phiPrevC[0]*rdx2SqVol[1]-2.0*phiPrevC[0]*rdx2SqVol[0]))/(2.0*rdx2SqVol[1]+2.0*rdx2SqVol[0]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]+(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiPrevC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiPrevC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiPrevC[0])*rdx2SqVol[1]*bcVals[6])*omega+((-2.0*phiPrevC[0]*rdx2SqVol[1])-2.0*phiPrevC[0]*rdx2SqVol[0])*bcVals[7]+4.0*phiPrevC[0]*rdx2SqVol[1]*bcVals[6]))/((2.0*rdx2SqVol[1]+2.0*rdx2SqVol[0])*bcVals[7]-4.0*rdx2SqVol[1]*bcVals[6]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*((6.0*bcVals[1]*rdx2SqVol[1]*bcVals[8]+6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiPrevC[0])*bcVals[1]*rdx2SqVol[1]+((phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiPrevC[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiPrevC[0])*rdx2SqVol[0])*omega-2.0*phiPrevC[0]*bcVals[1]*rdx2SqVol[1]-2.0*phiPrevC[0]*rdx2SqVol[0]*bcVals[1]+4.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0]))/(2.0*bcVals[1]*rdx2SqVol[1]+2.0*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*rdx2SqVol[0]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*((6.0*bcVals[1]*rdx2SqVol[1]*bcVals[8]+(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiPrevC[0])*bcVals[1]*rdx2SqVol[1]+((phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiPrevC[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiPrevC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiPrevC[0])*bcVals[1]*rdx2SqVol[1]*bcVals[6])*omega+((-2.0*phiPrevC[0]*bcVals[1]*rdx2SqVol[1])-2.0*phiPrevC[0]*rdx2SqVol[0]*bcVals[1]+4.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0])*bcVals[7]+4.0*phiPrevC[0]*bcVals[1]*rdx2SqVol[1]*bcVals[6]))/((2.0*bcVals[1]*rdx2SqVol[1]+2.0*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*rdx2SqVol[0])*bcVals[7]-4.0*bcVals[1]*rdx2SqVol[1]*bcVals[6]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = ((3.0*rdx2SqVol[1]*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+(phiUxLy[0]-2.0*phiUx[0]-4.0*phiPrevC[0]+2.0*phiLy[0])*rdx2SqVol[1]+6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-4.0*phiPrevC[0]-1.0*phiLy[0])*rdx2SqVol[0])*omega+4.0*phiPrevC[0]*rdx2SqVol[1]+4.0*phiPrevC[0]*rdx2SqVol[0])/(4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0]); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-6.0*rdx2SqVol0R2)*bcVals[2]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]-4.0*phiPrevC[0]+4.0*phiLy[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiUx[1]+12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]-20.0*phiPrevC[0]+2.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-7.0*phiPrevC[0]-2.0*phiLy[0])*rdx2SqVol0R2)*omega+4.0*phiPrevC[0]*rdx2SqVol1R2+20.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+7.0*phiPrevC[0]*rdx2SqVol0R2)/(4.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2); 
  phiC[1] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((-36.0*rdx2SqVol[0]*rdx2SqVol[1])-18.0*rdx2SqVol0R2)*bcVals[2]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+((-2.0*phiUx[1])-4.0*phiPrevC[1]+2.0*phiUxLy[0]+4.0*phiLy[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUx[1]-20.0*rdx2SqVol[0]*phiPrevC[1]+12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-4.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiUx[1]-7.0*rdx2SqVol0R2*phiPrevC[1]-6.0*rdx2SqVol[0]*rhoC[0]+(phiLy[0]-1.0*phiUxLy[0])*rdx2SqVol0R2)*omega+4.0*phiPrevC[1]*rdx2SqVol1R2+20.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiPrevC[1])/(4.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-6.0*rdx2SqVol0R2)*bcVals[2]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]-4.0*phiPrevC[0]+4.0*phiLy[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiUx[1]+12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]-20.0*phiPrevC[0]+2.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-7.0*phiPrevC[0]-2.0*phiLy[0])*rdx2SqVol0R2)*bcVals[10]+((-24.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[2])+(4.0*phiUxLy[0]-8.0*phiUx[0]-16.0*phiPrevC[0]+8.0*phiLy[0])*rdx2SqVol1R2+(6.0*rdx2SqVol[0]*phiUx[1]+24.0*rhoC[0]+(4.0*phiUxLy[0]+16.0*phiUx[0]-16.0*phiPrevC[0]-4.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[9])*omega+(4.0*phiPrevC[0]*rdx2SqVol1R2+20.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+7.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[10]+(16.0*phiPrevC[0]*rdx2SqVol1R2+16.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9])/((4.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 
  phiC[1] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(((-36.0*rdx2SqVol[0]*rdx2SqVol[1])-18.0*rdx2SqVol0R2)*bcVals[2]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+((-2.0*phiUx[1])-4.0*phiPrevC[1]+2.0*phiUxLy[0]+4.0*phiLy[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUx[1]-20.0*rdx2SqVol[0]*phiPrevC[1]+12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-4.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiUx[1]-7.0*rdx2SqVol0R2*phiPrevC[1]-6.0*rdx2SqVol[0]*rhoC[0]+(phiLy[0]-1.0*phiUxLy[0])*rdx2SqVol0R2)*bcVals[10]+(((-8.0*phiUx[1])-16.0*phiPrevC[1])*rdx2SqVol1R2+((-8.0*rdx2SqVol[0]*phiUx[1])-16.0*rdx2SqVol[0]*phiPrevC[1])*rdx2SqVol[1])*bcVals[9])*omega+(4.0*phiPrevC[1]*rdx2SqVol1R2+20.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[10]+(16.0*phiPrevC[1]*rdx2SqVol1R2+16.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1])*bcVals[9])/((4.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = (((3.0*bcVals[1]*rdx2SqVol[1]+2.0*bcVals[0]*rdx2SqVol[0])*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+(phiUxLy[0]-2.0*phiUx[0]-4.0*phiPrevC[0]+2.0*phiLy[0])*bcVals[1]*rdx2SqVol[1]+(6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-4.0*phiPrevC[0]-1.0*phiLy[0])*rdx2SqVol[0])*bcVals[1]+4.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0])*omega+4.0*phiPrevC[0]*bcVals[1]*rdx2SqVol[1]+4.0*phiPrevC[0]*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0])/(4.0*bcVals[1]*rdx2SqVol[1]+4.0*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*rdx2SqVol[0]); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);

  phiC[0] = (((12.0*bcVals1R2*rdx2SqVol1R2+(12.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-6.0*rdx2SqVol[0]*bcVals1R2)*rdx2SqVol[1])*bcVals[11]+((-24.0*rdx2SqVol[0]*bcVals[1]*rdx2SqVol[1])-6.0*rdx2SqVol0R2*bcVals[1]+12.0*bcVals[0]*rdx2SqVol0R2)*bcVals[2]+(12.0*bcVals1R2*rdx2SqVol[1]-6.0*rdx2SqVol[0]*bcVals1R2+12.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]-4.0*phiPrevC[0]+4.0*phiLy[0])*bcVals1R2*rdx2SqVol1R2+((9.0*rdx2SqVol[0]*bcVals1R2-6.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]-20.0*phiPrevC[0]+2.0*phiLy[0])*rdx2SqVol[0])*bcVals1R2+((-4.0*bcVals[0]*phiUxLy[0])+10.0*bcVals[0]*phiUx[0]+32.0*bcVals[0]*phiPrevC[0]-8.0*bcVals[0]*phiLy[0])*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+(12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-7.0*phiPrevC[0]-2.0*phiLy[0])*rdx2SqVol0R2)*bcVals1R2+(((-4.0*bcVals[0]*phiUxLy[0])-14.0*bcVals[0]*phiUx[0]+20.0*bcVals[0]*phiPrevC[0]+4.0*bcVals[0]*phiLy[0])*rdx2SqVol0R2-24.0*bcVals[0]*rdx2SqVol[0]*rhoC[0])*bcVals[1]-12.0*bcVals0R2*phiPrevC[0]*rdx2SqVol0R2)*omega+4.0*phiPrevC[0]*bcVals1R2*rdx2SqVol1R2+(20.0*phiPrevC[0]*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*phiPrevC[0]*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*phiPrevC[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*phiPrevC[0]*rdx2SqVol0R2)/(4.0*bcVals1R2*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2); 
  phiC[1] = (((24.0*bcVals1R2*rdx2SqVol1R2+(24.0*rdx2SqVol[0]*bcVals1R2-24.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[11]+((-36.0*rdx2SqVol[0]*bcVals[1]*rdx2SqVol[1])-18.0*rdx2SqVol0R2*bcVals[1]+12.0*bcVals[0]*rdx2SqVol0R2)*bcVals[2]+(24.0*bcVals1R2*rdx2SqVol[1]+24.0*rdx2SqVol[0]*bcVals1R2-24.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rhoC[1]+((-2.0*bcVals1R2*phiUx[1])-4.0*bcVals1R2*phiPrevC[1]+(2.0*phiUxLy[0]+4.0*phiLy[0])*bcVals1R2)*rdx2SqVol1R2+((5.0*rdx2SqVol[0]*bcVals1R2+6.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-20.0*rdx2SqVol[0]*bcVals1R2)*phiPrevC[1]+(12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-4.0*phiLy[0])*rdx2SqVol[0])*bcVals1R2+(2.0*bcVals[0]*phiUxLy[0]-8.0*bcVals[0]*phiUx[0]+4.0*bcVals[0]*phiLy[0])*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+(7.0*rdx2SqVol0R2*bcVals1R2-6.0*bcVals[0]*rdx2SqVol0R2*bcVals[1])*phiUx[1]+((-7.0*rdx2SqVol0R2*bcVals1R2)+20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]-12.0*bcVals0R2*rdx2SqVol0R2)*phiPrevC[1]+((phiLy[0]-1.0*phiUxLy[0])*rdx2SqVol0R2-6.0*rdx2SqVol[0]*rhoC[0])*bcVals1R2+(12.0*bcVals[0]*rdx2SqVol[0]*rhoC[0]+(2.0*bcVals[0]*phiUxLy[0]+4.0*bcVals[0]*phiUx[0]-2.0*bcVals[0]*phiLy[0])*rdx2SqVol0R2)*bcVals[1])*omega+4.0*bcVals1R2*phiPrevC[1]*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiPrevC[1]*rdx2SqVol[1]+(7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2)*phiPrevC[1])/(4.0*bcVals1R2*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);

  phiC[0] = (((12.0*bcVals1R2*rdx2SqVol1R2+(12.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-6.0*rdx2SqVol[0]*bcVals1R2)*rdx2SqVol[1])*bcVals[11]+(((-24.0*rdx2SqVol[0]*bcVals[1]*rdx2SqVol[1])-6.0*rdx2SqVol0R2*bcVals[1]+12.0*bcVals[0]*rdx2SqVol0R2)*bcVals[2]+(12.0*bcVals1R2*rdx2SqVol[1]-6.0*rdx2SqVol[0]*bcVals1R2+12.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]-4.0*phiPrevC[0]+4.0*phiLy[0])*bcVals1R2*rdx2SqVol1R2+((9.0*rdx2SqVol[0]*bcVals1R2-6.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]-20.0*phiPrevC[0]+2.0*phiLy[0])*rdx2SqVol[0])*bcVals1R2+((-4.0*bcVals[0]*phiUxLy[0])+10.0*bcVals[0]*phiUx[0]+32.0*bcVals[0]*phiPrevC[0]-8.0*bcVals[0]*phiLy[0])*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+(12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-7.0*phiPrevC[0]-2.0*phiLy[0])*rdx2SqVol0R2)*bcVals1R2+(((-4.0*bcVals[0]*phiUxLy[0])-14.0*bcVals[0]*phiUx[0]+20.0*bcVals[0]*phiPrevC[0]+4.0*bcVals[0]*phiLy[0])*rdx2SqVol0R2-24.0*bcVals[0]*rdx2SqVol[0]*rhoC[0])*bcVals[1]-12.0*bcVals0R2*phiPrevC[0]*rdx2SqVol0R2)*bcVals[10]+((-24.0*rdx2SqVol[0]*bcVals[1]*rdx2SqVol[1]*bcVals[2])+(4.0*phiUxLy[0]-8.0*phiUx[0]-16.0*phiPrevC[0]+8.0*phiLy[0])*bcVals1R2*rdx2SqVol1R2+((6.0*rdx2SqVol[0]*bcVals1R2-4.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(24.0*rhoC[0]+(4.0*phiUxLy[0]+16.0*phiUx[0]-16.0*phiPrevC[0]-4.0*phiLy[0])*rdx2SqVol[0])*bcVals1R2+16.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[9])*omega+(4.0*phiPrevC[0]*bcVals1R2*rdx2SqVol1R2+(20.0*phiPrevC[0]*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*phiPrevC[0]*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*phiPrevC[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*phiPrevC[0]*rdx2SqVol0R2)*bcVals[10]+(16.0*phiPrevC[0]*bcVals1R2*rdx2SqVol1R2+(16.0*phiPrevC[0]*rdx2SqVol[0]*bcVals1R2-16.0*bcVals[0]*phiPrevC[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[9])/((4.0*bcVals1R2*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2)*bcVals[10]+(16.0*bcVals1R2*rdx2SqVol1R2+(16.0*rdx2SqVol[0]*bcVals1R2-16.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[9]); 
  phiC[1] = (((24.0*bcVals1R2*rdx2SqVol1R2+(24.0*rdx2SqVol[0]*bcVals1R2-24.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[11]+(((-36.0*rdx2SqVol[0]*bcVals[1]*rdx2SqVol[1])-18.0*rdx2SqVol0R2*bcVals[1]+12.0*bcVals[0]*rdx2SqVol0R2)*bcVals[2]+(24.0*bcVals1R2*rdx2SqVol[1]+24.0*rdx2SqVol[0]*bcVals1R2-24.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rhoC[1]+((-2.0*bcVals1R2*phiUx[1])-4.0*bcVals1R2*phiPrevC[1]+(2.0*phiUxLy[0]+4.0*phiLy[0])*bcVals1R2)*rdx2SqVol1R2+((5.0*rdx2SqVol[0]*bcVals1R2+6.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-20.0*rdx2SqVol[0]*bcVals1R2)*phiPrevC[1]+(12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-4.0*phiLy[0])*rdx2SqVol[0])*bcVals1R2+(2.0*bcVals[0]*phiUxLy[0]-8.0*bcVals[0]*phiUx[0]+4.0*bcVals[0]*phiLy[0])*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+(7.0*rdx2SqVol0R2*bcVals1R2-6.0*bcVals[0]*rdx2SqVol0R2*bcVals[1])*phiUx[1]+((-7.0*rdx2SqVol0R2*bcVals1R2)+20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]-12.0*bcVals0R2*rdx2SqVol0R2)*phiPrevC[1]+((phiLy[0]-1.0*phiUxLy[0])*rdx2SqVol0R2-6.0*rdx2SqVol[0]*rhoC[0])*bcVals1R2+(12.0*bcVals[0]*rdx2SqVol[0]*rhoC[0]+(2.0*bcVals[0]*phiUxLy[0]+4.0*bcVals[0]*phiUx[0]-2.0*bcVals[0]*phiLy[0])*rdx2SqVol0R2)*bcVals[1])*bcVals[10]+(((-8.0*bcVals1R2*phiUx[1])-16.0*bcVals1R2*phiPrevC[1])*rdx2SqVol1R2+((8.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-8.0*rdx2SqVol[0]*bcVals1R2)*phiUx[1]+(16.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-16.0*rdx2SqVol[0]*bcVals1R2)*phiPrevC[1])*rdx2SqVol[1])*bcVals[9])*omega+(4.0*bcVals1R2*phiPrevC[1]*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiPrevC[1]*rdx2SqVol[1]+(7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2)*phiPrevC[1])*bcVals[10]+(16.0*bcVals1R2*phiPrevC[1]*rdx2SqVol1R2+(16.0*rdx2SqVol[0]*bcVals1R2-16.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiPrevC[1]*rdx2SqVol[1])*bcVals[9])/((4.0*bcVals1R2*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2)*bcVals[10]+(16.0*bcVals1R2*rdx2SqVol1R2+(16.0*rdx2SqVol[0]*bcVals1R2-16.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[9]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]-3.0*rdx2SqVol[0]*bcVals[5]+((-4.0*phiUy[0])+4.0*phiPrevC[0]-1.0*phiLxUy[0]+phiLx[0])*rdx2SqVol[1]-6.0*rhoC[0]+(2.0*phiUy[0]+4.0*phiPrevC[0]-1.0*phiLxUy[0]-2.0*phiLx[0])*rdx2SqVol[0])*omega-4.0*phiPrevC[0]*rdx2SqVol[1]-4.0*phiPrevC[0]*rdx2SqVol[0]))/(4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0]); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]+((-3.0*rdx2SqVol[0]*bcVals[5])+((-4.0*phiUy[0])+4.0*phiPrevC[0]-1.0*phiLxUy[0]+phiLx[0])*rdx2SqVol[1]-6.0*rhoC[0]+(2.0*phiUy[0]+4.0*phiPrevC[0]-1.0*phiLxUy[0]-2.0*phiLx[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*rdx2SqVol[1]*bcVals[5])-4.0*phiPrevC[0]*rdx2SqVol[1])*bcVals[6])*omega+((-4.0*phiPrevC[0]*rdx2SqVol[1])-4.0*phiPrevC[0]*rdx2SqVol[0])*bcVals[7]+4.0*phiPrevC[0]*rdx2SqVol[1]*bcVals[6]))/((4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0])*bcVals[7]-4.0*rdx2SqVol[1]*bcVals[6]); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = -(1.0*(((6.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[8]+(6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+(6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])+7.0*phiPrevC[0]-2.0*phiLxUy[0]+2.0*phiLx[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-12.0*rhoC[0]+((-5.0*phiUy[0])+20.0*phiPrevC[0]-4.0*phiLxUy[0]-2.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]+4.0*phiPrevC[0]-2.0*phiLxUy[0]-4.0*phiLx[0])*rdx2SqVol0R2)*omega-7.0*phiPrevC[0]*rdx2SqVol1R2-20.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]-4.0*phiPrevC[0]*rdx2SqVol0R2))/(7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2); 
  phiC[1] = -(1.0*(((18.0*rdx2SqVol1R2+36.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[8]+((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-24.0*rdx2SqVol0R2)*bcVals[5]+((-24.0*rdx2SqVol[1])-24.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[1])+7.0*phiPrevC[1]+phiLxUy[0]-1.0*phiLx[0])*rdx2SqVol1R2+((-5.0*rdx2SqVol[0]*phiUy[1])+20.0*rdx2SqVol[0]*phiPrevC[1]+6.0*rhoC[0]+((-18.0*phiUy[0])-1.0*phiLxUy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiUy[1]+4.0*rdx2SqVol0R2*phiPrevC[1]-12.0*rdx2SqVol[0]*rhoC[0]+((-2.0*phiLxUy[0])-4.0*phiLx[0])*rdx2SqVol0R2)*omega-7.0*phiPrevC[1]*rdx2SqVol1R2-20.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]-4.0*rdx2SqVol0R2*phiPrevC[1]))/(7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*((((6.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[7]-12.0*rdx2SqVol1R2*bcVals[6])*bcVals[8]+((6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+(6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])+7.0*phiPrevC[0]-2.0*phiLxUy[0]+2.0*phiLx[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-12.0*rhoC[0]+((-5.0*phiUy[0])+20.0*phiPrevC[0]-4.0*phiLxUy[0]-2.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]+4.0*phiPrevC[0]-2.0*phiLxUy[0]-4.0*phiLx[0])*rdx2SqVol0R2)*bcVals7R2+((-12.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[5])-12.0*rdx2SqVol[1]*rhoC[1]+(14.0*phiUy[0]-20.0*phiPrevC[0]+4.0*phiLxUy[0]-4.0*phiLx[0])*rdx2SqVol1R2+(6.0*rdx2SqVol[0]*phiUy[1]+24.0*rhoC[0]+((-10.0*phiUy[0])-32.0*phiPrevC[0]+4.0*phiLxUy[0]+8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[6]*bcVals[7]+12.0*phiPrevC[0]*rdx2SqVol1R2*bcVals6R2)*omega+((-7.0*phiPrevC[0]*rdx2SqVol1R2)-20.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]-4.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals7R2+(20.0*phiPrevC[0]*rdx2SqVol1R2+32.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[6]*bcVals[7]-12.0*phiPrevC[0]*rdx2SqVol1R2*bcVals6R2))/((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals7R2+((-20.0*rdx2SqVol1R2)-32.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[6]*bcVals[7]+12.0*rdx2SqVol1R2*bcVals6R2); 
  phiC[1] = -(1.0*((((18.0*rdx2SqVol1R2+36.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[7]-12.0*rdx2SqVol1R2*bcVals[6])*bcVals[8]+(((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-24.0*rdx2SqVol0R2)*bcVals[5]+((-24.0*rdx2SqVol[1])-24.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[1])+7.0*phiPrevC[1]+phiLxUy[0]-1.0*phiLx[0])*rdx2SqVol1R2+((-5.0*rdx2SqVol[0]*phiUy[1])+20.0*rdx2SqVol[0]*phiPrevC[1]+6.0*rhoC[0]+((-18.0*phiUy[0])-1.0*phiLxUy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiUy[1]+4.0*rdx2SqVol0R2*phiPrevC[1]-12.0*rdx2SqVol[0]*rhoC[0]+((-2.0*phiLxUy[0])-4.0*phiLx[0])*rdx2SqVol0R2)*bcVals7R2+(24.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[5]+24.0*rdx2SqVol[1]*rhoC[1]+(6.0*phiUy[1]-20.0*phiPrevC[1]-4.0*phiUy[0]-2.0*phiLxUy[0]+2.0*phiLx[0])*rdx2SqVol1R2+((-6.0*rdx2SqVol[0]*phiUy[1])-32.0*rdx2SqVol[0]*phiPrevC[1]-12.0*rhoC[0]+(8.0*phiUy[0]-2.0*phiLxUy[0]-4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[6]*bcVals[7]+12.0*phiPrevC[1]*rdx2SqVol1R2*bcVals6R2)*omega+((-7.0*phiPrevC[1]*rdx2SqVol1R2)-20.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]-4.0*rdx2SqVol0R2*phiPrevC[1])*bcVals7R2+(20.0*phiPrevC[1]*rdx2SqVol1R2+32.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1])*bcVals[6]*bcVals[7]-12.0*phiPrevC[1]*rdx2SqVol1R2*bcVals6R2))/((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals7R2+((-20.0*rdx2SqVol1R2)-32.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[6]*bcVals[7]+12.0*rdx2SqVol1R2*bcVals6R2); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = -(1.0*((((6.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]+24.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[3])*bcVals[8]+(6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+((6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])+7.0*phiPrevC[0]-2.0*phiLxUy[0]+2.0*phiLx[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-12.0*rhoC[0]+((-5.0*phiUy[0])+20.0*phiPrevC[0]-4.0*phiLxUy[0]-2.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]+4.0*phiPrevC[0]-2.0*phiLxUy[0]-4.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+((((-16.0*phiUy[0])+16.0*phiPrevC[0]-4.0*phiLxUy[0]+4.0*phiLx[0])*rdx2SqVol[0]-6.0*rdx2SqVol[0]*phiUy[1])*rdx2SqVol[1]-24.0*rdx2SqVol[0]*rhoC[0]+(8.0*phiUy[0]+16.0*phiPrevC[0]-4.0*phiLxUy[0]-8.0*phiLx[0])*rdx2SqVol0R2)*bcVals[3])*omega+((-7.0*phiPrevC[0]*rdx2SqVol1R2)-20.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]-4.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[4]+((-16.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1])-16.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[3]))/((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 
  phiC[1] = -(1.0*(((18.0*rdx2SqVol1R2+36.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]*bcVals[8]+((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-24.0*rdx2SqVol0R2)*bcVals[5]+(((-24.0*rdx2SqVol[1])-24.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[1])+7.0*phiPrevC[1]+phiLxUy[0]-1.0*phiLx[0])*rdx2SqVol1R2+((-5.0*rdx2SqVol[0]*phiUy[1])+20.0*rdx2SqVol[0]*phiPrevC[1]+6.0*rhoC[0]+((-18.0*phiUy[0])-1.0*phiLxUy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiUy[1]+4.0*rdx2SqVol0R2*phiPrevC[1]-12.0*rdx2SqVol[0]*rhoC[0]+((-2.0*phiLxUy[0])-4.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+((8.0*rdx2SqVol[0]*phiUy[1]+16.0*rdx2SqVol[0]*phiPrevC[1])*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiUy[1]+16.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[3])*omega+((-7.0*phiPrevC[1]*rdx2SqVol1R2)-20.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]-4.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[4]+((-16.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1])-16.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[3]))/((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*(((((6.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]+24.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[3])*bcVals[7]-12.0*rdx2SqVol1R2*bcVals[4]*bcVals[6])*bcVals[8]+((6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+((6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])+7.0*phiPrevC[0]-2.0*phiLxUy[0]+2.0*phiLx[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-12.0*rhoC[0]+((-5.0*phiUy[0])+20.0*phiPrevC[0]-4.0*phiLxUy[0]-2.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]+4.0*phiPrevC[0]-2.0*phiLxUy[0]-4.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+((((-16.0*phiUy[0])+16.0*phiPrevC[0]-4.0*phiLxUy[0]+4.0*phiLx[0])*rdx2SqVol[0]-6.0*rdx2SqVol[0]*phiUy[1])*rdx2SqVol[1]-24.0*rdx2SqVol[0]*rhoC[0]+(8.0*phiUy[0]+16.0*phiPrevC[0]-4.0*phiLxUy[0]-8.0*phiLx[0])*rdx2SqVol0R2)*bcVals[3])*bcVals7R2+((-12.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[5])+((-12.0*rdx2SqVol[1]*rhoC[1])+(14.0*phiUy[0]-20.0*phiPrevC[0]+4.0*phiLxUy[0]-4.0*phiLx[0])*rdx2SqVol1R2+(6.0*rdx2SqVol[0]*phiUy[1]+24.0*rhoC[0]+((-10.0*phiUy[0])-32.0*phiPrevC[0]+4.0*phiLxUy[0]+8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[4]+(4.0*rdx2SqVol[0]*phiUy[1]-16.0*phiPrevC[0]*rdx2SqVol[0])*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]+12.0*phiPrevC[0]*rdx2SqVol1R2*bcVals[4]*bcVals6R2)*omega+(((-7.0*phiPrevC[0]*rdx2SqVol1R2)-20.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]-4.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[4]+((-16.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1])-16.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[3])*bcVals7R2+((20.0*phiPrevC[0]*rdx2SqVol1R2+32.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]+16.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]-12.0*phiPrevC[0]*rdx2SqVol1R2*bcVals[4]*bcVals6R2))/(((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3])*bcVals7R2+(((-20.0*rdx2SqVol1R2)-32.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]-16.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]+12.0*rdx2SqVol1R2*bcVals[4]*bcVals6R2); 
  phiC[1] = -(1.0*((((18.0*rdx2SqVol1R2+36.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]*bcVals[7]-12.0*rdx2SqVol1R2*bcVals[4]*bcVals[6])*bcVals[8]+(((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-24.0*rdx2SqVol0R2)*bcVals[5]+(((-24.0*rdx2SqVol[1])-24.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[1])+7.0*phiPrevC[1]+phiLxUy[0]-1.0*phiLx[0])*rdx2SqVol1R2+((-5.0*rdx2SqVol[0]*phiUy[1])+20.0*rdx2SqVol[0]*phiPrevC[1]+6.0*rhoC[0]+((-18.0*phiUy[0])-1.0*phiLxUy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiUy[1]+4.0*rdx2SqVol0R2*phiPrevC[1]-12.0*rdx2SqVol[0]*rhoC[0]+((-2.0*phiLxUy[0])-4.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+((8.0*rdx2SqVol[0]*phiUy[1]+16.0*rdx2SqVol[0]*phiPrevC[1])*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiUy[1]+16.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[3])*bcVals7R2+(24.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[5]+(24.0*rdx2SqVol[1]*rhoC[1]+(6.0*phiUy[1]-20.0*phiPrevC[1]-4.0*phiUy[0]-2.0*phiLxUy[0]+2.0*phiLx[0])*rdx2SqVol1R2+((-6.0*rdx2SqVol[0]*phiUy[1])-32.0*rdx2SqVol[0]*phiPrevC[1]-12.0*rhoC[0]+(8.0*phiUy[0]-2.0*phiLxUy[0]-4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[4]+((-8.0*rdx2SqVol[0]*phiUy[1])-16.0*rdx2SqVol[0]*phiPrevC[1])*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]+12.0*phiPrevC[1]*rdx2SqVol1R2*bcVals[4]*bcVals6R2)*omega+(((-7.0*phiPrevC[1]*rdx2SqVol1R2)-20.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]-4.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[4]+((-16.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1])-16.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[3])*bcVals7R2+((20.0*phiPrevC[1]*rdx2SqVol1R2+32.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1])*bcVals[4]+16.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]-12.0*phiPrevC[1]*rdx2SqVol1R2*bcVals[4]*bcVals6R2))/(((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3])*bcVals7R2+(((-20.0*rdx2SqVol1R2)-32.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]-16.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]+12.0*rdx2SqVol1R2*bcVals[4]*bcVals6R2); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 


  phiC[0] = ((((5.0*dxC[1]+4.0*dxC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*dxC[1]-2.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+(((-2.0*dxC[1])-1.0*dxC[0])*rdx2SqVol[1]+4.0*rdx2SqVol[0]*dxC[1]+5.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+((dxC[1]+dxC[0])*phiLy[1]+(dxC[1]+dxC[0])*phiLx[1]+((-8.0*phiPrevC[0])+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0])*dxC[1]-8.0*dxC[0]*phiPrevC[0]+4.0*dxC[0]*phiLy[0]+dxC[0]*phiLxLy[0]-2.0*dxC[0]*phiLx[0])*rdx2SqVol[1]+(rdx2SqVol[0]*dxC[1]+dxC[0]*rdx2SqVol[0])*phiLy[1]+(rdx2SqVol[0]*dxC[1]+dxC[0]*rdx2SqVol[0])*phiLx[1]+(6.0*rhoC[0]+((-8.0*phiPrevC[0])-2.0*phiLy[0]+phiLxLy[0]+4.0*phiLx[0])*rdx2SqVol[0])*dxC[1]+6.0*dxC[0]*rhoC[0]+((-8.0*dxC[0]*phiPrevC[0])-2.0*dxC[0]*phiLy[0]+dxC[0]*phiLxLy[0]+4.0*dxC[0]*phiLx[0])*rdx2SqVol[0])*omega+(8.0*phiPrevC[0]*dxC[1]+8.0*dxC[0]*phiPrevC[0])*rdx2SqVol[1]+8.0*phiPrevC[0]*rdx2SqVol[0]*dxC[1]+8.0*dxC[0]*phiPrevC[0]*rdx2SqVol[0])/((8.0*dxC[1]+8.0*dxC[0])*rdx2SqVol[1]+8.0*rdx2SqVol[0]*dxC[1]+8.0*dxC[0]*rdx2SqVol[0]); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((-2.0*rdx2SqVol1R2)+14.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[5]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[2]+(2.0*phiLy[1]-8.0*phiPrevC[0]+8.0*phiLy[0]+2.0*phiLxLy[0]-2.0*phiLx[0])*rdx2SqVol1R2+(4.0*rdx2SqVol[0]*phiLy[1]+9.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+((-40.0*phiPrevC[0])+4.0*phiLy[0]+4.0*phiLxLy[0]+5.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiLy[1]+12.0*rdx2SqVol[0]*rhoC[0]+((-14.0*phiPrevC[0])-4.0*phiLy[0]+2.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol0R2)*omega+8.0*phiPrevC[0]*rdx2SqVol1R2+40.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+14.0*phiPrevC[0]*rdx2SqVol0R2)/(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2); 
  phiC[1] = (bcVals[5]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 
  phiC[2] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((-2.0*rdx2SqVol1R2)+23.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[5]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[2]+((-8.0*rdx2SqVol1R2)-40.0*rdx2SqVol[0]*rdx2SqVol[1]-14.0*rdx2SqVol0R2)*phiPrevC[2]+(2.0*phiLy[1]-2.0*phiLx[1]+8.0*phiLy[0]+2.0*phiLxLy[0])*rdx2SqVol1R2+(rdx2SqVol[0]*phiLy[1]+5.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+((-8.0*phiLy[0])+phiLxLy[0]+18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-1.0*rdx2SqVol0R2*phiLy[1]+7.0*rdx2SqVol0R2*phiLx[1]-6.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLy[0]-1.0*phiLxLy[0])*rdx2SqVol0R2)*omega+(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*phiPrevC[2])/(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(((-2.0*rdx2SqVol1R2)+14.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[5]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[2]+(2.0*phiLy[1]-8.0*phiPrevC[0]+8.0*phiLy[0]+2.0*phiLxLy[0]-2.0*phiLx[0])*rdx2SqVol1R2+(4.0*rdx2SqVol[0]*phiLy[1]+9.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+((-40.0*phiPrevC[0])+4.0*phiLy[0]+4.0*phiLxLy[0]+5.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiLy[1]+12.0*rdx2SqVol[0]*rhoC[0]+((-14.0*phiPrevC[0])-4.0*phiLy[0]+2.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol0R2)*bcVals[10]+((12.0*rdx2SqVol[0]*rdx2SqVol[1]-6.0*rdx2SqVol1R2)*bcVals[5]+(2.0*phiLy[1]+2.0*phiLx[1]-16.0*phiPrevC[0]+8.0*phiLy[0]+2.0*phiLxLy[0]-4.0*phiLx[0])*rdx2SqVol1R2+(2.0*rdx2SqVol[0]*phiLy[1]+2.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+((-16.0*phiPrevC[0])-4.0*phiLy[0]+2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[9])*omega+(8.0*phiPrevC[0]*rdx2SqVol1R2+40.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+14.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[10]+(16.0*phiPrevC[0]*rdx2SqVol1R2+16.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9])/((8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 
  phiC[1] = (bcVals[5]-1.0*phiPrevC[1])*omega+phiPrevC[1]; 
  phiC[2] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(((-2.0*rdx2SqVol1R2)+23.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[5]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[2]+((-8.0*rdx2SqVol1R2)-40.0*rdx2SqVol[0]*rdx2SqVol[1]-14.0*rdx2SqVol0R2)*phiPrevC[2]+(2.0*phiLy[1]-2.0*phiLx[1]+8.0*phiLy[0]+2.0*phiLxLy[0])*rdx2SqVol1R2+(rdx2SqVol[0]*phiLy[1]+5.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+((-8.0*phiLy[0])+phiLxLy[0]+18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-1.0*rdx2SqVol0R2*phiLy[1]+7.0*rdx2SqVol0R2*phiLx[1]-6.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLy[0]-1.0*phiLxLy[0])*rdx2SqVol0R2)*bcVals[10]+(((-8.0*rdx2SqVol1R2)-8.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[5]+((-16.0*rdx2SqVol1R2)-16.0*rdx2SqVol[0]*rdx2SqVol[1])*phiPrevC[2])*bcVals[9])*omega+(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*phiPrevC[2]*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*phiPrevC[2]*bcVals[9])/((8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((7.0*rdx2SqVol1R2+14.0*rdx2SqVol[0]*rdx2SqVol[1]-2.0*rdx2SqVol0R2)*bcVals[11]+(12.0*rdx2SqVol0R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[5]+(12.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*rhoC[1]+(2.0*phiLx[1]-14.0*phiPrevC[0]+7.0*phiLy[0]+2.0*phiLxLy[0]-4.0*phiLx[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiLy[1]+4.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+((-40.0*phiPrevC[0])+5.0*phiLy[0]+4.0*phiLxLy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiLx[1]+12.0*rdx2SqVol[0]*rhoC[0]+((-8.0*phiPrevC[0])-2.0*phiLy[0]+2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*omega+14.0*phiPrevC[0]*rdx2SqVol1R2+40.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+8.0*phiPrevC[0]*rdx2SqVol0R2)/(14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2); 
  phiC[1] = (((7.0*rdx2SqVol1R2+23.0*rdx2SqVol[0]*rdx2SqVol[1]-2.0*rdx2SqVol0R2)*bcVals[11]+(24.0*rdx2SqVol[0]*rdx2SqVol[1]+24.0*rdx2SqVol0R2)*bcVals[5]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+((-14.0*phiPrevC[1])+7.0*phiLy[1]-1.0*phiLx[1]-1.0*phiLxLy[0]+2.0*phiLx[0])*rdx2SqVol1R2+((-40.0*rdx2SqVol[0]*phiPrevC[1])+5.0*rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]-6.0*rhoC[0]+(18.0*phiLy[0]+phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-8.0*rdx2SqVol0R2*phiPrevC[1]-2.0*rdx2SqVol0R2*phiLy[1]+2.0*rdx2SqVol0R2*phiLx[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*omega+14.0*phiPrevC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiPrevC[1])/(14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol0R3 = std::pow(rdx2SqVol[0],3);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double rdx2SqVol1R3 = std::pow(rdx2SqVol[1],3);

  phiC[0] = (((42.0*rdx2SqVol1R3+822.0*rdx2SqVol[0]*rdx2SqVol1R2-84.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[11]+((-84.0*rdx2SqVol[0]*rdx2SqVol1R2)+822.0*rdx2SqVol0R2*rdx2SqVol[1]+42.0*rdx2SqVol0R3)*bcVals[5]+((-42.0*rdx2SqVol1R2)+564.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*rhoC[3]+(84.0*rdx2SqVol1R2+258.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*rhoC[2]+((-42.0*rdx2SqVol1R2)+258.0*rdx2SqVol[0]*rdx2SqVol[1]+84.0*rdx2SqVol0R2)*rhoC[1]+((-49.0*phiPrevC[0])+49.0*phiLy[0]+14.0*phiLxLy[0]-14.0*phiLx[0])*rdx2SqVol1R3+(189.0*rdx2SqVol[0]*phiLy[1]+81.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+((-651.0*phiPrevC[0])+336.0*phiLy[0]+96.0*phiLxLy[0]-51.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R2+(81.0*rdx2SqVol0R2*phiLy[1]+189.0*rdx2SqVol0R2*phiLx[1]+492.0*rdx2SqVol[0]*rhoC[0]+((-651.0*phiPrevC[0])-51.0*phiLy[0]+96.0*phiLxLy[0]+336.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol[1]+84.0*rdx2SqVol0R2*rhoC[0]+((-49.0*phiPrevC[0])-14.0*phiLy[0]+14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R3)*omega+49.0*phiPrevC[0]*rdx2SqVol1R3+651.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R2+651.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol[1]+49.0*phiPrevC[0]*rdx2SqVol0R3)/(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3); 
  phiC[1] = (((126.0*rdx2SqVol1R3+1080.0*rdx2SqVol[0]*rdx2SqVol1R2-126.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[11]+(336.0*rdx2SqVol[0]*rdx2SqVol1R2+1500.0*rdx2SqVol0R2*rdx2SqVol[1]+84.0*rdx2SqVol0R3)*bcVals[5]+(168.0*rdx2SqVol1R2+516.0*rdx2SqVol[0]*rdx2SqVol[1]-84.0*rdx2SqVol0R2)*rhoC[3]+((-42.0*rdx2SqVol1R2)+564.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*rhoC[2]+(168.0*rdx2SqVol1R2+984.0*rdx2SqVol[0]*rdx2SqVol[1]+168.0*rdx2SqVol0R2)*rhoC[1]+((-49.0*phiPrevC[1])+49.0*phiLy[1]-7.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol1R3+((-651.0*rdx2SqVol[0]*phiPrevC[1])+336.0*rdx2SqVol[0]*phiLy[1]-72.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(378.0*phiLy[0]+36.0*phiLxLy[0]-27.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R2+((-651.0*rdx2SqVol0R2*phiPrevC[1])-51.0*rdx2SqVol0R2*phiLy[1]+252.0*rdx2SqVol0R2*phiLx[1]+258.0*rdx2SqVol[0]*rhoC[0]+(162.0*phiLy[0]+57.0*phiLxLy[0]+231.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol[1]-49.0*rdx2SqVol0R3*phiPrevC[1]-14.0*rdx2SqVol0R3*phiLy[1]+84.0*rdx2SqVol0R2*rhoC[0]+(14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R3)*omega+49.0*phiPrevC[1]*rdx2SqVol1R3+651.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol[1]+49.0*rdx2SqVol0R3*phiPrevC[1])/(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3); 
  phiC[2] = (((84.0*rdx2SqVol1R3+1500.0*rdx2SqVol[0]*rdx2SqVol1R2+336.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[11]+((-126.0*rdx2SqVol[0]*rdx2SqVol1R2)+1080.0*rdx2SqVol0R2*rdx2SqVol[1]+126.0*rdx2SqVol0R3)*bcVals[5]+((-84.0*rdx2SqVol1R2)+516.0*rdx2SqVol[0]*rdx2SqVol[1]+168.0*rdx2SqVol0R2)*rhoC[3]+(168.0*rdx2SqVol1R2+984.0*rdx2SqVol[0]*rdx2SqVol[1]+168.0*rdx2SqVol0R2)*rhoC[2]+((-49.0*rdx2SqVol1R3)-651.0*rdx2SqVol[0]*rdx2SqVol1R2-651.0*rdx2SqVol0R2*rdx2SqVol[1]-49.0*rdx2SqVol0R3)*phiPrevC[2]+((-42.0*rdx2SqVol1R2)+564.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*rhoC[1]+((-14.0*phiLx[1])+49.0*phiLy[0]+14.0*phiLxLy[0])*rdx2SqVol1R3+(252.0*rdx2SqVol[0]*phiLy[1]-51.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(231.0*phiLy[0]+57.0*phiLxLy[0]+162.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R2+((-72.0*rdx2SqVol0R2*phiLy[1])+336.0*rdx2SqVol0R2*phiLx[1]+258.0*rdx2SqVol[0]*rhoC[0]+((-27.0*phiLy[0])+36.0*phiLxLy[0]+378.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol[1]+49.0*rdx2SqVol0R3*phiLx[1]-42.0*rdx2SqVol0R2*rhoC[0]+(7.0*phiLy[0]-7.0*phiLxLy[0])*rdx2SqVol0R3)*omega+(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3)*phiPrevC[2])/(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3); 
  phiC[3] = (((252.0*rdx2SqVol1R3+2484.0*rdx2SqVol[0]*rdx2SqVol1R2+504.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[11]+(504.0*rdx2SqVol[0]*rdx2SqVol1R2+2484.0*rdx2SqVol0R2*rdx2SqVol[1]+252.0*rdx2SqVol0R3)*bcVals[5]+(336.0*rdx2SqVol1R2+1968.0*rdx2SqVol[0]*rdx2SqVol[1]+336.0*rdx2SqVol0R2)*rhoC[3]+((-49.0*rdx2SqVol1R3)-651.0*rdx2SqVol[0]*rdx2SqVol1R2-651.0*rdx2SqVol0R2*rdx2SqVol[1]-49.0*rdx2SqVol0R3)*phiPrevC[3]+((-84.0*rdx2SqVol1R2)+516.0*rdx2SqVol[0]*rdx2SqVol[1]+168.0*rdx2SqVol0R2)*rhoC[2]+(168.0*rdx2SqVol1R2+516.0*rdx2SqVol[0]*rdx2SqVol[1]-84.0*rdx2SqVol0R2)*rhoC[1]+(49.0*phiLy[1]+7.0*phiLx[1]-7.0*phiLxLy[0])*rdx2SqVol1R3+(231.0*rdx2SqVol[0]*phiLy[1]-27.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(504.0*phiLy[0]+87.0*phiLxLy[0]-144.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R2+((-27.0*rdx2SqVol0R2*phiLy[1])+231.0*rdx2SqVol0R2*phiLx[1]+564.0*rdx2SqVol[0]*rhoC[0]+((-144.0*phiLy[0])+87.0*phiLxLy[0]+504.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol[1]+7.0*rdx2SqVol0R3*phiLy[1]+49.0*rdx2SqVol0R3*phiLx[1]-42.0*rdx2SqVol0R2*rhoC[0]-7.0*phiLxLy[0]*rdx2SqVol0R3)*omega+(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3)*phiPrevC[3])/(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol0R3 = std::pow(rdx2SqVol[0],3);
  const double rdx2SqVol0R4 = std::pow(rdx2SqVol[0],4);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double rdx2SqVol1R3 = std::pow(rdx2SqVol[1],3);
  const double rdx2SqVol1R4 = std::pow(rdx2SqVol[1],4);
  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = ((((42.0*rdx2SqVol1R4+864.0*rdx2SqVol[0]*rdx2SqVol1R3+738.0*rdx2SqVol0R2*rdx2SqVol1R2-84.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[10]+(168.0*rdx2SqVol1R4+336.0*rdx2SqVol[0]*rdx2SqVol1R3-48.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[9])*bcVals[11]+(((-84.0*rdx2SqVol[0]*rdx2SqVol1R3)+738.0*rdx2SqVol0R2*rdx2SqVol1R2+864.0*rdx2SqVol0R3*rdx2SqVol[1]+42.0*rdx2SqVol0R4)*bcVals[5]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[3]+(84.0*rdx2SqVol1R3+342.0*rdx2SqVol[0]*rdx2SqVol1R2+216.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+((-42.0*rdx2SqVol1R3)+216.0*rdx2SqVol[0]*rdx2SqVol1R2+342.0*rdx2SqVol0R2*rdx2SqVol[1]+84.0*rdx2SqVol0R3)*rhoC[1]+((-49.0*phiPrevC[0])+49.0*phiLy[0]+14.0*phiLxLy[0]-14.0*phiLx[0])*rdx2SqVol1R4+(189.0*rdx2SqVol[0]*phiLy[1]+81.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+((-700.0*phiPrevC[0])+385.0*phiLy[0]+110.0*phiLxLy[0]-65.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(270.0*rdx2SqVol0R2*phiLy[1]+270.0*rdx2SqVol0R2*phiLx[1]+576.0*rdx2SqVol[0]*rhoC[0]+((-1302.0*phiPrevC[0])+285.0*phiLy[0]+192.0*phiLxLy[0]+285.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(81.0*rdx2SqVol0R3*phiLy[1]+189.0*rdx2SqVol0R3*phiLx[1]+576.0*rdx2SqVol0R2*rhoC[0]+((-700.0*phiPrevC[0])-65.0*phiLy[0]+110.0*phiLxLy[0]+385.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+84.0*rdx2SqVol0R3*rhoC[0]+((-49.0*phiPrevC[0])-14.0*phiLy[0]+14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R4)*bcVals10R2+(((-372.0*rdx2SqVol[0]*rdx2SqVol1R3)+552.0*rdx2SqVol0R2*rdx2SqVol1R2+708.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[5]+((-168.0*rdx2SqVol1R3)+312.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+(336.0*rdx2SqVol1R3+24.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-204.0*rdx2SqVol1R3)+240.0*rdx2SqVol[0]*rdx2SqVol1R2+660.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(12.0*phiLx[1]-280.0*phiPrevC[0]+238.0*phiLy[0]+68.0*phiLxLy[0]-80.0*phiLx[0])*rdx2SqVol1R4+(402.0*rdx2SqVol[0]*phiLy[1]+396.0*rdx2SqVol[0]*phiLx[1]+408.0*rhoC[0]+((-2592.0*phiPrevC[0])+750.0*phiLy[0]+288.0*phiLxLy[0]-108.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(516.0*rdx2SqVol0R2*phiLy[1]+360.0*rdx2SqVol0R2*phiLx[1]+1320.0*rdx2SqVol[0]*rhoC[0]+((-2760.0*phiPrevC[0])+174.0*phiLy[0]+336.0*phiLxLy[0]+636.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(6.0*rdx2SqVol0R3*phiLy[1]+84.0*rdx2SqVol0R3*phiLx[1]+696.0*rdx2SqVol0R2*rhoC[0]+((-448.0*phiPrevC[0])-122.0*phiLy[0]+116.0*phiLxLy[0]+448.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[9]*bcVals[10]+((288.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol[0]*rdx2SqVol1R3)*bcVals[5]+(288.0*rdx2SqVol[0]*rdx2SqVol1R2-144.0*rdx2SqVol1R3)*rhoC[1]+(48.0*phiLx[1]-336.0*phiPrevC[0]+168.0*phiLy[0]+48.0*phiLxLy[0]-96.0*phiLx[0])*rdx2SqVol1R4+(216.0*rdx2SqVol[0]*phiLy[1]+96.0*rdx2SqVol[0]*phiLx[1]+288.0*rhoC[0]+((-960.0*phiPrevC[0])+120.0*phiLy[0]+96.0*phiLxLy[0]+96.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(48.0*rdx2SqVol0R2*phiLx[1]+288.0*rdx2SqVol[0]*rhoC[0]+((-192.0*phiPrevC[0])-48.0*phiLy[0]+48.0*phiLxLy[0]+192.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2)*bcVals9R2)*omega+(49.0*phiPrevC[0]*rdx2SqVol1R4+700.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+700.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1]+49.0*phiPrevC[0]*rdx2SqVol0R4)*bcVals10R2+(280.0*phiPrevC[0]*rdx2SqVol1R4+2592.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+448.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*phiPrevC[0]*rdx2SqVol1R4+960.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+192.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals9R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals9R2); 
  phiC[1] = ((((126.0*rdx2SqVol1R4+1206.0*rdx2SqVol[0]*rdx2SqVol1R3+954.0*rdx2SqVol0R2*rdx2SqVol1R2-126.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[10]+(168.0*rdx2SqVol1R4+552.0*rdx2SqVol[0]*rdx2SqVol1R3-48.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[9])*bcVals[11]+((336.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+1584.0*rdx2SqVol0R3*rdx2SqVol[1]+84.0*rdx2SqVol0R4)*bcVals[5]+(168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[3]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[1]+((-49.0*phiPrevC[1])+49.0*phiLy[1]-7.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol1R4+((-700.0*rdx2SqVol[0]*phiPrevC[1])+385.0*rdx2SqVol[0]*phiLy[1]-72.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(378.0*phiLy[0]+29.0*phiLxLy[0]-20.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-1302.0*rdx2SqVol0R2*phiPrevC[1])+285.0*rdx2SqVol0R2*phiLy[1]+180.0*rdx2SqVol0R2*phiLx[1]+216.0*rdx2SqVol[0]*rhoC[0]+(540.0*phiLy[0]+93.0*phiLxLy[0]+204.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-700.0*rdx2SqVol0R3*phiPrevC[1])-65.0*rdx2SqVol0R3*phiLy[1]+252.0*rdx2SqVol0R3*phiLx[1]+342.0*rdx2SqVol0R2*rhoC[0]+(162.0*phiLy[0]+71.0*phiLxLy[0]+280.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-49.0*rdx2SqVol0R4*phiPrevC[1]-14.0*rdx2SqVol0R4*phiLy[1]+84.0*rdx2SqVol0R3*rhoC[0]+(14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R4)*bcVals10R2+((984.0*rdx2SqVol[0]*rdx2SqVol1R3+2688.0*rdx2SqVol0R2*rdx2SqVol1R2+1272.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[5]+(336.0*rdx2SqVol1R3-192.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+((-168.0*rdx2SqVol1R3)+744.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(648.0*rdx2SqVol1R3+2880.0*rdx2SqVol[0]*rdx2SqVol1R2+1368.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+((-280.0*phiPrevC[1])+182.0*phiLy[1]-6.0*phiLx[1]-28.0*phiLy[0]-34.0*phiLxLy[0]+40.0*phiLx[0])*rdx2SqVol1R4+((-2592.0*rdx2SqVol[0]*phiPrevC[1])+858.0*rdx2SqVol[0]*phiLy[1]-174.0*rdx2SqVol[0]*phiLx[1]-204.0*rhoC[0]+(816.0*phiLy[0]+6.0*phiLxLy[0]-120.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-2760.0*rdx2SqVol0R2*phiPrevC[1])+126.0*rdx2SqVol0R2*phiLy[1]+390.0*rdx2SqVol0R2*phiLx[1]+240.0*rdx2SqVol[0]*rhoC[0]+(1068.0*phiLy[0]+150.0*phiLxLy[0]+72.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-448.0*rdx2SqVol0R3*phiPrevC[1])-118.0*rdx2SqVol0R3*phiLy[1]+126.0*rdx2SqVol0R3*phiLx[1]+660.0*rdx2SqVol0R2*rhoC[0]+(8.0*phiLy[0]+110.0*phiLxLy[0]+448.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[9]*bcVals[10]+((576.0*rdx2SqVol[0]*rdx2SqVol1R3+576.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[5]+(576.0*rdx2SqVol1R3+576.0*rdx2SqVol[0]*rdx2SqVol1R2)*rhoC[1]+((-336.0*phiPrevC[1])+168.0*phiLy[1]-24.0*phiLx[1]-24.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol1R4+((-960.0*rdx2SqVol[0]*phiPrevC[1])+120.0*rdx2SqVol[0]*phiLy[1]+24.0*rdx2SqVol[0]*phiLx[1]-144.0*rhoC[0]+(432.0*phiLy[0]+24.0*phiLxLy[0]-192.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-192.0*rdx2SqVol0R2*phiPrevC[1])-48.0*rdx2SqVol0R2*phiLy[1]+48.0*rdx2SqVol0R2*phiLx[1]+288.0*rdx2SqVol[0]*rhoC[0]+(48.0*phiLxLy[0]+192.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2)*bcVals9R2)*omega+(49.0*phiPrevC[1]*rdx2SqVol1R4+700.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+700.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiPrevC[1])*bcVals10R2+(280.0*phiPrevC[1]*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+448.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*phiPrevC[1]*rdx2SqVol1R4+960.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2)*bcVals9R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals9R2); 
  phiC[2] = ((((84.0*rdx2SqVol1R4+1584.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+336.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[9])*bcVals[11]+(((-126.0*rdx2SqVol[0]*rdx2SqVol1R3)+954.0*rdx2SqVol0R2*rdx2SqVol1R2+1206.0*rdx2SqVol0R3*rdx2SqVol[1]+126.0*rdx2SqVol0R4)*bcVals[5]+((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[3]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiPrevC[2]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[1]+((-14.0*phiLx[1])+49.0*phiLy[0]+14.0*phiLxLy[0])*rdx2SqVol1R4+(252.0*rdx2SqVol[0]*phiLy[1]-65.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(280.0*phiLy[0]+71.0*phiLxLy[0]+162.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(180.0*rdx2SqVol0R2*phiLy[1]+285.0*rdx2SqVol0R2*phiLx[1]+342.0*rdx2SqVol[0]*rhoC[0]+(204.0*phiLy[0]+93.0*phiLxLy[0]+540.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-72.0*rdx2SqVol0R3*phiLy[1])+385.0*rdx2SqVol0R3*phiLx[1]+216.0*rdx2SqVol0R2*rhoC[0]+((-20.0*phiLy[0])+29.0*phiLxLy[0]+378.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]+(7.0*phiLy[0]-7.0*phiLxLy[0])*rdx2SqVol0R4)*bcVals10R2+(((-504.0*rdx2SqVol[0]*rdx2SqVol1R3)-216.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[5]+((-336.0*rdx2SqVol1R3)-960.0*rdx2SqVol[0]*rdx2SqVol1R2-192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+(672.0*rdx2SqVol1R3+1920.0*rdx2SqVol[0]*rdx2SqVol1R2+384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-280.0*rdx2SqVol1R4)-2592.0*rdx2SqVol[0]*rdx2SqVol1R3-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[2]+((-168.0*rdx2SqVol1R3)+744.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+((-56.0*phiLx[1])+196.0*phiLy[0]+56.0*phiLxLy[0])*rdx2SqVol1R4+(336.0*rdx2SqVol[0]*phiLy[1]-36.0*rdx2SqVol[0]*phiLx[1]+336.0*rhoC[0]+(60.0*phiLxLy[0]+648.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-120.0*rdx2SqVol0R2*phiLy[1])+564.0*rdx2SqVol0R2*phiLx[1]+24.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLy[0]-12.0*phiLxLy[0]+432.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-24.0*rdx2SqVol0R3*phiLy[1])+112.0*rdx2SqVol0R3*phiLx[1]-96.0*rdx2SqVol0R2*rhoC[0]+(40.0*phiLy[0]-16.0*phiLxLy[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[9]*bcVals[10]+((-336.0*rdx2SqVol1R4)-960.0*rdx2SqVol[0]*rdx2SqVol1R3-192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiPrevC[2]*bcVals9R2)*omega+(49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiPrevC[2]*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[2]*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiPrevC[2]*bcVals9R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals9R2); 
  phiC[3] = ((((252.0*rdx2SqVol1R4+2736.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+504.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[9])*bcVals[11]+((504.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+2736.0*rdx2SqVol0R3*rdx2SqVol[1]+252.0*rdx2SqVol0R4)*bcVals[5]+(336.0*rdx2SqVol1R3+2304.0*rdx2SqVol[0]*rdx2SqVol1R2+2304.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[3]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiPrevC[3]+((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[1]+7.0*phiLx[1]-7.0*phiLxLy[0])*rdx2SqVol1R4+(280.0*rdx2SqVol[0]*phiLy[1]-20.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(504.0*phiLy[0]+80.0*phiLxLy[0]-144.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(204.0*rdx2SqVol0R2*phiLy[1]+204.0*rdx2SqVol0R2*phiLx[1]+522.0*rdx2SqVol[0]*rhoC[0]+(360.0*phiLy[0]+174.0*phiLxLy[0]+360.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-20.0*rdx2SqVol0R3*phiLy[1])+280.0*rdx2SqVol0R3*phiLx[1]+522.0*rdx2SqVol0R2*rhoC[0]+((-144.0*phiLy[0])+80.0*phiLxLy[0]+504.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+7.0*rdx2SqVol0R4*phiLy[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]-7.0*phiLxLy[0]*rdx2SqVol0R4)*bcVals10R2+((1008.0*rdx2SqVol[0]*rdx2SqVol1R3+1728.0*rdx2SqVol0R2*rdx2SqVol1R2+288.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[5]+(672.0*rdx2SqVol1R3+1920.0*rdx2SqVol[0]*rdx2SqVol1R2+384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+((-280.0*rdx2SqVol1R4)-2592.0*rdx2SqVol[0]*rdx2SqVol1R3-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[3]+((-336.0*rdx2SqVol1R3)-960.0*rdx2SqVol[0]*rdx2SqVol1R2-192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(336.0*rdx2SqVol1R3-192.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(84.0*phiLy[1]+28.0*phiLx[1]-56.0*phiLy[0]-28.0*phiLxLy[0])*rdx2SqVol1R4+((-96.0*rdx2SqVol[0]*phiLy[1])+72.0*rdx2SqVol[0]*phiLx[1]-168.0*rhoC[0]+(288.0*phiLy[0]+24.0*phiLxLy[0]-432.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(60.0*rdx2SqVol0R2*phiLy[1]-228.0*rdx2SqVol0R2*phiLx[1]+312.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLxLy[0]-120.0*phiLy[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(24.0*rdx2SqVol0R3*phiLy[1]-56.0*rdx2SqVol0R3*phiLx[1]+48.0*rdx2SqVol0R2*rhoC[0]+(8.0*phiLxLy[0]-32.0*phiLy[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[9]*bcVals[10]+((-336.0*rdx2SqVol1R4)-960.0*rdx2SqVol[0]*rdx2SqVol1R3-192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiPrevC[3]*bcVals9R2)*omega+(49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiPrevC[3]*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[3]*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiPrevC[3]*bcVals9R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals9R2); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = ((((7.0*rdx2SqVol1R2+14.0*rdx2SqVol[0]*rdx2SqVol[1]-2.0*rdx2SqVol0R2)*bcVals[4]+(12.0*rdx2SqVol[0]*rdx2SqVol[1]-6.0*rdx2SqVol0R2)*bcVals[3])*bcVals[11]+(12.0*rdx2SqVol0R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[5]+((12.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*rhoC[1]+(2.0*phiLx[1]-14.0*phiPrevC[0]+7.0*phiLy[0]+2.0*phiLxLy[0]-4.0*phiLx[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiLy[1]+4.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+((-40.0*phiPrevC[0])+5.0*phiLy[0]+4.0*phiLxLy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiLx[1]+12.0*rdx2SqVol[0]*rhoC[0]+((-8.0*phiPrevC[0])-2.0*phiLy[0]+2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+((2.0*rdx2SqVol[0]*phiLy[1]+2.0*rdx2SqVol[0]*phiLx[1]+((-16.0*phiPrevC[0])+8.0*phiLy[0]+2.0*phiLxLy[0]-4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiLy[1]+2.0*rdx2SqVol0R2*phiLx[1]+12.0*rdx2SqVol[0]*rhoC[0]+((-16.0*phiPrevC[0])-4.0*phiLy[0]+2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*bcVals[3])*omega+(14.0*phiPrevC[0]*rdx2SqVol1R2+40.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+8.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[4]+(16.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol[1]+16.0*phiPrevC[0]*rdx2SqVol0R2)*bcVals[3])/((14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 
  phiC[1] = ((((7.0*rdx2SqVol1R2+23.0*rdx2SqVol[0]*rdx2SqVol[1]-2.0*rdx2SqVol0R2)*bcVals[4]+((-8.0*rdx2SqVol[0]*rdx2SqVol[1])-8.0*rdx2SqVol0R2)*bcVals[3])*bcVals[11]+(24.0*rdx2SqVol[0]*rdx2SqVol[1]+24.0*rdx2SqVol0R2)*bcVals[5]+((24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+((-14.0*phiPrevC[1])+7.0*phiLy[1]-1.0*phiLx[1]-1.0*phiLxLy[0]+2.0*phiLx[0])*rdx2SqVol1R2+((-40.0*rdx2SqVol[0]*phiPrevC[1])+5.0*rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]-6.0*rhoC[0]+(18.0*phiLy[0]+phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-8.0*rdx2SqVol0R2*phiPrevC[1]-2.0*rdx2SqVol0R2*phiLy[1]+2.0*rdx2SqVol0R2*phiLx[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+((-16.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1])-16.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[3])*omega+(14.0*phiPrevC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[4]+(16.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol[1]+16.0*rdx2SqVol0R2*phiPrevC[1])*bcVals[3])/((14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 
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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol0R3 = std::pow(rdx2SqVol[0],3);
  const double rdx2SqVol0R4 = std::pow(rdx2SqVol[0],4);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double rdx2SqVol1R3 = std::pow(rdx2SqVol[1],3);
  const double rdx2SqVol1R4 = std::pow(rdx2SqVol[1],4);
  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);

  phiC[0] = ((((42.0*rdx2SqVol1R4+864.0*rdx2SqVol[0]*rdx2SqVol1R3+738.0*rdx2SqVol0R2*rdx2SqVol1R2-84.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(708.0*rdx2SqVol[0]*rdx2SqVol1R3+552.0*rdx2SqVol0R2*rdx2SqVol1R2-372.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(288.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[11]+(((-84.0*rdx2SqVol[0]*rdx2SqVol1R3)+738.0*rdx2SqVol0R2*rdx2SqVol1R2+864.0*rdx2SqVol0R3*rdx2SqVol[1]+42.0*rdx2SqVol0R4)*bcVals[4]+((-48.0*rdx2SqVol0R2*rdx2SqVol1R2)+336.0*rdx2SqVol0R3*rdx2SqVol[1]+168.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+(((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[3]+(84.0*rdx2SqVol1R3+342.0*rdx2SqVol[0]*rdx2SqVol1R2+216.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+((-42.0*rdx2SqVol1R3)+216.0*rdx2SqVol[0]*rdx2SqVol1R2+342.0*rdx2SqVol0R2*rdx2SqVol[1]+84.0*rdx2SqVol0R3)*rhoC[1]+((-49.0*phiPrevC[0])+49.0*phiLy[0]+14.0*phiLxLy[0]-14.0*phiLx[0])*rdx2SqVol1R4+(189.0*rdx2SqVol[0]*phiLy[1]+81.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+((-700.0*phiPrevC[0])+385.0*phiLy[0]+110.0*phiLxLy[0]-65.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(270.0*rdx2SqVol0R2*phiLy[1]+270.0*rdx2SqVol0R2*phiLx[1]+576.0*rdx2SqVol[0]*rhoC[0]+((-1302.0*phiPrevC[0])+285.0*phiLy[0]+192.0*phiLxLy[0]+285.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(81.0*rdx2SqVol0R3*phiLy[1]+189.0*rdx2SqVol0R3*phiLx[1]+576.0*rdx2SqVol0R2*rhoC[0]+((-700.0*phiPrevC[0])-65.0*phiLy[0]+110.0*phiLxLy[0]+385.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+84.0*rdx2SqVol0R3*rhoC[0]+((-49.0*phiPrevC[0])-14.0*phiLy[0]+14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R4)*bcVals4R2+((48.0*rdx2SqVol[0]*rdx2SqVol1R2+312.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((660.0*rdx2SqVol[0]*rdx2SqVol1R2+240.0*rdx2SqVol0R2*rdx2SqVol[1]-204.0*rdx2SqVol0R3)*rhoC[2]+((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)+24.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[1]+(84.0*rdx2SqVol[0]*phiLy[1]+6.0*rdx2SqVol[0]*phiLx[1]+((-448.0*phiPrevC[0])+448.0*phiLy[0]+116.0*phiLxLy[0]-122.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(360.0*rdx2SqVol0R2*phiLy[1]+516.0*rdx2SqVol0R2*phiLx[1]+696.0*rdx2SqVol[0]*rhoC[0]+((-2760.0*phiPrevC[0])+636.0*phiLy[0]+336.0*phiLxLy[0]+174.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(396.0*rdx2SqVol0R3*phiLy[1]+402.0*rdx2SqVol0R3*phiLx[1]+1320.0*rdx2SqVol0R2*rhoC[0]+((-2592.0*phiPrevC[0])-108.0*phiLy[0]+288.0*phiLxLy[0]+750.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+12.0*rdx2SqVol0R4*phiLy[1]+408.0*rdx2SqVol0R3*rhoC[0]+((-280.0*phiPrevC[0])-80.0*phiLy[0]+68.0*phiLxLy[0]+238.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((288.0*rdx2SqVol0R2*rdx2SqVol[1]-144.0*rdx2SqVol0R3)*rhoC[2]+(48.0*rdx2SqVol0R2*phiLy[1]+((-192.0*phiPrevC[0])+192.0*phiLy[0]+48.0*phiLxLy[0]-48.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(96.0*rdx2SqVol0R3*phiLy[1]+216.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+((-960.0*phiPrevC[0])+96.0*phiLy[0]+96.0*phiLxLy[0]+120.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+48.0*rdx2SqVol0R4*phiLy[1]+288.0*rdx2SqVol0R3*rhoC[0]+((-336.0*phiPrevC[0])-96.0*phiLy[0]+48.0*phiLxLy[0]+168.0*phiLx[0])*rdx2SqVol0R4)*bcVals3R2)*omega+(49.0*phiPrevC[0]*rdx2SqVol1R4+700.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+700.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1]+49.0*phiPrevC[0]*rdx2SqVol0R4)*bcVals4R2+(448.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1]+280.0*phiPrevC[0]*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+960.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1]+336.0*phiPrevC[0]*rdx2SqVol0R4)*bcVals3R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2); 
  phiC[1] = ((((126.0*rdx2SqVol1R4+1206.0*rdx2SqVol[0]*rdx2SqVol1R3+954.0*rdx2SqVol0R2*rdx2SqVol1R2-126.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+((-144.0*rdx2SqVol[0]*rdx2SqVol1R3)-216.0*rdx2SqVol0R2*rdx2SqVol1R2-504.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4])*bcVals[11]+((336.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+1584.0*rdx2SqVol0R3*rdx2SqVol[1]+84.0*rdx2SqVol0R4)*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+((168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[3]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[1]+((-49.0*phiPrevC[1])+49.0*phiLy[1]-7.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol1R4+((-700.0*rdx2SqVol[0]*phiPrevC[1])+385.0*rdx2SqVol[0]*phiLy[1]-72.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(378.0*phiLy[0]+29.0*phiLxLy[0]-20.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-1302.0*rdx2SqVol0R2*phiPrevC[1])+285.0*rdx2SqVol0R2*phiLy[1]+180.0*rdx2SqVol0R2*phiLx[1]+216.0*rdx2SqVol[0]*rhoC[0]+(540.0*phiLy[0]+93.0*phiLxLy[0]+204.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-700.0*rdx2SqVol0R3*phiPrevC[1])-65.0*rdx2SqVol0R3*phiLy[1]+252.0*rdx2SqVol0R3*phiLx[1]+342.0*rdx2SqVol0R2*rhoC[0]+(162.0*phiLy[0]+71.0*phiLxLy[0]+280.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-49.0*rdx2SqVol0R4*phiPrevC[1]-14.0*rdx2SqVol0R4*phiLy[1]+84.0*rdx2SqVol0R3*rhoC[0]+(14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R4)*bcVals4R2+(((-192.0*rdx2SqVol[0]*rdx2SqVol1R2)-960.0*rdx2SqVol0R2*rdx2SqVol[1]-336.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((48.0*rdx2SqVol[0]*rdx2SqVol1R2+744.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*rhoC[2]+(384.0*rdx2SqVol[0]*rdx2SqVol1R2+1920.0*rdx2SqVol0R2*rdx2SqVol[1]+672.0*rdx2SqVol0R3)*rhoC[1]+((-448.0*rdx2SqVol[0]*phiPrevC[1])+112.0*rdx2SqVol[0]*phiLy[1]-24.0*rdx2SqVol[0]*phiLx[1]+(40.0*phiLx[0]-16.0*phiLxLy[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-2760.0*rdx2SqVol0R2*phiPrevC[1])+564.0*rdx2SqVol0R2*phiLy[1]-120.0*rdx2SqVol0R2*phiLx[1]-96.0*rdx2SqVol[0]*rhoC[0]+(432.0*phiLy[0]-12.0*phiLxLy[0]+60.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-2592.0*rdx2SqVol0R3*phiPrevC[1])-36.0*rdx2SqVol0R3*phiLy[1]+336.0*rdx2SqVol0R3*phiLx[1]+24.0*rdx2SqVol0R2*rhoC[0]+(648.0*phiLy[0]+60.0*phiLxLy[0])*rdx2SqVol0R3)*rdx2SqVol[1]-280.0*rdx2SqVol0R4*phiPrevC[1]-56.0*rdx2SqVol0R4*phiLy[1]+336.0*rdx2SqVol0R3*rhoC[0]+(56.0*phiLxLy[0]+196.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((-192.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1]-336.0*rdx2SqVol0R4*phiPrevC[1])*bcVals3R2)*omega+(49.0*phiPrevC[1]*rdx2SqVol1R4+700.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+700.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiPrevC[1])*bcVals4R2+(448.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1]+280.0*rdx2SqVol0R4*phiPrevC[1])*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+960.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1]+336.0*rdx2SqVol0R4*phiPrevC[1])*bcVals3R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2); 
  phiC[2] = ((((84.0*rdx2SqVol1R4+1584.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+336.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1272.0*rdx2SqVol[0]*rdx2SqVol1R3+2688.0*rdx2SqVol0R2*rdx2SqVol1R2+984.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(576.0*rdx2SqVol0R2*rdx2SqVol1R2+576.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[11]+(((-126.0*rdx2SqVol[0]*rdx2SqVol1R3)+954.0*rdx2SqVol0R2*rdx2SqVol1R2+1206.0*rdx2SqVol0R3*rdx2SqVol[1]+126.0*rdx2SqVol0R4)*bcVals[4]+((-48.0*rdx2SqVol0R2*rdx2SqVol1R2)+552.0*rdx2SqVol0R3*rdx2SqVol[1]+168.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+(((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[3]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiPrevC[2]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[1]+((-14.0*phiLx[1])+49.0*phiLy[0]+14.0*phiLxLy[0])*rdx2SqVol1R4+(252.0*rdx2SqVol[0]*phiLy[1]-65.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(280.0*phiLy[0]+71.0*phiLxLy[0]+162.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(180.0*rdx2SqVol0R2*phiLy[1]+285.0*rdx2SqVol0R2*phiLx[1]+342.0*rdx2SqVol[0]*rhoC[0]+(204.0*phiLy[0]+93.0*phiLxLy[0]+540.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-72.0*rdx2SqVol0R3*phiLy[1])+385.0*rdx2SqVol0R3*phiLx[1]+216.0*rdx2SqVol0R2*rhoC[0]+((-20.0*phiLy[0])+29.0*phiLxLy[0]+378.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]+(7.0*phiLy[0]-7.0*phiLxLy[0])*rdx2SqVol0R4)*bcVals4R2+(((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)-192.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((1368.0*rdx2SqVol[0]*rdx2SqVol1R2+2880.0*rdx2SqVol0R2*rdx2SqVol[1]+648.0*rdx2SqVol0R3)*rhoC[2]+((-448.0*rdx2SqVol[0]*rdx2SqVol1R3)-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-2592.0*rdx2SqVol0R3*rdx2SqVol[1]-280.0*rdx2SqVol0R4)*phiPrevC[2]+(48.0*rdx2SqVol[0]*rdx2SqVol1R2+744.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*rhoC[1]+(126.0*rdx2SqVol[0]*phiLy[1]-118.0*rdx2SqVol[0]*phiLx[1]+(448.0*phiLy[0]+110.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(390.0*rdx2SqVol0R2*phiLy[1]+126.0*rdx2SqVol0R2*phiLx[1]+660.0*rdx2SqVol[0]*rhoC[0]+(72.0*phiLy[0]+150.0*phiLxLy[0]+1068.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-174.0*rdx2SqVol0R3*phiLy[1])+858.0*rdx2SqVol0R3*phiLx[1]+240.0*rdx2SqVol0R2*rhoC[0]+((-120.0*phiLy[0])+6.0*phiLxLy[0]+816.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-6.0*rdx2SqVol0R4*phiLy[1]+182.0*rdx2SqVol0R4*phiLx[1]-204.0*rdx2SqVol0R3*rhoC[0]+(40.0*phiLy[0]-34.0*phiLxLy[0]-28.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((576.0*rdx2SqVol0R2*rdx2SqVol[1]+576.0*rdx2SqVol0R3)*rhoC[2]+((-192.0*rdx2SqVol0R2*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*rdx2SqVol[1]-336.0*rdx2SqVol0R4)*phiPrevC[2]+(48.0*rdx2SqVol0R2*phiLy[1]-48.0*rdx2SqVol0R2*phiLx[1]+(192.0*phiLy[0]+48.0*phiLxLy[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(24.0*rdx2SqVol0R3*phiLy[1]+120.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+((-192.0*phiLy[0])+24.0*phiLxLy[0]+432.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-24.0*rdx2SqVol0R4*phiLy[1]+168.0*rdx2SqVol0R4*phiLx[1]-144.0*rdx2SqVol0R3*rhoC[0]+(48.0*phiLy[0]-24.0*phiLxLy[0])*rdx2SqVol0R4)*bcVals3R2)*omega+(49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiPrevC[2]*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*phiPrevC[2]*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*phiPrevC[2]*bcVals3R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2); 
  phiC[3] = ((((252.0*rdx2SqVol1R4+2736.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+504.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(288.0*rdx2SqVol[0]*rdx2SqVol1R3+1728.0*rdx2SqVol0R2*rdx2SqVol1R2+1008.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4])*bcVals[11]+((504.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+2736.0*rdx2SqVol0R3*rdx2SqVol[1]+252.0*rdx2SqVol0R4)*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+((336.0*rdx2SqVol1R3+2304.0*rdx2SqVol[0]*rdx2SqVol1R2+2304.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[3]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiPrevC[3]+((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[1]+7.0*phiLx[1]-7.0*phiLxLy[0])*rdx2SqVol1R4+(280.0*rdx2SqVol[0]*phiLy[1]-20.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(504.0*phiLy[0]+80.0*phiLxLy[0]-144.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(204.0*rdx2SqVol0R2*phiLy[1]+204.0*rdx2SqVol0R2*phiLx[1]+522.0*rdx2SqVol[0]*rhoC[0]+(360.0*phiLy[0]+174.0*phiLxLy[0]+360.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-20.0*rdx2SqVol0R3*phiLy[1])+280.0*rdx2SqVol0R3*phiLx[1]+522.0*rdx2SqVol0R2*rhoC[0]+((-144.0*phiLy[0])+80.0*phiLxLy[0]+504.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+7.0*rdx2SqVol0R4*phiLy[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]-7.0*phiLxLy[0]*rdx2SqVol0R4)*bcVals4R2+((384.0*rdx2SqVol[0]*rdx2SqVol1R2+1920.0*rdx2SqVol0R2*rdx2SqVol[1]+672.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((-448.0*rdx2SqVol[0]*rdx2SqVol1R3)-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-2592.0*rdx2SqVol0R3*rdx2SqVol[1]-280.0*rdx2SqVol0R4)*bcVals[3]*phiPrevC[3]+(((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)-192.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[2]+((-192.0*rdx2SqVol[0]*rdx2SqVol1R2)-960.0*rdx2SqVol0R2*rdx2SqVol[1]-336.0*rdx2SqVol0R3)*rhoC[1]+((-56.0*rdx2SqVol[0]*phiLy[1])+24.0*rdx2SqVol[0]*phiLx[1]+(8.0*phiLxLy[0]-32.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-228.0*rdx2SqVol0R2*phiLy[1])+60.0*rdx2SqVol0R2*phiLx[1]+48.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLxLy[0]-120.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(72.0*rdx2SqVol0R3*phiLy[1]-96.0*rdx2SqVol0R3*phiLx[1]+312.0*rdx2SqVol0R2*rhoC[0]+((-432.0*phiLy[0])+24.0*phiLxLy[0]+288.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+28.0*rdx2SqVol0R4*phiLy[1]+84.0*rdx2SqVol0R4*phiLx[1]-168.0*rdx2SqVol0R3*rhoC[0]+((-28.0*phiLxLy[0])-56.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((-192.0*rdx2SqVol0R2*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*rdx2SqVol[1]-336.0*rdx2SqVol0R4)*bcVals3R2*phiPrevC[3])*omega+(49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiPrevC[3]*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*phiPrevC[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2*phiPrevC[3])/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2); 

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
  double *phiUx = phiPrev[2]; 
  double *phiLy = phiPrev[3]; 
  double *phiUy = phiPrev[4]; 
  double *phiLxLy = phiPrev[5]; 
  double *phiLxUy = phiPrev[6]; 
  double *phiUxLy = phiPrev[7]; 
  double *phiUxUy = phiPrev[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol0R3 = std::pow(rdx2SqVol[0],3);
  const double rdx2SqVol0R4 = std::pow(rdx2SqVol[0],4);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double rdx2SqVol1R3 = std::pow(rdx2SqVol[1],3);
  const double rdx2SqVol1R4 = std::pow(rdx2SqVol[1],4);
  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);
  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = (((((42.0*rdx2SqVol1R4+864.0*rdx2SqVol[0]*rdx2SqVol1R3+738.0*rdx2SqVol0R2*rdx2SqVol1R2-84.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(708.0*rdx2SqVol[0]*rdx2SqVol1R3+552.0*rdx2SqVol0R2*rdx2SqVol1R2-372.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(288.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[10]+((168.0*rdx2SqVol1R4+336.0*rdx2SqVol[0]*rdx2SqVol1R3-48.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(288.0*rdx2SqVol[0]*rdx2SqVol1R3-144.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals[9])*bcVals[11]+((((-84.0*rdx2SqVol[0]*rdx2SqVol1R3)+738.0*rdx2SqVol0R2*rdx2SqVol1R2+864.0*rdx2SqVol0R3*rdx2SqVol[1]+42.0*rdx2SqVol0R4)*bcVals[4]+((-48.0*rdx2SqVol0R2*rdx2SqVol1R2)+336.0*rdx2SqVol0R3*rdx2SqVol[1]+168.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+(((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[3]+(84.0*rdx2SqVol1R3+342.0*rdx2SqVol[0]*rdx2SqVol1R2+216.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+((-42.0*rdx2SqVol1R3)+216.0*rdx2SqVol[0]*rdx2SqVol1R2+342.0*rdx2SqVol0R2*rdx2SqVol[1]+84.0*rdx2SqVol0R3)*rhoC[1]+((-49.0*phiPrevC[0])+49.0*phiLy[0]+14.0*phiLxLy[0]-14.0*phiLx[0])*rdx2SqVol1R4+(189.0*rdx2SqVol[0]*phiLy[1]+81.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+((-700.0*phiPrevC[0])+385.0*phiLy[0]+110.0*phiLxLy[0]-65.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(270.0*rdx2SqVol0R2*phiLy[1]+270.0*rdx2SqVol0R2*phiLx[1]+576.0*rdx2SqVol[0]*rhoC[0]+((-1302.0*phiPrevC[0])+285.0*phiLy[0]+192.0*phiLxLy[0]+285.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(81.0*rdx2SqVol0R3*phiLy[1]+189.0*rdx2SqVol0R3*phiLx[1]+576.0*rdx2SqVol0R2*rhoC[0]+((-700.0*phiPrevC[0])-65.0*phiLy[0]+110.0*phiLxLy[0]+385.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+84.0*rdx2SqVol0R3*rhoC[0]+((-49.0*phiPrevC[0])-14.0*phiLy[0]+14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R4)*bcVals4R2+((48.0*rdx2SqVol[0]*rdx2SqVol1R2+312.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((660.0*rdx2SqVol[0]*rdx2SqVol1R2+240.0*rdx2SqVol0R2*rdx2SqVol[1]-204.0*rdx2SqVol0R3)*rhoC[2]+((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)+24.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[1]+(84.0*rdx2SqVol[0]*phiLy[1]+6.0*rdx2SqVol[0]*phiLx[1]+((-448.0*phiPrevC[0])+448.0*phiLy[0]+116.0*phiLxLy[0]-122.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(360.0*rdx2SqVol0R2*phiLy[1]+516.0*rdx2SqVol0R2*phiLx[1]+696.0*rdx2SqVol[0]*rhoC[0]+((-2760.0*phiPrevC[0])+636.0*phiLy[0]+336.0*phiLxLy[0]+174.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(396.0*rdx2SqVol0R3*phiLy[1]+402.0*rdx2SqVol0R3*phiLx[1]+1320.0*rdx2SqVol0R2*rhoC[0]+((-2592.0*phiPrevC[0])-108.0*phiLy[0]+288.0*phiLxLy[0]+750.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+12.0*rdx2SqVol0R4*phiLy[1]+408.0*rdx2SqVol0R3*rhoC[0]+((-280.0*phiPrevC[0])-80.0*phiLy[0]+68.0*phiLxLy[0]+238.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((288.0*rdx2SqVol0R2*rdx2SqVol[1]-144.0*rdx2SqVol0R3)*rhoC[2]+(48.0*rdx2SqVol0R2*phiLy[1]+((-192.0*phiPrevC[0])+192.0*phiLy[0]+48.0*phiLxLy[0]-48.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(96.0*rdx2SqVol0R3*phiLy[1]+216.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+((-960.0*phiPrevC[0])+96.0*phiLy[0]+96.0*phiLxLy[0]+120.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+48.0*rdx2SqVol0R4*phiLy[1]+288.0*rdx2SqVol0R3*rhoC[0]+((-336.0*phiPrevC[0])-96.0*phiLy[0]+48.0*phiLxLy[0]+168.0*phiLx[0])*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((((-372.0*rdx2SqVol[0]*rdx2SqVol1R3)+552.0*rdx2SqVol0R2*rdx2SqVol1R2+708.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[4]+(288.0*rdx2SqVol0R3*rdx2SqVol[1]-144.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3])*bcVals[5]+(((-168.0*rdx2SqVol1R3)+312.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+(336.0*rdx2SqVol1R3+24.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-204.0*rdx2SqVol1R3)+240.0*rdx2SqVol[0]*rdx2SqVol1R2+660.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(12.0*phiLx[1]-280.0*phiPrevC[0]+238.0*phiLy[0]+68.0*phiLxLy[0]-80.0*phiLx[0])*rdx2SqVol1R4+(402.0*rdx2SqVol[0]*phiLy[1]+396.0*rdx2SqVol[0]*phiLx[1]+408.0*rhoC[0]+((-2592.0*phiPrevC[0])+750.0*phiLy[0]+288.0*phiLxLy[0]-108.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(516.0*rdx2SqVol0R2*phiLy[1]+360.0*rdx2SqVol0R2*phiLx[1]+1320.0*rdx2SqVol[0]*rhoC[0]+((-2760.0*phiPrevC[0])+174.0*phiLy[0]+336.0*phiLxLy[0]+636.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(6.0*rdx2SqVol0R3*phiLy[1]+84.0*rdx2SqVol0R3*phiLx[1]+696.0*rdx2SqVol0R2*rhoC[0]+((-448.0*phiPrevC[0])-122.0*phiLy[0]+116.0*phiLxLy[0]+448.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals4R2+((288.0*rdx2SqVol[0]*rdx2SqVol1R2-144.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(288.0*rdx2SqVol0R2*rdx2SqVol[1]-144.0*rdx2SqVol[0]*rdx2SqVol1R2)*rhoC[1]+(120.0*rdx2SqVol[0]*phiLy[1]+120.0*rdx2SqVol[0]*phiLx[1]+((-1104.0*phiPrevC[0])+648.0*phiLy[0]+168.0*phiLxLy[0]-288.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(456.0*rdx2SqVol0R2*phiLy[1]+456.0*rdx2SqVol0R2*phiLx[1]+1008.0*rdx2SqVol[0]*rhoC[0]+((-3072.0*phiPrevC[0])+360.0*phiLy[0]+336.0*phiLxLy[0]+360.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(120.0*rdx2SqVol0R3*phiLy[1]+120.0*rdx2SqVol0R3*phiLx[1]+1008.0*rdx2SqVol0R2*rhoC[0]+((-1104.0*phiPrevC[0])-288.0*phiLy[0]+168.0*phiLxLy[0]+648.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[3]*bcVals[4]+((48.0*rdx2SqVol0R2*phiLy[1]+48.0*rdx2SqVol0R2*phiLx[1]+((-384.0*phiPrevC[0])+192.0*phiLy[0]+48.0*phiLxLy[0]-96.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(48.0*rdx2SqVol0R3*phiLy[1]+48.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+((-384.0*phiPrevC[0])-96.0*phiLy[0]+48.0*phiLxLy[0]+192.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((288.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol[0]*rdx2SqVol1R3)*bcVals[4]*bcVals[5]+((288.0*rdx2SqVol[0]*rdx2SqVol1R2-144.0*rdx2SqVol1R3)*rhoC[1]+(48.0*phiLx[1]-336.0*phiPrevC[0]+168.0*phiLy[0]+48.0*phiLxLy[0]-96.0*phiLx[0])*rdx2SqVol1R4+(216.0*rdx2SqVol[0]*phiLy[1]+96.0*rdx2SqVol[0]*phiLx[1]+288.0*rhoC[0]+((-960.0*phiPrevC[0])+120.0*phiLy[0]+96.0*phiLxLy[0]+96.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(48.0*rdx2SqVol0R2*phiLx[1]+288.0*rdx2SqVol[0]*rhoC[0]+((-192.0*phiPrevC[0])-48.0*phiLy[0]+48.0*phiLxLy[0]+192.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2)*bcVals4R2+((48.0*rdx2SqVol[0]*phiLy[1]+48.0*rdx2SqVol[0]*phiLx[1]+((-384.0*phiPrevC[0])+192.0*phiLy[0]+48.0*phiLxLy[0]-96.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(48.0*rdx2SqVol0R2*phiLy[1]+48.0*rdx2SqVol0R2*phiLx[1]+288.0*rdx2SqVol[0]*rhoC[0]+((-384.0*phiPrevC[0])-96.0*phiLy[0]+48.0*phiLxLy[0]+192.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2)*omega+((49.0*phiPrevC[0]*rdx2SqVol1R4+700.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+700.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1]+49.0*phiPrevC[0]*rdx2SqVol0R4)*bcVals4R2+(448.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1]+280.0*phiPrevC[0]*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+960.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1]+336.0*phiPrevC[0]*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((280.0*phiPrevC[0]*rdx2SqVol1R4+2592.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+448.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1104.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2+384.0*phiPrevC[0]*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*phiPrevC[0]*rdx2SqVol1R4+960.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+192.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*phiPrevC[0]*rdx2SqVol[0]*rdx2SqVol1R3+384.0*phiPrevC[0]*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2)/(((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2); 
  phiC[1] = (((((126.0*rdx2SqVol1R4+1206.0*rdx2SqVol[0]*rdx2SqVol1R3+954.0*rdx2SqVol0R2*rdx2SqVol1R2-126.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+((-144.0*rdx2SqVol[0]*rdx2SqVol1R3)-216.0*rdx2SqVol0R2*rdx2SqVol1R2-504.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4])*bcVals[10]+((168.0*rdx2SqVol1R4+552.0*rdx2SqVol[0]*rdx2SqVol1R3-48.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+((-192.0*rdx2SqVol[0]*rdx2SqVol1R3)-192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals[9])*bcVals[11]+(((336.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+1584.0*rdx2SqVol0R3*rdx2SqVol[1]+84.0*rdx2SqVol0R4)*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+((168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[3]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[1]+((-49.0*phiPrevC[1])+49.0*phiLy[1]-7.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol1R4+((-700.0*rdx2SqVol[0]*phiPrevC[1])+385.0*rdx2SqVol[0]*phiLy[1]-72.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(378.0*phiLy[0]+29.0*phiLxLy[0]-20.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-1302.0*rdx2SqVol0R2*phiPrevC[1])+285.0*rdx2SqVol0R2*phiLy[1]+180.0*rdx2SqVol0R2*phiLx[1]+216.0*rdx2SqVol[0]*rhoC[0]+(540.0*phiLy[0]+93.0*phiLxLy[0]+204.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-700.0*rdx2SqVol0R3*phiPrevC[1])-65.0*rdx2SqVol0R3*phiLy[1]+252.0*rdx2SqVol0R3*phiLx[1]+342.0*rdx2SqVol0R2*rhoC[0]+(162.0*phiLy[0]+71.0*phiLxLy[0]+280.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-49.0*rdx2SqVol0R4*phiPrevC[1]-14.0*rdx2SqVol0R4*phiLy[1]+84.0*rdx2SqVol0R3*rhoC[0]+(14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R4)*bcVals4R2+(((-192.0*rdx2SqVol[0]*rdx2SqVol1R2)-960.0*rdx2SqVol0R2*rdx2SqVol[1]-336.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((48.0*rdx2SqVol[0]*rdx2SqVol1R2+744.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*rhoC[2]+(384.0*rdx2SqVol[0]*rdx2SqVol1R2+1920.0*rdx2SqVol0R2*rdx2SqVol[1]+672.0*rdx2SqVol0R3)*rhoC[1]+((-448.0*rdx2SqVol[0]*phiPrevC[1])+112.0*rdx2SqVol[0]*phiLy[1]-24.0*rdx2SqVol[0]*phiLx[1]+(40.0*phiLx[0]-16.0*phiLxLy[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-2760.0*rdx2SqVol0R2*phiPrevC[1])+564.0*rdx2SqVol0R2*phiLy[1]-120.0*rdx2SqVol0R2*phiLx[1]-96.0*rdx2SqVol[0]*rhoC[0]+(432.0*phiLy[0]-12.0*phiLxLy[0]+60.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-2592.0*rdx2SqVol0R3*phiPrevC[1])-36.0*rdx2SqVol0R3*phiLy[1]+336.0*rdx2SqVol0R3*phiLx[1]+24.0*rdx2SqVol0R2*rhoC[0]+(648.0*phiLy[0]+60.0*phiLxLy[0])*rdx2SqVol0R3)*rdx2SqVol[1]-280.0*rdx2SqVol0R4*phiPrevC[1]-56.0*rdx2SqVol0R4*phiLy[1]+336.0*rdx2SqVol0R3*rhoC[0]+(56.0*phiLxLy[0]+196.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((-192.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1]-336.0*rdx2SqVol0R4*phiPrevC[1])*bcVals3R2)*bcVals10R2+(((984.0*rdx2SqVol[0]*rdx2SqVol1R3+2688.0*rdx2SqVol0R2*rdx2SqVol1R2+1272.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3])*bcVals[5]+((336.0*rdx2SqVol1R3-192.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+((-168.0*rdx2SqVol1R3)+744.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(648.0*rdx2SqVol1R3+2880.0*rdx2SqVol[0]*rdx2SqVol1R2+1368.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+((-280.0*phiPrevC[1])+182.0*phiLy[1]-6.0*phiLx[1]-28.0*phiLy[0]-34.0*phiLxLy[0]+40.0*phiLx[0])*rdx2SqVol1R4+((-2592.0*rdx2SqVol[0]*phiPrevC[1])+858.0*rdx2SqVol[0]*phiLy[1]-174.0*rdx2SqVol[0]*phiLx[1]-204.0*rhoC[0]+(816.0*phiLy[0]+6.0*phiLxLy[0]-120.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-2760.0*rdx2SqVol0R2*phiPrevC[1])+126.0*rdx2SqVol0R2*phiLy[1]+390.0*rdx2SqVol0R2*phiLx[1]+240.0*rdx2SqVol[0]*rhoC[0]+(1068.0*phiLy[0]+150.0*phiLxLy[0]+72.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-448.0*rdx2SqVol0R3*phiPrevC[1])-118.0*rdx2SqVol0R3*phiLy[1]+126.0*rdx2SqVol0R3*phiLx[1]+660.0*rdx2SqVol0R2*rhoC[0]+(8.0*phiLy[0]+110.0*phiLxLy[0]+448.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals4R2+(((-384.0*rdx2SqVol[0]*rdx2SqVol1R2)-384.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[3]*rhoC[3]+((192.0*rdx2SqVol[0]*rdx2SqVol1R2+192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(768.0*rdx2SqVol[0]*rdx2SqVol1R2+768.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+((-1104.0*rdx2SqVol[0]*phiPrevC[1])+232.0*rdx2SqVol[0]*phiLy[1]-56.0*rdx2SqVol[0]*phiLx[1]+(32.0*phiLy[0]-24.0*phiLxLy[0]+80.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-3072.0*rdx2SqVol0R2*phiPrevC[1])+152.0*rdx2SqVol0R2*phiLy[1]+56.0*rdx2SqVol0R2*phiLx[1]-144.0*rdx2SqVol[0]*rhoC[0]+(496.0*phiLy[0]+24.0*phiLxLy[0]-128.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-1104.0*rdx2SqVol0R3*phiPrevC[1])-80.0*rdx2SqVol0R3*phiLy[1]+112.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+(32.0*phiLy[0]+48.0*phiLxLy[0]+224.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[3])*bcVals[4]+((-384.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2)-384.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((576.0*rdx2SqVol[0]*rdx2SqVol1R3+576.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[4]*bcVals[5]+((576.0*rdx2SqVol1R3+576.0*rdx2SqVol[0]*rdx2SqVol1R2)*rhoC[1]+((-336.0*phiPrevC[1])+168.0*phiLy[1]-24.0*phiLx[1]-24.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol1R4+((-960.0*rdx2SqVol[0]*phiPrevC[1])+120.0*rdx2SqVol[0]*phiLy[1]+24.0*rdx2SqVol[0]*phiLx[1]-144.0*rhoC[0]+(432.0*phiLy[0]+24.0*phiLxLy[0]-192.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-192.0*rdx2SqVol0R2*phiPrevC[1])-48.0*rdx2SqVol0R2*phiLy[1]+48.0*rdx2SqVol0R2*phiLx[1]+288.0*rdx2SqVol[0]*rhoC[0]+(48.0*phiLxLy[0]+192.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2)*bcVals4R2+((-384.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3)-384.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2)*omega+((49.0*phiPrevC[1]*rdx2SqVol1R4+700.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+700.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiPrevC[1])*bcVals4R2+(448.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1]+280.0*rdx2SqVol0R4*phiPrevC[1])*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+960.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1]+336.0*rdx2SqVol0R4*phiPrevC[1])*bcVals3R2)*bcVals10R2+((280.0*phiPrevC[1]*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+448.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1])*bcVals4R2+(1104.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2+384.0*rdx2SqVol0R3*phiPrevC[1]*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*phiPrevC[1]*rdx2SqVol1R4+960.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*phiPrevC[1]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*phiPrevC[1]*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2)/(((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2); 
  phiC[2] = (((((84.0*rdx2SqVol1R4+1584.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+336.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1272.0*rdx2SqVol[0]*rdx2SqVol1R3+2688.0*rdx2SqVol0R2*rdx2SqVol1R2+984.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(576.0*rdx2SqVol0R2*rdx2SqVol1R2+576.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals[9])*bcVals[11]+((((-126.0*rdx2SqVol[0]*rdx2SqVol1R3)+954.0*rdx2SqVol0R2*rdx2SqVol1R2+1206.0*rdx2SqVol0R3*rdx2SqVol[1]+126.0*rdx2SqVol0R4)*bcVals[4]+((-48.0*rdx2SqVol0R2*rdx2SqVol1R2)+552.0*rdx2SqVol0R3*rdx2SqVol[1]+168.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+(((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[3]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiPrevC[2]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[1]+((-14.0*phiLx[1])+49.0*phiLy[0]+14.0*phiLxLy[0])*rdx2SqVol1R4+(252.0*rdx2SqVol[0]*phiLy[1]-65.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(280.0*phiLy[0]+71.0*phiLxLy[0]+162.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(180.0*rdx2SqVol0R2*phiLy[1]+285.0*rdx2SqVol0R2*phiLx[1]+342.0*rdx2SqVol[0]*rhoC[0]+(204.0*phiLy[0]+93.0*phiLxLy[0]+540.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-72.0*rdx2SqVol0R3*phiLy[1])+385.0*rdx2SqVol0R3*phiLx[1]+216.0*rdx2SqVol0R2*rhoC[0]+((-20.0*phiLy[0])+29.0*phiLxLy[0]+378.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]+(7.0*phiLy[0]-7.0*phiLxLy[0])*rdx2SqVol0R4)*bcVals4R2+(((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)-192.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((1368.0*rdx2SqVol[0]*rdx2SqVol1R2+2880.0*rdx2SqVol0R2*rdx2SqVol[1]+648.0*rdx2SqVol0R3)*rhoC[2]+((-448.0*rdx2SqVol[0]*rdx2SqVol1R3)-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-2592.0*rdx2SqVol0R3*rdx2SqVol[1]-280.0*rdx2SqVol0R4)*phiPrevC[2]+(48.0*rdx2SqVol[0]*rdx2SqVol1R2+744.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*rhoC[1]+(126.0*rdx2SqVol[0]*phiLy[1]-118.0*rdx2SqVol[0]*phiLx[1]+(448.0*phiLy[0]+110.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(390.0*rdx2SqVol0R2*phiLy[1]+126.0*rdx2SqVol0R2*phiLx[1]+660.0*rdx2SqVol[0]*rhoC[0]+(72.0*phiLy[0]+150.0*phiLxLy[0]+1068.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-174.0*rdx2SqVol0R3*phiLy[1])+858.0*rdx2SqVol0R3*phiLx[1]+240.0*rdx2SqVol0R2*rhoC[0]+((-120.0*phiLy[0])+6.0*phiLxLy[0]+816.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-6.0*rdx2SqVol0R4*phiLy[1]+182.0*rdx2SqVol0R4*phiLx[1]-204.0*rdx2SqVol0R3*rhoC[0]+(40.0*phiLy[0]-34.0*phiLxLy[0]-28.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((576.0*rdx2SqVol0R2*rdx2SqVol[1]+576.0*rdx2SqVol0R3)*rhoC[2]+((-192.0*rdx2SqVol0R2*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*rdx2SqVol[1]-336.0*rdx2SqVol0R4)*phiPrevC[2]+(48.0*rdx2SqVol0R2*phiLy[1]-48.0*rdx2SqVol0R2*phiLx[1]+(192.0*phiLy[0]+48.0*phiLxLy[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(24.0*rdx2SqVol0R3*phiLy[1]+120.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+((-192.0*phiLy[0])+24.0*phiLxLy[0]+432.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-24.0*rdx2SqVol0R4*phiLy[1]+168.0*rdx2SqVol0R4*phiLx[1]-144.0*rdx2SqVol0R3*rhoC[0]+(48.0*phiLy[0]-24.0*phiLxLy[0])*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((((-504.0*rdx2SqVol[0]*rdx2SqVol1R3)-216.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[4]+((-192.0*rdx2SqVol0R2*rdx2SqVol1R2)-192.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3])*bcVals[5]+(((-336.0*rdx2SqVol1R3)-960.0*rdx2SqVol[0]*rdx2SqVol1R2-192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+(672.0*rdx2SqVol1R3+1920.0*rdx2SqVol[0]*rdx2SqVol1R2+384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-280.0*rdx2SqVol1R4)-2592.0*rdx2SqVol[0]*rdx2SqVol1R3-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[2]+((-168.0*rdx2SqVol1R3)+744.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+((-56.0*phiLx[1])+196.0*phiLy[0]+56.0*phiLxLy[0])*rdx2SqVol1R4+(336.0*rdx2SqVol[0]*phiLy[1]-36.0*rdx2SqVol[0]*phiLx[1]+336.0*rhoC[0]+(60.0*phiLxLy[0]+648.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-120.0*rdx2SqVol0R2*phiLy[1])+564.0*rdx2SqVol0R2*phiLx[1]+24.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLy[0]-12.0*phiLxLy[0]+432.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-24.0*rdx2SqVol0R3*phiLy[1])+112.0*rdx2SqVol0R3*phiLx[1]-96.0*rdx2SqVol0R2*rhoC[0]+(40.0*phiLy[0]-16.0*phiLxLy[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals4R2+(((-384.0*rdx2SqVol[0]*rdx2SqVol1R2)-384.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[3]*rhoC[3]+((768.0*rdx2SqVol[0]*rdx2SqVol1R2+768.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-1104.0*rdx2SqVol[0]*rdx2SqVol1R3)-3072.0*rdx2SqVol0R2*rdx2SqVol1R2-1104.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[2]+(192.0*rdx2SqVol[0]*rdx2SqVol1R2+192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(112.0*rdx2SqVol[0]*phiLy[1]-80.0*rdx2SqVol[0]*phiLx[1]+(224.0*phiLy[0]+48.0*phiLxLy[0]+32.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(56.0*rdx2SqVol0R2*phiLy[1]+152.0*rdx2SqVol0R2*phiLx[1]+288.0*rdx2SqVol[0]*rhoC[0]+((-128.0*phiLy[0])+24.0*phiLxLy[0]+496.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-56.0*rdx2SqVol0R3*phiLy[1])+232.0*rdx2SqVol0R3*phiLx[1]-144.0*rdx2SqVol0R2*rhoC[0]+(80.0*phiLy[0]-24.0*phiLxLy[0]+32.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[3])*bcVals[4]+((-384.0*rdx2SqVol0R2*rdx2SqVol1R2)-384.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[2]*bcVals3R2)*bcVals[9]*bcVals[10]+(((-336.0*rdx2SqVol1R4)-960.0*rdx2SqVol[0]*rdx2SqVol1R3-192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiPrevC[2]*bcVals4R2+((-384.0*rdx2SqVol[0]*rdx2SqVol1R3)-384.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiPrevC[2]*bcVals[3]*bcVals[4])*bcVals9R2)*omega+((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiPrevC[2]*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*phiPrevC[2]*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*phiPrevC[2]*bcVals3R2)*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[2]*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[2]*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[2]*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiPrevC[2]*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiPrevC[2]*bcVals[3]*bcVals[4])*bcVals9R2)/(((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2); 
  phiC[3] = (((((252.0*rdx2SqVol1R4+2736.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+504.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(288.0*rdx2SqVol[0]*rdx2SqVol1R3+1728.0*rdx2SqVol0R2*rdx2SqVol1R2+1008.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4])*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals[9])*bcVals[11]+(((504.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+2736.0*rdx2SqVol0R3*rdx2SqVol[1]+252.0*rdx2SqVol0R4)*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+((336.0*rdx2SqVol1R3+2304.0*rdx2SqVol[0]*rdx2SqVol1R2+2304.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[3]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiPrevC[3]+((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[1]+7.0*phiLx[1]-7.0*phiLxLy[0])*rdx2SqVol1R4+(280.0*rdx2SqVol[0]*phiLy[1]-20.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(504.0*phiLy[0]+80.0*phiLxLy[0]-144.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(204.0*rdx2SqVol0R2*phiLy[1]+204.0*rdx2SqVol0R2*phiLx[1]+522.0*rdx2SqVol[0]*rhoC[0]+(360.0*phiLy[0]+174.0*phiLxLy[0]+360.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-20.0*rdx2SqVol0R3*phiLy[1])+280.0*rdx2SqVol0R3*phiLx[1]+522.0*rdx2SqVol0R2*rhoC[0]+((-144.0*phiLy[0])+80.0*phiLxLy[0]+504.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+7.0*rdx2SqVol0R4*phiLy[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]-7.0*phiLxLy[0]*rdx2SqVol0R4)*bcVals4R2+((384.0*rdx2SqVol[0]*rdx2SqVol1R2+1920.0*rdx2SqVol0R2*rdx2SqVol[1]+672.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((-448.0*rdx2SqVol[0]*rdx2SqVol1R3)-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-2592.0*rdx2SqVol0R3*rdx2SqVol[1]-280.0*rdx2SqVol0R4)*bcVals[3]*phiPrevC[3]+(((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)-192.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[2]+((-192.0*rdx2SqVol[0]*rdx2SqVol1R2)-960.0*rdx2SqVol0R2*rdx2SqVol[1]-336.0*rdx2SqVol0R3)*rhoC[1]+((-56.0*rdx2SqVol[0]*phiLy[1])+24.0*rdx2SqVol[0]*phiLx[1]+(8.0*phiLxLy[0]-32.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-228.0*rdx2SqVol0R2*phiLy[1])+60.0*rdx2SqVol0R2*phiLx[1]+48.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLxLy[0]-120.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(72.0*rdx2SqVol0R3*phiLy[1]-96.0*rdx2SqVol0R3*phiLx[1]+312.0*rdx2SqVol0R2*rhoC[0]+((-432.0*phiLy[0])+24.0*phiLxLy[0]+288.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+28.0*rdx2SqVol0R4*phiLy[1]+84.0*rdx2SqVol0R4*phiLx[1]-168.0*rdx2SqVol0R3*rhoC[0]+((-28.0*phiLxLy[0])-56.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((-192.0*rdx2SqVol0R2*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*rdx2SqVol[1]-336.0*rdx2SqVol0R4)*bcVals3R2*phiPrevC[3])*bcVals10R2+(((1008.0*rdx2SqVol[0]*rdx2SqVol1R3+1728.0*rdx2SqVol0R2*rdx2SqVol1R2+288.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3])*bcVals[5]+((672.0*rdx2SqVol1R3+1920.0*rdx2SqVol[0]*rdx2SqVol1R2+384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+((-280.0*rdx2SqVol1R4)-2592.0*rdx2SqVol[0]*rdx2SqVol1R3-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[3]+((-336.0*rdx2SqVol1R3)-960.0*rdx2SqVol[0]*rdx2SqVol1R2-192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(336.0*rdx2SqVol1R3-192.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(84.0*phiLy[1]+28.0*phiLx[1]-56.0*phiLy[0]-28.0*phiLxLy[0])*rdx2SqVol1R4+((-96.0*rdx2SqVol[0]*phiLy[1])+72.0*rdx2SqVol[0]*phiLx[1]-168.0*rhoC[0]+(288.0*phiLy[0]+24.0*phiLxLy[0]-432.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(60.0*rdx2SqVol0R2*phiLy[1]-228.0*rdx2SqVol0R2*phiLx[1]+312.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLxLy[0]-120.0*phiLy[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(24.0*rdx2SqVol0R3*phiLy[1]-56.0*rdx2SqVol0R3*phiLx[1]+48.0*rdx2SqVol0R2*rhoC[0]+(8.0*phiLxLy[0]-32.0*phiLy[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals4R2+((768.0*rdx2SqVol[0]*rdx2SqVol1R2+768.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[3]*rhoC[3]+((-1104.0*rdx2SqVol[0]*rdx2SqVol1R3)-3072.0*rdx2SqVol0R2*rdx2SqVol1R2-1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*phiPrevC[3]+(((-384.0*rdx2SqVol[0]*rdx2SqVol1R2)-384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-384.0*rdx2SqVol[0]*rdx2SqVol1R2)-384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+((-128.0*rdx2SqVol[0]*phiLy[1])+64.0*rdx2SqVol[0]*phiLx[1]+((-64.0*phiLy[0])-64.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-64.0*rdx2SqVol0R2*phiLy[1])-64.0*rdx2SqVol0R2*phiLx[1]+((-128.0*phiLy[0])-128.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(64.0*rdx2SqVol0R3*phiLy[1]-128.0*rdx2SqVol0R3*phiLx[1]+((-64.0*phiLy[0])-64.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[3])*bcVals[4]+((-384.0*rdx2SqVol0R2*rdx2SqVol1R2)-384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2*phiPrevC[3])*bcVals[9]*bcVals[10]+(((-336.0*rdx2SqVol1R4)-960.0*rdx2SqVol[0]*rdx2SqVol1R3-192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiPrevC[3]*bcVals4R2+((-384.0*rdx2SqVol[0]*rdx2SqVol1R3)-384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*phiPrevC[3]*bcVals[4])*bcVals9R2)*omega+((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiPrevC[3]*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*phiPrevC[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2*phiPrevC[3])*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiPrevC[3]*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*phiPrevC[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2*phiPrevC[3])*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiPrevC[3]*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*phiPrevC[3]*bcVals[4])*bcVals9R2)/(((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = (((4.0*phiUy[0]+phiUxUy[0]+phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+6.0*rhoC[0]+((-2.0*phiUy[0])+phiUxUy[0]+phiUxLy[0]+4.0*phiUx[0]-2.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0])*omega+8.0*phiC[0]*rdx2SqVol[1]+8.0*phiC[0]*rdx2SqVol[0])/(8.0*rdx2SqVol[1]+8.0*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]-1.0*phiUxLy[0]+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(phiUy[0]-1.0*phiUxUy[0]-1.0*phiUxLy[0]-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0])*omega-4.0*phiC[0]*rdx2SqVol[1]-4.0*phiC[0]*rdx2SqVol[0]))/(4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]-1.0*phiUxLy[0]+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*bcVals[1]*rdx2SqVol[1]+((phiUy[0]-1.0*phiUxUy[0]-1.0*phiUxLy[0]-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiC[0])*rdx2SqVol[0])*omega-4.0*phiC[0]*bcVals[1]*rdx2SqVol[1]-4.0*phiC[0]*rdx2SqVol[0]*bcVals[1]+4.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/(4.0*bcVals[1]*rdx2SqVol[1]+4.0*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*(((rdx2SqVol[1]-5.0*rdx2SqVol[0])*bcVals[5]+((-1.0*phiLy[1])-4.0*phiUy[0]-4.0*phiLy[0]-1.0*phiLxUy[0]-1.0*phiLxLy[0]+2.0*phiLx[0]+8.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-6.0*rhoC[0]+(2.0*phiUy[0]+2.0*phiLy[0]-1.0*phiLxUy[0]-1.0*phiLxLy[0]-4.0*phiLx[0]+8.0*phiC[0])*rdx2SqVol[0])*omega-8.0*phiC[0]*rdx2SqVol[1]-8.0*phiC[0]*rdx2SqVol[0]))/(8.0*rdx2SqVol[1]+8.0*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = -(1.0*(((6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+(6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])-7.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]+4.0*phiLx[0]+14.0*phiC[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-9.0*rdx2SqVol[0]*phiLy[1]-12.0*rhoC[0]+((-5.0*phiUy[0])-5.0*phiLy[0]-4.0*phiLxUy[0]-4.0*phiLxLy[0]-4.0*phiLx[0]+40.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]+2.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]-8.0*phiLx[0]+8.0*phiC[0])*rdx2SqVol0R2)*omega-14.0*phiC[0]*rdx2SqVol1R2-40.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]-8.0*phiC[0]*rdx2SqVol0R2))/(14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2); 
  phiC[1] = (((24.0*rdx2SqVol[0]*rdx2SqVol[1]+24.0*rdx2SqVol0R2)*bcVals[5]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+(7.0*phiUy[1]+7.0*phiLy[1]-14.0*phiC[1]-1.0*phiLxUy[0]-1.0*phiLxLy[0]+2.0*phiLx[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUy[1]+5.0*rdx2SqVol[0]*phiLy[1]-40.0*rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+(18.0*phiUy[0]+18.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-2.0*rdx2SqVol0R2*phiUy[1]-2.0*rdx2SqVol0R2*phiLy[1]-8.0*rdx2SqVol0R2*phiC[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLxUy[0]+2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*omega+14.0*phiC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiC[1])/(14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = -(1.0*(((6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+((6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])-7.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]+4.0*phiLx[0]+14.0*phiC[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-9.0*rdx2SqVol[0]*phiLy[1]-12.0*rhoC[0]+((-5.0*phiUy[0])-5.0*phiLy[0]-4.0*phiLxUy[0]-4.0*phiLxLy[0]-4.0*phiLx[0]+40.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]+2.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]-8.0*phiLx[0]+8.0*phiC[0])*rdx2SqVol0R2)*bcVals[4]+(((-4.0*rdx2SqVol[0]*phiUy[1])-2.0*rdx2SqVol[0]*phiLy[1]+((-8.0*phiUy[0])-8.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]+4.0*phiLx[0]+16.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiUy[1]-2.0*rdx2SqVol0R2*phiLy[1]-12.0*rdx2SqVol[0]*rhoC[0]+(4.0*phiUy[0]+4.0*phiLy[0]-2.0*phiLxUy[0]-2.0*phiLxLy[0]-8.0*phiLx[0]+16.0*phiC[0])*rdx2SqVol0R2)*bcVals[3])*omega+((-14.0*phiC[0]*rdx2SqVol1R2)-40.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]-8.0*phiC[0]*rdx2SqVol0R2)*bcVals[4]+((-16.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1])-16.0*phiC[0]*rdx2SqVol0R2)*bcVals[3]))/((14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 
  phiC[1] = (((24.0*rdx2SqVol[0]*rdx2SqVol[1]+24.0*rdx2SqVol0R2)*bcVals[5]+((24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+(7.0*phiUy[1]+7.0*phiLy[1]-14.0*phiC[1]-1.0*phiLxUy[0]-1.0*phiLxLy[0]+2.0*phiLx[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUy[1]+5.0*rdx2SqVol[0]*phiLy[1]-40.0*rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+(18.0*phiUy[0]+18.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-2.0*rdx2SqVol0R2*phiUy[1]-2.0*rdx2SqVol0R2*phiLy[1]-8.0*rdx2SqVol0R2*phiC[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLxUy[0]+2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+(((-8.0*rdx2SqVol[0]*phiUy[1])-16.0*rdx2SqVol[0]*phiC[1])*rdx2SqVol[1]-8.0*rdx2SqVol0R2*phiUy[1]-16.0*rdx2SqVol0R2*phiC[1])*bcVals[3])*omega+(14.0*phiC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiC[1])*bcVals[4]+(16.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+16.0*rdx2SqVol0R2*phiC[1])*bcVals[3])/((14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]+((-4.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*omega-4.0*phiC[0]*rdx2SqVol[1]-4.0*phiC[0]*rdx2SqVol[0]))/(4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]+(((-4.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6])*omega+((-4.0*phiC[0]*rdx2SqVol[1])-4.0*phiC[0]*rdx2SqVol[0])*bcVals[7]+4.0*phiC[0]*rdx2SqVol[1]*bcVals[6]))/((4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0])*bcVals[7]-4.0*rdx2SqVol[1]*bcVals[6]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = (((5.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*bcVals[11]+(phiLx[1]+phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLx[1]+6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-2.0*phiLy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0])*omega+8.0*phiC[0]*rdx2SqVol[1]+8.0*phiC[0]*rdx2SqVol[0])/(8.0*rdx2SqVol[1]+8.0*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]+8.0*phiLy[0]+2.0*phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiUx[1]+9.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]+4.0*phiLy[0]+4.0*phiLxLy[0]+5.0*phiLx[0]-40.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-4.0*phiLy[0]+2.0*phiLxLy[0]+7.0*phiLx[0]-14.0*phiC[0])*rdx2SqVol0R2)*omega+8.0*phiC[0]*rdx2SqVol1R2+40.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+14.0*phiC[0]*rdx2SqVol0R2)/(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2); 
  phiC[1] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+((-2.0*phiUx[1])-2.0*phiLx[1]-8.0*phiC[1]+2.0*phiUxLy[0]+8.0*phiLy[0]+2.0*phiLxLy[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUx[1]+5.0*rdx2SqVol[0]*phiLx[1]-40.0*rdx2SqVol[0]*phiC[1]+12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-8.0*phiLy[0]+phiLxLy[0]+18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiUx[1]+7.0*rdx2SqVol0R2*phiLx[1]-14.0*rdx2SqVol0R2*phiC[1]-6.0*rdx2SqVol[0]*rhoC[0]+((-1.0*phiUxLy[0])+2.0*phiLy[0]-1.0*phiLxLy[0])*rdx2SqVol0R2)*omega+8.0*phiC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+14.0*rdx2SqVol0R2*phiC[1])/(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]+8.0*phiLy[0]+2.0*phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiUx[1]+9.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]+4.0*phiLy[0]+4.0*phiLxLy[0]+5.0*phiLx[0]-40.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-4.0*phiLy[0]+2.0*phiLxLy[0]+7.0*phiLx[0]-14.0*phiC[0])*rdx2SqVol0R2)*bcVals[10]+(((-2.0*phiUx[1])+2.0*phiLx[1]+2.0*phiUxLy[0]-4.0*phiUx[0]+8.0*phiLy[0]+2.0*phiLxLy[0]-4.0*phiLx[0]-16.0*phiC[0])*rdx2SqVol1R2+(4.0*rdx2SqVol[0]*phiUx[1]+2.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(2.0*phiUxLy[0]+8.0*phiUx[0]-4.0*phiLy[0]+2.0*phiLxLy[0]+8.0*phiLx[0]-16.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[9])*omega+(8.0*phiC[0]*rdx2SqVol1R2+40.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+14.0*phiC[0]*rdx2SqVol0R2)*bcVals[10]+(16.0*phiC[0]*rdx2SqVol1R2+16.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9])/((8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 
  phiC[1] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+((-2.0*phiUx[1])-2.0*phiLx[1]-8.0*phiC[1]+2.0*phiUxLy[0]+8.0*phiLy[0]+2.0*phiLxLy[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUx[1]+5.0*rdx2SqVol[0]*phiLx[1]-40.0*rdx2SqVol[0]*phiC[1]+12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-8.0*phiLy[0]+phiLxLy[0]+18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiUx[1]+7.0*rdx2SqVol0R2*phiLx[1]-14.0*rdx2SqVol0R2*phiC[1]-6.0*rdx2SqVol[0]*rhoC[0]+((-1.0*phiUxLy[0])+2.0*phiLy[0]-1.0*phiLxLy[0])*rdx2SqVol0R2)*bcVals[10]+(((-8.0*phiUx[1])-16.0*phiC[1])*rdx2SqVol1R2+((-8.0*rdx2SqVol[0]*phiUx[1])-16.0*rdx2SqVol[0]*phiC[1])*rdx2SqVol[1])*bcVals[9])*omega+(8.0*phiC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+14.0*rdx2SqVol0R2*phiC[1])*bcVals[10]+(16.0*phiC[1]*rdx2SqVol1R2+16.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1])*bcVals[9])/((8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]+6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0])*omega-2.0*phiC[0]*rdx2SqVol[1]-2.0*phiC[0]*rdx2SqVol[0]))/(2.0*rdx2SqVol[1]+2.0*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]+(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6])*omega+((-2.0*phiC[0]*rdx2SqVol[1])-2.0*phiC[0]*rdx2SqVol[0])*bcVals[7]+4.0*phiC[0]*rdx2SqVol[1]*bcVals[6]))/((2.0*rdx2SqVol[1]+2.0*rdx2SqVol[0])*bcVals[7]-4.0*rdx2SqVol[1]*bcVals[6]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*((6.0*bcVals[1]*rdx2SqVol[1]*bcVals[8]+6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*bcVals[1]*rdx2SqVol[1]+((phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiC[0])*rdx2SqVol[0])*omega-2.0*phiC[0]*bcVals[1]*rdx2SqVol[1]-2.0*phiC[0]*rdx2SqVol[0]*bcVals[1]+4.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/(2.0*bcVals[1]*rdx2SqVol[1]+2.0*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*((6.0*bcVals[1]*rdx2SqVol[1]*bcVals[8]+(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*bcVals[1]*rdx2SqVol[1]+((phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiC[0])*bcVals[1]*rdx2SqVol[1]*bcVals[6])*omega+((-2.0*phiC[0]*bcVals[1]*rdx2SqVol[1])-2.0*phiC[0]*rdx2SqVol[0]*bcVals[1]+4.0*bcVals[0]*phiC[0]*rdx2SqVol[0])*bcVals[7]+4.0*phiC[0]*bcVals[1]*rdx2SqVol[1]*bcVals[6]))/((2.0*bcVals[1]*rdx2SqVol[1]+2.0*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*rdx2SqVol[0])*bcVals[7]-4.0*bcVals[1]*rdx2SqVol[1]*bcVals[6]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = ((3.0*rdx2SqVol[1]*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+(phiUxLy[0]-2.0*phiUx[0]+2.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol[1]+6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-1.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol[0])*omega+4.0*phiC[0]*rdx2SqVol[1]+4.0*phiC[0]*rdx2SqVol[0])/(4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-6.0*rdx2SqVol0R2)*bcVals[2]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiUx[1]+12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]+2.0*phiLy[0]-20.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-2.0*phiLy[0]-7.0*phiC[0])*rdx2SqVol0R2)*omega+4.0*phiC[0]*rdx2SqVol1R2+20.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+7.0*phiC[0]*rdx2SqVol0R2)/(4.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2); 
  phiC[1] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((-36.0*rdx2SqVol[0]*rdx2SqVol[1])-18.0*rdx2SqVol0R2)*bcVals[2]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+((-2.0*phiUx[1])-4.0*phiC[1]+2.0*phiUxLy[0]+4.0*phiLy[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUx[1]-20.0*rdx2SqVol[0]*phiC[1]+12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-4.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiUx[1]-7.0*rdx2SqVol0R2*phiC[1]-6.0*rdx2SqVol[0]*rhoC[0]+(phiLy[0]-1.0*phiUxLy[0])*rdx2SqVol0R2)*omega+4.0*phiC[1]*rdx2SqVol1R2+20.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiC[1])/(4.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-6.0*rdx2SqVol0R2)*bcVals[2]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiUx[1]+12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]+2.0*phiLy[0]-20.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-2.0*phiLy[0]-7.0*phiC[0])*rdx2SqVol0R2)*bcVals[10]+((-24.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[2])+(4.0*phiUxLy[0]-8.0*phiUx[0]+8.0*phiLy[0]-16.0*phiC[0])*rdx2SqVol1R2+(6.0*rdx2SqVol[0]*phiUx[1]+24.0*rhoC[0]+(4.0*phiUxLy[0]+16.0*phiUx[0]-4.0*phiLy[0]-16.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[9])*omega+(4.0*phiC[0]*rdx2SqVol1R2+20.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+7.0*phiC[0]*rdx2SqVol0R2)*bcVals[10]+(16.0*phiC[0]*rdx2SqVol1R2+16.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9])/((4.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 
  phiC[1] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(((-36.0*rdx2SqVol[0]*rdx2SqVol[1])-18.0*rdx2SqVol0R2)*bcVals[2]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+((-2.0*phiUx[1])-4.0*phiC[1]+2.0*phiUxLy[0]+4.0*phiLy[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiUx[1]-20.0*rdx2SqVol[0]*phiC[1]+12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-4.0*phiLy[0])*rdx2SqVol[0])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiUx[1]-7.0*rdx2SqVol0R2*phiC[1]-6.0*rdx2SqVol[0]*rhoC[0]+(phiLy[0]-1.0*phiUxLy[0])*rdx2SqVol0R2)*bcVals[10]+(((-8.0*phiUx[1])-16.0*phiC[1])*rdx2SqVol1R2+((-8.0*rdx2SqVol[0]*phiUx[1])-16.0*rdx2SqVol[0]*phiC[1])*rdx2SqVol[1])*bcVals[9])*omega+(4.0*phiC[1]*rdx2SqVol1R2+20.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+7.0*rdx2SqVol0R2*phiC[1])*bcVals[10]+(16.0*phiC[1]*rdx2SqVol1R2+16.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1])*bcVals[9])/((4.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = (((3.0*bcVals[1]*rdx2SqVol[1]+2.0*bcVals[0]*rdx2SqVol[0])*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+(phiUxLy[0]-2.0*phiUx[0]+2.0*phiLy[0]-4.0*phiC[0])*bcVals[1]*rdx2SqVol[1]+(6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-1.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol[0])*bcVals[1]+4.0*bcVals[0]*phiC[0]*rdx2SqVol[0])*omega+4.0*phiC[0]*bcVals[1]*rdx2SqVol[1]+4.0*phiC[0]*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*phiC[0]*rdx2SqVol[0])/(4.0*bcVals[1]*rdx2SqVol[1]+4.0*rdx2SqVol[0]*bcVals[1]-4.0*bcVals[0]*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);

  phiC[0] = (((12.0*bcVals1R2*rdx2SqVol1R2+(12.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-6.0*rdx2SqVol[0]*bcVals1R2)*rdx2SqVol[1])*bcVals[11]+((-24.0*rdx2SqVol[0]*bcVals[1]*rdx2SqVol[1])-6.0*rdx2SqVol0R2*bcVals[1]+12.0*bcVals[0]*rdx2SqVol0R2)*bcVals[2]+(12.0*bcVals1R2*rdx2SqVol[1]-6.0*rdx2SqVol[0]*bcVals1R2+12.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]-4.0*phiC[0])*bcVals1R2*rdx2SqVol1R2+((9.0*rdx2SqVol[0]*bcVals1R2-6.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]+2.0*phiLy[0]-20.0*phiC[0])*rdx2SqVol[0])*bcVals1R2+((-4.0*bcVals[0]*phiUxLy[0])+10.0*bcVals[0]*phiUx[0]-8.0*bcVals[0]*phiLy[0]+32.0*bcVals[0]*phiC[0])*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+(12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-2.0*phiLy[0]-7.0*phiC[0])*rdx2SqVol0R2)*bcVals1R2+(((-4.0*bcVals[0]*phiUxLy[0])-14.0*bcVals[0]*phiUx[0]+4.0*bcVals[0]*phiLy[0]+20.0*bcVals[0]*phiC[0])*rdx2SqVol0R2-24.0*bcVals[0]*rdx2SqVol[0]*rhoC[0])*bcVals[1]-12.0*bcVals0R2*phiC[0]*rdx2SqVol0R2)*omega+4.0*phiC[0]*bcVals1R2*rdx2SqVol1R2+(20.0*phiC[0]*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*phiC[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*phiC[0]*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*phiC[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*phiC[0]*rdx2SqVol0R2)/(4.0*bcVals1R2*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2); 
  phiC[1] = (((24.0*bcVals1R2*rdx2SqVol1R2+(24.0*rdx2SqVol[0]*bcVals1R2-24.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[11]+((-36.0*rdx2SqVol[0]*bcVals[1]*rdx2SqVol[1])-18.0*rdx2SqVol0R2*bcVals[1]+12.0*bcVals[0]*rdx2SqVol0R2)*bcVals[2]+(24.0*bcVals1R2*rdx2SqVol[1]+24.0*rdx2SqVol[0]*bcVals1R2-24.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rhoC[1]+((-2.0*bcVals1R2*phiUx[1])-4.0*bcVals1R2*phiC[1]+(2.0*phiUxLy[0]+4.0*phiLy[0])*bcVals1R2)*rdx2SqVol1R2+((5.0*rdx2SqVol[0]*bcVals1R2+6.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-20.0*rdx2SqVol[0]*bcVals1R2)*phiC[1]+(12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-4.0*phiLy[0])*rdx2SqVol[0])*bcVals1R2+(2.0*bcVals[0]*phiUxLy[0]-8.0*bcVals[0]*phiUx[0]+4.0*bcVals[0]*phiLy[0])*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+(7.0*rdx2SqVol0R2*bcVals1R2-6.0*bcVals[0]*rdx2SqVol0R2*bcVals[1])*phiUx[1]+((-7.0*rdx2SqVol0R2*bcVals1R2)+20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]-12.0*bcVals0R2*rdx2SqVol0R2)*phiC[1]+((phiLy[0]-1.0*phiUxLy[0])*rdx2SqVol0R2-6.0*rdx2SqVol[0]*rhoC[0])*bcVals1R2+(12.0*bcVals[0]*rdx2SqVol[0]*rhoC[0]+(2.0*bcVals[0]*phiUxLy[0]+4.0*bcVals[0]*phiUx[0]-2.0*bcVals[0]*phiLy[0])*rdx2SqVol0R2)*bcVals[1])*omega+4.0*bcVals1R2*phiC[1]*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiC[1]*rdx2SqVol[1]+(7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2)*phiC[1])/(4.0*bcVals1R2*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double bcVals0R2 = std::pow(bcVals[0],2);
  const double bcVals1R2 = std::pow(bcVals[1],2);

  phiC[0] = (((12.0*bcVals1R2*rdx2SqVol1R2+(12.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-6.0*rdx2SqVol[0]*bcVals1R2)*rdx2SqVol[1])*bcVals[11]+(((-24.0*rdx2SqVol[0]*bcVals[1]*rdx2SqVol[1])-6.0*rdx2SqVol0R2*bcVals[1]+12.0*bcVals[0]*rdx2SqVol0R2)*bcVals[2]+(12.0*bcVals1R2*rdx2SqVol[1]-6.0*rdx2SqVol[0]*bcVals1R2+12.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rhoC[1]+(2.0*phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]-4.0*phiC[0])*bcVals1R2*rdx2SqVol1R2+((9.0*rdx2SqVol[0]*bcVals1R2-6.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(12.0*rhoC[0]+(4.0*phiUxLy[0]+5.0*phiUx[0]+2.0*phiLy[0]-20.0*phiC[0])*rdx2SqVol[0])*bcVals1R2+((-4.0*bcVals[0]*phiUxLy[0])+10.0*bcVals[0]*phiUx[0]-8.0*bcVals[0]*phiLy[0]+32.0*bcVals[0]*phiC[0])*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+(12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUxLy[0]+7.0*phiUx[0]-2.0*phiLy[0]-7.0*phiC[0])*rdx2SqVol0R2)*bcVals1R2+(((-4.0*bcVals[0]*phiUxLy[0])-14.0*bcVals[0]*phiUx[0]+4.0*bcVals[0]*phiLy[0]+20.0*bcVals[0]*phiC[0])*rdx2SqVol0R2-24.0*bcVals[0]*rdx2SqVol[0]*rhoC[0])*bcVals[1]-12.0*bcVals0R2*phiC[0]*rdx2SqVol0R2)*bcVals[10]+((-24.0*rdx2SqVol[0]*bcVals[1]*rdx2SqVol[1]*bcVals[2])+(4.0*phiUxLy[0]-8.0*phiUx[0]+8.0*phiLy[0]-16.0*phiC[0])*bcVals1R2*rdx2SqVol1R2+((6.0*rdx2SqVol[0]*bcVals1R2-4.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(24.0*rhoC[0]+(4.0*phiUxLy[0]+16.0*phiUx[0]-4.0*phiLy[0]-16.0*phiC[0])*rdx2SqVol[0])*bcVals1R2+16.0*bcVals[0]*phiC[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[9])*omega+(4.0*phiC[0]*bcVals1R2*rdx2SqVol1R2+(20.0*phiC[0]*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*phiC[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*phiC[0]*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*phiC[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*phiC[0]*rdx2SqVol0R2)*bcVals[10]+(16.0*phiC[0]*bcVals1R2*rdx2SqVol1R2+(16.0*phiC[0]*rdx2SqVol[0]*bcVals1R2-16.0*bcVals[0]*phiC[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[9])/((4.0*bcVals1R2*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2)*bcVals[10]+(16.0*bcVals1R2*rdx2SqVol1R2+(16.0*rdx2SqVol[0]*bcVals1R2-16.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[9]); 
  phiC[1] = (((24.0*bcVals1R2*rdx2SqVol1R2+(24.0*rdx2SqVol[0]*bcVals1R2-24.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[11]+(((-36.0*rdx2SqVol[0]*bcVals[1]*rdx2SqVol[1])-18.0*rdx2SqVol0R2*bcVals[1]+12.0*bcVals[0]*rdx2SqVol0R2)*bcVals[2]+(24.0*bcVals1R2*rdx2SqVol[1]+24.0*rdx2SqVol[0]*bcVals1R2-24.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rhoC[1]+((-2.0*bcVals1R2*phiUx[1])-4.0*bcVals1R2*phiC[1]+(2.0*phiUxLy[0]+4.0*phiLy[0])*bcVals1R2)*rdx2SqVol1R2+((5.0*rdx2SqVol[0]*bcVals1R2+6.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiUx[1]+(32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-20.0*rdx2SqVol[0]*bcVals1R2)*phiC[1]+(12.0*rhoC[0]+(phiUxLy[0]+18.0*phiUx[0]-4.0*phiLy[0])*rdx2SqVol[0])*bcVals1R2+(2.0*bcVals[0]*phiUxLy[0]-8.0*bcVals[0]*phiUx[0]+4.0*bcVals[0]*phiLy[0])*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+(7.0*rdx2SqVol0R2*bcVals1R2-6.0*bcVals[0]*rdx2SqVol0R2*bcVals[1])*phiUx[1]+((-7.0*rdx2SqVol0R2*bcVals1R2)+20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]-12.0*bcVals0R2*rdx2SqVol0R2)*phiC[1]+((phiLy[0]-1.0*phiUxLy[0])*rdx2SqVol0R2-6.0*rdx2SqVol[0]*rhoC[0])*bcVals1R2+(12.0*bcVals[0]*rdx2SqVol[0]*rhoC[0]+(2.0*bcVals[0]*phiUxLy[0]+4.0*bcVals[0]*phiUx[0]-2.0*bcVals[0]*phiLy[0])*rdx2SqVol0R2)*bcVals[1])*bcVals[10]+(((-8.0*bcVals1R2*phiUx[1])-16.0*bcVals1R2*phiC[1])*rdx2SqVol1R2+((8.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-8.0*rdx2SqVol[0]*bcVals1R2)*phiUx[1]+(16.0*bcVals[0]*rdx2SqVol[0]*bcVals[1]-16.0*rdx2SqVol[0]*bcVals1R2)*phiC[1])*rdx2SqVol[1])*bcVals[9])*omega+(4.0*bcVals1R2*phiC[1]*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiC[1]*rdx2SqVol[1]+(7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2)*phiC[1])*bcVals[10]+(16.0*bcVals1R2*phiC[1]*rdx2SqVol1R2+(16.0*rdx2SqVol[0]*bcVals1R2-16.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*phiC[1]*rdx2SqVol[1])*bcVals[9])/((4.0*bcVals1R2*rdx2SqVol1R2+(20.0*rdx2SqVol[0]*bcVals1R2-32.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1]+7.0*rdx2SqVol0R2*bcVals1R2-20.0*bcVals[0]*rdx2SqVol0R2*bcVals[1]+12.0*bcVals0R2*rdx2SqVol0R2)*bcVals[10]+(16.0*bcVals1R2*rdx2SqVol1R2+(16.0*rdx2SqVol[0]*bcVals1R2-16.0*bcVals[0]*rdx2SqVol[0]*bcVals[1])*rdx2SqVol[1])*bcVals[9]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]-3.0*rdx2SqVol[0]*bcVals[5]+((-4.0*phiUy[0])-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*omega-4.0*phiC[0]*rdx2SqVol[1]-4.0*phiC[0]*rdx2SqVol[0]))/(4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = -(1.0*((6.0*rdx2SqVol[1]*bcVals[8]+((-3.0*rdx2SqVol[0]*bcVals[5])+((-4.0*phiUy[0])-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*rdx2SqVol[1]*bcVals[5])-4.0*phiC[0]*rdx2SqVol[1])*bcVals[6])*omega+((-4.0*phiC[0]*rdx2SqVol[1])-4.0*phiC[0]*rdx2SqVol[0])*bcVals[7]+4.0*phiC[0]*rdx2SqVol[1]*bcVals[6]))/((4.0*rdx2SqVol[1]+4.0*rdx2SqVol[0])*bcVals[7]-4.0*rdx2SqVol[1]*bcVals[6]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = -(1.0*(((6.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[8]+(6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+(6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])-2.0*phiLxUy[0]+2.0*phiLx[0]+7.0*phiC[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-12.0*rhoC[0]+((-5.0*phiUy[0])-4.0*phiLxUy[0]-2.0*phiLx[0]+20.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]-2.0*phiLxUy[0]-4.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol0R2)*omega-7.0*phiC[0]*rdx2SqVol1R2-20.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]-4.0*phiC[0]*rdx2SqVol0R2))/(7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2); 
  phiC[1] = -(1.0*(((18.0*rdx2SqVol1R2+36.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[8]+((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-24.0*rdx2SqVol0R2)*bcVals[5]+((-24.0*rdx2SqVol[1])-24.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[1])+7.0*phiC[1]+phiLxUy[0]-1.0*phiLx[0])*rdx2SqVol1R2+((-5.0*rdx2SqVol[0]*phiUy[1])+20.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-18.0*phiUy[0])-1.0*phiLxUy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiUy[1]+4.0*rdx2SqVol0R2*phiC[1]-12.0*rdx2SqVol[0]*rhoC[0]+((-2.0*phiLxUy[0])-4.0*phiLx[0])*rdx2SqVol0R2)*omega-7.0*phiC[1]*rdx2SqVol1R2-20.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]-4.0*rdx2SqVol0R2*phiC[1]))/(7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*((((6.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[7]-12.0*rdx2SqVol1R2*bcVals[6])*bcVals[8]+((6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+(6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])-2.0*phiLxUy[0]+2.0*phiLx[0]+7.0*phiC[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-12.0*rhoC[0]+((-5.0*phiUy[0])-4.0*phiLxUy[0]-2.0*phiLx[0]+20.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]-2.0*phiLxUy[0]-4.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol0R2)*bcVals7R2+((-12.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[5])-12.0*rdx2SqVol[1]*rhoC[1]+(14.0*phiUy[0]+4.0*phiLxUy[0]-4.0*phiLx[0]-20.0*phiC[0])*rdx2SqVol1R2+(6.0*rdx2SqVol[0]*phiUy[1]+24.0*rhoC[0]+((-10.0*phiUy[0])+4.0*phiLxUy[0]+8.0*phiLx[0]-32.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[6]*bcVals[7]+12.0*phiC[0]*rdx2SqVol1R2*bcVals6R2)*omega+((-7.0*phiC[0]*rdx2SqVol1R2)-20.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]-4.0*phiC[0]*rdx2SqVol0R2)*bcVals7R2+(20.0*phiC[0]*rdx2SqVol1R2+32.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[6]*bcVals[7]-12.0*phiC[0]*rdx2SqVol1R2*bcVals6R2))/((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals7R2+((-20.0*rdx2SqVol1R2)-32.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[6]*bcVals[7]+12.0*rdx2SqVol1R2*bcVals6R2); 
  phiC[1] = -(1.0*((((18.0*rdx2SqVol1R2+36.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[7]-12.0*rdx2SqVol1R2*bcVals[6])*bcVals[8]+(((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-24.0*rdx2SqVol0R2)*bcVals[5]+((-24.0*rdx2SqVol[1])-24.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[1])+7.0*phiC[1]+phiLxUy[0]-1.0*phiLx[0])*rdx2SqVol1R2+((-5.0*rdx2SqVol[0]*phiUy[1])+20.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-18.0*phiUy[0])-1.0*phiLxUy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiUy[1]+4.0*rdx2SqVol0R2*phiC[1]-12.0*rdx2SqVol[0]*rhoC[0]+((-2.0*phiLxUy[0])-4.0*phiLx[0])*rdx2SqVol0R2)*bcVals7R2+(24.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[5]+24.0*rdx2SqVol[1]*rhoC[1]+(6.0*phiUy[1]-20.0*phiC[1]-4.0*phiUy[0]-2.0*phiLxUy[0]+2.0*phiLx[0])*rdx2SqVol1R2+((-6.0*rdx2SqVol[0]*phiUy[1])-32.0*rdx2SqVol[0]*phiC[1]-12.0*rhoC[0]+(8.0*phiUy[0]-2.0*phiLxUy[0]-4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[6]*bcVals[7]+12.0*phiC[1]*rdx2SqVol1R2*bcVals6R2)*omega+((-7.0*phiC[1]*rdx2SqVol1R2)-20.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]-4.0*rdx2SqVol0R2*phiC[1])*bcVals7R2+(20.0*phiC[1]*rdx2SqVol1R2+32.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1])*bcVals[6]*bcVals[7]-12.0*phiC[1]*rdx2SqVol1R2*bcVals6R2))/((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals7R2+((-20.0*rdx2SqVol1R2)-32.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[6]*bcVals[7]+12.0*rdx2SqVol1R2*bcVals6R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = -(1.0*((((6.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]+24.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[3])*bcVals[8]+(6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+((6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])-2.0*phiLxUy[0]+2.0*phiLx[0]+7.0*phiC[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-12.0*rhoC[0]+((-5.0*phiUy[0])-4.0*phiLxUy[0]-2.0*phiLx[0]+20.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]-2.0*phiLxUy[0]-4.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol0R2)*bcVals[4]+((((-16.0*phiUy[0])-4.0*phiLxUy[0]+4.0*phiLx[0]+16.0*phiC[0])*rdx2SqVol[0]-6.0*rdx2SqVol[0]*phiUy[1])*rdx2SqVol[1]-24.0*rdx2SqVol[0]*rhoC[0]+(8.0*phiUy[0]-4.0*phiLxUy[0]-8.0*phiLx[0]+16.0*phiC[0])*rdx2SqVol0R2)*bcVals[3])*omega+((-7.0*phiC[0]*rdx2SqVol1R2)-20.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]-4.0*phiC[0]*rdx2SqVol0R2)*bcVals[4]+((-16.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1])-16.0*phiC[0]*rdx2SqVol0R2)*bcVals[3]))/((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 
  phiC[1] = -(1.0*(((18.0*rdx2SqVol1R2+36.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]*bcVals[8]+((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-24.0*rdx2SqVol0R2)*bcVals[5]+(((-24.0*rdx2SqVol[1])-24.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[1])+7.0*phiC[1]+phiLxUy[0]-1.0*phiLx[0])*rdx2SqVol1R2+((-5.0*rdx2SqVol[0]*phiUy[1])+20.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-18.0*phiUy[0])-1.0*phiLxUy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiUy[1]+4.0*rdx2SqVol0R2*phiC[1]-12.0*rdx2SqVol[0]*rhoC[0]+((-2.0*phiLxUy[0])-4.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+((8.0*rdx2SqVol[0]*phiUy[1]+16.0*rdx2SqVol[0]*phiC[1])*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiUy[1]+16.0*rdx2SqVol0R2*phiC[1])*bcVals[3])*omega+((-7.0*phiC[1]*rdx2SqVol1R2)-20.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]-4.0*rdx2SqVol0R2*phiC[1])*bcVals[4]+((-16.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1])-16.0*rdx2SqVol0R2*phiC[1])*bcVals[3]))/((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double bcVals6R2 = std::pow(bcVals[6],2);
  const double bcVals7R2 = std::pow(bcVals[7],2);

  phiC[0] = -(1.0*(((((6.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]+24.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[3])*bcVals[7]-12.0*rdx2SqVol1R2*bcVals[4]*bcVals[6])*bcVals[8]+((6.0*rdx2SqVol[0]*rdx2SqVol[1]-12.0*rdx2SqVol0R2)*bcVals[5]+((6.0*rdx2SqVol[1]-12.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[0])-2.0*phiLxUy[0]+2.0*phiLx[0]+7.0*phiC[0])*rdx2SqVol1R2+((-9.0*rdx2SqVol[0]*phiUy[1])-12.0*rhoC[0]+((-5.0*phiUy[0])-4.0*phiLxUy[0]-2.0*phiLx[0]+20.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]-12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiUy[0]-2.0*phiLxUy[0]-4.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol0R2)*bcVals[4]+((((-16.0*phiUy[0])-4.0*phiLxUy[0]+4.0*phiLx[0]+16.0*phiC[0])*rdx2SqVol[0]-6.0*rdx2SqVol[0]*phiUy[1])*rdx2SqVol[1]-24.0*rdx2SqVol[0]*rhoC[0]+(8.0*phiUy[0]-4.0*phiLxUy[0]-8.0*phiLx[0]+16.0*phiC[0])*rdx2SqVol0R2)*bcVals[3])*bcVals7R2+((-12.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[5])+((-12.0*rdx2SqVol[1]*rhoC[1])+(14.0*phiUy[0]+4.0*phiLxUy[0]-4.0*phiLx[0]-20.0*phiC[0])*rdx2SqVol1R2+(6.0*rdx2SqVol[0]*phiUy[1]+24.0*rhoC[0]+((-10.0*phiUy[0])+4.0*phiLxUy[0]+8.0*phiLx[0]-32.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[4]+(4.0*rdx2SqVol[0]*phiUy[1]-16.0*phiC[0]*rdx2SqVol[0])*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]+12.0*phiC[0]*rdx2SqVol1R2*bcVals[4]*bcVals6R2)*omega+(((-7.0*phiC[0]*rdx2SqVol1R2)-20.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]-4.0*phiC[0]*rdx2SqVol0R2)*bcVals[4]+((-16.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1])-16.0*phiC[0]*rdx2SqVol0R2)*bcVals[3])*bcVals7R2+((20.0*phiC[0]*rdx2SqVol1R2+32.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]+16.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]-12.0*phiC[0]*rdx2SqVol1R2*bcVals[4]*bcVals6R2))/(((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3])*bcVals7R2+(((-20.0*rdx2SqVol1R2)-32.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]-16.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]+12.0*rdx2SqVol1R2*bcVals[4]*bcVals6R2); 
  phiC[1] = -(1.0*((((18.0*rdx2SqVol1R2+36.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]*bcVals[7]-12.0*rdx2SqVol1R2*bcVals[4]*bcVals[6])*bcVals[8]+(((-24.0*rdx2SqVol[0]*rdx2SqVol[1])-24.0*rdx2SqVol0R2)*bcVals[5]+(((-24.0*rdx2SqVol[1])-24.0*rdx2SqVol[0])*rhoC[1]+((-7.0*phiUy[1])+7.0*phiC[1]+phiLxUy[0]-1.0*phiLx[0])*rdx2SqVol1R2+((-5.0*rdx2SqVol[0]*phiUy[1])+20.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-18.0*phiUy[0])-1.0*phiLxUy[0]+4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiUy[1]+4.0*rdx2SqVol0R2*phiC[1]-12.0*rdx2SqVol[0]*rhoC[0]+((-2.0*phiLxUy[0])-4.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+((8.0*rdx2SqVol[0]*phiUy[1]+16.0*rdx2SqVol[0]*phiC[1])*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiUy[1]+16.0*rdx2SqVol0R2*phiC[1])*bcVals[3])*bcVals7R2+(24.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[5]+(24.0*rdx2SqVol[1]*rhoC[1]+(6.0*phiUy[1]-20.0*phiC[1]-4.0*phiUy[0]-2.0*phiLxUy[0]+2.0*phiLx[0])*rdx2SqVol1R2+((-6.0*rdx2SqVol[0]*phiUy[1])-32.0*rdx2SqVol[0]*phiC[1]-12.0*rhoC[0]+(8.0*phiUy[0]-2.0*phiLxUy[0]-4.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[4]+((-8.0*rdx2SqVol[0]*phiUy[1])-16.0*rdx2SqVol[0]*phiC[1])*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]+12.0*phiC[1]*rdx2SqVol1R2*bcVals[4]*bcVals6R2)*omega+(((-7.0*phiC[1]*rdx2SqVol1R2)-20.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]-4.0*rdx2SqVol0R2*phiC[1])*bcVals[4]+((-16.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1])-16.0*rdx2SqVol0R2*phiC[1])*bcVals[3])*bcVals7R2+((20.0*phiC[1]*rdx2SqVol1R2+32.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1])*bcVals[4]+16.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]-12.0*phiC[1]*rdx2SqVol1R2*bcVals[4]*bcVals6R2))/(((7.0*rdx2SqVol1R2+20.0*rdx2SqVol[0]*rdx2SqVol[1]+4.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3])*bcVals7R2+(((-20.0*rdx2SqVol1R2)-32.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[4]-16.0*rdx2SqVol[0]*rdx2SqVol[1]*bcVals[3])*bcVals[6]*bcVals[7]+12.0*rdx2SqVol1R2*bcVals[4]*bcVals6R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  phiC[0] = ((((5.0*dxC[1]+4.0*dxC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*dxC[1]-2.0*dxC[0]*rdx2SqVol[0])*bcVals[11]+(((-2.0*dxC[1])-1.0*dxC[0])*rdx2SqVol[1]+4.0*rdx2SqVol[0]*dxC[1]+5.0*dxC[0]*rdx2SqVol[0])*bcVals[5]+((dxC[1]+dxC[0])*phiLy[1]+(dxC[1]+dxC[0])*phiLx[1]+(4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*dxC[1]+4.0*dxC[0]*phiLy[0]+dxC[0]*phiLxLy[0]-2.0*dxC[0]*phiLx[0]-8.0*dxC[0]*phiC[0])*rdx2SqVol[1]+(rdx2SqVol[0]*dxC[1]+dxC[0]*rdx2SqVol[0])*phiLy[1]+(rdx2SqVol[0]*dxC[1]+dxC[0]*rdx2SqVol[0])*phiLx[1]+(6.0*rhoC[0]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0])*dxC[1]+6.0*dxC[0]*rhoC[0]+((-2.0*dxC[0]*phiLy[0])+dxC[0]*phiLxLy[0]+4.0*dxC[0]*phiLx[0]-8.0*dxC[0]*phiC[0])*rdx2SqVol[0])*omega+(8.0*phiC[0]*dxC[1]+8.0*dxC[0]*phiC[0])*rdx2SqVol[1]+8.0*phiC[0]*rdx2SqVol[0]*dxC[1]+8.0*dxC[0]*phiC[0]*rdx2SqVol[0])/((8.0*dxC[1]+8.0*dxC[0])*rdx2SqVol[1]+8.0*rdx2SqVol[0]*dxC[1]+8.0*dxC[0]*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((-2.0*rdx2SqVol1R2)+14.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[5]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[2]+(2.0*phiLy[1]+8.0*phiLy[0]+2.0*phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol1R2+(4.0*rdx2SqVol[0]*phiLy[1]+9.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(4.0*phiLy[0]+4.0*phiLxLy[0]+5.0*phiLx[0]-40.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiLy[1]+12.0*rdx2SqVol[0]*rhoC[0]+((-4.0*phiLy[0])+2.0*phiLxLy[0]+7.0*phiLx[0]-14.0*phiC[0])*rdx2SqVol0R2)*omega+8.0*phiC[0]*rdx2SqVol1R2+40.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+14.0*phiC[0]*rdx2SqVol0R2)/(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2); 
  phiC[1] = (bcVals[5]-1.0*phiC[1])*omega+phiC[1]; 
  phiC[2] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+((-2.0*rdx2SqVol1R2)+23.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[5]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[2]+((-8.0*rdx2SqVol1R2)-40.0*rdx2SqVol[0]*rdx2SqVol[1]-14.0*rdx2SqVol0R2)*phiC[2]+(2.0*phiLy[1]-2.0*phiLx[1]+8.0*phiLy[0]+2.0*phiLxLy[0])*rdx2SqVol1R2+(rdx2SqVol[0]*phiLy[1]+5.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+((-8.0*phiLy[0])+phiLxLy[0]+18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-1.0*rdx2SqVol0R2*phiLy[1]+7.0*rdx2SqVol0R2*phiLx[1]-6.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLy[0]-1.0*phiLxLy[0])*rdx2SqVol0R2)*omega+(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*phiC[2])/(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((12.0*rdx2SqVol1R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(((-2.0*rdx2SqVol1R2)+14.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[5]+(12.0*rdx2SqVol[1]-6.0*rdx2SqVol[0])*rhoC[2]+(2.0*phiLy[1]+8.0*phiLy[0]+2.0*phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol1R2+(4.0*rdx2SqVol[0]*phiLy[1]+9.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(4.0*phiLy[0]+4.0*phiLxLy[0]+5.0*phiLx[0]-40.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiLy[1]+12.0*rdx2SqVol[0]*rhoC[0]+((-4.0*phiLy[0])+2.0*phiLxLy[0]+7.0*phiLx[0]-14.0*phiC[0])*rdx2SqVol0R2)*bcVals[10]+((12.0*rdx2SqVol[0]*rdx2SqVol[1]-6.0*rdx2SqVol1R2)*bcVals[5]+(2.0*phiLy[1]+2.0*phiLx[1]+8.0*phiLy[0]+2.0*phiLxLy[0]-4.0*phiLx[0]-16.0*phiC[0])*rdx2SqVol1R2+(2.0*rdx2SqVol[0]*phiLy[1]+2.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+((-4.0*phiLy[0])+2.0*phiLxLy[0]+8.0*phiLx[0]-16.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1])*bcVals[9])*omega+(8.0*phiC[0]*rdx2SqVol1R2+40.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+14.0*phiC[0]*rdx2SqVol0R2)*bcVals[10]+(16.0*phiC[0]*rdx2SqVol1R2+16.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9])/((8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 
  phiC[1] = (bcVals[5]-1.0*phiC[1])*omega+phiC[1]; 
  phiC[2] = (((24.0*rdx2SqVol1R2+24.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[11]+(((-2.0*rdx2SqVol1R2)+23.0*rdx2SqVol[0]*rdx2SqVol[1]+7.0*rdx2SqVol0R2)*bcVals[5]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[2]+((-8.0*rdx2SqVol1R2)-40.0*rdx2SqVol[0]*rdx2SqVol[1]-14.0*rdx2SqVol0R2)*phiC[2]+(2.0*phiLy[1]-2.0*phiLx[1]+8.0*phiLy[0]+2.0*phiLxLy[0])*rdx2SqVol1R2+(rdx2SqVol[0]*phiLy[1]+5.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+((-8.0*phiLy[0])+phiLxLy[0]+18.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-1.0*rdx2SqVol0R2*phiLy[1]+7.0*rdx2SqVol0R2*phiLx[1]-6.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLy[0]-1.0*phiLxLy[0])*rdx2SqVol0R2)*bcVals[10]+(((-8.0*rdx2SqVol1R2)-8.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[5]+((-16.0*rdx2SqVol1R2)-16.0*rdx2SqVol[0]*rdx2SqVol[1])*phiC[2])*bcVals[9])*omega+(8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*phiC[2]*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*phiC[2]*bcVals[9])/((8.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+14.0*rdx2SqVol0R2)*bcVals[10]+(16.0*rdx2SqVol1R2+16.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[9]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = (((7.0*rdx2SqVol1R2+14.0*rdx2SqVol[0]*rdx2SqVol[1]-2.0*rdx2SqVol0R2)*bcVals[11]+(12.0*rdx2SqVol0R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[5]+(12.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*rhoC[1]+(2.0*phiLx[1]+7.0*phiLy[0]+2.0*phiLxLy[0]-4.0*phiLx[0]-14.0*phiC[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiLy[1]+4.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(5.0*phiLy[0]+4.0*phiLxLy[0]+4.0*phiLx[0]-40.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiLx[1]+12.0*rdx2SqVol[0]*rhoC[0]+((-2.0*phiLy[0])+2.0*phiLxLy[0]+8.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol0R2)*omega+14.0*phiC[0]*rdx2SqVol1R2+40.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+8.0*phiC[0]*rdx2SqVol0R2)/(14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2); 
  phiC[1] = (((7.0*rdx2SqVol1R2+23.0*rdx2SqVol[0]*rdx2SqVol[1]-2.0*rdx2SqVol0R2)*bcVals[11]+(24.0*rdx2SqVol[0]*rdx2SqVol[1]+24.0*rdx2SqVol0R2)*bcVals[5]+(24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+(7.0*phiLy[1]-1.0*phiLx[1]-14.0*phiC[1]-1.0*phiLxLy[0]+2.0*phiLx[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]-40.0*rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+(18.0*phiLy[0]+phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-2.0*rdx2SqVol0R2*phiLy[1]+2.0*rdx2SqVol0R2*phiLx[1]-8.0*rdx2SqVol0R2*phiC[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*omega+14.0*phiC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiC[1])/(14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol0R3 = std::pow(rdx2SqVol[0],3);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double rdx2SqVol1R3 = std::pow(rdx2SqVol[1],3);

  phiC[0] = (((42.0*rdx2SqVol1R3+822.0*rdx2SqVol[0]*rdx2SqVol1R2-84.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[11]+((-84.0*rdx2SqVol[0]*rdx2SqVol1R2)+822.0*rdx2SqVol0R2*rdx2SqVol[1]+42.0*rdx2SqVol0R3)*bcVals[5]+((-42.0*rdx2SqVol1R2)+564.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*rhoC[3]+(84.0*rdx2SqVol1R2+258.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*rhoC[2]+((-42.0*rdx2SqVol1R2)+258.0*rdx2SqVol[0]*rdx2SqVol[1]+84.0*rdx2SqVol0R2)*rhoC[1]+(49.0*phiLy[0]+14.0*phiLxLy[0]-14.0*phiLx[0]-49.0*phiC[0])*rdx2SqVol1R3+(189.0*rdx2SqVol[0]*phiLy[1]+81.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(336.0*phiLy[0]+96.0*phiLxLy[0]-51.0*phiLx[0]-651.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R2+(81.0*rdx2SqVol0R2*phiLy[1]+189.0*rdx2SqVol0R2*phiLx[1]+492.0*rdx2SqVol[0]*rhoC[0]+((-51.0*phiLy[0])+96.0*phiLxLy[0]+336.0*phiLx[0]-651.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol[1]+84.0*rdx2SqVol0R2*rhoC[0]+((-14.0*phiLy[0])+14.0*phiLxLy[0]+49.0*phiLx[0]-49.0*phiC[0])*rdx2SqVol0R3)*omega+49.0*phiC[0]*rdx2SqVol1R3+651.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R2+651.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol[1]+49.0*phiC[0]*rdx2SqVol0R3)/(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3); 
  phiC[1] = (((126.0*rdx2SqVol1R3+1080.0*rdx2SqVol[0]*rdx2SqVol1R2-126.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[11]+(336.0*rdx2SqVol[0]*rdx2SqVol1R2+1500.0*rdx2SqVol0R2*rdx2SqVol[1]+84.0*rdx2SqVol0R3)*bcVals[5]+(168.0*rdx2SqVol1R2+516.0*rdx2SqVol[0]*rdx2SqVol[1]-84.0*rdx2SqVol0R2)*rhoC[3]+((-42.0*rdx2SqVol1R2)+564.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*rhoC[2]+(168.0*rdx2SqVol1R2+984.0*rdx2SqVol[0]*rdx2SqVol[1]+168.0*rdx2SqVol0R2)*rhoC[1]+(49.0*phiLy[1]-49.0*phiC[1]-7.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol1R3+(336.0*rdx2SqVol[0]*phiLy[1]-72.0*rdx2SqVol[0]*phiLx[1]-651.0*rdx2SqVol[0]*phiC[1]-42.0*rhoC[0]+(378.0*phiLy[0]+36.0*phiLxLy[0]-27.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R2+((-51.0*rdx2SqVol0R2*phiLy[1])+252.0*rdx2SqVol0R2*phiLx[1]-651.0*rdx2SqVol0R2*phiC[1]+258.0*rdx2SqVol[0]*rhoC[0]+(162.0*phiLy[0]+57.0*phiLxLy[0]+231.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol[1]-14.0*rdx2SqVol0R3*phiLy[1]-49.0*rdx2SqVol0R3*phiC[1]+84.0*rdx2SqVol0R2*rhoC[0]+(14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R3)*omega+49.0*phiC[1]*rdx2SqVol1R3+651.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol[1]+49.0*rdx2SqVol0R3*phiC[1])/(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3); 
  phiC[2] = (((84.0*rdx2SqVol1R3+1500.0*rdx2SqVol[0]*rdx2SqVol1R2+336.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[11]+((-126.0*rdx2SqVol[0]*rdx2SqVol1R2)+1080.0*rdx2SqVol0R2*rdx2SqVol[1]+126.0*rdx2SqVol0R3)*bcVals[5]+((-84.0*rdx2SqVol1R2)+516.0*rdx2SqVol[0]*rdx2SqVol[1]+168.0*rdx2SqVol0R2)*rhoC[3]+(168.0*rdx2SqVol1R2+984.0*rdx2SqVol[0]*rdx2SqVol[1]+168.0*rdx2SqVol0R2)*rhoC[2]+((-49.0*rdx2SqVol1R3)-651.0*rdx2SqVol[0]*rdx2SqVol1R2-651.0*rdx2SqVol0R2*rdx2SqVol[1]-49.0*rdx2SqVol0R3)*phiC[2]+((-42.0*rdx2SqVol1R2)+564.0*rdx2SqVol[0]*rdx2SqVol[1]-42.0*rdx2SqVol0R2)*rhoC[1]+((-14.0*phiLx[1])+49.0*phiLy[0]+14.0*phiLxLy[0])*rdx2SqVol1R3+(252.0*rdx2SqVol[0]*phiLy[1]-51.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(231.0*phiLy[0]+57.0*phiLxLy[0]+162.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R2+((-72.0*rdx2SqVol0R2*phiLy[1])+336.0*rdx2SqVol0R2*phiLx[1]+258.0*rdx2SqVol[0]*rhoC[0]+((-27.0*phiLy[0])+36.0*phiLxLy[0]+378.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol[1]+49.0*rdx2SqVol0R3*phiLx[1]-42.0*rdx2SqVol0R2*rhoC[0]+(7.0*phiLy[0]-7.0*phiLxLy[0])*rdx2SqVol0R3)*omega+(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3)*phiC[2])/(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3); 
  phiC[3] = (((252.0*rdx2SqVol1R3+2484.0*rdx2SqVol[0]*rdx2SqVol1R2+504.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[11]+(504.0*rdx2SqVol[0]*rdx2SqVol1R2+2484.0*rdx2SqVol0R2*rdx2SqVol[1]+252.0*rdx2SqVol0R3)*bcVals[5]+(336.0*rdx2SqVol1R2+1968.0*rdx2SqVol[0]*rdx2SqVol[1]+336.0*rdx2SqVol0R2)*rhoC[3]+((-49.0*rdx2SqVol1R3)-651.0*rdx2SqVol[0]*rdx2SqVol1R2-651.0*rdx2SqVol0R2*rdx2SqVol[1]-49.0*rdx2SqVol0R3)*phiC[3]+((-84.0*rdx2SqVol1R2)+516.0*rdx2SqVol[0]*rdx2SqVol[1]+168.0*rdx2SqVol0R2)*rhoC[2]+(168.0*rdx2SqVol1R2+516.0*rdx2SqVol[0]*rdx2SqVol[1]-84.0*rdx2SqVol0R2)*rhoC[1]+(49.0*phiLy[1]+7.0*phiLx[1]-7.0*phiLxLy[0])*rdx2SqVol1R3+(231.0*rdx2SqVol[0]*phiLy[1]-27.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(504.0*phiLy[0]+87.0*phiLxLy[0]-144.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R2+((-27.0*rdx2SqVol0R2*phiLy[1])+231.0*rdx2SqVol0R2*phiLx[1]+564.0*rdx2SqVol[0]*rhoC[0]+((-144.0*phiLy[0])+87.0*phiLxLy[0]+504.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol[1]+7.0*rdx2SqVol0R3*phiLy[1]+49.0*rdx2SqVol0R3*phiLx[1]-42.0*rdx2SqVol0R2*rhoC[0]-7.0*phiLxLy[0]*rdx2SqVol0R3)*omega+(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3)*phiC[3])/(49.0*rdx2SqVol1R3+651.0*rdx2SqVol[0]*rdx2SqVol1R2+651.0*rdx2SqVol0R2*rdx2SqVol[1]+49.0*rdx2SqVol0R3); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol0R3 = std::pow(rdx2SqVol[0],3);
  const double rdx2SqVol0R4 = std::pow(rdx2SqVol[0],4);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double rdx2SqVol1R3 = std::pow(rdx2SqVol[1],3);
  const double rdx2SqVol1R4 = std::pow(rdx2SqVol[1],4);
  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = ((((42.0*rdx2SqVol1R4+864.0*rdx2SqVol[0]*rdx2SqVol1R3+738.0*rdx2SqVol0R2*rdx2SqVol1R2-84.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[10]+(168.0*rdx2SqVol1R4+336.0*rdx2SqVol[0]*rdx2SqVol1R3-48.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[9])*bcVals[11]+(((-84.0*rdx2SqVol[0]*rdx2SqVol1R3)+738.0*rdx2SqVol0R2*rdx2SqVol1R2+864.0*rdx2SqVol0R3*rdx2SqVol[1]+42.0*rdx2SqVol0R4)*bcVals[5]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[3]+(84.0*rdx2SqVol1R3+342.0*rdx2SqVol[0]*rdx2SqVol1R2+216.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+((-42.0*rdx2SqVol1R3)+216.0*rdx2SqVol[0]*rdx2SqVol1R2+342.0*rdx2SqVol0R2*rdx2SqVol[1]+84.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[0]+14.0*phiLxLy[0]-14.0*phiLx[0]-49.0*phiC[0])*rdx2SqVol1R4+(189.0*rdx2SqVol[0]*phiLy[1]+81.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(385.0*phiLy[0]+110.0*phiLxLy[0]-65.0*phiLx[0]-700.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(270.0*rdx2SqVol0R2*phiLy[1]+270.0*rdx2SqVol0R2*phiLx[1]+576.0*rdx2SqVol[0]*rhoC[0]+(285.0*phiLy[0]+192.0*phiLxLy[0]+285.0*phiLx[0]-1302.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(81.0*rdx2SqVol0R3*phiLy[1]+189.0*rdx2SqVol0R3*phiLx[1]+576.0*rdx2SqVol0R2*rhoC[0]+((-65.0*phiLy[0])+110.0*phiLxLy[0]+385.0*phiLx[0]-700.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1]+84.0*rdx2SqVol0R3*rhoC[0]+((-14.0*phiLy[0])+14.0*phiLxLy[0]+49.0*phiLx[0]-49.0*phiC[0])*rdx2SqVol0R4)*bcVals10R2+(((-372.0*rdx2SqVol[0]*rdx2SqVol1R3)+552.0*rdx2SqVol0R2*rdx2SqVol1R2+708.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[5]+((-168.0*rdx2SqVol1R3)+312.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+(336.0*rdx2SqVol1R3+24.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-204.0*rdx2SqVol1R3)+240.0*rdx2SqVol[0]*rdx2SqVol1R2+660.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(12.0*phiLx[1]+238.0*phiLy[0]+68.0*phiLxLy[0]-80.0*phiLx[0]-280.0*phiC[0])*rdx2SqVol1R4+(402.0*rdx2SqVol[0]*phiLy[1]+396.0*rdx2SqVol[0]*phiLx[1]+408.0*rhoC[0]+(750.0*phiLy[0]+288.0*phiLxLy[0]-108.0*phiLx[0]-2592.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(516.0*rdx2SqVol0R2*phiLy[1]+360.0*rdx2SqVol0R2*phiLx[1]+1320.0*rdx2SqVol[0]*rhoC[0]+(174.0*phiLy[0]+336.0*phiLxLy[0]+636.0*phiLx[0]-2760.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(6.0*rdx2SqVol0R3*phiLy[1]+84.0*rdx2SqVol0R3*phiLx[1]+696.0*rdx2SqVol0R2*rhoC[0]+((-122.0*phiLy[0])+116.0*phiLxLy[0]+448.0*phiLx[0]-448.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[9]*bcVals[10]+((288.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol[0]*rdx2SqVol1R3)*bcVals[5]+(288.0*rdx2SqVol[0]*rdx2SqVol1R2-144.0*rdx2SqVol1R3)*rhoC[1]+(48.0*phiLx[1]+168.0*phiLy[0]+48.0*phiLxLy[0]-96.0*phiLx[0]-336.0*phiC[0])*rdx2SqVol1R4+(216.0*rdx2SqVol[0]*phiLy[1]+96.0*rdx2SqVol[0]*phiLx[1]+288.0*rhoC[0]+(120.0*phiLy[0]+96.0*phiLxLy[0]+96.0*phiLx[0]-960.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(48.0*rdx2SqVol0R2*phiLx[1]+288.0*rdx2SqVol[0]*rhoC[0]+((-48.0*phiLy[0])+48.0*phiLxLy[0]+192.0*phiLx[0]-192.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2)*bcVals9R2)*omega+(49.0*phiC[0]*rdx2SqVol1R4+700.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+700.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1]+49.0*phiC[0]*rdx2SqVol0R4)*bcVals10R2+(280.0*phiC[0]*rdx2SqVol1R4+2592.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+448.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*phiC[0]*rdx2SqVol1R4+960.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+192.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals9R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals9R2); 
  phiC[1] = ((((126.0*rdx2SqVol1R4+1206.0*rdx2SqVol[0]*rdx2SqVol1R3+954.0*rdx2SqVol0R2*rdx2SqVol1R2-126.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[10]+(168.0*rdx2SqVol1R4+552.0*rdx2SqVol[0]*rdx2SqVol1R3-48.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[9])*bcVals[11]+((336.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+1584.0*rdx2SqVol0R3*rdx2SqVol[1]+84.0*rdx2SqVol0R4)*bcVals[5]+(168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[3]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[1]-49.0*phiC[1]-7.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol1R4+(385.0*rdx2SqVol[0]*phiLy[1]-72.0*rdx2SqVol[0]*phiLx[1]-700.0*rdx2SqVol[0]*phiC[1]-42.0*rhoC[0]+(378.0*phiLy[0]+29.0*phiLxLy[0]-20.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(285.0*rdx2SqVol0R2*phiLy[1]+180.0*rdx2SqVol0R2*phiLx[1]-1302.0*rdx2SqVol0R2*phiC[1]+216.0*rdx2SqVol[0]*rhoC[0]+(540.0*phiLy[0]+93.0*phiLxLy[0]+204.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-65.0*rdx2SqVol0R3*phiLy[1])+252.0*rdx2SqVol0R3*phiLx[1]-700.0*rdx2SqVol0R3*phiC[1]+342.0*rdx2SqVol0R2*rhoC[0]+(162.0*phiLy[0]+71.0*phiLxLy[0]+280.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-14.0*rdx2SqVol0R4*phiLy[1]-49.0*rdx2SqVol0R4*phiC[1]+84.0*rdx2SqVol0R3*rhoC[0]+(14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R4)*bcVals10R2+((984.0*rdx2SqVol[0]*rdx2SqVol1R3+2688.0*rdx2SqVol0R2*rdx2SqVol1R2+1272.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[5]+(336.0*rdx2SqVol1R3-192.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+((-168.0*rdx2SqVol1R3)+744.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(648.0*rdx2SqVol1R3+2880.0*rdx2SqVol[0]*rdx2SqVol1R2+1368.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(182.0*phiLy[1]-6.0*phiLx[1]-280.0*phiC[1]-28.0*phiLy[0]-34.0*phiLxLy[0]+40.0*phiLx[0])*rdx2SqVol1R4+(858.0*rdx2SqVol[0]*phiLy[1]-174.0*rdx2SqVol[0]*phiLx[1]-2592.0*rdx2SqVol[0]*phiC[1]-204.0*rhoC[0]+(816.0*phiLy[0]+6.0*phiLxLy[0]-120.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(126.0*rdx2SqVol0R2*phiLy[1]+390.0*rdx2SqVol0R2*phiLx[1]-2760.0*rdx2SqVol0R2*phiC[1]+240.0*rdx2SqVol[0]*rhoC[0]+(1068.0*phiLy[0]+150.0*phiLxLy[0]+72.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-118.0*rdx2SqVol0R3*phiLy[1])+126.0*rdx2SqVol0R3*phiLx[1]-448.0*rdx2SqVol0R3*phiC[1]+660.0*rdx2SqVol0R2*rhoC[0]+(8.0*phiLy[0]+110.0*phiLxLy[0]+448.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[9]*bcVals[10]+((576.0*rdx2SqVol[0]*rdx2SqVol1R3+576.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[5]+(576.0*rdx2SqVol1R3+576.0*rdx2SqVol[0]*rdx2SqVol1R2)*rhoC[1]+(168.0*phiLy[1]-24.0*phiLx[1]-336.0*phiC[1]-24.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol1R4+(120.0*rdx2SqVol[0]*phiLy[1]+24.0*rdx2SqVol[0]*phiLx[1]-960.0*rdx2SqVol[0]*phiC[1]-144.0*rhoC[0]+(432.0*phiLy[0]+24.0*phiLxLy[0]-192.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-48.0*rdx2SqVol0R2*phiLy[1])+48.0*rdx2SqVol0R2*phiLx[1]-192.0*rdx2SqVol0R2*phiC[1]+288.0*rdx2SqVol[0]*rhoC[0]+(48.0*phiLxLy[0]+192.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2)*bcVals9R2)*omega+(49.0*phiC[1]*rdx2SqVol1R4+700.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+700.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiC[1])*bcVals10R2+(280.0*phiC[1]*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+448.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*phiC[1]*rdx2SqVol1R4+960.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2)*bcVals9R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals9R2); 
  phiC[2] = ((((84.0*rdx2SqVol1R4+1584.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+336.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[9])*bcVals[11]+(((-126.0*rdx2SqVol[0]*rdx2SqVol1R3)+954.0*rdx2SqVol0R2*rdx2SqVol1R2+1206.0*rdx2SqVol0R3*rdx2SqVol[1]+126.0*rdx2SqVol0R4)*bcVals[5]+((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[3]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiC[2]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[1]+((-14.0*phiLx[1])+49.0*phiLy[0]+14.0*phiLxLy[0])*rdx2SqVol1R4+(252.0*rdx2SqVol[0]*phiLy[1]-65.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(280.0*phiLy[0]+71.0*phiLxLy[0]+162.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(180.0*rdx2SqVol0R2*phiLy[1]+285.0*rdx2SqVol0R2*phiLx[1]+342.0*rdx2SqVol[0]*rhoC[0]+(204.0*phiLy[0]+93.0*phiLxLy[0]+540.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-72.0*rdx2SqVol0R3*phiLy[1])+385.0*rdx2SqVol0R3*phiLx[1]+216.0*rdx2SqVol0R2*rhoC[0]+((-20.0*phiLy[0])+29.0*phiLxLy[0]+378.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]+(7.0*phiLy[0]-7.0*phiLxLy[0])*rdx2SqVol0R4)*bcVals10R2+(((-504.0*rdx2SqVol[0]*rdx2SqVol1R3)-216.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[5]+((-336.0*rdx2SqVol1R3)-960.0*rdx2SqVol[0]*rdx2SqVol1R2-192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+(672.0*rdx2SqVol1R3+1920.0*rdx2SqVol[0]*rdx2SqVol1R2+384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-280.0*rdx2SqVol1R4)-2592.0*rdx2SqVol[0]*rdx2SqVol1R3-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[2]+((-168.0*rdx2SqVol1R3)+744.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+((-56.0*phiLx[1])+196.0*phiLy[0]+56.0*phiLxLy[0])*rdx2SqVol1R4+(336.0*rdx2SqVol[0]*phiLy[1]-36.0*rdx2SqVol[0]*phiLx[1]+336.0*rhoC[0]+(60.0*phiLxLy[0]+648.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-120.0*rdx2SqVol0R2*phiLy[1])+564.0*rdx2SqVol0R2*phiLx[1]+24.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLy[0]-12.0*phiLxLy[0]+432.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-24.0*rdx2SqVol0R3*phiLy[1])+112.0*rdx2SqVol0R3*phiLx[1]-96.0*rdx2SqVol0R2*rhoC[0]+(40.0*phiLy[0]-16.0*phiLxLy[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[9]*bcVals[10]+((-336.0*rdx2SqVol1R4)-960.0*rdx2SqVol[0]*rdx2SqVol1R3-192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiC[2]*bcVals9R2)*omega+(49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiC[2]*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[2]*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiC[2]*bcVals9R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals9R2); 
  phiC[3] = ((((252.0*rdx2SqVol1R4+2736.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+504.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[9])*bcVals[11]+((504.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+2736.0*rdx2SqVol0R3*rdx2SqVol[1]+252.0*rdx2SqVol0R4)*bcVals[5]+(336.0*rdx2SqVol1R3+2304.0*rdx2SqVol[0]*rdx2SqVol1R2+2304.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[3]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiC[3]+((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[1]+7.0*phiLx[1]-7.0*phiLxLy[0])*rdx2SqVol1R4+(280.0*rdx2SqVol[0]*phiLy[1]-20.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(504.0*phiLy[0]+80.0*phiLxLy[0]-144.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(204.0*rdx2SqVol0R2*phiLy[1]+204.0*rdx2SqVol0R2*phiLx[1]+522.0*rdx2SqVol[0]*rhoC[0]+(360.0*phiLy[0]+174.0*phiLxLy[0]+360.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-20.0*rdx2SqVol0R3*phiLy[1])+280.0*rdx2SqVol0R3*phiLx[1]+522.0*rdx2SqVol0R2*rhoC[0]+((-144.0*phiLy[0])+80.0*phiLxLy[0]+504.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+7.0*rdx2SqVol0R4*phiLy[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]-7.0*phiLxLy[0]*rdx2SqVol0R4)*bcVals10R2+((1008.0*rdx2SqVol[0]*rdx2SqVol1R3+1728.0*rdx2SqVol0R2*rdx2SqVol1R2+288.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[5]+(672.0*rdx2SqVol1R3+1920.0*rdx2SqVol[0]*rdx2SqVol1R2+384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+((-280.0*rdx2SqVol1R4)-2592.0*rdx2SqVol[0]*rdx2SqVol1R3-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[3]+((-336.0*rdx2SqVol1R3)-960.0*rdx2SqVol[0]*rdx2SqVol1R2-192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(336.0*rdx2SqVol1R3-192.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(84.0*phiLy[1]+28.0*phiLx[1]-56.0*phiLy[0]-28.0*phiLxLy[0])*rdx2SqVol1R4+((-96.0*rdx2SqVol[0]*phiLy[1])+72.0*rdx2SqVol[0]*phiLx[1]-168.0*rhoC[0]+(288.0*phiLy[0]+24.0*phiLxLy[0]-432.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(60.0*rdx2SqVol0R2*phiLy[1]-228.0*rdx2SqVol0R2*phiLx[1]+312.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLxLy[0]-120.0*phiLy[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(24.0*rdx2SqVol0R3*phiLy[1]-56.0*rdx2SqVol0R3*phiLx[1]+48.0*rdx2SqVol0R2*rhoC[0]+(8.0*phiLxLy[0]-32.0*phiLy[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[9]*bcVals[10]+((-336.0*rdx2SqVol1R4)-960.0*rdx2SqVol[0]*rdx2SqVol1R3-192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiC[3]*bcVals9R2)*omega+(49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiC[3]*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[3]*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiC[3]*bcVals9R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals10R2+(280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[9]*bcVals[10]+(336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals9R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);

  phiC[0] = ((((7.0*rdx2SqVol1R2+14.0*rdx2SqVol[0]*rdx2SqVol[1]-2.0*rdx2SqVol0R2)*bcVals[4]+(12.0*rdx2SqVol[0]*rdx2SqVol[1]-6.0*rdx2SqVol0R2)*bcVals[3])*bcVals[11]+(12.0*rdx2SqVol0R2-6.0*rdx2SqVol[0]*rdx2SqVol[1])*bcVals[5]+((12.0*rdx2SqVol[0]-6.0*rdx2SqVol[1])*rhoC[1]+(2.0*phiLx[1]+7.0*phiLy[0]+2.0*phiLxLy[0]-4.0*phiLx[0]-14.0*phiC[0])*rdx2SqVol1R2+(9.0*rdx2SqVol[0]*phiLy[1]+4.0*rdx2SqVol[0]*phiLx[1]+12.0*rhoC[0]+(5.0*phiLy[0]+4.0*phiLxLy[0]+4.0*phiLx[0]-40.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiLx[1]+12.0*rdx2SqVol[0]*rhoC[0]+((-2.0*phiLy[0])+2.0*phiLxLy[0]+8.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol0R2)*bcVals[4]+((2.0*rdx2SqVol[0]*phiLy[1]+2.0*rdx2SqVol[0]*phiLx[1]+(8.0*phiLy[0]+2.0*phiLxLy[0]-4.0*phiLx[0]-16.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol[1]+2.0*rdx2SqVol0R2*phiLy[1]+2.0*rdx2SqVol0R2*phiLx[1]+12.0*rdx2SqVol[0]*rhoC[0]+((-4.0*phiLy[0])+2.0*phiLxLy[0]+8.0*phiLx[0]-16.0*phiC[0])*rdx2SqVol0R2)*bcVals[3])*omega+(14.0*phiC[0]*rdx2SqVol1R2+40.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+8.0*phiC[0]*rdx2SqVol0R2)*bcVals[4]+(16.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol[1]+16.0*phiC[0]*rdx2SqVol0R2)*bcVals[3])/((14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 
  phiC[1] = ((((7.0*rdx2SqVol1R2+23.0*rdx2SqVol[0]*rdx2SqVol[1]-2.0*rdx2SqVol0R2)*bcVals[4]+((-8.0*rdx2SqVol[0]*rdx2SqVol[1])-8.0*rdx2SqVol0R2)*bcVals[3])*bcVals[11]+(24.0*rdx2SqVol[0]*rdx2SqVol[1]+24.0*rdx2SqVol0R2)*bcVals[5]+((24.0*rdx2SqVol[1]+24.0*rdx2SqVol[0])*rhoC[1]+(7.0*phiLy[1]-1.0*phiLx[1]-14.0*phiC[1]-1.0*phiLxLy[0]+2.0*phiLx[0])*rdx2SqVol1R2+(5.0*rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]-40.0*rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+(18.0*phiLy[0]+phiLxLy[0]-8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol[1]-2.0*rdx2SqVol0R2*phiLy[1]+2.0*rdx2SqVol0R2*phiLx[1]-8.0*rdx2SqVol0R2*phiC[1]+12.0*rdx2SqVol[0]*rhoC[0]+(2.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol0R2)*bcVals[4]+((-16.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1])-16.0*rdx2SqVol0R2*phiC[1])*bcVals[3])*omega+(14.0*phiC[1]*rdx2SqVol1R2+40.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+8.0*rdx2SqVol0R2*phiC[1])*bcVals[4]+(16.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol[1]+16.0*rdx2SqVol0R2*phiC[1])*bcVals[3])/((14.0*rdx2SqVol1R2+40.0*rdx2SqVol[0]*rdx2SqVol[1]+8.0*rdx2SqVol0R2)*bcVals[4]+(16.0*rdx2SqVol[0]*rdx2SqVol[1]+16.0*rdx2SqVol0R2)*bcVals[3]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol0R3 = std::pow(rdx2SqVol[0],3);
  const double rdx2SqVol0R4 = std::pow(rdx2SqVol[0],4);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double rdx2SqVol1R3 = std::pow(rdx2SqVol[1],3);
  const double rdx2SqVol1R4 = std::pow(rdx2SqVol[1],4);
  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);

  phiC[0] = ((((42.0*rdx2SqVol1R4+864.0*rdx2SqVol[0]*rdx2SqVol1R3+738.0*rdx2SqVol0R2*rdx2SqVol1R2-84.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(708.0*rdx2SqVol[0]*rdx2SqVol1R3+552.0*rdx2SqVol0R2*rdx2SqVol1R2-372.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(288.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[11]+(((-84.0*rdx2SqVol[0]*rdx2SqVol1R3)+738.0*rdx2SqVol0R2*rdx2SqVol1R2+864.0*rdx2SqVol0R3*rdx2SqVol[1]+42.0*rdx2SqVol0R4)*bcVals[4]+((-48.0*rdx2SqVol0R2*rdx2SqVol1R2)+336.0*rdx2SqVol0R3*rdx2SqVol[1]+168.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+(((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[3]+(84.0*rdx2SqVol1R3+342.0*rdx2SqVol[0]*rdx2SqVol1R2+216.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+((-42.0*rdx2SqVol1R3)+216.0*rdx2SqVol[0]*rdx2SqVol1R2+342.0*rdx2SqVol0R2*rdx2SqVol[1]+84.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[0]+14.0*phiLxLy[0]-14.0*phiLx[0]-49.0*phiC[0])*rdx2SqVol1R4+(189.0*rdx2SqVol[0]*phiLy[1]+81.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(385.0*phiLy[0]+110.0*phiLxLy[0]-65.0*phiLx[0]-700.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(270.0*rdx2SqVol0R2*phiLy[1]+270.0*rdx2SqVol0R2*phiLx[1]+576.0*rdx2SqVol[0]*rhoC[0]+(285.0*phiLy[0]+192.0*phiLxLy[0]+285.0*phiLx[0]-1302.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(81.0*rdx2SqVol0R3*phiLy[1]+189.0*rdx2SqVol0R3*phiLx[1]+576.0*rdx2SqVol0R2*rhoC[0]+((-65.0*phiLy[0])+110.0*phiLxLy[0]+385.0*phiLx[0]-700.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1]+84.0*rdx2SqVol0R3*rhoC[0]+((-14.0*phiLy[0])+14.0*phiLxLy[0]+49.0*phiLx[0]-49.0*phiC[0])*rdx2SqVol0R4)*bcVals4R2+((48.0*rdx2SqVol[0]*rdx2SqVol1R2+312.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((660.0*rdx2SqVol[0]*rdx2SqVol1R2+240.0*rdx2SqVol0R2*rdx2SqVol[1]-204.0*rdx2SqVol0R3)*rhoC[2]+((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)+24.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[1]+(84.0*rdx2SqVol[0]*phiLy[1]+6.0*rdx2SqVol[0]*phiLx[1]+(448.0*phiLy[0]+116.0*phiLxLy[0]-122.0*phiLx[0]-448.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(360.0*rdx2SqVol0R2*phiLy[1]+516.0*rdx2SqVol0R2*phiLx[1]+696.0*rdx2SqVol[0]*rhoC[0]+(636.0*phiLy[0]+336.0*phiLxLy[0]+174.0*phiLx[0]-2760.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(396.0*rdx2SqVol0R3*phiLy[1]+402.0*rdx2SqVol0R3*phiLx[1]+1320.0*rdx2SqVol0R2*rhoC[0]+((-108.0*phiLy[0])+288.0*phiLxLy[0]+750.0*phiLx[0]-2592.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1]+12.0*rdx2SqVol0R4*phiLy[1]+408.0*rdx2SqVol0R3*rhoC[0]+((-80.0*phiLy[0])+68.0*phiLxLy[0]+238.0*phiLx[0]-280.0*phiC[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((288.0*rdx2SqVol0R2*rdx2SqVol[1]-144.0*rdx2SqVol0R3)*rhoC[2]+(48.0*rdx2SqVol0R2*phiLy[1]+(192.0*phiLy[0]+48.0*phiLxLy[0]-48.0*phiLx[0]-192.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(96.0*rdx2SqVol0R3*phiLy[1]+216.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+(96.0*phiLy[0]+96.0*phiLxLy[0]+120.0*phiLx[0]-960.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1]+48.0*rdx2SqVol0R4*phiLy[1]+288.0*rdx2SqVol0R3*rhoC[0]+((-96.0*phiLy[0])+48.0*phiLxLy[0]+168.0*phiLx[0]-336.0*phiC[0])*rdx2SqVol0R4)*bcVals3R2)*omega+(49.0*phiC[0]*rdx2SqVol1R4+700.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+700.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1]+49.0*phiC[0]*rdx2SqVol0R4)*bcVals4R2+(448.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1]+280.0*phiC[0]*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+960.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1]+336.0*phiC[0]*rdx2SqVol0R4)*bcVals3R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2); 
  phiC[1] = ((((126.0*rdx2SqVol1R4+1206.0*rdx2SqVol[0]*rdx2SqVol1R3+954.0*rdx2SqVol0R2*rdx2SqVol1R2-126.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+((-144.0*rdx2SqVol[0]*rdx2SqVol1R3)-216.0*rdx2SqVol0R2*rdx2SqVol1R2-504.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4])*bcVals[11]+((336.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+1584.0*rdx2SqVol0R3*rdx2SqVol[1]+84.0*rdx2SqVol0R4)*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+((168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[3]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[1]-49.0*phiC[1]-7.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol1R4+(385.0*rdx2SqVol[0]*phiLy[1]-72.0*rdx2SqVol[0]*phiLx[1]-700.0*rdx2SqVol[0]*phiC[1]-42.0*rhoC[0]+(378.0*phiLy[0]+29.0*phiLxLy[0]-20.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(285.0*rdx2SqVol0R2*phiLy[1]+180.0*rdx2SqVol0R2*phiLx[1]-1302.0*rdx2SqVol0R2*phiC[1]+216.0*rdx2SqVol[0]*rhoC[0]+(540.0*phiLy[0]+93.0*phiLxLy[0]+204.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-65.0*rdx2SqVol0R3*phiLy[1])+252.0*rdx2SqVol0R3*phiLx[1]-700.0*rdx2SqVol0R3*phiC[1]+342.0*rdx2SqVol0R2*rhoC[0]+(162.0*phiLy[0]+71.0*phiLxLy[0]+280.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-14.0*rdx2SqVol0R4*phiLy[1]-49.0*rdx2SqVol0R4*phiC[1]+84.0*rdx2SqVol0R3*rhoC[0]+(14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R4)*bcVals4R2+(((-192.0*rdx2SqVol[0]*rdx2SqVol1R2)-960.0*rdx2SqVol0R2*rdx2SqVol[1]-336.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((48.0*rdx2SqVol[0]*rdx2SqVol1R2+744.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*rhoC[2]+(384.0*rdx2SqVol[0]*rdx2SqVol1R2+1920.0*rdx2SqVol0R2*rdx2SqVol[1]+672.0*rdx2SqVol0R3)*rhoC[1]+(112.0*rdx2SqVol[0]*phiLy[1]-24.0*rdx2SqVol[0]*phiLx[1]-448.0*rdx2SqVol[0]*phiC[1]+(40.0*phiLx[0]-16.0*phiLxLy[0])*rdx2SqVol[0])*rdx2SqVol1R3+(564.0*rdx2SqVol0R2*phiLy[1]-120.0*rdx2SqVol0R2*phiLx[1]-2760.0*rdx2SqVol0R2*phiC[1]-96.0*rdx2SqVol[0]*rhoC[0]+(432.0*phiLy[0]-12.0*phiLxLy[0]+60.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-36.0*rdx2SqVol0R3*phiLy[1])+336.0*rdx2SqVol0R3*phiLx[1]-2592.0*rdx2SqVol0R3*phiC[1]+24.0*rdx2SqVol0R2*rhoC[0]+(648.0*phiLy[0]+60.0*phiLxLy[0])*rdx2SqVol0R3)*rdx2SqVol[1]-56.0*rdx2SqVol0R4*phiLy[1]-280.0*rdx2SqVol0R4*phiC[1]+336.0*rdx2SqVol0R3*rhoC[0]+(56.0*phiLxLy[0]+196.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((-192.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1]-336.0*rdx2SqVol0R4*phiC[1])*bcVals3R2)*omega+(49.0*phiC[1]*rdx2SqVol1R4+700.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+700.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiC[1])*bcVals4R2+(448.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1]+280.0*rdx2SqVol0R4*phiC[1])*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+960.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1]+336.0*rdx2SqVol0R4*phiC[1])*bcVals3R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2); 
  phiC[2] = ((((84.0*rdx2SqVol1R4+1584.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+336.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1272.0*rdx2SqVol[0]*rdx2SqVol1R3+2688.0*rdx2SqVol0R2*rdx2SqVol1R2+984.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(576.0*rdx2SqVol0R2*rdx2SqVol1R2+576.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[11]+(((-126.0*rdx2SqVol[0]*rdx2SqVol1R3)+954.0*rdx2SqVol0R2*rdx2SqVol1R2+1206.0*rdx2SqVol0R3*rdx2SqVol[1]+126.0*rdx2SqVol0R4)*bcVals[4]+((-48.0*rdx2SqVol0R2*rdx2SqVol1R2)+552.0*rdx2SqVol0R3*rdx2SqVol[1]+168.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+(((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[3]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiC[2]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[1]+((-14.0*phiLx[1])+49.0*phiLy[0]+14.0*phiLxLy[0])*rdx2SqVol1R4+(252.0*rdx2SqVol[0]*phiLy[1]-65.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(280.0*phiLy[0]+71.0*phiLxLy[0]+162.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(180.0*rdx2SqVol0R2*phiLy[1]+285.0*rdx2SqVol0R2*phiLx[1]+342.0*rdx2SqVol[0]*rhoC[0]+(204.0*phiLy[0]+93.0*phiLxLy[0]+540.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-72.0*rdx2SqVol0R3*phiLy[1])+385.0*rdx2SqVol0R3*phiLx[1]+216.0*rdx2SqVol0R2*rhoC[0]+((-20.0*phiLy[0])+29.0*phiLxLy[0]+378.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]+(7.0*phiLy[0]-7.0*phiLxLy[0])*rdx2SqVol0R4)*bcVals4R2+(((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)-192.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((1368.0*rdx2SqVol[0]*rdx2SqVol1R2+2880.0*rdx2SqVol0R2*rdx2SqVol[1]+648.0*rdx2SqVol0R3)*rhoC[2]+((-448.0*rdx2SqVol[0]*rdx2SqVol1R3)-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-2592.0*rdx2SqVol0R3*rdx2SqVol[1]-280.0*rdx2SqVol0R4)*phiC[2]+(48.0*rdx2SqVol[0]*rdx2SqVol1R2+744.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*rhoC[1]+(126.0*rdx2SqVol[0]*phiLy[1]-118.0*rdx2SqVol[0]*phiLx[1]+(448.0*phiLy[0]+110.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(390.0*rdx2SqVol0R2*phiLy[1]+126.0*rdx2SqVol0R2*phiLx[1]+660.0*rdx2SqVol[0]*rhoC[0]+(72.0*phiLy[0]+150.0*phiLxLy[0]+1068.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-174.0*rdx2SqVol0R3*phiLy[1])+858.0*rdx2SqVol0R3*phiLx[1]+240.0*rdx2SqVol0R2*rhoC[0]+((-120.0*phiLy[0])+6.0*phiLxLy[0]+816.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-6.0*rdx2SqVol0R4*phiLy[1]+182.0*rdx2SqVol0R4*phiLx[1]-204.0*rdx2SqVol0R3*rhoC[0]+(40.0*phiLy[0]-34.0*phiLxLy[0]-28.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((576.0*rdx2SqVol0R2*rdx2SqVol[1]+576.0*rdx2SqVol0R3)*rhoC[2]+((-192.0*rdx2SqVol0R2*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*rdx2SqVol[1]-336.0*rdx2SqVol0R4)*phiC[2]+(48.0*rdx2SqVol0R2*phiLy[1]-48.0*rdx2SqVol0R2*phiLx[1]+(192.0*phiLy[0]+48.0*phiLxLy[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(24.0*rdx2SqVol0R3*phiLy[1]+120.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+((-192.0*phiLy[0])+24.0*phiLxLy[0]+432.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-24.0*rdx2SqVol0R4*phiLy[1]+168.0*rdx2SqVol0R4*phiLx[1]-144.0*rdx2SqVol0R3*rhoC[0]+(48.0*phiLy[0]-24.0*phiLxLy[0])*rdx2SqVol0R4)*bcVals3R2)*omega+(49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiC[2]*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*phiC[2]*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*phiC[2]*bcVals3R2)/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2); 
  phiC[3] = ((((252.0*rdx2SqVol1R4+2736.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+504.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(288.0*rdx2SqVol[0]*rdx2SqVol1R3+1728.0*rdx2SqVol0R2*rdx2SqVol1R2+1008.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4])*bcVals[11]+((504.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+2736.0*rdx2SqVol0R3*rdx2SqVol[1]+252.0*rdx2SqVol0R4)*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+((336.0*rdx2SqVol1R3+2304.0*rdx2SqVol[0]*rdx2SqVol1R2+2304.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[3]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiC[3]+((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[1]+7.0*phiLx[1]-7.0*phiLxLy[0])*rdx2SqVol1R4+(280.0*rdx2SqVol[0]*phiLy[1]-20.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(504.0*phiLy[0]+80.0*phiLxLy[0]-144.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(204.0*rdx2SqVol0R2*phiLy[1]+204.0*rdx2SqVol0R2*phiLx[1]+522.0*rdx2SqVol[0]*rhoC[0]+(360.0*phiLy[0]+174.0*phiLxLy[0]+360.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-20.0*rdx2SqVol0R3*phiLy[1])+280.0*rdx2SqVol0R3*phiLx[1]+522.0*rdx2SqVol0R2*rhoC[0]+((-144.0*phiLy[0])+80.0*phiLxLy[0]+504.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+7.0*rdx2SqVol0R4*phiLy[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]-7.0*phiLxLy[0]*rdx2SqVol0R4)*bcVals4R2+((384.0*rdx2SqVol[0]*rdx2SqVol1R2+1920.0*rdx2SqVol0R2*rdx2SqVol[1]+672.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((-448.0*rdx2SqVol[0]*rdx2SqVol1R3)-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-2592.0*rdx2SqVol0R3*rdx2SqVol[1]-280.0*rdx2SqVol0R4)*bcVals[3]*phiC[3]+(((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)-192.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[2]+((-192.0*rdx2SqVol[0]*rdx2SqVol1R2)-960.0*rdx2SqVol0R2*rdx2SqVol[1]-336.0*rdx2SqVol0R3)*rhoC[1]+((-56.0*rdx2SqVol[0]*phiLy[1])+24.0*rdx2SqVol[0]*phiLx[1]+(8.0*phiLxLy[0]-32.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-228.0*rdx2SqVol0R2*phiLy[1])+60.0*rdx2SqVol0R2*phiLx[1]+48.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLxLy[0]-120.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(72.0*rdx2SqVol0R3*phiLy[1]-96.0*rdx2SqVol0R3*phiLx[1]+312.0*rdx2SqVol0R2*rhoC[0]+((-432.0*phiLy[0])+24.0*phiLxLy[0]+288.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+28.0*rdx2SqVol0R4*phiLy[1]+84.0*rdx2SqVol0R4*phiLx[1]-168.0*rdx2SqVol0R3*rhoC[0]+((-28.0*phiLxLy[0])-56.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((-192.0*rdx2SqVol0R2*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*rdx2SqVol[1]-336.0*rdx2SqVol0R4)*bcVals3R2*phiC[3])*omega+(49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiC[3]*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*phiC[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2*phiC[3])/((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 

  const double rdx2SqVol0R2 = std::pow(rdx2SqVol[0],2);
  const double rdx2SqVol0R3 = std::pow(rdx2SqVol[0],3);
  const double rdx2SqVol0R4 = std::pow(rdx2SqVol[0],4);
  const double rdx2SqVol1R2 = std::pow(rdx2SqVol[1],2);
  const double rdx2SqVol1R3 = std::pow(rdx2SqVol[1],3);
  const double rdx2SqVol1R4 = std::pow(rdx2SqVol[1],4);
  const double bcVals3R2 = std::pow(bcVals[3],2);
  const double bcVals4R2 = std::pow(bcVals[4],2);
  const double bcVals9R2 = std::pow(bcVals[9],2);
  const double bcVals10R2 = std::pow(bcVals[10],2);

  phiC[0] = (((((42.0*rdx2SqVol1R4+864.0*rdx2SqVol[0]*rdx2SqVol1R3+738.0*rdx2SqVol0R2*rdx2SqVol1R2-84.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(708.0*rdx2SqVol[0]*rdx2SqVol1R3+552.0*rdx2SqVol0R2*rdx2SqVol1R2-372.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(288.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[10]+((168.0*rdx2SqVol1R4+336.0*rdx2SqVol[0]*rdx2SqVol1R3-48.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(288.0*rdx2SqVol[0]*rdx2SqVol1R3-144.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals[9])*bcVals[11]+((((-84.0*rdx2SqVol[0]*rdx2SqVol1R3)+738.0*rdx2SqVol0R2*rdx2SqVol1R2+864.0*rdx2SqVol0R3*rdx2SqVol[1]+42.0*rdx2SqVol0R4)*bcVals[4]+((-48.0*rdx2SqVol0R2*rdx2SqVol1R2)+336.0*rdx2SqVol0R3*rdx2SqVol[1]+168.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+(((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[3]+(84.0*rdx2SqVol1R3+342.0*rdx2SqVol[0]*rdx2SqVol1R2+216.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+((-42.0*rdx2SqVol1R3)+216.0*rdx2SqVol[0]*rdx2SqVol1R2+342.0*rdx2SqVol0R2*rdx2SqVol[1]+84.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[0]+14.0*phiLxLy[0]-14.0*phiLx[0]-49.0*phiC[0])*rdx2SqVol1R4+(189.0*rdx2SqVol[0]*phiLy[1]+81.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(385.0*phiLy[0]+110.0*phiLxLy[0]-65.0*phiLx[0]-700.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(270.0*rdx2SqVol0R2*phiLy[1]+270.0*rdx2SqVol0R2*phiLx[1]+576.0*rdx2SqVol[0]*rhoC[0]+(285.0*phiLy[0]+192.0*phiLxLy[0]+285.0*phiLx[0]-1302.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(81.0*rdx2SqVol0R3*phiLy[1]+189.0*rdx2SqVol0R3*phiLx[1]+576.0*rdx2SqVol0R2*rhoC[0]+((-65.0*phiLy[0])+110.0*phiLxLy[0]+385.0*phiLx[0]-700.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1]+84.0*rdx2SqVol0R3*rhoC[0]+((-14.0*phiLy[0])+14.0*phiLxLy[0]+49.0*phiLx[0]-49.0*phiC[0])*rdx2SqVol0R4)*bcVals4R2+((48.0*rdx2SqVol[0]*rdx2SqVol1R2+312.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((660.0*rdx2SqVol[0]*rdx2SqVol1R2+240.0*rdx2SqVol0R2*rdx2SqVol[1]-204.0*rdx2SqVol0R3)*rhoC[2]+((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)+24.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[1]+(84.0*rdx2SqVol[0]*phiLy[1]+6.0*rdx2SqVol[0]*phiLx[1]+(448.0*phiLy[0]+116.0*phiLxLy[0]-122.0*phiLx[0]-448.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(360.0*rdx2SqVol0R2*phiLy[1]+516.0*rdx2SqVol0R2*phiLx[1]+696.0*rdx2SqVol[0]*rhoC[0]+(636.0*phiLy[0]+336.0*phiLxLy[0]+174.0*phiLx[0]-2760.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(396.0*rdx2SqVol0R3*phiLy[1]+402.0*rdx2SqVol0R3*phiLx[1]+1320.0*rdx2SqVol0R2*rhoC[0]+((-108.0*phiLy[0])+288.0*phiLxLy[0]+750.0*phiLx[0]-2592.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1]+12.0*rdx2SqVol0R4*phiLy[1]+408.0*rdx2SqVol0R3*rhoC[0]+((-80.0*phiLy[0])+68.0*phiLxLy[0]+238.0*phiLx[0]-280.0*phiC[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((288.0*rdx2SqVol0R2*rdx2SqVol[1]-144.0*rdx2SqVol0R3)*rhoC[2]+(48.0*rdx2SqVol0R2*phiLy[1]+(192.0*phiLy[0]+48.0*phiLxLy[0]-48.0*phiLx[0]-192.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(96.0*rdx2SqVol0R3*phiLy[1]+216.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+(96.0*phiLy[0]+96.0*phiLxLy[0]+120.0*phiLx[0]-960.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1]+48.0*rdx2SqVol0R4*phiLy[1]+288.0*rdx2SqVol0R3*rhoC[0]+((-96.0*phiLy[0])+48.0*phiLxLy[0]+168.0*phiLx[0]-336.0*phiC[0])*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((((-372.0*rdx2SqVol[0]*rdx2SqVol1R3)+552.0*rdx2SqVol0R2*rdx2SqVol1R2+708.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[4]+(288.0*rdx2SqVol0R3*rdx2SqVol[1]-144.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3])*bcVals[5]+(((-168.0*rdx2SqVol1R3)+312.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+(336.0*rdx2SqVol1R3+24.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-204.0*rdx2SqVol1R3)+240.0*rdx2SqVol[0]*rdx2SqVol1R2+660.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(12.0*phiLx[1]+238.0*phiLy[0]+68.0*phiLxLy[0]-80.0*phiLx[0]-280.0*phiC[0])*rdx2SqVol1R4+(402.0*rdx2SqVol[0]*phiLy[1]+396.0*rdx2SqVol[0]*phiLx[1]+408.0*rhoC[0]+(750.0*phiLy[0]+288.0*phiLxLy[0]-108.0*phiLx[0]-2592.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(516.0*rdx2SqVol0R2*phiLy[1]+360.0*rdx2SqVol0R2*phiLx[1]+1320.0*rdx2SqVol[0]*rhoC[0]+(174.0*phiLy[0]+336.0*phiLxLy[0]+636.0*phiLx[0]-2760.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(6.0*rdx2SqVol0R3*phiLy[1]+84.0*rdx2SqVol0R3*phiLx[1]+696.0*rdx2SqVol0R2*rhoC[0]+((-122.0*phiLy[0])+116.0*phiLxLy[0]+448.0*phiLx[0]-448.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals4R2+((288.0*rdx2SqVol[0]*rdx2SqVol1R2-144.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(288.0*rdx2SqVol0R2*rdx2SqVol[1]-144.0*rdx2SqVol[0]*rdx2SqVol1R2)*rhoC[1]+(120.0*rdx2SqVol[0]*phiLy[1]+120.0*rdx2SqVol[0]*phiLx[1]+(648.0*phiLy[0]+168.0*phiLxLy[0]-288.0*phiLx[0]-1104.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(456.0*rdx2SqVol0R2*phiLy[1]+456.0*rdx2SqVol0R2*phiLx[1]+1008.0*rdx2SqVol[0]*rhoC[0]+(360.0*phiLy[0]+336.0*phiLxLy[0]+360.0*phiLx[0]-3072.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(120.0*rdx2SqVol0R3*phiLy[1]+120.0*rdx2SqVol0R3*phiLx[1]+1008.0*rdx2SqVol0R2*rhoC[0]+((-288.0*phiLy[0])+168.0*phiLxLy[0]+648.0*phiLx[0]-1104.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[3]*bcVals[4]+((48.0*rdx2SqVol0R2*phiLy[1]+48.0*rdx2SqVol0R2*phiLx[1]+(192.0*phiLy[0]+48.0*phiLxLy[0]-96.0*phiLx[0]-384.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(48.0*rdx2SqVol0R3*phiLy[1]+48.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+((-96.0*phiLy[0])+48.0*phiLxLy[0]+192.0*phiLx[0]-384.0*phiC[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((288.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol[0]*rdx2SqVol1R3)*bcVals[4]*bcVals[5]+((288.0*rdx2SqVol[0]*rdx2SqVol1R2-144.0*rdx2SqVol1R3)*rhoC[1]+(48.0*phiLx[1]+168.0*phiLy[0]+48.0*phiLxLy[0]-96.0*phiLx[0]-336.0*phiC[0])*rdx2SqVol1R4+(216.0*rdx2SqVol[0]*phiLy[1]+96.0*rdx2SqVol[0]*phiLx[1]+288.0*rhoC[0]+(120.0*phiLy[0]+96.0*phiLxLy[0]+96.0*phiLx[0]-960.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(48.0*rdx2SqVol0R2*phiLx[1]+288.0*rdx2SqVol[0]*rhoC[0]+((-48.0*phiLy[0])+48.0*phiLxLy[0]+192.0*phiLx[0]-192.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2)*bcVals4R2+((48.0*rdx2SqVol[0]*phiLy[1]+48.0*rdx2SqVol[0]*phiLx[1]+(192.0*phiLy[0]+48.0*phiLxLy[0]-96.0*phiLx[0]-384.0*phiC[0])*rdx2SqVol[0])*rdx2SqVol1R3+(48.0*rdx2SqVol0R2*phiLy[1]+48.0*rdx2SqVol0R2*phiLx[1]+288.0*rdx2SqVol[0]*rhoC[0]+((-96.0*phiLy[0])+48.0*phiLxLy[0]+192.0*phiLx[0]-384.0*phiC[0])*rdx2SqVol0R2)*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2)*omega+((49.0*phiC[0]*rdx2SqVol1R4+700.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+700.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1]+49.0*phiC[0]*rdx2SqVol0R4)*bcVals4R2+(448.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1]+280.0*phiC[0]*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+960.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1]+336.0*phiC[0]*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((280.0*phiC[0]*rdx2SqVol1R4+2592.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+448.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1104.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2+384.0*phiC[0]*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*phiC[0]*rdx2SqVol1R4+960.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+192.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*phiC[0]*rdx2SqVol[0]*rdx2SqVol1R3+384.0*phiC[0]*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2)/(((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2); 
  phiC[1] = (((((126.0*rdx2SqVol1R4+1206.0*rdx2SqVol[0]*rdx2SqVol1R3+954.0*rdx2SqVol0R2*rdx2SqVol1R2-126.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+((-144.0*rdx2SqVol[0]*rdx2SqVol1R3)-216.0*rdx2SqVol0R2*rdx2SqVol1R2-504.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4])*bcVals[10]+((168.0*rdx2SqVol1R4+552.0*rdx2SqVol[0]*rdx2SqVol1R3-48.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+((-192.0*rdx2SqVol[0]*rdx2SqVol1R3)-192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals[9])*bcVals[11]+(((336.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+1584.0*rdx2SqVol0R3*rdx2SqVol[1]+84.0*rdx2SqVol0R4)*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+((168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[3]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[1]-49.0*phiC[1]-7.0*phiLxLy[0]+7.0*phiLx[0])*rdx2SqVol1R4+(385.0*rdx2SqVol[0]*phiLy[1]-72.0*rdx2SqVol[0]*phiLx[1]-700.0*rdx2SqVol[0]*phiC[1]-42.0*rhoC[0]+(378.0*phiLy[0]+29.0*phiLxLy[0]-20.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(285.0*rdx2SqVol0R2*phiLy[1]+180.0*rdx2SqVol0R2*phiLx[1]-1302.0*rdx2SqVol0R2*phiC[1]+216.0*rdx2SqVol[0]*rhoC[0]+(540.0*phiLy[0]+93.0*phiLxLy[0]+204.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-65.0*rdx2SqVol0R3*phiLy[1])+252.0*rdx2SqVol0R3*phiLx[1]-700.0*rdx2SqVol0R3*phiC[1]+342.0*rdx2SqVol0R2*rhoC[0]+(162.0*phiLy[0]+71.0*phiLxLy[0]+280.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-14.0*rdx2SqVol0R4*phiLy[1]-49.0*rdx2SqVol0R4*phiC[1]+84.0*rdx2SqVol0R3*rhoC[0]+(14.0*phiLxLy[0]+49.0*phiLx[0])*rdx2SqVol0R4)*bcVals4R2+(((-192.0*rdx2SqVol[0]*rdx2SqVol1R2)-960.0*rdx2SqVol0R2*rdx2SqVol[1]-336.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((48.0*rdx2SqVol[0]*rdx2SqVol1R2+744.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*rhoC[2]+(384.0*rdx2SqVol[0]*rdx2SqVol1R2+1920.0*rdx2SqVol0R2*rdx2SqVol[1]+672.0*rdx2SqVol0R3)*rhoC[1]+(112.0*rdx2SqVol[0]*phiLy[1]-24.0*rdx2SqVol[0]*phiLx[1]-448.0*rdx2SqVol[0]*phiC[1]+(40.0*phiLx[0]-16.0*phiLxLy[0])*rdx2SqVol[0])*rdx2SqVol1R3+(564.0*rdx2SqVol0R2*phiLy[1]-120.0*rdx2SqVol0R2*phiLx[1]-2760.0*rdx2SqVol0R2*phiC[1]-96.0*rdx2SqVol[0]*rhoC[0]+(432.0*phiLy[0]-12.0*phiLxLy[0]+60.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-36.0*rdx2SqVol0R3*phiLy[1])+336.0*rdx2SqVol0R3*phiLx[1]-2592.0*rdx2SqVol0R3*phiC[1]+24.0*rdx2SqVol0R2*rhoC[0]+(648.0*phiLy[0]+60.0*phiLxLy[0])*rdx2SqVol0R3)*rdx2SqVol[1]-56.0*rdx2SqVol0R4*phiLy[1]-280.0*rdx2SqVol0R4*phiC[1]+336.0*rdx2SqVol0R3*rhoC[0]+(56.0*phiLxLy[0]+196.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((-192.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1]-336.0*rdx2SqVol0R4*phiC[1])*bcVals3R2)*bcVals10R2+(((984.0*rdx2SqVol[0]*rdx2SqVol1R3+2688.0*rdx2SqVol0R2*rdx2SqVol1R2+1272.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3])*bcVals[5]+((336.0*rdx2SqVol1R3-192.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+((-168.0*rdx2SqVol1R3)+744.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(648.0*rdx2SqVol1R3+2880.0*rdx2SqVol[0]*rdx2SqVol1R2+1368.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(182.0*phiLy[1]-6.0*phiLx[1]-280.0*phiC[1]-28.0*phiLy[0]-34.0*phiLxLy[0]+40.0*phiLx[0])*rdx2SqVol1R4+(858.0*rdx2SqVol[0]*phiLy[1]-174.0*rdx2SqVol[0]*phiLx[1]-2592.0*rdx2SqVol[0]*phiC[1]-204.0*rhoC[0]+(816.0*phiLy[0]+6.0*phiLxLy[0]-120.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(126.0*rdx2SqVol0R2*phiLy[1]+390.0*rdx2SqVol0R2*phiLx[1]-2760.0*rdx2SqVol0R2*phiC[1]+240.0*rdx2SqVol[0]*rhoC[0]+(1068.0*phiLy[0]+150.0*phiLxLy[0]+72.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-118.0*rdx2SqVol0R3*phiLy[1])+126.0*rdx2SqVol0R3*phiLx[1]-448.0*rdx2SqVol0R3*phiC[1]+660.0*rdx2SqVol0R2*rhoC[0]+(8.0*phiLy[0]+110.0*phiLxLy[0]+448.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals4R2+(((-384.0*rdx2SqVol[0]*rdx2SqVol1R2)-384.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[3]*rhoC[3]+((192.0*rdx2SqVol[0]*rdx2SqVol1R2+192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(768.0*rdx2SqVol[0]*rdx2SqVol1R2+768.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(232.0*rdx2SqVol[0]*phiLy[1]-56.0*rdx2SqVol[0]*phiLx[1]-1104.0*rdx2SqVol[0]*phiC[1]+(32.0*phiLy[0]-24.0*phiLxLy[0]+80.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(152.0*rdx2SqVol0R2*phiLy[1]+56.0*rdx2SqVol0R2*phiLx[1]-3072.0*rdx2SqVol0R2*phiC[1]-144.0*rdx2SqVol[0]*rhoC[0]+(496.0*phiLy[0]+24.0*phiLxLy[0]-128.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-80.0*rdx2SqVol0R3*phiLy[1])+112.0*rdx2SqVol0R3*phiLx[1]-1104.0*rdx2SqVol0R3*phiC[1]+288.0*rdx2SqVol0R2*rhoC[0]+(32.0*phiLy[0]+48.0*phiLxLy[0]+224.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[3])*bcVals[4]+((-384.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2)-384.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((576.0*rdx2SqVol[0]*rdx2SqVol1R3+576.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[4]*bcVals[5]+((576.0*rdx2SqVol1R3+576.0*rdx2SqVol[0]*rdx2SqVol1R2)*rhoC[1]+(168.0*phiLy[1]-24.0*phiLx[1]-336.0*phiC[1]-24.0*phiLxLy[0]+48.0*phiLx[0])*rdx2SqVol1R4+(120.0*rdx2SqVol[0]*phiLy[1]+24.0*rdx2SqVol[0]*phiLx[1]-960.0*rdx2SqVol[0]*phiC[1]-144.0*rhoC[0]+(432.0*phiLy[0]+24.0*phiLxLy[0]-192.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-48.0*rdx2SqVol0R2*phiLy[1])+48.0*rdx2SqVol0R2*phiLx[1]-192.0*rdx2SqVol0R2*phiC[1]+288.0*rdx2SqVol[0]*rhoC[0]+(48.0*phiLxLy[0]+192.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2)*bcVals4R2+((-384.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3)-384.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2)*omega+((49.0*phiC[1]*rdx2SqVol1R4+700.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+700.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiC[1])*bcVals4R2+(448.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1]+280.0*rdx2SqVol0R4*phiC[1])*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+960.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1]+336.0*rdx2SqVol0R4*phiC[1])*bcVals3R2)*bcVals10R2+((280.0*phiC[1]*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+448.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1])*bcVals4R2+(1104.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2+384.0*rdx2SqVol0R3*phiC[1]*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*phiC[1]*rdx2SqVol1R4+960.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*phiC[1]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*phiC[1]*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2)/(((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2); 
  phiC[2] = (((((84.0*rdx2SqVol1R4+1584.0*rdx2SqVol[0]*rdx2SqVol1R3+1836.0*rdx2SqVol0R2*rdx2SqVol1R2+336.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1272.0*rdx2SqVol[0]*rdx2SqVol1R3+2688.0*rdx2SqVol0R2*rdx2SqVol1R2+984.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(576.0*rdx2SqVol0R2*rdx2SqVol1R2+576.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals[9])*bcVals[11]+((((-126.0*rdx2SqVol[0]*rdx2SqVol1R3)+954.0*rdx2SqVol0R2*rdx2SqVol1R2+1206.0*rdx2SqVol0R3*rdx2SqVol[1]+126.0*rdx2SqVol0R4)*bcVals[4]+((-48.0*rdx2SqVol0R2*rdx2SqVol1R2)+552.0*rdx2SqVol0R3*rdx2SqVol[1]+168.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+(((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[3]+(168.0*rdx2SqVol1R3+1152.0*rdx2SqVol[0]*rdx2SqVol1R2+1152.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiC[2]+((-42.0*rdx2SqVol1R3)+522.0*rdx2SqVol[0]*rdx2SqVol1R2+522.0*rdx2SqVol0R2*rdx2SqVol[1]-42.0*rdx2SqVol0R3)*rhoC[1]+((-14.0*phiLx[1])+49.0*phiLy[0]+14.0*phiLxLy[0])*rdx2SqVol1R4+(252.0*rdx2SqVol[0]*phiLy[1]-65.0*rdx2SqVol[0]*phiLx[1]+84.0*rhoC[0]+(280.0*phiLy[0]+71.0*phiLxLy[0]+162.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(180.0*rdx2SqVol0R2*phiLy[1]+285.0*rdx2SqVol0R2*phiLx[1]+342.0*rdx2SqVol[0]*rhoC[0]+(204.0*phiLy[0]+93.0*phiLxLy[0]+540.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-72.0*rdx2SqVol0R3*phiLy[1])+385.0*rdx2SqVol0R3*phiLx[1]+216.0*rdx2SqVol0R2*rhoC[0]+((-20.0*phiLy[0])+29.0*phiLxLy[0]+378.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]+(7.0*phiLy[0]-7.0*phiLxLy[0])*rdx2SqVol0R4)*bcVals4R2+(((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)-192.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((1368.0*rdx2SqVol[0]*rdx2SqVol1R2+2880.0*rdx2SqVol0R2*rdx2SqVol[1]+648.0*rdx2SqVol0R3)*rhoC[2]+((-448.0*rdx2SqVol[0]*rdx2SqVol1R3)-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-2592.0*rdx2SqVol0R3*rdx2SqVol[1]-280.0*rdx2SqVol0R4)*phiC[2]+(48.0*rdx2SqVol[0]*rdx2SqVol1R2+744.0*rdx2SqVol0R2*rdx2SqVol[1]-168.0*rdx2SqVol0R3)*rhoC[1]+(126.0*rdx2SqVol[0]*phiLy[1]-118.0*rdx2SqVol[0]*phiLx[1]+(448.0*phiLy[0]+110.0*phiLxLy[0]+8.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(390.0*rdx2SqVol0R2*phiLy[1]+126.0*rdx2SqVol0R2*phiLx[1]+660.0*rdx2SqVol[0]*rhoC[0]+(72.0*phiLy[0]+150.0*phiLxLy[0]+1068.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-174.0*rdx2SqVol0R3*phiLy[1])+858.0*rdx2SqVol0R3*phiLx[1]+240.0*rdx2SqVol0R2*rhoC[0]+((-120.0*phiLy[0])+6.0*phiLxLy[0]+816.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-6.0*rdx2SqVol0R4*phiLy[1]+182.0*rdx2SqVol0R4*phiLx[1]-204.0*rdx2SqVol0R3*rhoC[0]+(40.0*phiLy[0]-34.0*phiLxLy[0]-28.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((576.0*rdx2SqVol0R2*rdx2SqVol[1]+576.0*rdx2SqVol0R3)*rhoC[2]+((-192.0*rdx2SqVol0R2*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*rdx2SqVol[1]-336.0*rdx2SqVol0R4)*phiC[2]+(48.0*rdx2SqVol0R2*phiLy[1]-48.0*rdx2SqVol0R2*phiLx[1]+(192.0*phiLy[0]+48.0*phiLxLy[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(24.0*rdx2SqVol0R3*phiLy[1]+120.0*rdx2SqVol0R3*phiLx[1]+288.0*rdx2SqVol0R2*rhoC[0]+((-192.0*phiLy[0])+24.0*phiLxLy[0]+432.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]-24.0*rdx2SqVol0R4*phiLy[1]+168.0*rdx2SqVol0R4*phiLx[1]-144.0*rdx2SqVol0R3*rhoC[0]+(48.0*phiLy[0]-24.0*phiLxLy[0])*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((((-504.0*rdx2SqVol[0]*rdx2SqVol1R3)-216.0*rdx2SqVol0R2*rdx2SqVol1R2-144.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[4]+((-192.0*rdx2SqVol0R2*rdx2SqVol1R2)-192.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3])*bcVals[5]+(((-336.0*rdx2SqVol1R3)-960.0*rdx2SqVol[0]*rdx2SqVol1R2-192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+(672.0*rdx2SqVol1R3+1920.0*rdx2SqVol[0]*rdx2SqVol1R2+384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-280.0*rdx2SqVol1R4)-2592.0*rdx2SqVol[0]*rdx2SqVol1R3-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[2]+((-168.0*rdx2SqVol1R3)+744.0*rdx2SqVol[0]*rdx2SqVol1R2+48.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+((-56.0*phiLx[1])+196.0*phiLy[0]+56.0*phiLxLy[0])*rdx2SqVol1R4+(336.0*rdx2SqVol[0]*phiLy[1]-36.0*rdx2SqVol[0]*phiLx[1]+336.0*rhoC[0]+(60.0*phiLxLy[0]+648.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-120.0*rdx2SqVol0R2*phiLy[1])+564.0*rdx2SqVol0R2*phiLx[1]+24.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLy[0]-12.0*phiLxLy[0]+432.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-24.0*rdx2SqVol0R3*phiLy[1])+112.0*rdx2SqVol0R3*phiLx[1]-96.0*rdx2SqVol0R2*rhoC[0]+(40.0*phiLy[0]-16.0*phiLxLy[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals4R2+(((-384.0*rdx2SqVol[0]*rdx2SqVol1R2)-384.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[3]*rhoC[3]+((768.0*rdx2SqVol[0]*rdx2SqVol1R2+768.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-1104.0*rdx2SqVol[0]*rdx2SqVol1R3)-3072.0*rdx2SqVol0R2*rdx2SqVol1R2-1104.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[2]+(192.0*rdx2SqVol[0]*rdx2SqVol1R2+192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(112.0*rdx2SqVol[0]*phiLy[1]-80.0*rdx2SqVol[0]*phiLx[1]+(224.0*phiLy[0]+48.0*phiLxLy[0]+32.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(56.0*rdx2SqVol0R2*phiLy[1]+152.0*rdx2SqVol0R2*phiLx[1]+288.0*rdx2SqVol[0]*rhoC[0]+((-128.0*phiLy[0])+24.0*phiLxLy[0]+496.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-56.0*rdx2SqVol0R3*phiLy[1])+232.0*rdx2SqVol0R3*phiLx[1]-144.0*rdx2SqVol0R2*rhoC[0]+(80.0*phiLy[0]-24.0*phiLxLy[0]+32.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[3])*bcVals[4]+((-384.0*rdx2SqVol0R2*rdx2SqVol1R2)-384.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[2]*bcVals3R2)*bcVals[9]*bcVals[10]+(((-336.0*rdx2SqVol1R4)-960.0*rdx2SqVol[0]*rdx2SqVol1R3-192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiC[2]*bcVals4R2+((-384.0*rdx2SqVol[0]*rdx2SqVol1R3)-384.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiC[2]*bcVals[3]*bcVals[4])*bcVals9R2)*omega+((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiC[2]*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*phiC[2]*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*phiC[2]*bcVals3R2)*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[2]*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[2]*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[2]*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiC[2]*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiC[2]*bcVals[3]*bcVals[4])*bcVals9R2)/(((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2); 
  phiC[3] = (((((252.0*rdx2SqVol1R4+2736.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+504.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(288.0*rdx2SqVol[0]*rdx2SqVol1R3+1728.0*rdx2SqVol0R2*rdx2SqVol1R2+1008.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4])*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals[9])*bcVals[11]+(((504.0*rdx2SqVol[0]*rdx2SqVol1R3+2988.0*rdx2SqVol0R2*rdx2SqVol1R2+2736.0*rdx2SqVol0R3*rdx2SqVol[1]+252.0*rdx2SqVol0R4)*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals[3])*bcVals[5]+((336.0*rdx2SqVol1R3+2304.0*rdx2SqVol[0]*rdx2SqVol1R2+2304.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[3]+((-49.0*rdx2SqVol1R4)-700.0*rdx2SqVol[0]*rdx2SqVol1R3-1302.0*rdx2SqVol0R2*rdx2SqVol1R2-700.0*rdx2SqVol0R3*rdx2SqVol[1]-49.0*rdx2SqVol0R4)*phiC[3]+((-84.0*rdx2SqVol1R3)+432.0*rdx2SqVol[0]*rdx2SqVol1R2+684.0*rdx2SqVol0R2*rdx2SqVol[1]+168.0*rdx2SqVol0R3)*rhoC[2]+(168.0*rdx2SqVol1R3+684.0*rdx2SqVol[0]*rdx2SqVol1R2+432.0*rdx2SqVol0R2*rdx2SqVol[1]-84.0*rdx2SqVol0R3)*rhoC[1]+(49.0*phiLy[1]+7.0*phiLx[1]-7.0*phiLxLy[0])*rdx2SqVol1R4+(280.0*rdx2SqVol[0]*phiLy[1]-20.0*rdx2SqVol[0]*phiLx[1]-42.0*rhoC[0]+(504.0*phiLy[0]+80.0*phiLxLy[0]-144.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(204.0*rdx2SqVol0R2*phiLy[1]+204.0*rdx2SqVol0R2*phiLx[1]+522.0*rdx2SqVol[0]*rhoC[0]+(360.0*phiLy[0]+174.0*phiLxLy[0]+360.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+((-20.0*rdx2SqVol0R3*phiLy[1])+280.0*rdx2SqVol0R3*phiLx[1]+522.0*rdx2SqVol0R2*rhoC[0]+((-144.0*phiLy[0])+80.0*phiLxLy[0]+504.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+7.0*rdx2SqVol0R4*phiLy[1]+49.0*rdx2SqVol0R4*phiLx[1]-42.0*rdx2SqVol0R3*rhoC[0]-7.0*phiLxLy[0]*rdx2SqVol0R4)*bcVals4R2+((384.0*rdx2SqVol[0]*rdx2SqVol1R2+1920.0*rdx2SqVol0R2*rdx2SqVol[1]+672.0*rdx2SqVol0R3)*bcVals[3]*rhoC[3]+((-448.0*rdx2SqVol[0]*rdx2SqVol1R3)-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-2592.0*rdx2SqVol0R3*rdx2SqVol[1]-280.0*rdx2SqVol0R4)*bcVals[3]*phiC[3]+(((-96.0*rdx2SqVol[0]*rdx2SqVol1R2)-192.0*rdx2SqVol0R2*rdx2SqVol[1]+336.0*rdx2SqVol0R3)*rhoC[2]+((-192.0*rdx2SqVol[0]*rdx2SqVol1R2)-960.0*rdx2SqVol0R2*rdx2SqVol[1]-336.0*rdx2SqVol0R3)*rhoC[1]+((-56.0*rdx2SqVol[0]*phiLy[1])+24.0*rdx2SqVol[0]*phiLx[1]+(8.0*phiLxLy[0]-32.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-228.0*rdx2SqVol0R2*phiLy[1])+60.0*rdx2SqVol0R2*phiLx[1]+48.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLxLy[0]-120.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(72.0*rdx2SqVol0R3*phiLy[1]-96.0*rdx2SqVol0R3*phiLx[1]+312.0*rdx2SqVol0R2*rhoC[0]+((-432.0*phiLy[0])+24.0*phiLxLy[0]+288.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1]+28.0*rdx2SqVol0R4*phiLy[1]+84.0*rdx2SqVol0R4*phiLx[1]-168.0*rdx2SqVol0R3*rhoC[0]+((-28.0*phiLxLy[0])-56.0*phiLx[0])*rdx2SqVol0R4)*bcVals[3])*bcVals[4]+((-192.0*rdx2SqVol0R2*rdx2SqVol1R2)-960.0*rdx2SqVol0R3*rdx2SqVol[1]-336.0*rdx2SqVol0R4)*bcVals3R2*phiC[3])*bcVals10R2+(((1008.0*rdx2SqVol[0]*rdx2SqVol1R3+1728.0*rdx2SqVol0R2*rdx2SqVol1R2+288.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3])*bcVals[5]+((672.0*rdx2SqVol1R3+1920.0*rdx2SqVol[0]*rdx2SqVol1R2+384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[3]+((-280.0*rdx2SqVol1R4)-2592.0*rdx2SqVol[0]*rdx2SqVol1R3-2760.0*rdx2SqVol0R2*rdx2SqVol1R2-448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[3]+((-336.0*rdx2SqVol1R3)-960.0*rdx2SqVol[0]*rdx2SqVol1R2-192.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+(336.0*rdx2SqVol1R3-192.0*rdx2SqVol[0]*rdx2SqVol1R2-96.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+(84.0*phiLy[1]+28.0*phiLx[1]-56.0*phiLy[0]-28.0*phiLxLy[0])*rdx2SqVol1R4+((-96.0*rdx2SqVol[0]*phiLy[1])+72.0*rdx2SqVol[0]*phiLx[1]-168.0*rhoC[0]+(288.0*phiLy[0]+24.0*phiLxLy[0]-432.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+(60.0*rdx2SqVol0R2*phiLy[1]-228.0*rdx2SqVol0R2*phiLx[1]+312.0*rdx2SqVol[0]*rhoC[0]+(60.0*phiLxLy[0]-120.0*phiLy[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(24.0*rdx2SqVol0R3*phiLy[1]-56.0*rdx2SqVol0R3*phiLx[1]+48.0*rdx2SqVol0R2*rhoC[0]+(8.0*phiLxLy[0]-32.0*phiLy[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals4R2+((768.0*rdx2SqVol[0]*rdx2SqVol1R2+768.0*rdx2SqVol0R2*rdx2SqVol[1])*bcVals[3]*rhoC[3]+((-1104.0*rdx2SqVol[0]*rdx2SqVol1R3)-3072.0*rdx2SqVol0R2*rdx2SqVol1R2-1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*phiC[3]+(((-384.0*rdx2SqVol[0]*rdx2SqVol1R2)-384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[2]+((-384.0*rdx2SqVol[0]*rdx2SqVol1R2)-384.0*rdx2SqVol0R2*rdx2SqVol[1])*rhoC[1]+((-128.0*rdx2SqVol[0]*phiLy[1])+64.0*rdx2SqVol[0]*phiLx[1]+((-64.0*phiLy[0])-64.0*phiLx[0])*rdx2SqVol[0])*rdx2SqVol1R3+((-64.0*rdx2SqVol0R2*phiLy[1])-64.0*rdx2SqVol0R2*phiLx[1]+((-128.0*phiLy[0])-128.0*phiLx[0])*rdx2SqVol0R2)*rdx2SqVol1R2+(64.0*rdx2SqVol0R3*phiLy[1]-128.0*rdx2SqVol0R3*phiLx[1]+((-64.0*phiLy[0])-64.0*phiLx[0])*rdx2SqVol0R3)*rdx2SqVol[1])*bcVals[3])*bcVals[4]+((-384.0*rdx2SqVol0R2*rdx2SqVol1R2)-384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2*phiC[3])*bcVals[9]*bcVals[10]+(((-336.0*rdx2SqVol1R4)-960.0*rdx2SqVol[0]*rdx2SqVol1R3-192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiC[3]*bcVals4R2+((-384.0*rdx2SqVol[0]*rdx2SqVol1R3)-384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*phiC[3]*bcVals[4])*bcVals9R2)*omega+((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*phiC[3]*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*phiC[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2*phiC[3])*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*phiC[3]*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*phiC[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2*phiC[3])*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*phiC[3]*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*phiC[3]*bcVals[4])*bcVals9R2)/(((49.0*rdx2SqVol1R4+700.0*rdx2SqVol[0]*rdx2SqVol1R3+1302.0*rdx2SqVol0R2*rdx2SqVol1R2+700.0*rdx2SqVol0R3*rdx2SqVol[1]+49.0*rdx2SqVol0R4)*bcVals4R2+(448.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+2592.0*rdx2SqVol0R3*rdx2SqVol[1]+280.0*rdx2SqVol0R4)*bcVals[3]*bcVals[4]+(192.0*rdx2SqVol0R2*rdx2SqVol1R2+960.0*rdx2SqVol0R3*rdx2SqVol[1]+336.0*rdx2SqVol0R4)*bcVals3R2)*bcVals10R2+((280.0*rdx2SqVol1R4+2592.0*rdx2SqVol[0]*rdx2SqVol1R3+2760.0*rdx2SqVol0R2*rdx2SqVol1R2+448.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals4R2+(1104.0*rdx2SqVol[0]*rdx2SqVol1R3+3072.0*rdx2SqVol0R2*rdx2SqVol1R2+1104.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals[3]*bcVals[4]+(384.0*rdx2SqVol0R2*rdx2SqVol1R2+384.0*rdx2SqVol0R3*rdx2SqVol[1])*bcVals3R2)*bcVals[9]*bcVals[10]+((336.0*rdx2SqVol1R4+960.0*rdx2SqVol[0]*rdx2SqVol1R3+192.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals4R2+(384.0*rdx2SqVol[0]*rdx2SqVol1R3+384.0*rdx2SqVol0R2*rdx2SqVol1R2)*bcVals[3]*bcVals[4])*bcVals9R2); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((4.0*phiUy[0]+phiUxUy[0]+phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+6.0*rhoC[0]+((-2.0*phiUy[0])+phiUxUy[0]+phiUxLy[0]+4.0*phiUx[0]-2.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]-1.0*phiUxLy[0]+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(phiUy[0]-1.0*phiUxUy[0]-1.0*phiUxLy[0]-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]-1.0*phiUxLy[0]+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*bcVals[1]*rdx2SqVol[1]+((phiUy[0]-1.0*phiUxUy[0]-1.0*phiUxLy[0]-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiC[0])*rdx2SqVol[0]))/bcVals[1]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*bcVals[5]+(phiLy[1]-2.0*phiC[1]+4.0*phiUy[0]+4.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiUy[0])-2.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((phiUy[1]+phiLy[1]-2.0*phiC[1]+4.0*phiUy[0]+4.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+rdx2SqVol[0]*phiLy[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiUy[0])-2.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+6.0*rhoC[1]+(2.0*phiUy[1]+2.0*phiLy[1]-4.0*phiC[1]+phiUy[0]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiUy[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((phiUy[1]+phiLy[1]-2.0*phiC[1]+4.0*phiUy[0]+4.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+rdx2SqVol[0]*phiLy[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiUy[0])-2.0*phiLy[0]+phiLxUy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+(6.0*rhoC[1]+(2.0*phiUy[1]+2.0*phiLy[1]-4.0*phiC[1]+phiUy[0]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiUy[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[4]+((-2.0*rdx2SqVol[0]*phiUy[1])-4.0*rdx2SqVol[0]*phiC[1])*bcVals[3]))/bcVals[4]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+((-4.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+(((-4.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*bcVals[11]+(phiLx[1]+4.0*phiC[1]+phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLx[1]-2.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-2.0*phiLy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((phiUx[1]+phiLx[1]+4.0*phiC[1]+phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUx[1]+rdx2SqVol[0]*phiLx[1]-2.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-2.0*phiLy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+6.0*rhoC[1]+((-1.0*phiUx[1])-1.0*phiLx[1]-4.0*phiC[1]+phiUx[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiUx[1]+2.0*rdx2SqVol[0]*phiLx[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiUx[0]+phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((phiUx[1]+phiLx[1]+4.0*phiC[1]+phiUxLy[0]-2.0*phiUx[0]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUx[1]+rdx2SqVol[0]*phiLx[1]-2.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-2.0*phiLy[0]+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+(6.0*rhoC[1]+((-1.0*phiUx[1])-1.0*phiLx[1]-4.0*phiC[1]+phiUx[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiUx[1]+2.0*rdx2SqVol[0]*phiLx[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiUx[0]+phiLx[0]-2.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((-2.0*phiUx[1])-4.0*phiC[1])*rdx2SqVol[1]*bcVals[9]))/bcVals[10]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*rdx2SqVol[1]-6.0*rhoC[0]+(phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -(0.1666666666666667*(6.0*bcVals[1]*rdx2SqVol[1]*bcVals[8]+6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*bcVals[1]*rdx2SqVol[1]+((phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiC[0])*rdx2SqVol[0]))/bcVals[1]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -(0.1666666666666667*(6.0*bcVals[1]*rdx2SqVol[1]*bcVals[8]+(6.0*rdx2SqVol[0]*bcVals[2]+((-2.0*phiUy[0])-1.0*phiUxUy[0]+phiUx[0]+2.0*phiC[0])*bcVals[1]*rdx2SqVol[1]+((phiUy[0]-1.0*phiUxUy[0]-2.0*phiUx[0]+2.0*phiC[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]+((-2.0*bcVals[0]*phiUy[0])-4.0*bcVals[0]*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiUx[0])-4.0*phiC[0])*bcVals[1]*rdx2SqVol[1]*bcVals[6]))/(bcVals[1]*bcVals[7]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+(2.0*phiC[1]+phiUxLy[0]-2.0*phiUx[0]+2.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-1.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-1.0*phiUx[1])-2.0*phiC[1]-1.0*phiUxLy[0]+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUx[1]+rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+((-1.0*phiUxLy[0])-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+6.0*rhoC[1]+((-1.0*phiUx[1])-2.0*phiC[1]+phiUx[0]+2.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiUx[1]-2.0*rdx2SqVol[0]*phiC[1]+(phiUx[0]-1.0*phiC[0])*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-1.0*phiUx[1])-2.0*phiC[1]-1.0*phiUxLy[0]+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUx[1]+rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+((-1.0*phiUxLy[0])-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+((-6.0*rdx2SqVol[0]*bcVals[2])+6.0*rhoC[1]+((-1.0*phiUx[1])-2.0*phiC[1]+phiUx[0]+2.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiUx[1]-2.0*rdx2SqVol[0]*phiC[1]+(phiUx[0]-1.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((-2.0*phiUx[1])-4.0*phiC[1])*rdx2SqVol[1]*bcVals[9]))/bcVals[10]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = (0.1666666666666667*((bcVals[1]*rdx2SqVol[1]+rdx2SqVol[0]*bcVals[1])*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+(2.0*bcVals[1]*phiC[1]+(phiUxLy[0]-2.0*phiUx[0]+2.0*phiLy[0]-4.0*phiC[0])*bcVals[1])*rdx2SqVol[1]+(2.0*bcVals[0]*rdx2SqVol[0]-1.0*rdx2SqVol[0]*bcVals[1])*phiC[1]+(6.0*rhoC[0]+(phiUxLy[0]+4.0*phiUx[0]-1.0*phiLy[0]-4.0*phiC[0])*rdx2SqVol[0])*bcVals[1]+4.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/bcVals[1]; 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-1.0*bcVals[1]*phiUx[1])-2.0*bcVals[1]*phiC[1]+((-1.0*phiUxLy[0])+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*bcVals[1])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*bcVals[1]*phiUx[1]+(rdx2SqVol[0]*bcVals[1]-2.0*bcVals[0]*rdx2SqVol[0])*phiC[1]+(((-1.0*phiUxLy[0])-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]-4.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/bcVals[1]; 
  resOut[1] = (0.1666666666666667*(6.0*bcVals[1]*rdx2SqVol[1]*bcVals[11]-6.0*rdx2SqVol[0]*bcVals[2]+6.0*bcVals[1]*rhoC[1]+((-1.0*bcVals[1]*phiUx[1])-2.0*bcVals[1]*phiC[1]+(phiUx[0]+2.0*phiC[0])*bcVals[1])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*bcVals[1]*phiUx[1]+(4.0*bcVals[0]*rdx2SqVol[0]-2.0*rdx2SqVol[0]*bcVals[1])*phiC[1]+(phiUx[0]-1.0*phiC[0])*rdx2SqVol[0]*bcVals[1]+2.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/bcVals[1]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[2]+((-1.0*bcVals[1]*phiUx[1])-2.0*bcVals[1]*phiC[1]+((-1.0*phiUxLy[0])+2.0*phiUx[0]-2.0*phiLy[0]+4.0*phiC[0])*bcVals[1])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*bcVals[1]*phiUx[1]+(rdx2SqVol[0]*bcVals[1]-2.0*bcVals[0]*rdx2SqVol[0])*phiC[1]+(((-1.0*phiUxLy[0])-4.0*phiUx[0]+phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]-6.0*rhoC[0])*bcVals[1]-4.0*bcVals[0]*phiC[0]*rdx2SqVol[0]))/bcVals[1]; 
  resOut[1] = (0.1666666666666667*(6.0*bcVals[1]*rdx2SqVol[1]*bcVals[11]+((-6.0*rdx2SqVol[0]*bcVals[2])+6.0*bcVals[1]*rhoC[1]+((-1.0*bcVals[1]*phiUx[1])-2.0*bcVals[1]*phiC[1]+(phiUx[0]+2.0*phiC[0])*bcVals[1])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*bcVals[1]*phiUx[1]+(4.0*bcVals[0]*rdx2SqVol[0]-2.0*rdx2SqVol[0]*bcVals[1])*phiC[1]+(phiUx[0]-1.0*phiC[0])*rdx2SqVol[0]*bcVals[1]+2.0*bcVals[0]*phiC[0]*rdx2SqVol[0])*bcVals[10]+((-2.0*bcVals[1]*phiUx[1])-4.0*bcVals[1]*phiC[1])*rdx2SqVol[1]*bcVals[9]))/(bcVals[1]*bcVals[10]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+((-1.0*rdx2SqVol[1])-1.0*rdx2SqVol[0])*bcVals[5]+(phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-2.0*rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+(((-1.0*rdx2SqVol[1])-1.0*rdx2SqVol[0])*bcVals[5]+(phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-2.0*rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiC[1])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+((-1.0*phiUy[1])+phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-2.0*rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]-6.0*rdx2SqVol[0]*bcVals[5]-6.0*rhoC[1]+((-2.0*phiUy[1])+2.0*phiC[1]-1.0*phiUy[0]+phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+2.0*rdx2SqVol[0]*phiC[1]+((-1.0*phiUy[0])-2.0*phiC[0])*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+(((-1.0*phiUy[1])+phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-2.0*rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiC[1])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 
  resOut[1] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+((-6.0*rdx2SqVol[0]*bcVals[5])-6.0*rhoC[1]+((-2.0*phiUy[1])+2.0*phiC[1]-1.0*phiUy[0]+phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+2.0*rdx2SqVol[0]*phiC[1]+((-1.0*phiUy[0])-2.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-4.0*phiC[1])-2.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+((-1.0*phiUy[1])+phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-2.0*rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[4]*bcVals[8]-6.0*rdx2SqVol[0]*bcVals[5]+((-6.0*rhoC[1])+((-2.0*phiUy[1])+2.0*phiC[1]-1.0*phiUy[0]+phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+2.0*rdx2SqVol[0]*phiC[1]+((-1.0*phiUy[0])-2.0*phiC[0])*rdx2SqVol[0])*bcVals[4]+(2.0*rdx2SqVol[0]*phiUy[1]+4.0*rdx2SqVol[0]*phiC[1])*bcVals[3]))/bcVals[4]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[8]+(((-1.0*phiUy[1])+phiC[1]-4.0*phiUy[0]-1.0*phiLxUy[0]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiUy[1]-2.0*rdx2SqVol[0]*phiC[1]-6.0*rhoC[0]+(2.0*phiUy[0]-1.0*phiLxUy[0]-2.0*phiLx[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[7]+((-2.0*phiC[1])-4.0*phiC[0])*rdx2SqVol[1]*bcVals[6]))/bcVals[7]; 
  resOut[1] = -(0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[4]*bcVals[8]+((-6.0*rdx2SqVol[0]*bcVals[5])+((-6.0*rhoC[1])+((-2.0*phiUy[1])+2.0*phiC[1]-1.0*phiUy[0]+phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiUy[1]+2.0*rdx2SqVol[0]*phiC[1]+((-1.0*phiUy[0])-2.0*phiC[0])*rdx2SqVol[0])*bcVals[4]+(2.0*rdx2SqVol[0]*phiUy[1]+4.0*rdx2SqVol[0]*phiC[1])*bcVals[3])*bcVals[7]+((-4.0*phiC[1])-2.0*phiC[0])*rdx2SqVol[1]*bcVals[4]*bcVals[6]))/(bcVals[4]*bcVals[7]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[5]; 
  resOut[2] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+6.0*rhoC[2]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = phiC[1]-1.0*bcVals[5]; 
  resOut[2] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+((2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+6.0*rhoC[2]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((-2.0*rdx2SqVol[1]*phiC[3])-4.0*rdx2SqVol[1]*phiC[2])*bcVals[9]))/bcVals[10]; 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+(2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+6.0*rhoC[1]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+(2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+6.0*rhoC[1]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  resOut[2] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+6.0*rhoC[2]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  resOut[3] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+6.0*rdx2SqVol[0]*bcVals[5]+6.0*rhoC[3]+((-2.0*rdx2SqVol[1])-2.0*rdx2SqVol[0])*phiC[3]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[2]+(2.0*phiC[1]+phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0]); 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = 0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+(2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+6.0*rhoC[1]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0]); 
  resOut[2] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+((2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+6.0*rhoC[2]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((-2.0*rdx2SqVol[1]*phiC[3])-4.0*rdx2SqVol[1]*phiC[2])*bcVals[9]))/bcVals[10]; 
  resOut[3] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+(6.0*rdx2SqVol[0]*bcVals[5]+6.0*rhoC[3]+((-2.0*rdx2SqVol[1])-2.0*rdx2SqVol[0])*phiC[3]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[2]+(2.0*phiC[1]+phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0])*bcVals[10]+((-4.0*rdx2SqVol[1]*phiC[3])-2.0*rdx2SqVol[1]*phiC[2])*bcVals[9]))/bcVals[10]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+((2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+6.0*rhoC[1]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[4]-2.0*rdx2SqVol[0]*bcVals[3]*phiC[3]-4.0*rdx2SqVol[0]*phiC[1]*bcVals[3]))/bcVals[4]; 
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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+((2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+6.0*rhoC[1]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[4]-2.0*rdx2SqVol[0]*bcVals[3]*phiC[3]-4.0*rdx2SqVol[0]*phiC[1]*bcVals[3]))/bcVals[4]; 
  resOut[2] = 0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+6.0*rhoC[2]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0]); 
  resOut[3] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[4]*bcVals[11]+6.0*rdx2SqVol[0]*bcVals[5]+(6.0*rhoC[3]+((-2.0*rdx2SqVol[1])-2.0*rdx2SqVol[0])*phiC[3]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[2]+(2.0*phiC[1]+phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0])*bcVals[4]-4.0*rdx2SqVol[0]*bcVals[3]*phiC[3]-2.0*rdx2SqVol[0]*phiC[1]*bcVals[3]))/bcVals[4]; 

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
  double *phiUx = phi[2]; 
  double *phiLy = phi[3]; 
  double *phiUy = phi[4]; 
  double *phiLxLy = phi[5]; 
  double *phiLxUy = phi[6]; 
  double *phiUxLy = phi[7]; 
  double *phiUxUy = phi[8]; 


  resOut[0] = 0.1666666666666667*((rdx2SqVol[1]+rdx2SqVol[0])*phiC[3]+(4.0*rdx2SqVol[1]-2.0*rdx2SqVol[0])*phiC[2]+(phiLy[1]+phiLx[1]-2.0*phiC[1]+4.0*phiLy[0]+phiLxLy[0]-2.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[1]+rdx2SqVol[0]*phiLy[1]+rdx2SqVol[0]*phiLx[1]+4.0*rdx2SqVol[0]*phiC[1]+6.0*rhoC[0]+((-2.0*phiLy[0])+phiLxLy[0]+4.0*phiLx[0]-8.0*phiC[0])*rdx2SqVol[0]); 
  resOut[1] = (0.1666666666666667*(6.0*rdx2SqVol[0]*bcVals[5]+((2.0*rdx2SqVol[1]-1.0*rdx2SqVol[0])*phiC[3]+(rdx2SqVol[1]+rdx2SqVol[0])*phiC[2]+6.0*rhoC[1]+(2.0*phiLy[1]-4.0*phiC[1]+phiLy[0]-2.0*phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiLy[1]-4.0*rdx2SqVol[0]*phiC[1]+(phiLy[0]+4.0*phiC[0])*rdx2SqVol[0])*bcVals[4]-2.0*rdx2SqVol[0]*bcVals[3]*phiC[3]-4.0*rdx2SqVol[0]*phiC[1]*bcVals[3]))/bcVals[4]; 
  resOut[2] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[11]+((2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[3]+6.0*rhoC[2]+((-4.0*rdx2SqVol[1])-4.0*rdx2SqVol[0])*phiC[2]+((-1.0*phiLx[1])+phiC[1]+phiLx[0]+4.0*phiC[0])*rdx2SqVol[1]+2.0*rdx2SqVol[0]*phiLx[1]+rdx2SqVol[0]*phiC[1]+(phiLx[0]-2.0*phiC[0])*rdx2SqVol[0])*bcVals[10]+((-2.0*rdx2SqVol[1]*phiC[3])-4.0*rdx2SqVol[1]*phiC[2])*bcVals[9]))/bcVals[10]; 
  resOut[3] = (0.1666666666666667*(6.0*rdx2SqVol[1]*bcVals[4]*bcVals[11]+(6.0*rdx2SqVol[0]*bcVals[5]+(6.0*rhoC[3]+((-2.0*rdx2SqVol[1])-2.0*rdx2SqVol[0])*phiC[3]+(2.0*rdx2SqVol[0]-1.0*rdx2SqVol[1])*phiC[2]+(2.0*phiC[1]+phiC[0])*rdx2SqVol[1]-1.0*rdx2SqVol[0]*phiC[1]+phiC[0]*rdx2SqVol[0])*bcVals[4]-4.0*rdx2SqVol[0]*bcVals[3]*phiC[3]-2.0*rdx2SqVol[0]*phiC[1]*bcVals[3])*bcVals[10]+((-4.0*rdx2SqVol[1]*phiC[3])-2.0*rdx2SqVol[1]*phiC[2])*bcVals[4]*bcVals[9]))/(bcVals[4]*bcVals[10]); 

}

void MGpoissonFEML2norm2xSer_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 
  double *femFldUy = femFld[2]; 
  double *femFldUxUy = femFld[3]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldUx0R2 = std::pow(femFldUx[0],2);
  const double femFldUy0R2 = std::pow(femFldUy[0],2);
  const double femFldUxUy0R2 = std::pow(femFldUxUy[0],2);

  normOut[0] += 0.1111111111111111*(4.0*femFldUy0R2+(4.0*femFldUxUy[0]+2.0*femFldUx[0]+4.0*femFldC[0])*femFldUy[0]+4.0*femFldUxUy0R2+(4.0*femFldUx[0]+2.0*femFldC[0])*femFldUxUy[0]+4.0*femFldUx0R2+4.0*femFldC[0]*femFldUx[0]+4.0*femFldC0R2)*volFac; 
}

void MGpoissonFEML2norm2xSer_LxNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 
  double *femFldUy = femFld[2]; 
  double *femFldUxUy = femFld[3]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldUx0R2 = std::pow(femFldUx[0],2);
  const double femFldUy0R2 = std::pow(femFldUy[0],2);
  const double femFldUxUy0R2 = std::pow(femFldUxUy[0],2);

  normOut[0] += 0.1111111111111111*(4.0*femFldUy0R2+(4.0*femFldUxUy[0]+2.0*femFldUx[0]+4.0*femFldC[0])*femFldUy[0]+4.0*femFldUxUy0R2+(4.0*femFldUx[0]+2.0*femFldC[0])*femFldUxUy[0]+4.0*femFldUx0R2+4.0*femFldC[0]*femFldUx[0]+4.0*femFldC0R2)*volFac; 
}

void MGpoissonFEML2norm2xSer_UxNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUy = femFld[2]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldC1R2 = std::pow(femFldC[1],2);
  const double femFldUy0R2 = std::pow(femFldUy[0],2);
  const double femFldUy1R2 = std::pow(femFldUy[1],2);

  normOut[0] += 0.1111111111111111*(4.0*femFldUy1R2+(4.0*femFldC[1]+4.0*femFldUy[0]+2.0*femFldC[0])*femFldUy[1]+4.0*femFldC1R2+(2.0*femFldUy[0]+4.0*femFldC[0])*femFldC[1]+4.0*femFldUy0R2+4.0*femFldC[0]*femFldUy[0]+4.0*femFldC0R2)*volFac; 
}

void MGpoissonFEML2norm2xSer_LyNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 
  double *femFldUy = femFld[2]; 
  double *femFldUxUy = femFld[3]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldUx0R2 = std::pow(femFldUx[0],2);
  const double femFldUy0R2 = std::pow(femFldUy[0],2);
  const double femFldUxUy0R2 = std::pow(femFldUxUy[0],2);

  normOut[0] += 0.1111111111111111*(4.0*femFldUy0R2+(4.0*femFldUxUy[0]+2.0*femFldUx[0]+4.0*femFldC[0])*femFldUy[0]+4.0*femFldUxUy0R2+(4.0*femFldUx[0]+2.0*femFldC[0])*femFldUxUy[0]+4.0*femFldUx0R2+4.0*femFldC[0]*femFldUx[0]+4.0*femFldC0R2)*volFac; 
}

void MGpoissonFEML2norm2xSer_UyNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldUx0R2 = std::pow(femFldUx[0],2);
  const double femFldC1R2 = std::pow(femFldC[1],2);
  const double femFldUx1R2 = std::pow(femFldUx[1],2);

  normOut[0] += 0.1111111111111111*(4.0*femFldUx1R2+(4.0*femFldC[1]+4.0*femFldUx[0]+2.0*femFldC[0])*femFldUx[1]+4.0*femFldC1R2+(2.0*femFldUx[0]+4.0*femFldC[0])*femFldC[1]+4.0*femFldUx0R2+4.0*femFldC[0]*femFldUx[0]+4.0*femFldC0R2)*volFac; 
}

void MGpoissonFEML2norm2xSer_LxNonPeriodicUyNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldUx0R2 = std::pow(femFldUx[0],2);
  const double femFldC1R2 = std::pow(femFldC[1],2);
  const double femFldUx1R2 = std::pow(femFldUx[1],2);

  normOut[0] += 0.1111111111111111*(4.0*femFldUx1R2+(4.0*femFldC[1]+4.0*femFldUx[0]+2.0*femFldC[0])*femFldUx[1]+4.0*femFldC1R2+(2.0*femFldUx[0]+4.0*femFldC[0])*femFldC[1]+4.0*femFldUx0R2+4.0*femFldC[0]*femFldUx[0]+4.0*femFldC0R2)*volFac; 
}

void MGpoissonFEML2norm2xSer_UxNonPeriodicLyNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUy = femFld[2]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldC1R2 = std::pow(femFldC[1],2);
  const double femFldUy0R2 = std::pow(femFldUy[0],2);
  const double femFldUy1R2 = std::pow(femFldUy[1],2);

  normOut[0] += 0.1111111111111111*(4.0*femFldUy1R2+(4.0*femFldC[1]+4.0*femFldUy[0]+2.0*femFldC[0])*femFldUy[1]+4.0*femFldC1R2+(2.0*femFldUy[0]+4.0*femFldC[0])*femFldC[1]+4.0*femFldUy0R2+4.0*femFldC[0]*femFldUy[0]+4.0*femFldC0R2)*volFac; 
}

void MGpoissonFEML2norm2xSer_UxNonPeriodicUyNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 

  const double femFldC0R2 = std::pow(femFldC[0],2);
  const double femFldC1R2 = std::pow(femFldC[1],2);
  const double femFldC2R2 = std::pow(femFldC[2],2);
  const double femFldC3R2 = std::pow(femFldC[3],2);

  normOut[0] += 0.1111111111111111*(4.0*femFldC3R2+(4.0*femFldC[2]+4.0*femFldC[1]+2.0*femFldC[0])*femFldC[3]+4.0*femFldC2R2+(2.0*femFldC[1]+4.0*femFldC[0])*femFldC[2]+4.0*femFldC1R2+4.0*femFldC[0]*femFldC[1]+4.0*femFldC0R2)*volFac; 
}

void MGpoissonFEMM0norm2xSer_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 
  double *femFldUy = femFld[2]; 
  double *femFldUxUy = femFld[3]; 


  normOut[0] += (femFldUy[0]+femFldUxUy[0]+femFldUx[0]+femFldC[0])*volFac; 
}

void MGpoissonFEMM0norm2xSer_LxNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 
  double *femFldUy = femFld[2]; 
  double *femFldUxUy = femFld[3]; 


  normOut[0] += (femFldUy[0]+femFldUxUy[0]+femFldUx[0]+femFldC[0])*volFac; 
}

void MGpoissonFEMM0norm2xSer_UxNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUy = femFld[2]; 


  normOut[0] += (femFldUy[1]+femFldC[1]+femFldUy[0]+femFldC[0])*volFac; 
}

void MGpoissonFEMM0norm2xSer_LyNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 
  double *femFldUy = femFld[2]; 
  double *femFldUxUy = femFld[3]; 


  normOut[0] += (femFldUy[0]+femFldUxUy[0]+femFldUx[0]+femFldC[0])*volFac; 
}

void MGpoissonFEMM0norm2xSer_UyNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 


  normOut[0] += (femFldUx[1]+femFldC[1]+femFldUx[0]+femFldC[0])*volFac; 
}

void MGpoissonFEMM0norm2xSer_LxNonPeriodicUyNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUx = femFld[1]; 


  normOut[0] += (femFldUx[1]+femFldC[1]+femFldUx[0]+femFldC[0])*volFac; 
}

void MGpoissonFEMM0norm2xSer_UxNonPeriodicLyNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 
  double *femFldUy = femFld[2]; 


  normOut[0] += (femFldUy[1]+femFldC[1]+femFldUy[0]+femFldC[0])*volFac; 
}

void MGpoissonFEMM0norm2xSer_UxNonPeriodicUyNonPeriodic_P1(const double *dxC, double **femFld, double *normOut) 
{ 
  // femFld:  FEM field in neighboring cells.
  // normOut: norm.

  double volFac = 0.25*dxC[0]*dxC[1]; 

  double *femFldC = femFld[0]; 


  normOut[0] += (femFldC[3]+femFldC[2]+femFldC[1]+femFldC[0])*volFac; 
}

void MGpoissonFEMaccuConst2xSer_P1(const double constIn, double *femFld) 
{ 
  // constIn: constant to accumulate.
  // femFld:  FEM field to accumulate.

  femFld[0] += constIn; 
}

void MGpoissonFEMaccuConst2xSer_LxNonPeriodic_P1(const double constIn, double *femFld) 
{ 
  // constIn: constant to accumulate.
  // femFld:  FEM field to accumulate.

  femFld[0] += constIn; 
}

void MGpoissonFEMaccuConst2xSer_UxNonPeriodic_P1(const double constIn, double *femFld) 
{ 
  // constIn: constant to accumulate.
  // femFld:  FEM field to accumulate.

  femFld[0] += constIn; 
  femFld[1] += constIn; 
}

void MGpoissonFEMaccuConst2xSer_LyNonPeriodic_P1(const double constIn, double *femFld) 
{ 
  // constIn: constant to accumulate.
  // femFld:  FEM field to accumulate.

  femFld[0] += constIn; 
}

void MGpoissonFEMaccuConst2xSer_UyNonPeriodic_P1(const double constIn, double *femFld) 
{ 
  // constIn: constant to accumulate.
  // femFld:  FEM field to accumulate.

  femFld[0] += constIn; 
  femFld[1] += constIn; 
}

void MGpoissonFEMaccuConst2xSer_LxNonPeriodicUyNonPeriodic_P1(const double constIn, double *femFld) 
{ 
  // constIn: constant to accumulate.
  // femFld:  FEM field to accumulate.

  femFld[0] += constIn; 
  femFld[1] += constIn; 
}

void MGpoissonFEMaccuConst2xSer_UxNonPeriodicLyNonPeriodic_P1(const double constIn, double *femFld) 
{ 
  // constIn: constant to accumulate.
  // femFld:  FEM field to accumulate.

  femFld[0] += constIn; 
  femFld[1] += constIn; 
}

void MGpoissonFEMaccuConst2xSer_UxNonPeriodicUyNonPeriodic_P1(const double constIn, double *femFld) 
{ 
  // constIn: constant to accumulate.
  // femFld:  FEM field to accumulate.

  femFld[0] += constIn; 
  femFld[1] += constIn; 
  femFld[2] += constIn; 
  femFld[3] += constIn; 
}

