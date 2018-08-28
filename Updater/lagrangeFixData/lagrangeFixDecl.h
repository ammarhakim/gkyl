#ifndef LAGRANGE_FIX_H 
#define LAGRANGE_FIX_H 

extern "C" { 

double gkyl_ipow(double base, int exp);

void lagrangeFixSer1x1v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer1x1v2p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer1x1v3p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer1x2v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer1x2v2p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer1x2v3p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer1x3v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer1x3v2p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer1x3v3p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer2x2v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer2x2v2p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer2x2v3p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer2x3v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer2x3v2p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer2x3v3p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x1v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x1v2p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x1v3p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x2v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x2v2p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x2v3p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x3v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x3v2p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x3v3p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax2x2v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax2x2v2p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax2x2v3p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax2x3v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax2x3v2p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax2x3v3p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

} 
#endif