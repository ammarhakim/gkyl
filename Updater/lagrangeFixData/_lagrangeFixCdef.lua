local ffi  = require "ffi"
ffi.cdef [[

void lagrangeFixSer1x1v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixSer1x2v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x1v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x2v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

]]