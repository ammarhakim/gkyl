local ffi  = require "ffi"
ffi.cdef [[

void GkLagrangeFixSer1x2v1p(double *dm0, double *dm1, double *dm2, double *B, double mass, double *L, double *Nv, double *vc, double *f);

]]