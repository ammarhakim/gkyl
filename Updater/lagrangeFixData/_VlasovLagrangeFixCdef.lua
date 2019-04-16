local ffi  = require "ffi"
ffi.cdef [[

void VlasovLagrangeFixSer1x1v1p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixSer1x1v2p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixSer1x2v1p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixSer1x2v2p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixSer2x2v1p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixSer2x2v2p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixMax1x1v1p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixMax1x1v2p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixMax1x2v1p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixMax1x2v2p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixMax2x2v1p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

void VlasovLagrangeFixMax2x2v2p(double *dm0, double *dm1, double *dm2, double *lo, double *L, double *Nv, double *vc, double *f);

]]