-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the Vlasov solver
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"


ffi.cdef [[
void VlasovVolStream1x1vMaxP1(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x1vMaxP2(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x1vMaxP3(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x1vMaxP4(const double *w, const double *dxv, const double *f, double *out);

void VlasovVolStream1x2vMaxP1(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x2vMaxP2(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x2vMaxP3(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x2vMaxP4(const double *w, const double *dxv, const double *f, double *out);

void VlasovVolStream1x3vMaxP1(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x3vMaxP2(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x3vMaxP3(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x3vMaxP4(const double *w, const double *dxv, const double *f, double *out);

void VlasovVolStream2x2vMaxP1(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x2vMaxP2(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x2vMaxP3(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x2vMaxP4(const double *w, const double *dxv, const double *f, double *out);

void VlasovVolStream2x3vMaxP1(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x3vMaxP2(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x3vMaxP3(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x3vMaxP4(const double *w, const double *dxv, const double *f, double *out);


 
void VlasovVolStream1x1vSerP1(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x1vSerP2(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x1vSerP3(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x1vSerP4(const double *w, const double *dxv, const double *f, double *out);

void VlasovVolStream1x2vSerP1(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x2vSerP2(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x2vSerP3(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x2vSerP4(const double *w, const double *dxv, const double *f, double *out);

void VlasovVolStream1x3vSerP1(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x3vSerP2(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x3vSerP3(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream1x3vSerP4(const double *w, const double *dxv, const double *f, double *out);

void VlasovVolStream2x2vSerP1(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x2vSerP2(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x2vSerP3(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x2vSerP4(const double *w, const double *dxv, const double *f, double *out);

void VlasovVolStream2x3vSerP1(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x3vSerP2(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x3vSerP3(const double *w, const double *dxv, const double *f, double *out);
void VlasovVolStream2x3vSerP4(const double *w, const double *dxv, const double *f, double *out);
]]
