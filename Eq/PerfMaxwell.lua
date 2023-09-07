-- Gkyl ------------------------------------------------------------------------
--
-- Perfectly Hyperbolic Maxwell equations. See
-- http://ammar-hakim.org/sj/maxwell-eigensystem.html and references
-- on that page for details.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local EqBase = require "Eq.EqBase"
local Proto = require "Lib.Proto"
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"

ffi.cdef [[ 
/**
 * Create a new Maxwell equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param lightSpeed Speed of light
 * @param elcErrorSpeedFactor Factor multiplying lightSpeed for div E correction
 * @param mgnErrorSpeedFactor Factor multiplying lightSpeed for div B correction
 * @return Pointer to Maxwell equation object
 */
struct gkyl_dg_eqn* gkyl_dg_maxwell_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor, bool use_gpu);
]]

-- Perfectly hyperbolic Maxwell equations
local PerfMaxwell = Proto(EqBase)

-- ctor
function PerfMaxwell:init(tbl)
   self._c = assert(tbl.lightSpeed, "Eq.PerfMaxwell: Must specify light speed (lightSpeed)")
   self._ce = tbl.elcErrorSpeedFactor and tbl.elcErrorSpeedFactor or 0.0
   self._cb = tbl.mgnErrorSpeedFactor and tbl.mgnErrorSpeedFactor or 0.0

   self._basis = assert(tbl.basis, "Eq.PerfMaxwell: Must specify basis functions to use using 'basis'")

   -- store pointers to C equation object
   self._zero = ffi.gc(ffiC.gkyl_dg_maxwell_new(self._basis._zero, self._c, self._ce, self._cb, GKYL_USE_GPU or 0),
                       ffiC.gkyl_dg_eqn_release)
end

return PerfMaxwell
