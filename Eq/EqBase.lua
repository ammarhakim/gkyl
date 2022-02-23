-- Gkyl ------------------------------------------------------------------------
--
-- A base object for equation object
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local Proto = require "Lib.Proto"
local ffi = require "ffi"
require "Lib.ZeroUtil"

ffi.cdef [[ 
  typedef struct GkylEquation_t GkylEquation_t; 

struct gkyl_dg_eqn;

// Function pointer type for volume kernel
typedef double (*vol_termf_t)(const struct gkyl_dg_eqn *eqn,
  const double*  xc, const double*  dx, const int*  idx,
  const double* qIn, double* qRhsOut);

// Function pointer type for surface kernel
typedef void (*surf_termf_t)(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR,
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* qRhsOut);

// Function pointer type for surface kernel
typedef void (*boundary_surf_termf_t)(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* qRhsOut);

struct gkyl_dg_eqn {
  int num_equations; // number of equations in system
  vol_termf_t vol_term; // volume term kernel
  surf_termf_t surf_term; // surface term kernel
  boundary_surf_termf_t boundary_surf_term; // boundary surface term kernel

  uint32_t flags;
  struct gkyl_ref_count ref_count; // reference count
  struct  gkyl_dg_eqn *on_dev; // pointer to itself or device data
};
]]

local EqBase = Proto()

-- for FV scheme
function EqBase:numEquations() return 1 end
function EqBase:numWaves() return 1 end
function EqBase:flux(dir, qIn, fOut) end
function EqBase:speeds(dir, qIn, sOut) end
function EqBase:maxAbsSpeed(dir, qIn) return 0.0 end
function EqBase:isPositive(q) return true end
function EqBase:rp(dir, delta, ql, qr, waves, s) end
function EqBase:qFluctuations(dir, ql, qr, waves, s, amdq, apdq) end
-- for DG scheme
function EqBase:setAuxFields(auxFields) end
function EqBase:volTerm(w, dx, idx, q, out) end
function EqBase:surfTerm(dir, cfll, cflr, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr) end
function EqBase:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr) end

return EqBase
