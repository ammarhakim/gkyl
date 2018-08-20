-- Gkyl ------------------------------------------------------------------------
--
-- A base object for equation object
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local Proto = require "Lib.Proto"

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
function EqBase:surfTerm(dir, cfl, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr) end
function EqBase:boundarySurfTerm(dir, wl, wr, dxl, dxr, maxs, idxl, idxr, ql, qr, outl, outr) end

return EqBase
