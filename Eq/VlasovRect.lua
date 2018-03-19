-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov equation on a rectangular mesh.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local BoundaryCondition = require "Updater.BoundaryCondition"
local Proto = require "Lib.Proto"
local ffi = require "ffi"

-- Vlasov equation on a rectangular mesh
local VlasovRect = Proto()

-- ctor
function VlasovRect:init(tbl)

end

-- Methods
function VlasovRect:numEquations() return 1 end
function VlasovRect:numWaves() return 1 end
function VlasovRect:isPositive(q)
   if q[0] > 0.0 return true else return false
end

-- flux in direction dir
function VlasovRect:flux(dir, qIn, fOut)
   assert(false, "VlasovRect:flux: NYI!")
end

-- Riemann problem for Vlasov equation
function VlasovRect:rp(dir, delta, ql, qr, waves, s)
   assert(false, "VlasovRect:rp: NYI!")
end

-- Compute q-fluctuations
function VlasovRect:qFluctuations(dir, ql, qr, waves, s, amdq, apdq)
   assert(false, "VlasovRect:qFluctuations: NYI!")
end

-- Maximum wave speed
function VlasovRect:maxSpeed(dir, w, dx, q)
   return 0.0
end       
-- Volume integral term for use in DG scheme
function VlasovRect:volTerm(w, dx, idx, q, out)
   -- TODO
end
-- Surface integral term for use in DG scheme
function VlasovRect:surfTerm(dir, w, dx, maxs, idxl, idxr, ql, qr, outl, outr)
   -- TODO
end

function VlasovRect:setAuxFields(auxFields)
end

return VlasovRect
