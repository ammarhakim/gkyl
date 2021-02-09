-- Gkyl ------------------------------------------------------------------------
--
-- Updater to perform ten-moment gradient-based closure based on 
-- \partial_k q_{ijk} = \partial_k (-v_t/|k_s| \chi \partial[_i T_{jk}]),
-- where [_i A_{jk}] is the symmetrization operator
-- NOTE: This updater computes only the RHS, i.e., \partial_k q_{ijk} 
-- with a second order central difference
-- ADDITIONAL NOTE: gradient-based closure only affects pressure tensor, NOT 
-- full stress tensor -> updater expects to receive primitive variables (just P_{ij})
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

-- system libraries
local ffi = require "ffi"

-- Define C types for storing private data for use in updater
ffi.cdef [[
void gkylTenMomentHeatFlux(const double alpha, const double* dT1, const double* dT2, const double* dT3, const double* f, double* q);
void gkylTenMomentAccumulateGradClosure(const double* divQ1, const double* divQ2, const double* divQ3, double* f);
void gkylTenMomentGradT(const int dir, const double* dxv, const double* fL, const double* fR, double* dT);
void gkylTenMomentDivQX(const double* dxv, const double* qL, const double* qR, double* divQ);
void gkylTenMomentDivQY(const double* dxv, const double* qL, const double* qR, double* divQ);
void gkylTenMomentDivQZ(const double* dxv, const double* qL, const double* qR, double* divQ);
]]

-- Function to generate symmetrized heat flux tensor from gradient of temperature tensor
local function heatFlux(alpha, dT1, dT2, dT3, fPtr, qPtr)
   ffi.C.gkylTenMomentHeatFlux(alpha, dT1, dT2, dT3, fPtr, qPtr)
end

-- Function to accumulate the resulting divergence of the heat flux onto the stress tensor
local function accumulateGradClosure(divQ1, divQ2, divQ3, fPtr)
   ffi.C.gkylTenMomentAccumulateGradClosure(divQ1, divQ2, divQ3, fPtr)
end

-- Function to compute the gradient of the temperature tensor in a direction 'dir'
local function gradT(dir, dx, fLPtr, fRPtr, dTPtr)
   ffi.C.gkylTenMomentGradT(dir, dx, fLPtr, fRPtr, dTPtr)
end

-- Functions to compute the divergence of the heat flux tensor (in a given direction)
local function divQX(dx, qLPtr, qRPtr, divQPtr)
   ffi.C.gkylTenMomentDivQX(dx, qLPtr, qRPtr, divQPtr)
end
local function divQY(dx, qLPtr, qRPtr, divQPtr)
   ffi.C.gkylTenMomentDivQY(dx, qLPtr, qRPtr, divQPtr)
end
local function divQZ(dx, qLPtr, qRPtr, divQPtr)
   ffi.C.gkylTenMomentDivQZ(dx, qLPtr, qRPtr, divQPtr)
end
-- Ten-moment source updater object
local TenMomentGrad = Proto(UpdaterBase)

function TenMomentGrad:init(tbl)
   TenMomentGrad.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid, "Updater.TenMomentGrad: Must provide grid object using 'onGrid'")
   self._ndim = self._onGrid:ndim()
   -- Setup grad(T) and div(Q) functions,
   -- NOTE: div(Q) will be different in each direction and depend on dimensionality
   self._gradT = gradT
   if self._ndim == 1 then
      self._divQX = divQX
   elseif self._ndim == 2 then
      self._divQX = divQX
      self._divQY = divQY
   else 
      self._divQX = divQX
      self._divQY = divQY
      self._divQZ = divQZ
   end
   -- Setup heat flux and accumulate methods, which are dimensionally independent
   self._heatFlux = heatFlux
   self._accumulateGradClosure = accumulateGradClosure

   self._alpha = assert(tbl.alpha, "Updater.TenMomentGrad: Must provide coefficient alpha for strength of thermal conductivity")

   -- Six small helper arrays to hold the gradients of the temperature tensor and divergence of the heat flux tensor in each direction
   self._dT1, self._dT2, self._dT3 =  Lin.Vec(6), Lin.Vec(6), Lin.Vec(6)
   self._divQ1, self._divQ2, self._divQ3 = Lin.Vec(6), Lin.Vec(6), Lin.Vec(6)
end

-- advance method
function TenMomentGrad:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local ndim = self._ndim

   -- additional coefficient for the strength of "thermal conductivity"
   -- alpha has units of length to make units of heat flux work 
   local alpha = self._alpha

   -- NOTE: the inFld is a temporary array used to store the heat flux
   -- the outFld is the incremented source term
   -- Because this update is in-place, the outFld also contains the necessary input information (e.g. pressure at previous time-step)
   local f =  assert(inFld[1], "TenMomentGrad.advance: Must specify an input field")
   local q = assert(inFld[2], "TenMomentGrad.advance: Must specify an input field")
   local diff = assert(outFld[1], "TenMomentGrad.advance: Must specify an output field")

   local localRange = f:localRange()
   local qIdxr, fIdxr, diffIdxr = q:genIndexer(), f:genIndexer(), diff:genIndexer() -- indexer functions into fields

   -- to store grid info
   local dx = Lin.Vec(ndim) -- cell shape on left/center/right
   local idxl, idxc, idxr = Lin.IntVec(ndim), Lin.IntVec(ndim), Lin.IntVec(ndim) -- index on left/center/right

   -- pointers for (re)use in update
   local qL, qC, qR = q:get(1), q:get(1), q:get(1)
   local fL, fC, fR = f:get(1), f:get(1), f:get(1)
   local diffC = diff:get(1)
   -- clear the inFld that will be used to store the heat fluxes before beginning loops
   q:clear(0.0)
   -- zero out helper arrays before start of loop
   for i = 1,6 do
      self._dT1[i], self._dT2[i], self._dT3[i] = 0.0, 0.0, 0.0
      self._divQ1[i], self._divQ2[i], self._divQ3[i] = 0.0, 0.0, 0.0
   end
   -- first loop over the whole domain and compute the heat flux, q_{ijk} = (-v_t/|k_s| \chi \partial[_i T_{jk}])
   for idx in localRange:rowMajorIter() do
      idx:copyInto(idxl)
      idx:copyInto(idxc)
      idx:copyInto(idxr)
      -- get grid information and set pointers for current cell
      grid:setIndex(idxc)
      grid:getDx(dx)

      q:fill(qIdxr(idxc), qC)
      f:fill(fIdxr(idxc), fC)

      -- loop over directions for the update 
      for dir = 1, ndim do
         -- cell left/right of cell 'i' in direction dir
         idxl[dir] = idxl[dir] - 1
         idxr[dir] = idxr[dir] + 1

         grid:setIndex(idxl)
         grid:setIndex(idxr)

         f:fill(fIdxr(idxl), fL)
         f:fill(fIdxr(idxr), fR)
         -- Store the gradient of the temperature tensor in each direction
         -- NOTE: Functions expect to be passed PRIMITIVE variables
         if dir == 1 then
            self._gradT(dir-1, dx:data(), fL:data(), fR:data(), self._dT1:data()) 
         elseif dir == 2 then
            self._gradT(dir-1, dx:data(), fL:data(), fR:data(), self._dT2:data()) 
         else
            self._gradT(dir-1, dx:data(), fL:data(), fR:data(), self._dT3:data()) 
         end 
      end
      -- Use the pointer to the inFld to store the evaluation of the heat flux
      -- NOTE: includes the additional corresponding factors of density and vt so the units of heat-flux are correct
      -- These factors are computed from the outFld, which AT THIS POINT still contains the previous time-step fluid data
      self._heatFlux(alpha, self._dT1:data(), self._dT2:data(), self._dT3:data(), fC:data(), qC:data())
   end

   -- loop over the grid again to evaluate div(Q) = \partial_k q_{ijk} and increment the result to advance the solution in time
   -- NOTE: In this second loop, outFld will be *incremented* to the new time-step
   for idx in localRange:rowMajorIter() do
      idx:copyInto(idxl)
      idx:copyInto(idxc)
      idx:copyInto(idxr)
      -- get grid information and set pointers for current cell
      grid:setIndex(idxc)
      grid:getDx(dx)
      q:fill(qIdxr(idxc), qC)
      diff:fill(diffIdxr(idxc), diffC)

      -- loop over directions for the update 
      for dir = 1, ndim do
         -- cell left/right of cell 'i' in direction dir
         idxl[dir] = idxl[dir] - 1
         idxr[dir] = idxr[dir] + 1

         grid:setIndex(idxl)
         grid:setIndex(idxr)

         q:fill(qIdxr(idxl), qL)
         q:fill(qIdxr(idxr), qR)
         -- Store the divergence of the heat flux tensor in each direction
         if dir == 1 then
            self._divQX(dx:data(), qL:data(), qR:data(), self._divQ1:data()) 
         elseif dir == 2 then
            self._divQY(dx:data(), qL:data(), qR:data(), self._divQ2:data())  
         else
            self._divQZ(dx:data(), qL:data(), qR:data(), self._divQ3:data())  
         end 
      end
      -- Accumulate the components of the divergence
      self._accumulateGradClosure(self._divQ1:data(), self._divQ2:data(), self._divQ3:data(), diffC:data())
   end
   return true, GKYL_MAX_DOUBLE
end

return TenMomentGrad