-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute RHS or forward Euler update for Vlasov equation
-- with Discontinuous Galerkin scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Base = require "Updater.Base"
local Lin = require "Lib.Linalg"
local VlasovModDecl = require "Updater.vlasovData.VlasovModDecl"

-- Vlasov DG solver updater object
local VlasovDisCont = {}

function VlasovDisCont:new(tbl)
   local self = setmetatable({}, VlasovDisCont)
   Base.setup(self, tbl) -- setup base object

   -- read data from input file
   self._onGrid = assert(
      tbl.onGrid, "Updater.VlasovDisCont: Must provide grid object using 'onGrid'")

   self._phaseBasis = assert(
      tbl.phaseBasis, "Updater.VlasovDisCont: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(
      tbl.confBasis, "Updater.VlasovDisCont: Must specify configuration-space basis functions to use using 'confBasis'")
   self._charge = assert(tbl.charge, "Updater.VlasovDisCont: must specify charge using 'charge' ")
   self._mass = assert(tbl.mass, "Updater.VlasovDisCont: must specify mass using 'mass' ")

   assert(self._onGrid:ndim() == self._phaseBasis:ndim(), "Dimensions of basis and grid must match")

   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   -- CFL number
   self._cfl = assert(tbl.cfl, "Updater.VlasovDisCont: Must specify CFL number using 'cfl'")
   self._cflm = tbl.cflm and tbl.cflm or 1.1*self._cfl -- no larger than this

   -- functions to perform streaming updates
   self._volStreamUpdate = VlasovModDecl.selectVolStream(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   self._surfStreamUpdate = VlasovModDecl.selectSurfStream(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(VlasovDisCont, { __call = function (self, o) return self.new(self, o) end })

-- advance method
local function advance(self, tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local fIn = assert(inFld[1], "VlasovDisCont.advance: Must specify an input dist-function")
   local fOut = assert(outFld[1], "VlasovDisCont.advance: Must specify an output dist-function")
   
   local localRange = fOut:localRange()
   local fInIdxr, fOutIdxr = fIn:genIndexer(), fOut:genIndexer()
 
   -- pointers for (re)use in update
   local fInPtr, fOutPtr = fIn:get(1), fOut:get(1)
   local fInL, fInR = fIn:get(1), fIn:get(1)
   local fOutL, fOutR = fOut:get(1), fOut:get(1)

   local pdim, cdim, vdim = self._pdim, self._cdim, self._vdim
   local dx, xc = Lin.Vec(pdim), Lin.Vec(pdim)
   local cfl, cflm = self._cfl, self._cflm
   local cfla = 0.0 -- actual CFL number used

   -- for ease in indexing coordinate arrays
   local X, Y, Z, VX, VY, VZ = 1, 1, 1, 1, 1, 1
   if cdim == 1 then
      X, VX, VY, VZ = 1, 2, 3, 4
   elseif cdim == 2 then
      X, Y, VX, VY, VZ = 1, 2, 3, 4, 5
   elseif cdim == 2 then
      X, Y, Z, VX, VY, VZ = 1, 2, 3, 4, 5, 6
   end
   
   fOut:clear(0.0) -- compute increments
   -- accumulate contributions from volume integrals
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)
      for d = 1, pdim do dx[d] = grid:dx(d) end
      grid:cellCenter(xc)

      fIn:fill(fInIdxr(idx), fInPtr)
      fOut:fill(fOutIdxr(idx), fOutPtr)
      self._volStreamUpdate(xc:data(), dx:data(), fInPtr:data(), fOutPtr:data())
   end

   -- accumulate contributions from surface integrals in streaming direction
   for dir = 1, cdim do
      -- lower/upper bounds in direction 'dir': these are edge indices (one more edge than cell)
      local dirLoIdx, dirUpIdx = localRange:lower(dir), localRange:upper(dir)+1
      local perpRange = localRange:shorten(dir) -- range orthogonal to 'dir'

      -- outer loop is over directions orthogonal to 'dir' and inner
      -- loop is over 1D slice in `dir`.
      for idx in perpRange:colMajorIter() do
	 grid:setIndex(idx)
	 for d = 1, pdim do dx[d] = grid:dx(d) end
	 grid:cellCenter(xc)

	 -- compute local CFL number
	 local vel = math.abs(xc[dir+cdim]) + 0.5*dx[dir+cdim] -- ptcl velocity at cell edge
	 cfla = math.max(cfla, vel*dt/dx[dir])
	 
	 local idxp, idxm = idx:copy(), idx:copy()
   	 for i = dirLoIdx, dirUpIdx do -- this loop is over edges
	    idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'

	    fIn:fill(fInIdxr(idxm), fInL); fIn:fill(fInIdxr(idxp), fInR)
	    fOut:fill(fOutIdxr(idxm), fOutL); fOut:fill(fOutIdxr(idxp), fOutR)
	    self._surfStreamUpdate[dir](xc:data(), dx:data(), fInL:data(), fInR:data(), fOutL:data(), fOutR:data())
	 end
      end
   end

   -- return if time-step was too large
   if cfla > cflm then return false, dt*cfl/cfla end

   -- NEED TO ADD fIn after scaling fOut by dt

   return true, dt*cfl/cfla
end

-- Methods in updater
VlasovDisCont.__index = { advance = Base.advanceFuncWrap(advance) }

return VlasovDisCont
