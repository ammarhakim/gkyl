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
local ffi = require "ffi"
local Time = require "Lib.Time"
local VlasovModDecl = require "Updater.vlasovData.VlasovModDecl"
local xsys = require "xsys"

-- for incrementing in updater
ffi.cdef [[ void vlasovIncr(unsigned n, const double *aIn, double a, double *aOut); ]]

-- compute accelOut = q/m*forceIn
local function calcAccelFromForce(qbym, forceIn, accelOut)
   ffi.C.vlasovIncr(#accelOut, forceIn:data(), qbym, accelOut:data())
end

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

   local charge = assert(tbl.charge, "Updater.VlasovDisCont: must specify charge using 'charge' ")
   local mass = assert(tbl.mass, "Updater.VlasovDisCont: must specify mass using 'mass' ")

   self._qbym = charge/mass -- only q/m ratio is ever needed
   
   assert(self._onGrid:ndim() == self._phaseBasis:ndim(), "Dimensions of basis and grid must match")

   self._pdim = self._phaseBasis:ndim()
   self._cdim = self._confBasis:ndim()
   self._vdim = self._pdim-self._cdim

   -- by default, compute forward Euler: if onlyIncrement is true,
   -- then only increments are computed. NOTE: The increments are NOT
   -- multiplied by dt.
   self._onlyIncrement = xsys.pickBool(tbl.onlyIncrement, false)

   -- CFL number
   self._cfl = assert(tbl.cfl, "Updater.VlasovDisCont: Must specify CFL number using 'cfl'")
   self._cflm = tbl.cflm and tbl.cflm or 1.1*self._cfl -- no larger than this

   self._tmStreamUpdate = 0.0 -- time in streaming terms
   self._tmForceUpdate = 0.0 -- time in field terms
   self._tmIncrement = 0.0 -- time 

   -- functions to perform streaming updates
   self._volStreamUpdate = VlasovModDecl.selectVolStream(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   self._surfStreamUpdate = VlasovModDecl.selectSurfStream(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())

   -- check if we have a electric and magnetic field
   local hasElcField = xsys.pickBool(tbl.hasElectricField, true)
   local hasMagField = xsys.pickBool(tbl.hasMagneticField, true)

   self._hasForceTerm = false -- flag to indicate if we have any force terms at all
   if hasElcField or hasMagField then
      self._hasForceTerm = true
   end
   
   -- functions to perform force updates from electric field
   self._volForceUpdate = VlasovModDecl.selectVolElc(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())
   self._surfForceUpdate = VlasovModDecl.selectSurfElc(self._phaseBasis:id(), self._cdim, self._vdim, self._phaseBasis:polyOrder())

   if hasMagField then -- more complicated functions if there is a magnetic field
      assert(false, "VlasovDisCont: hasMagField NYI!")
   end
   
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
   elseif cdim == 3 then
      X, Y, Z, VX, VY, VZ = 1, 2, 3, 4, 5, 6
   end
   
   fOut:clear(0.0) -- compute increments

   local tmStreamStart = Time.clock()
   -- accumulate contributions from streaming direction
   for dir = 1, cdim do
      -- This flag is needed as the volume integral already contains
      -- contributions from all streaming directions. Hence, we can
      -- only accumulate the volume contribution once, skipping it for
      -- other directions
      local firstDir = true
      
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

	    if firstDir then
	       -- accumulate contribution from volume streaming terms
	       -- only if first time around
	       self._volStreamUpdate(xc:data(), dx:data(), fInR:data(), fOutR:data())
	    end
	    -- accumulate contribution from surface streaming terms
	    self._surfStreamUpdate[dir](xc:data(), dx:data(), fInL:data(), fInR:data(), fOutL:data(), fOutR:data())
	 end
      end
      firstDir = false
   end
   self._tmStreamUpdate = self._tmStreamUpdate + Time.clock()-tmStreamStart

   if self._hasForceTerm then
      local emIn = assert(inFld[2], "VlasovDisCont.advance: Must specify an input EM field")
      local emPtr = emIn:get(1)
      local emIdxr = emIn:genIndexer()
      local emAccel = Lin.Vec(emIn:numComponents())
      
      local tmForceStart = Time.clock()
      -- accumulate contributions from force directions
      for dir = cdim+1, pdim do
	 -- See comment for streaming terms above to understand role
	 -- of this flag
	 local firstDir = true
	 
	 -- lower/upper bounds in direction 'dir': these are edge indices
	 local dirLoIdx, dirUpIdx = localRange:lower(dir), localRange:upper(dir)+1
	 local perpRange = localRange:shorten(dir) -- range orthogonal to 'dir'

	 -- outer loop is over directions orthogonal to 'dir' and inner
	 -- loop is over 1D slice in `dir`.
	 for idx in perpRange:colMajorIter() do
	    grid:setIndex(idx)
	    for d = 1, cdim do dx[d] = grid:dx(d) end
	    grid:cellCenter(xc)

	    emIn:fill(emIdxr(idx), emPtr) -- get pointer to EM field
	    calcAccelFromForce(self._qbym, emPtr, emAccel) -- compute acceleration

	    -- NEED TO COMPUTE CFL!!!
	    
	    local idxp, idxm = idx:copy(), idx:copy()
	    for i = dirLoIdx, dirUpIdx do -- this loop is over edges
	       idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'

	       fIn:fill(fInIdxr(idxm), fInL); fIn:fill(fInIdxr(idxp), fInR)
	       fOut:fill(fOutIdxr(idxm), fOutL); fOut:fill(fOutIdxr(idxp), fOutR)

	       if firstDir then
		  -- accumulate contribution from volume force terms
		  self._volForceUpdate(xc:data(), dx:data(), emAccel:data(), fInR:data(), fOutR:data())
	       end
	       if i > dirLoIdx and i < dirUpIdx then
		  -- accumulate contribution from surface force terms
		  -- (skiping contribution from boundary faces as
		  -- assuming zero particle flux)
		  self._surfForceUpdate[dir-cdim](
		     xc:data(), dx:data(), emAccel:data(), fInL:data(), fInR:data(), fOutL:data(), fOutR:data())
	       end
	    end
	 end
	 firstDir = false
      end
      self._tmForceUpdate = self._tmForceUpdate + Time.clock()-tmForceStart
   end

   -- return failure if time-step was too large
   if cfla > cflm then return false, dt*cfl/cfla end

   local tmIncrementStart = Time.clock()
   -- accumulate full solution if not computing increments
   if not self._onlyIncrement then
      fOut:scale(dt); fOut:accumulate(1.0, fIn) -- fOut = fIn + dt*fOut
   end
   self._tmIncrement = self._tmIncrement + Time.clock()-tmIncrementStart

   return true, dt*cfl/cfla
end

-- Methods in updater
VlasovDisCont.__index = {
   advance = Base.advanceFuncWrap(advance),
   streamTime = function(self)
      return self._tmStreamUpdate
   end,
   forceTime = function(self)
      return self._tmForceUpdate
   end,
   incrementTime = function(self)
      return self._tmIncrement
   end
}

return VlasovDisCont
