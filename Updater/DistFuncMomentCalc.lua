-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments of distribution function on a
-- rectangular (but potentially non-uniform) grid.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Base = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
      "new, copy, fill, sizeof, typeof, metatype")

-- makes dispatching into appropriate module easier
local maxModNames = {
   {"MaxOrderMom1x1v", "MaxOrderMom1x2v", "MaxOrderMom1x3v"},
   {nil, "MaxOrderMom2x2v", "MaxOrderMom2x3v"},
   {nil, nil, "MaxOrderMom3x3v"},
}
local serModNames = {
   {"SerendipMom1x1v", "SerendipMom1x2v", "SerendipMom1x3v"},
   {nil, "SerendipMom2x2v", "SerendipMom2x3v"},
   {nil, nil, "SerendipMom3x3v"},
}

local maxPrModNames = {
   {"MaxOrderPressureMom1x1v", "MaxOrderPressureMom1x2v", "MaxOrderPressureMom1x3v"},
   {nil, "MaxOrderPressureMom2x2v", "MaxOrderPressureMom2x3v"},
   {nil, nil, "MaxOrderPressureMom3x3v"},
}
local serPrModNames = {
   {"SerendipPressureMom1x1v", "SerendipPressureMom1x2v", "SerendipPressureMom1x3v"},
   {nil, "SerendipPressureMom2x2v", "SerendipPressureMom2x3v"},
   {nil, nil, "SerendipPressureMom3x3v"},
}

-- select appropriate module
local function selectMod(cDim, vDim, bId, polyOrder)
   local nm = nil
   if bId == "serendipity" then
      nm = serModNames[cDim][vDim]
   elseif bId == "maximal-order" then
      nm = maxModNames[cDim][vDim]
   end
   return require("Updater._projectOnBasisData." .. nm)
end
-- select appropriate module for energy and pressure-tensor (these are
-- stored in different files)
local function selectPrMod(cDim, vDim, bId, polyOrder)
   local nm = nil
   if bId == "serendipity" then
      nm = serPrModNames[cDim][vDim]
   elseif bId == "maximal-order" then
      nm = maxPrModNames[cDim][vDim]
   end
   return require("Updater._projectOnBasisData." .. nm)
end

-- pick function to compute various moments
local function pickNumDensityFunc(cDim, vDim, bId, polyOrder)
   local m = selectMod(cDim, vDim, bId, polyOrder)
   return m.numDensity[polyOrder]
end
local function pickMomentumFunc(cDim, vDim, bId, polyOrder)
   local m = selectMod(cDim, vDim, bId, polyOrder)
   return m.momentum[polyOrder]
end
local function pickPressureFunc(cDim, vDim, bId, polyOrder)
   local m = selectPrMod(cDim, vDim, bId, polyOrder)
   return m.pressureTensor[polyOrder]
end
local function pickEnergyFunc(cDim, vDim, bId, polyOrder)
   local m = selectPrMod(cDim, vDim, bId, polyOrder)
   return m.energy[polyOrder]
end

-- set field values to zero
local function scaleField(field, scale)
   local localRange, indexer = field:localRange(), field:genIndexer()
   for idx in localRange:colMajorIter() do
      local fitr = field:get(indexer(idx))
      for k = 1, field:numComponents() do
	 fitr[k] = fitr[k]*scale
      end
   end
end

-- templates for various functions
local calcGridQuantTempl = xsys.template([[
return function (grid, xcp, dv, w)
  local vol = 1.0
|for d = 1, VDIM do
   dv[${d}] = grid:dx(${CDIM}+${d})
   w[${d}] = xcp[${CDIM}+${d}]
   vol = vol*dv[${d}]
|end
   return vol
end
]])

-- Moments updater object
local DistFuncMomentCalc = {}

function DistFuncMomentCalc:new(tbl)
   local self = setmetatable({}, DistFuncMomentCalc)
   Base.setup(self, tbl) -- setup base object

   self._onGrid = assert(
      tbl.onGrid, "Updater.ProjectOnBasis: Must provide grid object using 'onGrid'")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.ProjectOnBasis: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis, "Updater.ProjectOnBasis: Must provide configuration-space basis object using 'confBasis'")

   -- dimension of spaces
   self._pDim = phaseBasis:ndim() 
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- ensure sanity
   assert(phaseBasis:polyOrder() == confBasis:polyOrder(),
	  "Polynomial orders of phase-space and config-space basis must match")
   assert(phaseBasis:id() == confBasis:id(),
	  "Type of phase-space and config-space basis must match")
   
   local mom = assert(
      tbl.moment, "Updater.ProjectOnBasis: Must provide moment to compute using 'moment'.")

   local id, polyOrder = phaseBasis:id(), phaseBasis:polyOrder()
   
   self._momCalcFun = nil -- function to compute requested moments
   if mom == "number-density" then
      self._momCalcFun = pickNumDensityFunc(self._cDim, self._vDim, id, polyOrder)
   elseif mom == "momentum" then
      self._momCalcFun = pickMomentumFunc(self._cDim, self._vDim, id, polyOrder)
   elseif mom == "energy" then
      self._momCalcFun = pickEnergyFunc(self._cDim, self._vDim, id, polyOrder)
   elseif mom == "pressure-tensor" then
      self._momCalcFun = pickPressureFunc(self._cDim, self._vDim, id, polyOrder)
   elseif mom == "heat-flux-tensor" then
      assert(true, "NYI!")
   elseif mom == "heat-flux-vector" then
      assert(true, "NYI!")
   else
      assert(false, "Did not recognize moment type " .. mom)
   end

   -- construct various functions from template representations
   self._calcGridQuant = loadstring(calcGridQuantTempl {VDIM = self._vDim, CDIM = self._cDim} )()
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(DistFuncMomentCalc, { __call = function (self, o) return self.new(self, o) end })

-- advance method
local function advance(self, tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local distf = assert(inFld[1], "DistFuncMomentCalc:advance: Must supply input distribution function")
   local mom = assert(outFld[1], "DistFuncMomentCalc:advance: Must supply output moment field")

   -- sanity checks
   assert(distf:ndim() == self._pDim,
	  "DistFuncMomentCalc:advance: Phase-space dimensions do not match input distribution function ")
   assert(mom:ndim() == self._cDim,
	  "DistFuncMomentCalc:advance: Config-space dimensions do not match output moment field")

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local dv = Lin.Vec(vDim) -- cell shape
   local xcp = Lin.Vec(pDim) -- phase-space cell center
   local w = Lin.Vec(vDim) -- cell center
   local momCell = Lin.Vec(mom:numComponents()) -- moments in a cell
   
   local localExtRange = distf:localExtRange()
   local distfIndexer = distf:genIndexer()
   local momIndexer = mom:genIndexer()

   local distfItr, momItr = distf:get(0), mom:get(0)

   scaleField(mom, 0.0) -- zero out moments
   local refVol = 2^vDim -- volume of velocity space reference cell
   
   -- loop, computing moments in each cell
   for idx in localExtRange:colMajorIter() do
      -- get cell shape, center coordinates
      grid:setIndex(idx); grid:cellCenter(xcp)
      -- compute velocity cell spacing, cell-center and volume
      -- local vol = self._calcGridQuant(grid, xcp, dv, w)
      local vol = 1.0
      for d = 1, vDim do
      	 dv[d] = grid:dx(cDim+d)
      	 w[d] = xcp[cDim+d]
      	 vol = vol*dv[d]
      end

      -- compute moment in cell
      distf:fill(distfIndexer(idx), distfItr)
      self._momCalcFun(distfItr, momCell, dv, w)

      -- accumulate it with proper normalization factor (Vv/2^dv)
      mom:fill(momIndexer(idx), momItr)
      for k = 1, mom:numComponents() do
      	 momItr[k] = momItr[k] + vol/refVol*momCell[k]
      end
   end
   
   return true, GKYL_MAX_DOUBLE
end

-- Methods in updater
DistFuncMomentCalc.__index = { advance = Base.advanceFuncWrap(advance) }

return {
   DistFuncMomentCalc = DistFuncMomentCalc
}
