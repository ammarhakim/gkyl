-- Gkyl ------------------------------------------------------------------------
--
-- Adiabatic species.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto        = require "Lib.Proto"
local FluidSpecies = require "App.Species.FluidSpecies"
local BasicBC      = require "App.BCs.AdiabaticBasic"

local AdiabaticSpecies = Proto(FluidSpecies)

-- ............. Backwards compatible treatment of BCs .....................--
-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB = 1
local SP_BC_COPY   = 6
AdiabaticSpecies.bcAbsorb = SP_BC_ABSORB   -- Absorb all particles.
AdiabaticSpecies.bcCopy   = SP_BC_COPY     -- Copy stuff.

function AdiabaticSpecies:makeBcApp(bcIn)
   local bcOut
   if bcIn == SP_BC_COPY then
      bcOut = BasicBC{kind="copy", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_ABSORB then
      bcOut = BasicBC{kind="absorb", diagnostics={}, saveFlux=false}
   end
   return bcOut
end

-- ............. End of backwards compatibility for BCs .....................--

function AdiabaticSpecies:fullInit(appTbl)
   AdiabaticSpecies.super.fullInit(self, appTbl)

   self.nMoments = 1

   self._temp = self.tbl.temp

   self.initFunc = self.tbl.init

   -- Adiabatic species does not require fluctuationBCs (MF 2021/04/10: is this always true?).
   self.fluctuationBCs = false

   assert(self.evolve==false, "AdiabaticSpecies: cannot evolve an adiabatic species")
end

function AdiabaticSpecies:createSolver(field, externalField)

   -- Compute density in center of domain.
   local gridCenter = {}
   for d = 1, self.ndim do gridCenter[d] = (self.grid:upper(d) + self.grid:lower(d))/2 end
   self._dens0 = self.initFunc(0., gridCenter)

   self.qneutFac = self:allocMoment()

   -- Set up jacobian for general geometry
   -- and use it to scale initial density
   -- (consistent with scaling distribution function with jacobian).
   if externalField then
      self.jacobGeoFunc = externalField.jacobGeoFunc

      local initFuncWithoutJacobian = self.initFunc
      self.initFunc = function (t, xn)
         local J = self.jacobGeoFunc(t,xn)
         local f = initFuncWithoutJacobian(t,xn)
         return J*f
      end
   end

   if self.deltaF then
      self.calcCouplingMomentsFunc = function(tCurr, rkIdx)
         self.couplingMoments:clear(0.0)
      end
   else
      self.calcCouplingMomentsFunc = function(tCurr, rkIdx)
         local fIn = self:rkStepperFields()[rkIdx]
         self.couplingMoments:copy(fIn)
      end
   end

   self.suggestDtFunc = function() return FluidSpecies["suggestDtDontEvolve"](self) end
   self.applyBcFunc   = function(tCurr, momIn) return FluidSpecies["applyBcDontEvolve"](self, tCurr, momIn) end

   -- Empty methods needed in case positivity is used (for the kinetic species).
   self.checkPositivity      = function(tCurr, idx) end
   self.posRescaler          = {advance=function(tCurr, inFlds, outFlds, computeDiagnostics, zeroOut) end}
   self.posRescalerDiffAdv   = function(tCurr, rkIdx, computeDiagnostics, zeroOut) end
   self.posRescalerDiffWrite = function(tm, fr) end
end

function AdiabaticSpecies:calcCouplingMomentsEvolve(tCurr, rkIdx)
   self.calcCouplingMomentsFunc(tCurr, rkIdx)
end

function AdiabaticSpecies:fluidMoments()
   return { self.couplingMoments }
end

-- For interfacing with GkField.
function AdiabaticSpecies:getNumDensity(rkIdx)
   if rkIdx == nil then return self.couplingMoments 
   else 
      self.couplingMoments:copy(self:rkStepperFields()[rkIdx])
      return self.couplingMoments
   end
end

function AdiabaticSpecies:temp() return self._temp end

function AdiabaticSpecies:dens0() return self._dens0 end

function AdiabaticSpecies:getQneutFacLin()
   -- Return the factor on potential in in charge neutrality equation,
   -- assuming a linearized polarization.
   return self._dens0*self.charge^2/self._temp
end
function AdiabaticSpecies:getQneutFacNotLin()
   -- Return the factor on potential in in charge neutrality equation,
   -- not assuming a linearized polarization.
   self.qneutFac:combine(self.charge^2/self._temp, self.couplingMoments)
   return self.qneutFac
end

return AdiabaticSpecies
