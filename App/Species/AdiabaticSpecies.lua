-- Gkyl ------------------------------------------------------------------------
--
-- Adiabatic species.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto        = require "Lib.Proto"
local FluidSpecies = require "App.Species.FluidSpecies"

local AdiabaticSpecies = Proto(FluidSpecies)

-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB = 1
local SP_BC_COPY   = 6
AdiabaticSpecies.bcAbsorb = SP_BC_ABSORB      -- Absorb all particles.
AdiabaticSpecies.bcCopy   = SP_BC_COPY        -- Copy stuff.

function AdiabaticSpecies:fullInit(appTbl)
   AdiabaticSpecies.super.fullInit(self, appTbl)

   self.nMoments = 1

   self._temp = self.tbl.temp

   self.initFunc = self.tbl.init

   -- Adiabatic species does not require fluctuationBCs (MF 2021/04/10: is this always true?).
   self.fluctuationBCs = false

   assert(self.evolve==false, "AdiabaticSpecies: cannot evolve an adiabatic species")
end

function AdiabaticSpecies:createSolver(hasE, hasB, externalField)

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

end

function AdiabaticSpecies:bcAbsorbFunc(dir, tm, idxIn, fIn, fOut)
   -- The idea is that by setting the plasma quantities to zero in the
   -- ghost cell nothing is transported into the domain, and whatever is transported
   -- out is lost. We can't set them to exactly zero or else the sound speed
   -- and drift velocity would diverge, so we set them to something small.
   local numB = self.basis:numBasis()
   for i = 1, numB do fOut[i] = 1.e-10*fIn[i] end   -- Density. 
end

function AdiabaticSpecies:appendBoundaryConditions(dir, edge, bcType)
   -- Need to wrap member functions so that self is passed.
   local function bcAbsorbFunc(...)  return self:bcAbsorbFunc(...) end
   local function bcCopyFunc(...)    return self:bcCopyFunc(...) end

   if bcType == SP_BC_ABSORB then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, { bcAbsorbFunc }, "pointwise"))
   elseif bcType == SP_BC_COPY then
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, { bcCopyFunc }, "pointwise"))
   else
      assert(false, "AdiabaticSpecies: Unsupported BC type!")
   end
end

-- Nothing to calculate, just copy.
function AdiabaticSpecies:calcCouplingMoments(tCurr, rkIdx)
   if self.deltaF then
      self.couplingMoments:clear(0.0)
   else
      local fIn = self:rkStepperFields()[rkIdx]
      self.couplingMoments:copy(fIn)
   end
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

function AdiabaticSpecies:temp()
   return self._temp
end

function AdiabaticSpecies:dens0()
   return self._dens0
end

-- This is factor on potential in qneut equation.
function AdiabaticSpecies:getQneutFac(linearized)
   if linearized == false then
      self.qneutFac:combine(self.charge^2/self._temp, self.couplingMoments)
      return self.qneutFac
   else
      return self._dens0*self.charge^2/self._temp
   end
end

return AdiabaticSpecies
