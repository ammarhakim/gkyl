-- Gkyl ------------------------------------------------------------------------
--
-- Adiabatic species
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local FluidSpecies = require "App.Species.FluidSpecies"

local AdiabaticSpecies = Proto(FluidSpecies)

function AdiabaticSpecies:fullInit(appTbl)
   AdiabaticSpecies.super.fullInit(self, appTbl)

   self.nMoments = 1

   self._temp = self.tbl.temp

   assert(self.evolve==false, "AdiabaticSpecies: cannot evolve an adiabatic species")
end

--function AdiabaticSpecies:initDist()
--   self.moments[1]:combine(1/self.charge, kineticChargeDens)
--end

function AdiabaticSpecies:createSolver(hasE, hasB, funcField)

   -- compute density in center of domain
   local gridCenter = {}
   for d = 1, self.ndim do
      gridCenter[d] = (self.grid:upper(d) + self.grid:lower(d))/2
   end
   self._dens0 = self.initFunc(0, gridCenter)

   self.qneutFac = self:allocMoment()

   -- set up jacobian for general geometry
   -- and use it to scale initial density
   -- (consistent with scaling distribution function with jacobian)
   if funcField then
      self.jacobGeoFunc = funcField.jacobGeoFunc

      local initFuncWithoutJacobian = self.initFunc
      self.initFunc = function (t, xn)
         local J = self.jacobGeoFunc(t,xn)
         local f = initFuncWithoutJacobian(t,xn)
         return J*f
      end
   end

end

-- nothing to calculate, just copy
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

-- for interfacing with GkField
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

-- this is factor on potential in qneut equation
function AdiabaticSpecies:getQneutFac(linearized)
   if linearized == false then
      self.qneutFac:combine(self.charge^2/self._temp, self.couplingMoments)
      return self.qneutFac
   else
      return self._dens0*self.charge^2/self._temp
   end
end

return AdiabaticSpecies
