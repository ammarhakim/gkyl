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
local Updater      = require "Updater"

local AdiabaticSpecies = Proto(FluidSpecies)

-- ............. Backwards compatible treatment of BCs .....................--
-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB = 1
local SP_BC_COPY   = 6
AdiabaticSpecies.bcAbsorb = SP_BC_ABSORB   -- Absorb all particles.
AdiabaticSpecies.bcCopy   = SP_BC_COPY     -- Copy stuff.

function AdiabaticSpecies:makeBcApp(bcIn)
   local bcOut
   print("AdiabaticSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
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

   -- Compute the factor multiplying the potential in quasineutrality equation.
   self.QneutFac = self:allocMoment()
   if field.linearizedPolarization then
      if field.uniformPolarization then
         -- Compute density in center of domain.
         local gridCenter = {}
         for d = 1, self.ndim do gridCenter[d] = self.grid:mid(d) end
         local dens0 = self.initFunc(0., gridCenter)
         self.QneutFac:combine(dens0*(self.charge^2)/self._temp, self.unitField)
      else
         local evOnNodes = Updater.EvalOnNodes {
            onGrid = self.grid,   evaluate = self.initFunc,
            basis  = self.basis,  onGhosts = false, --true,
         }
         evOnNodes:advance(0., {}, {self.QneutFac})

         -- Apply open BCs:
         local function makeOpenBcUpdater(dir, edge)
            local bcOpen = function(dir, tm, idxIn, fIn, fOut)
               self.confBasis:flipSign(dir, fIn, fOut)   -- Requires skinLoop = "pointwise".
            end
            return Updater.Bc {
               onGrid   = self.grid,  edge = edge,
               dir      = dir,        boundaryConditions = {bcOpen},
               skinLoop = "pointwise",
            }
         end
         openBCupdaters = {}
         for dir = 1, self.ndim do
            if not lume.any(self.grid:getPeriodicDirs(), function(t) return t==dir end) then
               openBCupdaters["lower"] = makeOpenBcUpdater(dir, "lower")
               openBCupdaters["upper"] = makeOpenBcUpdater(dir, "upper")
            end
         end
         for _, bc in pairs(openBCupdaters) do bc:advance(0., {}, {self.QneutFac}) end
         self.QneutFac:sync(true)

         self.QneutFac:scale((self.charge^2)/self._temp)
      end

      -- Scale by the Jacobian.
      if externalField.geo then
         local weakMultiply = Updater.CartFieldBinOp {
            onGrid    = self.grid,   operation = "Multiply",
            weakBasis = self.basis,  onGhosts  = true,
         }

         local jacob = externalField.geo.jacobGeo
         if jacob then weakMultiply:advance(0, {self.QneutFac, jacob}, {self.QneutFac}) end
      end
   
      self.getQneutFacFunc = function() return self.QneutFac end
   else
      self.getQneutFacFunc = function() 
         self.QneutFac:combine(self.charge^2/self._temp, self.couplingMoments)
         return self.QneutFac
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

function AdiabaticSpecies:getQneutFac()
   -- Return the factor on potential in charge neutrality equation.
   return self.getQneutFacFunc()
end

return AdiabaticSpecies
