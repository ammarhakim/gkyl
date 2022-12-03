-- Gkyl ------------------------------------------------------------------------
--
-- Euler Species
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Basis          = require "Basis"
local DataStruct     = require "DataStruct"
local Grid           = require "Grid"
local FluidSpecies = require "App.Species.FluidSpecies"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"
local Projection     = require "App.Projection"
local Proto          = require "Lib.Proto"
local Source         = require "App.Sources.VmSource"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local DiagsApp       = require "App.Diagnostics.SpeciesDiagnostics"
local BCsBase        = require "App.BCs.BCsBase"
local xsys           = require "xsys"
local lume           = require "Lib.lume"

local EulerSpecies = Proto(FluidSpecies)

function EulerSpecies:createGrid(...)
   EulerSpecies.super.createGrid(self, ...)

   self.nMoments = self.ndim + 1 --TODO: this should be nMoments + 2 when euler is used
end

function EulerSpecies:alloc(nRkDup)
   -- Allocate statevec [rho, rhoux, rhouy, rhouz] for isoeuler, [rho, rhoux, rhouy, pres] for euler
   self.nMoments = 4 --TODO: add euler (i.e. nMoments = 5 when doing euler)
   EulerSpecies.super.alloc(self, nRkDup)

end

function EulerSpecies:fullInit(appTbl)
   EulerSpecies.super.fullInit(self, appTbl)

   local tbl = self.tbl
end

function EulerSpecies:createSolver(field, externalField)


   EulerSpecies.super.createSolver(self)--, field)

   -- Alloc aux var that contains u_i  = rho u_i / rho from weak division
   self.uvar = self.super.allocVectorMoment(self,3)

   -- Need to overwrite weak Division operator from Fluid Species
   self.weakDivide = Updater.CartFieldBinOp {
      operation = "Divide",  weakBasis = self.basis,
      onRange   = self.uvar:localExtRange(),  onGhosts = true,
   }

   self.confBasis = self.basis
   -- Create updater to advance solution by one time-step.
   self.solver = Updater.FluidDG {
      onGrid     = self.grid,
      confBasis  = self.confBasis,
      confRange  = self.uvar:localRange(),
      eqnId = "GKYL_EQN_ISO_EULER" --TODO: dont hard code this "GKYL_EQN_ISO_EULER"
    }

end

-- Nothing to calculate, just copy.
function EulerSpecies:calcCouplingMomentsEvolve(tCurr, rkIdx)
   local fIn = self:rkStepperFields()[rkIdx]
end

function EulerSpecies:setActiveRKidx(rkIdx) self.activeRKidx = rkIdx end

function EulerSpecies:advance(tCurr, population, emIn, inIdx, outIdx)
   self:setActiveRKidx(inIdx)
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]
   fRhsOut:clear(0.0)

   -- Do weak division --TODO: figure out how to separte these rho and {rhoux,rhouy,rhouz}
   --TODO: might need to pass through emIn (think this is used for aux fields)
   -- local rho = fIn[1] --WARNING: This is not how to separate state vec
   -- local rhouivec = {fIn[2],fIn[3],fIn[4]}
   -- self.weakDivide:advance(tm, {rhouivec,rho}, {self.uvar})

   -- TODO: add explicit diffusion if isoeuler

   -- update system
   self.solver:advance(tCurr, {fIn, self.uvar}, {self.moments, self.cflRateByCell}) --(Something like this from Gyrofluid)

   for _, bc in pairs(self.nonPeriodicBCs) do
      bc:storeBoundaryFlux(tCurr, outIdx, fRhsOut)   -- Save boundary fluxes.
   end
end

function EulerSpecies:projToSource(proj)
   -- For backwards compatibility: in the case where the source is specified
   -- as a projection object in the input file, this function turns that
   -- projection object into a Source object.
   local tbl = proj.tbl
   local pow = tbl.power
   return Source { profile = proj, power = pow }
end

return EulerSpecies
