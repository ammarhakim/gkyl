-- Gkyl ------------------------------------------------------------------------
--
-- Basic boundary condition for a Gyrokinetic species, i.e. those that can be
-- applied with Updater/Bc.lua using just a function (e.g. bcAbsorb, bcOpen)
-- and no additional setup.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase    = require "App.BCs.BCsBase"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"
local Mpi        = require "Comm.Mpi"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"
local Range      = require "Lib.Range"
local Lin        = require "Lib.Linalg"
local Grid       = require "Grid"
local DiagsApp   = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags    = require "App.Diagnostics.GkDiagnostics"
local xsys       = require "xsys"

local GkBasicBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function GkBasicBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkBasicBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   if tbl.kind=="function" or tbl.bcFunction then
      assert(type(tbl.bcFunction)=="function", "GkBasicBC: bcFunction must be a function.")
      self.bcKind   = "function"
      self.bcFuncIn = assert(tbl.bcFunction, "GkBasicBC: must specify the BC function in 'bcFunc' when using 'function' BC kind.")
      self.feedback = xsys.pickBool(tbl.feedback, false) 
   else
      self.bcKind = assert(tbl.kind, "GkBasicBC: must specify the type of BC in 'kind'.")

      self.phiWallFunc = tbl.phiWall
      if self.phiWallFunc then assert(type(self.phiWallFunc)=="function", "GkBasicBC: phiWall must be a function (t, xn).") end
   end
   self.evolve = xsys.pickBool(tbl.evolve, self.feedback or false) 

   self.saveFlux = tbl.saveFlux or false
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
         self.saveFlux       = true
      end
   end
   self.createBoundaryTools = function(mySpecies, field, externalField) return self:createBoundaryToolsGK(mySpecies, field, externalField)
end

function GkBasicBC:setName(nm) self.name = self.speciesName.."_"..nm end

function GkBasicBC:bcOpen(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "pointwise".
   self.basis:flipSign(dir, fIn:data(), fOut:data())
end

function GkBasicBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()
   self:createBoundaryTools(mySpecies,field,externalField)
   self.bcBuffer = self.allocDistf() -- Buffer used by BasicBc updater.

   -- Sheath BCs use phi and phiWall. Set here so non-sheath BCs call these empty functions.
   self.setPhiWall = {advance = function(tCurr,inFlds,OutFlds) end}
   self.getPhi     = function(fieldIn, inIdx) return nil end 

   local bcFunc, skinType
   if self.bcKind == "copy" or self.bcKind == "absorb" or self.bcKind == "reflect" then
      self.bcSolver = Updater.BasicBc{
         onGrid  = self.grid,   edge   = self.bcEdge,  
         cdim    = self.cdim,   basis  = self.basis,
         dir     = self.bcDir,  bcType = self.bcKind,
         onField = mySpecies:rkStepperFields()[1],
      }
      self.bcSolverAdvance = function(tm, inFlds, outFlds)
         self.bcSolver:advance(tm, {inFlds[1]}, outFlds)
      end
   elseif self.bcKind == "sheath" then
      local charge, mass = mySpecies.charge, mySpecies.mass

      self.getPhi = function(fieldIn, inIdx) return fieldIn:rkStepperFields()[inIdx].phi end 

      -- Create field and function for calculating wall potential according to user-provided function.
      self.phiWallFld = DataStruct.Field {
         onGrid        = field.grid,
         numComponents = field.basis:numBasis(),
         ghost         = {1,1},
         metaData      = {polyOrder = field.basis:polyOrder(),
                          basisType = field.basis:id()},
         syncPeriodicDirs = false,
      }
      self.phiWallFld:clear(0.)
      if self.phiWallFunc then
         self.setPhiWall = Updater.EvalOnNodes {
            onGrid = field.grid,   evaluate = self.phiWallFunc,
            basis  = field.basis,  onGhosts = true,
         }
         self.setPhiWall:advance(0.0, {}, {self.phiWallFld})
      end
      self.phiWallFld:sync(false)

      self.bcSolver = Updater.GkSheathBc{
         onGrid  = self.grid,   edge    = self.bcEdge,  
         cdim    = self.cdim,   basis   = self.basis,
         dir     = self.bcDir,  phiWall = self.phiWallFld,
         charge  = charge,      mass    = mass,
         onField = mySpecies:rkStepperFields()[1],
      }
      self.bcSolverAdvance = function(tm, inFlds, outFlds)
         self.bcSolver:advance(tm, {inFlds[2]}, outFlds)
      end
   else
      -- g2, to be deleted.
      if self.bcKind == "open" then
         bcFunc   = function(...) return self:bcOpen(...) end
         skinType = "pointwise"
      elseif self.bcKind == "function" then
         bcFunc   = function(...) return self:bcCopy(...) end
         skinType = "pointwise"
      else
         assert(false, "GkBasicBC: BC kind not recognized.")
      end

      local vdir = nil
      if self.bcDir==self.cdim then vdir = self.cdim+1 end

      self.bcSolver = Updater.Bc{
         onGrid   = self.grid,   edge               = self.bcEdge,  
         cdim     = self.cdim,   boundaryConditions = {bcFunc},   
         dir      = self.bcDir,  evaluate           = self.bcFuncIn,
         vdir     = vdir,        evolveFn           = self.evolve,
         skinLoop = skinType,    feedback           = self.feedback,
         basis    = self.basis,  confBasis          = self.confBasis,
         advanceArgs = {{mySpecies:rkStepperFields()[1]}, {mySpecies:rkStepperFields()[1]}},
      }
      self.bcSolverAdvance = function(tm, inFlds, outFlds)
         self.bcSolver:advance(tm, {}, outFlds)
      end
   end
end

function GkBasicBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function GkBasicBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function GkBasicBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function GkBasicBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function GkBasicBC:createDiagnostics(mySpecies, field)
   -- Create BC diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = GkDiags()}
      self.diagnostics:fullInit(mySpecies, field, self)
      -- Presently boundary diagnostics are boundary flux diagnostics. Append 'flux' to the diagnostic's
      -- name so files are named accordingly. Re-design this when non-flux diagnostics are implemented
      self.diagnostics.name = self.diagnostics.name..'_flux'
   end
   return self.diagnostics
end

-- These are needed to recycle the GkDiagnostics with GkBasicBC.
function GkBasicBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                             self.boundaryFluxRate, self.boundaryFluxRate} end
function GkBasicBC:getFlucF() return self.boundaryFluxRate end

function GkBasicBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)

   self.setPhiWall:advance(tCurr, {}, {self.phiWallFld}) -- Compute wall potential if needed (i.e. sheath BC).
   self.phi = self.getPhi(field, inIdx)              -- If needed get the current plasma potential (for sheath BC).

   local fIn = mySpecies:rkStepperFields()[outIdx] 

   self.bcSolverAdvance(tCurr, {self.bcBuffer, self.phi}, {fIn})
end

function GkBasicBC:getBoundaryFluxFields() return self.boundaryFluxFields end

-- ................... Classes meant as aliases to simplify input files ...................... --
local GkAbsorbBC = Proto(GkBasicBC)
function GkAbsorbBC:fullInit(mySpecies)
   self.tbl.kind  = "absorb"
   GkAbsorbBC.super.fullInit(self, mySpecies)
end

local GkReflectBC = Proto(GkBasicBC)
function GkReflectBC:fullInit(mySpecies)
   self.tbl.kind  = "reflect"
   GkReflectBC.super.fullInit(self, mySpecies)
end

local GkCopyBC = Proto(GkBasicBC)
function GkCopyBC:fullInit(mySpecies)
   self.tbl.kind  = "copy"
   GkCopyBC.super.fullInit(self, mySpecies)
end

local GkOpenBC = Proto(GkBasicBC)
function GkOpenBC:fullInit(mySpecies)
   self.tbl.kind  = "copy"
   GkOpenBC.super.fullInit(self, mySpecies)
end

local GkZeroFluxBC = Proto()
function GkZeroFluxBC:init(tbl)
   self.tbl      = tbl
   self.tbl.kind = "zeroFlux"
end

local GkSheathBC = Proto(GkBasicBC)
function GkSheathBC:fullInit(mySpecies)
   self.tbl.kind  = "sheath"
   GkSheathBC.super.fullInit(self, mySpecies)
end
-- ................... End of GkBasicBC alias classes .................... --

return {GkBasic    = GkBasicBC,
        GkAbsorb   = GkAbsorbBC,
        GkCopy     = GkCopyBC,
        GkOpen     = GkOpenBC,
        GkReflect  = GkReflectBC,
        GkZeroFlux = GkZeroFluxBC,
        GkSheath   = GkSheathBC}
