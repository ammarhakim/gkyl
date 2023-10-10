-- Gkyl ------------------------------------------------------------------------
--
-- Basic boundary condition for a Vlasov species, i.e. those that can be
-- applied with Updater/Bc.lua using just a function (e.g. bcAbsorb, bcOpen)
-- and no additional setup.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase     = require "App.BCs.BCsBase"
local DataStruct  = require "DataStruct"
local Updater     = require "Updater"
local Mpi         = require "Comm.Mpi"
local Proto       = require "Lib.Proto"
local Time        = require "Lib.Time"
local Range       = require "Lib.Range"
local Lin         = require "Lib.Linalg"
local Grid        = require "Grid"
local DiagsApp    = require "App.Diagnostics.SpeciesDiagnostics"
local VlasovDiags = require "App.Diagnostics.VlasovDiagnostics"
local xsys        = require "xsys"

local VlasovBasicBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function VlasovBasicBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VlasovBasicBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   if tbl.kind=="function" or tbl.bcFunction then
      assert(type(tbl.bcFunction)=="function", "VlasovBasicBC: bcFunction must be a function.")
      self.bcKind   = "function"
      self.bcFuncIn = assert(tbl.bcFunction, "VlasovBasicBC: must specify the BC function in 'bcFunc' when using 'function' BC kind.")
      self.feedback = xsys.pickBool(tbl.feedback, false) 
      self.evolve   = xsys.pickBool(tbl.evolve, self.feedback) 
   else
      self.bcKind = assert(tbl.kind, "VlasovBasicBC: must specify the type of BC in 'kind'.")
   end

   self.saveFlux = tbl.saveFlux or false
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
         self.saveFlux       = true
      end
   end
end

function VlasovBasicBC:setName(nm) self.name = self.speciesName.."_"..nm end

function VlasovBasicBC:bcOpen(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "pointwise".
   self.basis:flipSign(dir, fIn:data(), fOut:data())
end

function VlasovBasicBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()

   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   local globalGhostRange = self.bcEdge=="lower" and distf:localGhostRangeLower()[self.bcDir]
                                                  or distf:localGhostRangeUpper()[self.bcDir]
   self:createBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())


   -- Need to define methods to allocate fields defined on boundary grid (used by diagnostics).
   self.allocCartField = function(self, grid, nComp, ghosts, metaData)
      local f = DataStruct.Field {
         onGrid        = grid,   ghost    = ghosts,
         numComponents = nComp,  metaData = metaData,
      }
      f:clear(0.0)
      return f
   end
   local allocDistf = function()
      return self:allocCartField(self.boundaryGrid, self.basis:numBasis(), {0,0}, distf:getMetaData())
   end

   self.bcBuffer = allocDistf() -- Buffer used by BasicBc updater.

   local bcFunc, skinType
   if self.bcKind == "copy" or self.bcKind == "absorb" or self.bcKind == "reflect" then
      self.bcSolver = Updater.BasicBc{
         onGrid  = self.grid,   edge   = self.bcEdge,  
         cdim    = self.cdim,   basis  = self.basis,
         dir     = self.bcDir,  bcType = self.bcKind,
         onField = mySpecies:rkStepperFields()[1],
      }
   else
      -- g2, to be deleted.
      if self.bcKind == "open" then
         bcFunc   = function(...) return self:bcOpen(...) end
         skinType = "pointwise"
      elseif self.bcKind == "function" then
         bcFunc   = function(...) return self:bcCopy(...) end
         skinType = "pointwise"
      else
         assert(false, "VlasovBasicBC: BC kind not recognized.")
      end

      local vdir = self.bcDir+self.cdim

      self.bcSolver = Updater.Bc{
         onGrid   = self.grid,   edge               = self.bcEdge,  
         cdim     = self.cdim,   boundaryConditions = {bcFunc},   
         dir      = self.bcDir,  evaluate           = self.bcFuncIn,
         vdir     = vdir,        evolveFn           = self.evolve,
         skinLoop = skinType,    feedback           = self.feedback,
         basis    = self.basis,  confBasis          = self.confBasis,
         advanceArgs = {{mySpecies:rkStepperFields()[1]}, {mySpecies:rkStepperFields()[1]}},
      }
   end

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then

      -- Part of global ghost range this rank owns.
      self.myGlobalGhostRange = self.bcEdge=="lower" and distf:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                      or distf:localGlobalGhostRangeIntersectUpper()[self.bcDir]

      -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
      self:createConfBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())

      self.allocMoment = function(self)
         return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
      end
      self.allocIntThreeMoments = function(self)
         return self:allocCartField(self.confBoundaryGrid, self.vdim+2, {0,0}, numDensity:getMetaData())
      end
      self.allocVectorMoment = function(self, dim)
         return self:allocCartField(self.confBoundaryGrid, dim*self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
      end
      self.allocIntMoment = function(self, comp)
         local metaData = {charge = self.charge,  mass = self.mass,}
         local ncomp = comp or 1
         local gridWriteRank = self.confBoundaryGrid:commSet().writeRank
         local f = DataStruct.DynVector{numComponents = ncomp,     writeRank = gridWriteRank<0 and gridWriteRank or 0,
                                        metaData      = metaData,  comm      = self.confBoundaryGrid:commSet().comm,}
         return f
      end

      -- Allocate fields needed.
      self.boundaryFluxFields = {}  -- Fluxes through the boundary, into ghost region, from each RK stage.
      self.distfInIdxr        = distf:genIndexer()
      for i = 1, #mySpecies:rkStepperFields() do
         self.boundaryFluxFields[i] = allocDistf()
      end
      self.boundaryFluxRate      = allocDistf()
      self.boundaryFluxFieldPrev = allocDistf()

      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      self.confBoundaryField = self:allocMoment()
      -- Range spanning ghost cells.
      self.myGlobalConfGhostRange = self.bcEdge=="lower" and numDensity:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                          or numDensity:localGlobalGhostRangeIntersectUpper()[self.bcDir]

      local jacobGeo = externalField.geo and externalField.geo.jacobGeo
      if jacobGeo then
         self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeo:getMetaData())
         self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo, self.confBoundaryField))
      end
      local jacobGeoInv = externalField.geo and externalField.geo.jacobGeoInv
      if jacobGeoInv then
         self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeoInv:getMetaData())
         self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv, self.confBoundaryField))
      end

      self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut)
         self.boundaryFluxFields[rkIdx]:copyRangeToRange(qOut, self.boundaryFluxFields[rkIdx]:localRange(), self.myGlobalGhostRange)
      end
      self.copyBoundaryFluxFieldFunc = function(inIdx, outIdx)
         self.boundaryFluxFields[outIdx]:copy(self.boundaryFluxFields[inIdx])
      end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...)
         local args  = {...} -- Package up rest of args as table.
         local nFlds = #args/2
         self.boundaryFluxFields[outIdx]:combine(a, self.boundaryFluxFields[aIdx])
         for i = 1, nFlds do -- Accumulate rest of the fields.
            self.boundaryFluxFields[outIdx]:accumulate(args[2*i-1], self.boundaryFluxFields[args[2*i]])
         end
      end

      -- Number density calculator. Needed regardless of diagnostics (for recycling BCs).
      self.numDensityCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
         phaseBasis = self.basis,         moment     = "M0",
      }

      -- Integrated number density calculator. Needed regardless of diagnostics (for steady state sources).
      self.integNumDensityCalc = Updater.DistFuncMomentDG {
         onGrid     = self.boundaryGrid,   confBasis  = self.confBasis,
         phaseBasis = self.basis,          moment     = "M0",
         isIntegrated = true,              model_id   = "GKYL_MODEL_DEFAULT", 
      }

      if not self.anyDiagnostics then
         self.calcBoundaryFluxRateFunc = function(dtIn) end
      else
         self.calcBoundaryFluxRateFunc = function(dtIn)
            -- Compute boundary flux rate ~ (fGhost_new - fGhost_old)/dt.
            self.boundaryFluxRate:combine( 1.0/dtIn, self.boundaryFluxFields[1],
                                          -1.0/dtIn, self.boundaryFluxFieldPrev)
            self.boundaryFluxFieldPrev:copy(self.boundaryFluxFields[1])
         end
         -- Set up weak multiplication and division operators (for diagnostics).
         self.confWeakMultiply = Updater.CartFieldBinOp {
            weakBasis = self.confBasis,  operation = "Multiply",
            onGhosts  = true,
         }
         self.confWeakDivide = Updater.CartFieldBinOp {
            weakBasis = self.confBasis,  operation = "Divide",
            onRange   = self.confBoundaryField:localRange(),  onGhosts = false,
         }
         self.confWeakDotProduct = Updater.CartFieldBinOp {
            weakBasis = self.confBasis,  operation = "DotProduct",
            onGhosts  = true,
         }
         -- Volume integral operator (for diagnostics).
         self.volIntegral = {
            scalar = Updater.CartFieldIntegratedQuantCalc {
               onGrid = self.confBoundaryGrid,  numComponents = 1,
               basis  = self.confBasis,         quantity      = "V",
            },
            vector = Updater.CartFieldIntegratedQuantCalc {
               onGrid = self.confBoundaryGrid,  numComponents = self.vdim,
               basis  = self.confBasis,         quantity      = "V",
            },
         }
         -- Moment calculators (for diagnostics).
         self.momDensityCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
            phaseBasis = self.basis,         moment     = "M1i",
         }
         self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis  = self.confBasis,
            phaseBasis = self.basis,         moment     = "M2",
         }
         self.M2ijCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         moment    = "M2ij",
         }
         self.M3iCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         moment    = "M3i",
         }
         self.divideByJacobGeo = self.jacobGeoInv
            and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeoInv}, {fldOut}) end
            or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
         self.multiplyByJacobGeo = self.jacobGeo
            and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeo}, {fldOut}) end
            or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
      end
   else
      self.storeBoundaryFluxFunc        = function(tCurr, rkIdx, qOut) end
      self.copyBoundaryFluxFieldFunc    = function(inIdx, outIdx) end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...) end
      self.calcBoundaryFluxRateFunc     = function(dtIn) end
   end
end

function VlasovBasicBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function VlasovBasicBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function VlasovBasicBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function VlasovBasicBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function VlasovBasicBC:createDiagnostics(mySpecies, field)
   -- Create BC diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = VlasovDiags()}
      self.diagnostics:fullInit(mySpecies, field, self)
      -- Presently boundary diagnostics are boundary flux diagnostics. Append 'flux' to the diagnostic's
      -- name so files are named accordingly. Re-design this when non-flux diagnostics are implemented
      self.diagnostics.name = self.diagnostics.name..'_flux'
   end
   return self.diagnostics
end

-- These are needed to recycle the VlasovDiagnostics with VlasovBasicBC.
function VlasovBasicBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                                 self.boundaryFluxRate, self.boundaryFluxRate} end
function VlasovBasicBC:getFlucF() return self.boundaryFluxRate end

function VlasovBasicBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   local fIn = mySpecies:rkStepperFields()[outIdx]

   self.bcSolver:advance(tCurr, {self.bcBuffer}, {fIn})
end

function VlasovBasicBC:getBoundaryFluxFields() return self.boundaryFluxFields end

-- ................... Classes meant as aliases to simplify input files ...................... --
local VlasovAbsorbBC = Proto(VlasovBasicBC)
function VlasovAbsorbBC:fullInit(mySpecies)
   self.tbl.kind  = "absorb"
   VlasovAbsorbBC.super.fullInit(self, mySpecies)
end

local VlasovReflectBC = Proto(VlasovBasicBC)
function VlasovReflectBC:fullInit(mySpecies)
   self.tbl.kind  = "reflect"
   VlasovReflectBC.super.fullInit(self, mySpecies)
end

local VlasovCopyBC = Proto(VlasovBasicBC)
function VlasovCopyBC:fullInit(mySpecies)
   self.tbl.kind  = "copy"
   VlasovCopyBC.super.fullInit(self, mySpecies)
end

local VlasovOpenBC = Proto(VlasovBasicBC)
function VlasovOpenBC:fullInit(mySpecies)
   self.tbl.kind  = "copy"
   VlasovOpenBC.super.fullInit(self, mySpecies)
end

local VlasovZeroFluxBC = Proto(BCsBase)
function VlasovZeroFluxBC:init(tbl)
   self.tbl      = tbl
   self.tbl.kind = "zeroFlux"
end
-- ................... End of VlasovBasicBC alias classes .................... --

return {VlasovBasic    = VlasovBasicBC,
        VlasovAbsorb   = VlasovAbsorbBC,
        VlasovCopy     = VlasovCopyBC,
        VlasovOpen     = VlasovOpenBC,
        VlasovReflect  = VlasovReflectBC,
        VlasovZeroFlux = VlasovZeroFluxBC}
