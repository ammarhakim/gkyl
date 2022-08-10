-- Gkyl ------------------------------------------------------------------------
--
-- Basic boundary condition for a Vlasov species, i.e. those that can be
-- applied with Updater/Bc.lua using just a function (e.g. bcAbsorb, bcOpen)
-- and no additional setup.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase        = require "App.BCs.BCsBase"
local DataStruct     = require "DataStruct"
local Updater        = require "Updater"
local Mpi            = require "Comm.Mpi"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Range          = require "Lib.Range"
local Lin            = require "Lib.Linalg"
local LinearDecomp   = require "Lib.LinearDecomp"
local CartDecomp     = require "Lib.CartDecomp"
local Grid           = require "Grid"
local DiagsApp       = require "App.Diagnostics.SpeciesDiagnostics"
local VlasovDiags    = require "App.Diagnostics.VlasovDiagnostics"
local xsys           = require "xsys"

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
   local vdir = self.bcDir+self.cdim

   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   self:createBoundaryGrid()

   local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()
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
      return self:allocCartField(self.boundaryGrid, self.basis:numBasis(),
                                 {distf:lowerGhost(),distf:upperGhost()}, distf:getMetaData())
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

      -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
      self:createConfBoundaryGrid()

      self.allocMoment = function(self)
         return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                    {numDensity:lowerGhost(),numDensity:upperGhost()}, numDensity:getMetaData())
      end
      self.allocVectorMoment = function(self, dim)
         return self:allocCartField(self.confBoundaryGrid, dim*self.basis:numBasis(),
                                    {numDensity:lowerGhost(),numDensity:upperGhost()}, numDensity:getMetaData())
      end

      -- Allocate fields needed.
      self.boundaryFluxFields = {}  -- Fluxes through the boundary, into ghost region, from each RK stage.
      self.boundaryPtr        = {}
      self.distfInIdxr        = distf:genIndexer()
      for i = 1, #mySpecies:rkStepperFields() do
         self.boundaryFluxFields[i] = allocDistf()
         self.boundaryPtr[i]        = self.boundaryFluxFields[i]:get(1)
      end
      self.boundaryFluxRate      = allocDistf()
      self.boundaryFluxFieldPrev = allocDistf()
      self.boundaryIdxr          = self.boundaryFluxFields[1]:genIndexer()

      self.idxOut = Lin.IntVec(self.grid:ndim())

      -- Create the range needed to loop over ghosts.
      local global, globalExt, localExtRange = distf:globalRange(), distf:globalExtRange(), distf:localExtRange()
      self.ghostRange = localExtRange:intersect(self:getGhostRange(global, globalExt))
      -- Decompose ghost region into threads.
      self.ghostRangeDecomp = LinearDecomp.LinearDecompRange{range=self.ghostRange, numSplit=self.grid:numSharedProcs()}
      self.tId              = self.grid:subGridSharedId() -- Local thread ID.

      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      self.confBoundaryField    = self:allocMoment()
      self.confBoundaryFieldPtr = self.confBoundaryField:get(1)
      self.confBoundaryIdxr     = self.confBoundaryField:genIndexer()
      local confGlobal        = numDensity:globalRange()
      local confGlobalExt     = numDensity:globalExtRange()
      local confLocalExtRange = numDensity:localExtRange()
      self.confGhostRange = confLocalExtRange:intersect(self:getGhostRange(confGlobal, confGlobalExt)) -- Range spanning ghost cells.
      -- Decompose ghost region into threads.
      self.confGhostRangeDecomp = LinearDecomp.LinearDecompRange {range=self.confGhostRange, numSplit=self.grid:numSharedProcs()}

      local jacobGeo = externalField.geo and externalField.geo.jacobGeo
      if jacobGeo then
         self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                             {jacobGeo:lowerGhost(),jacobGeo:upperGhost()}, jacobGeo:getMetaData())
         self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo))
      end
      local jacobGeoInv = externalField.geo and externalField.geo.jacobGeoInv
      if jacobGeoInv then
         self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                                {jacobGeoInv:lowerGhost(),jacobGeoInv:upperGhost()}, jacobGeoInv:getMetaData())
         self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv))
      end

      self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut)
         local ptrOut = qOut:get(1)
         for idx in self.ghostRangeDecomp:rowMajorIter(self.tId) do
            idx:copyInto(self.idxOut)
            qOut:fill(self.distfInIdxr(idx), ptrOut)

            -- Before operating on ghosts, store ghost values for later flux diagnostics
            self.idxOut[self.bcDir] = 1
            self.boundaryFluxFields[rkIdx]:fill(self.boundaryIdxr(self.idxOut), self.boundaryPtr[rkIdx])
            for c = 1, qOut:numComponents() do self.boundaryPtr[rkIdx][c] = ptrOut[c] end
         end
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
            onGrid    = self.confBoundaryGrid,  operation = "Multiply",
            weakBasis = self.confBasis,         onGhosts  = true,
         }
         self.confWeakDivide = Updater.CartFieldBinOp {
            onGrid    = self.confBoundaryGrid,  operation = "Divide",
            weakBasis = self.confBasis,         onGhosts  = true,
         }
         self.confWeakDotProduct = Updater.CartFieldBinOp {
            onGrid    = self.confBoundaryGrid,  operation = "DotProduct",
            weakBasis = self.confBasis,         onGhosts  = true,
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

local VlasovZeroFluxBC = Proto()
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
