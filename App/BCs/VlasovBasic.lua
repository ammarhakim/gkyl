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

   self.diagnostics = tbl.diagnostics or {}
   self.saveFlux    = tbl.saveFlux or false
   if tbl.diagnostics then
     if #tbl.diagnostics>0 then self.saveFlux = true end
   end
end

function VlasovBasicBC:setName(nm) self.name = self.speciesName.."_"..nm end

function VlasovBasicBC:bcAbsorb(dir, tm, idxIn, fIn, fOut)
   -- Note that for bcAbsorb there is no operation on fIn,
   -- so skinLoop (which determines indexing of fIn) does not matter
   for i = 1, self.basis:numBasis() do fOut[i] = 0.0 end
end
function VlasovBasicBC:bcOpen(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "pointwise".
   self.basis:flipSign(dir, fIn, fOut)
end
function VlasovBasicBC:bcCopy(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "pointwise".
   for i = 1, self.basis:numBasis() do fOut[i] = fIn[i] end
end
function VlasovBasicBC:bcReflect(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "flip".
   self.basis:flipSign(dir, fIn, fOut)
   self.basis:flipSign(dir+self.cdim, fOut, fOut)
end

function VlasovBasicBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   local bcFunc, skinType
   if self.bcKind == "copy" then
      bcFunc   = function(...) return self:bcCopy(...) end
      skinType = "pointwise"
   elseif self.bcKind == "absorb" then
      bcFunc   = function(...) return self:bcAbsorb(...) end
      skinType = "pointwise"
   elseif self.bcKind == "open" then
      bcFunc   = function(...) return self:bcOpen(...) end
      skinType = "pointwise"
   elseif self.bcKind == "reflect" then
      bcFunc   = function(...) return self:bcReflect(...) end
      skinType = "flip"
   elseif self.bcKind == "function" then
      bcFunc   = function(...) return self:bcCopy(...) end
      skinType = "pointwise"
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

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then
      -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
      self:createBoundaryGrid()

      -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
      self:createConfBoundaryGrid()

      local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()
      -- Need to define methods to allocate fields defined on boundary grid (used by diagnostics).
      self.allocCartField = function(self, grid, nComp, ghosts, metaData)
         local f = DataStruct.Field {
            onGrid        = grid,
            numComponents = nComp,
            ghost         = ghosts,
            metaData      = metaData,
         }
         f:clear(0.0)
         return f
      end
      local allocDistf = function()
         return self:allocCartField(self.boundaryGrid,self.basis:numBasis(),
                                    {distf:lowerGhost(),distf:upperGhost()},distf:getMetaData())
      end
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
      self.ghostRng = localExtRange:intersect(self:getGhostRange(global, globalExt))
      -- Decompose ghost region into threads.
      self.ghostRangeDecomp = LinearDecomp.LinearDecompRange{range=self.ghostRng, numSplit=self.grid:numSharedProcs()}
      self.tId              = self.grid:subGridSharedId() -- Local thread ID.

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
      self.calcBoundaryFluxRateFunc = function(dtIn)
         -- Compute boundary flux rate ~ (fGhost_new - fGhost_old)/dt.
         self.boundaryFluxRate:combine( 1.0/dtIn, self.boundaryFluxFields[1],
                                       -1.0/dtIn, self.boundaryFluxFieldPrev)
         self.boundaryFluxFieldPrev:copy(self.boundaryFluxFields[1])
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
   end
   return self.diagnostics
end

-- These are needed to recycle the VlasovDiagnostics with VlasovBasicBC.
function VlasovBasicBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                                 self.boundaryFluxRate, self.boundaryFluxRate} end
function VlasovBasicBC:getFlucF() return self.boundaryFluxRate end

function VlasovBasicBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   local fIn = mySpecies:rkStepperFields()[outIdx]
   self.bcSolver:advance(tCurr, {fIn}, {fIn})
end

function VlasovBasicBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return VlasovBasicBC
