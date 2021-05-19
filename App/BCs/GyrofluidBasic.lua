-- Gkyl ------------------------------------------------------------------------
--
-- Basic boundary condition for a gyrofluid species, i.e. those that can be
-- applied with Updater/Bc.lua.
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
local GyrofluidDiags = require "App.Diagnostics.GyrofluidDiagnostics"

local GyrofluidBasicBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function GyrofluidBasicBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GyrofluidBasicBC:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.bcKind      = assert(tbl.kind, "GyrofluidBasicBC: must specify the type of BC in 'kind'.")
   self.diagnostics = tbl.diagnostics or {}
   self.saveFlux    = tbl.saveFlux or false
end

function GyrofluidBasicBC:setName(nm) self.name = self.speciesName.."_"..nm end

function GyrofluidBasicBC:bcCopy(dir, tm, idxIn, fIn, fOut)
   for i = 1, self.nMoments*self.basis:numBasis() do fOut[i] = fIn[i] end
end

function GyrofluidBasicBC:bcAbsorb(dir, tm, idxIn, fIn, fOut)
   -- The idea is that by setting the plasma quantities to zero in the
   -- ghost cell nothing is transported into the domain, and whatever is transported
   -- out is lost. We can't set them to exactly zero or else the sound speed
   -- and drift velocity would diverge, so we set them to something small.
   local numB = self.basis:numBasis()
   for i = 1, numB do fOut[0*numB+i] = 1.e-10*fIn[0*numB+i] end   -- Mass density.
   for i = 1, numB do fOut[1*numB+i] = 0. end                     -- Momentum density.
   for i = 1, numB do fOut[2*numB+i] = 1.e-10*fIn[2*numB+i] end   -- Energy density.
   for i = 1, numB do fOut[3*numB+i] = 1.e-10*fIn[3*numB+i] end   -- Perpendicular pressure (divided by B).
end

function GyrofluidBasicBC:getGhostRange(global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if self.bcEdge == "lower" then
      uv[self.bcDir] = global:lower(self.bcDir)-1   -- For ghost cells on "left".
   else
      lv[self.bcDir] = global:upper(self.bcDir)+1   -- For ghost cells on "right".
   end
   return Range.Range(lv, uv)
end

function GyrofluidBasicBC:createSolver(mySpecies)
   local bcFunc, skinType
   if self.bcKind == "copy" then
      bcFunc   = function(...) return self:bcCopy(...) end
      skinType = "pointwise"
   elseif self.bcKind == "absorb" then
      bcFunc   = function(...) return self:bcAbsorb(...) end
      skinType = "pointwise"
   end

   self.bcSolver = Updater.Bc {
      onGrid             = self.grid,
      cdim               = self.grid:ndim(),
      dir                = self.bcDir,
      edge               = self.bcEdge,
      boundaryConditions = {bcFunc},
      skinLoop           = skinType,
   }

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then
      -- Create reduced boundary grid with 1 cell in dimension of self._dir.
      if self.grid:isShared() then
         assert(false, "GyrofluidBasicBC: shared memory implementation of boundary flux diagnostics not ready.")
      else
         local reducedLower, reducedUpper, reducedNumCells, reducedCuts = {}, {}, {}, {}
         for d = 1, self.grid:ndim() do
            if d==self.bcDir then
               table.insert(reducedLower, -self.grid:dx(d)/2.)
               table.insert(reducedUpper,  self.grid:dx(d)/2.)
               table.insert(reducedNumCells, 1)
               table.insert(reducedCuts, 1)
            else
               table.insert(reducedLower, self.grid:lower(d))
               table.insert(reducedUpper, self.grid:upper(d))
               table.insert(reducedNumCells, self.grid:numCells(d))
               table.insert(reducedCuts, self.grid:cuts(d))
            end
         end
         local commSet  = self.grid:commSet()
         local worldComm, nodeComm = commSet.comm, commSet.nodeComm
         local nodeRank = Mpi.Comm_rank(nodeComm)
         local dirRank  = nodeRank
         local cuts     = {}
         for d=1,3 do cuts[d] = self.grid:cuts(d) or 1 end
         local writeRank = -1
         if self.bcDir == 1 then
            dirRank = nodeRank % (cuts[1]*cuts[2]) % cuts[1]
         elseif self.bcDir == 2 then
            dirRank = math.floor(nodeRank/cuts[1]) % cuts[2]
         elseif self.bcDir == 3 then
            dirRank = math.floor(nodeRank/cuts[1]/cuts[2])
         end
         self._splitComm = Mpi.Comm_split(worldComm, dirRank, nodeRank)
         -- Set up which ranks to write from.
         if self.bcEdge == "lower" and dirRank == 0 then
            writeRank = nodeRank
         elseif self.bcEdge == "upper" and dirRank == self.grid:cuts(self.bcDir)-1 then
            writeRank = nodeRank
         end
         self.writeRank = writeRank
         local reducedDecomp = CartDecomp.CartProd {
            comm      = self._splitComm,  cuts      = reducedCuts,
            writeRank = writeRank,        useShared = self.grid:isShared(),
         }
         self.boundaryGrid = Grid.RectCart {
            lower = reducedLower,  cells = reducedNumCells,
            upper = reducedUpper,  decomposition = reducedDecomp,
         }
      end

      local moms = mySpecies:getMoments()
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
      self.allocMoment = function(self)
         return self:allocCartField(self.boundaryGrid, self.basis:numBasis(), 
                                    {moms:lowerGhost(),moms:upperGhost()}, moms:getMetaData())
      end
      self.allocVectorMoment = function(self, dim)
         return self:allocCartField(self.boundaryGrid, dim*self.basis:numBasis(), 
                                    {moms:lowerGhost(),moms:upperGhost()}, moms:getMetaData())
      end

      -- Allocate fields needed.
      self.boundaryFluxFields = {}  -- Fluxes through the boundary, into ghost region, from each RK stage.
      self.boundaryPtr        = {}
      self.momInIdxr          = moms:genIndexer()
      for i = 1, #mySpecies:rkStepperFields() do
         self.boundaryFluxFields[i] = self:allocMoment()
         self.boundaryPtr[i]        = self.boundaryFluxFields[i]:get(1)
      end
      self.boundaryFluxRate      = self:allocMoment()
      self.boundaryFluxFieldPrev = self:allocMoment()
      self.boundaryIdxr          = self.boundaryFluxFields[1]:genIndexer()

      self.idxOut = Lin.IntVec(self.grid:ndim())

      -- Create the range needed to loop over ghosts.
      local global, globalExt, localExtRange = moms:globalRange(), moms:globalExtRange(), moms:localExtRange()
      self.ghostRng = localExtRange:intersect(self:getGhostRange(global, globalExt))
      -- Decompose ghost region into threads.
      self.ghostRangeDecomp = LinearDecomp.LinearDecompRange{range=self.ghostRng, numSplit=self.grid:numSharedProcs()}
      self.tId              = self.grid:subGridSharedId() -- Local thread ID.

      self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut)
         local ptrOut = qOut:get(1)
         for idx in self.ghostRangeDecomp:rowMajorIter(self.tId) do
            idx:copyInto(self.idxOut)
            qOut:fill(self.momInIdxr(idx), ptrOut)

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
      self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut) end
      self.copyBoundaryFluxFieldFunc = function(inIdx, outIdx) end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...) end
      self.calcBoundaryFluxRateFunc = function(dtIn) end
   end
end

function GyrofluidBasicBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end
function GyrofluidBasicBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function GyrofluidBasicBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function GyrofluidBasicBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function GyrofluidBasicBC:createDiagnostics(mySpecies)
   -- Create BC diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = GyrofluidDiags()}
      self.diagnostics:fullInit(mySpecies, self)
   end
   return self.diagnostics
end

function GyrofluidBasicBC:getNoJacMoments() return self.boundaryFluxRate end  -- Used by diagnostics.

function GyrofluidBasicBC:advance(tCurr, species, inFlds)
   self.bcSolver:advance(tCurr, {}, inFlds)
end

function GyrofluidBasicBC:getBoundaryFluxFields()
   return self.boundaryFluxFields
end

return GyrofluidBasicBC
