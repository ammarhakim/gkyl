-- Gkyl ------------------------------------------------------------------------
--
-- Electron reflection BC based on the Bronold & Fehske model.
--   F.X. Bronold, H. Fehske, Absorption of an electron by a dielectric wall, Phys. Rev. Lett. 115 (22) (2015) 225001.
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

local BnFReflectionBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function BnFReflectionBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function BnFReflectionBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.diagnostics = tbl.diagnostics or {}
   self.saveFlux    = tbl.saveFlux or false

   self.electronAffinity = tbl.electronAffinity
   self.effectiveMass    = tbl.effectiveMass   
   self.elemCharge       = tbl.elemCharge      
   self.electronMass     = tbl.electronMass    
   self.ignore           = tbl.ignore

   if tbl.diagnostics then
     if #tbl.diagnostics>0 then self.saveFlux = true end
   end
end

function BnFReflectionBC:setName(nm) self.name = self.speciesName.."_"..nm end

function BnFReflectionBC:bcBnFReflectionFunc(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "flip".
   local numBasis = self.basis:numBasis()
   local velIdx = Lin.IntVec(self.ndim)
   velIdx[1] = 1
   for d = 1, self.vdim do
      velIdx[d + 1] = idxIn[self.cdim + d]
   end
   local exIdxr = self.externalBCFunction:genIndexer()
   local externalBCFunction = self.externalBCFunction:get(exIdxr(velIdx))
   if velIdx[1] ~= 0 and velIdx[1] ~= self.grid:numCells(2) + 1 then
      for i = 1, numBasis do
        fOut[i] = 0
         for j = 1, numBasis do
            fOut[i] = fOut[i] + fIn[j]*externalBCFunction[(i - 1)*numBasis + j]
         end
      end
   end
   return fOut
end

function BnFReflectionBC:createSolver(mySpecies)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   local lower, upper, cells = {}, {}, {}
   local GridConstructor     = Grid.RectCart
   local coordinateMap       = {}
   if mySpecies.coordinateMap then
      for d = 1, self.cdim do table.insert(coordinateMap, function (z) return z end) end
      for d = 1, self.vdim do table.insert(coordinateMap, mySpecies.coordinateMap[d]) end
      GridConstructor = Grid.NonUniformRectCart
   end
   for d = 1, self.cdim do lower[d], upper[d], cells[d] = 0, 1, 1 end
   for d = 1, self.vdim do
      lower[d + self.cdim] = mySpecies.lower[d]
      upper[d + self.cdim] = mySpecies.upper[d]
      cells[d + self.cdim] = self.grid:numCells(d + self.cdim)
   end
   local grid = GridConstructor {
      lower = lower,  periodicDirs = {},
      upper = upper,  mappings     = coordinateMap,
      cells = cells,
   }
   self.externalBCFunction = DataStruct.Field {
      onGrid        = grid,
      numComponents = self.basis:numBasis()*self.basis:numBasis(),
      metaData      = {polyOrder = self.basis:polyOrder(),  charge = self.charge,
                       basisType = self.basis:id(),         mass   = self.mass,  },
   }
   -- Calculate Bronold and Fehske reflection function coefficients
   local evaluateBronold = Updater.EvaluateBronoldFehskeBC {
      onGrid = grid,        electronAffinity = self.electronAffinity,
      basis  = self.basis,  effectiveMass    = self.effectiveMass,
      cdim   = self.cdim,   elemCharge       = self.elemCharge,
      vdim   = self.vdim,   me               = self.electronMass,
      ignore = self.ignore,
   }
   evaluateBronold:advance(0.0, {}, {self.externalBCFunction})

   local bcFunc   = function(...) return self:bcBnFReflectionFunc(...) end
   local skinType = "flip"

   local vdir = self.bcDir+self.cdim
   self.bcSolver = Updater.Bc{
      onGrid = self.grid,   edge               = self.bcEdge,  
      cdim   = self.cdim,   skinLoop           = skinType,
      dir    = self.bcDir,  boundaryConditions = {bcFunc},   
      vdir   = vdir,      
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

function BnFReflectionBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function BnFReflectionBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function BnFReflectionBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function BnFReflectionBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function BnFReflectionBC:createDiagnostics(mySpecies, field)
   -- Create BC diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = VlasovDiags()}
      self.diagnostics:fullInit(mySpecies, field, self)
   end
   return self.diagnostics
end

-- These are needed to recycle the VlasovDiagnostics with BnFReflectionBC.
function BnFReflectionBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                                 self.boundaryFluxRate, self.boundaryFluxRate} end
function BnFReflectionBC:getFlucF() return self.boundaryFluxRate end

function BnFReflectionBC:advance(tCurr, species, inFlds)
   self.bcSolver:advance(tCurr, {}, inFlds)
end

function BnFReflectionBC:getBoundaryFluxFields()
   return self.boundaryFluxFields
end

return BnFReflectionBC
