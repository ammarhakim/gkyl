-- Gkyl ------------------------------------------------------------------------
--
-- Apply boundary conditions specified as BoundaryCondition objects.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local CartDecomp = require "Lib.CartDecomp"
local CartFieldBinOp = require "Updater.CartFieldBinOp"
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"
local DataStruct = require "DataStruct"
local DistFuncMomentCalc = require "Updater.DistFuncMomentCalc"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Mpi = require "Comm.Mpi"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase  = require "Updater.Base"

-- Boundary condition updater.
local Bc = Proto(UpdaterBase)

local dirlabel = {"X", "Y", "Z"}

function Bc:init(tbl)
   Bc.super.init(self, tbl) -- Setup base object.

   self._isFirst = true -- Will be reset first time _advance() is called.

   self._grid = assert(tbl.onGrid, "Updater.Bc: Must specify grid to use with 'onGrid''")
   self._dir  = assert(tbl.dir, "Updater.Bc: Must specify direction to apply BCs with 'dir'")
   self._dirlabel = dirlabel[self._dir]

   self._edge = assert(
      tbl.edge, "Updater.Bc: Must specify edge to apply BCs with 'edge'. Must be one of 'upper', or 'lower'.")
   if self._edge ~= "lower" and self._edge ~= "upper" then
      error("Updater.Bc: 'edge' must be one of 'lower' or 'upper'. Was " .. self._edge .. " instead")
   end
   self._bcList = assert(
      tbl.boundaryConditions, "Updater.Bc: Must specify boundary conditions to apply with 'boundaryConditions'")

   self._skinLoop = tbl.skinLoop and tbl.skinLoop or "pointwise"
   if self._skinLoop == "flip" then
      self._cdim = assert(
	 tbl.cdim,
	 "Updater.Bc: Must specify configuration space dimensions to apply with 'cdim'")
      self._vdir = assert(
	 tbl.vdir,
	 "Updater.Bc: Must specify velocity direction to flip with 'vdir'")
   end

   self._ghostRangeDecomp = nil -- Will be constructed on first call to advance.


   -- For diagnostics: create reduced boundary grid with 1 cell in dimension of self._dir.
   if self._grid:isShared() then 
      -- Shared memory implementation needs more work...
      self._boundaryGrid = self._grid
   else
      local reducedLower, reducedUpper, reducedNumCells, reducedCuts = {}, {}, {}, {}
      for d = 1, self._grid:ndim() do
         if d==self._dir then
            table.insert(reducedLower, -self._grid:dx(d)/2)
            table.insert(reducedUpper, self._grid:dx(d)/2)
            table.insert(reducedNumCells, 1)
            table.insert(reducedCuts, 1)
         else
            table.insert(reducedLower, self._grid:lower(d))
            table.insert(reducedUpper, self._grid:upper(d))
            table.insert(reducedNumCells, self._grid:numCells(d))
            table.insert(reducedCuts, self._grid:cuts(d))
         end
      end
      local commSet = self._grid:commSet()
      local worldComm = commSet.comm
      local nodeComm = commSet.nodeComm
      local nodeRank = Mpi.Comm_rank(nodeComm)
      local dirRank = nodeRank
      local cutsX = self._grid:cuts(1) or 1
      local cutsY = self._grid:cuts(2) or 1
      local cutsZ = self._grid:cuts(3) or 1
      local writeRank = -1
      if self._dir == 1 then 
         dirRank = nodeRank % (cutsX*cutsY) % cutsX
      elseif self._dir == 2 then 
         dirRank = math.floor(nodeRank / cutsX) % cutsY
      elseif self._dir == 3 then
         dirRank = math.floor(nodeRank/cutsX/cutsY)
      end
      self._splitComm = Mpi.Comm_split(worldComm, dirRank, nodeRank)
      -- set up which ranks to write from 
      if self._edge == "lower" and dirRank == 0 then 
         writeRank = nodeRank
      elseif self._edge == "upper" and dirRank == self._grid:cuts(self._dir)-1 then
         writeRank = nodeRank
      end
      self.writeRank = writeRank
      
      local reducedDecomp = CartDecomp.CartProd {
         comm = self._splitComm,
         writeRank = writeRank,
         cuts = reducedCuts,
         useShared = self._grid:isShared(),
      }

      self._boundaryGrid = Grid.RectCart {
         lower = reducedLower,
         upper = reducedUpper,
         cells = reducedNumCells,
         decomposition = reducedDecomp,
      }
   end

   -- A function can be specied in the ghost cell layer to be used as
   -- a boundary condition. Additionally, a feedback function of fluid
   -- moments at the boundary can be set up.
   self._evaluateFn = tbl.evaluate
   if self._evaluateFn then
      self._basis = assert(tbl.basis, "Bc.init: Evaluate is currently implemented only for DG; 'basis' must be specified.")
      self._evolveFn = xsys.pickBool(tbl.evolveFn, false)
      self._feedback = xsys.pickBool(tbl.feedback, false)
      if self._feedback then
         assert(tbl.cdim == 1, "Bc.init: Feedback boundary condition is implemented only for 1X simulations.")
         local confBasis = assert(tbl.confBasis, "Bc.init: Must specify 'confBasis' when 'feedback' is true.")
         local numConfDims = confBasis:ndim()
         self.numConfBasis = confBasis:numBasis()
         local node = Lin.Vec(numConfDims)
         if self._edge == "lower" then node[1] = -1.0 else node[1] = 1.0 end
         self._confBasisEdge = Lin.Vec(self.numConfBasis)
         confBasis:evalBasis(node, self._confBasisEdge)
      end
   end
end

function Bc:getGhostRange(global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if self._edge == "lower" then
      uv[self._dir] = global:lower(self._dir)-1 -- For ghost cells on "left".
   else
      lv[self._dir] = global:upper(self._dir)+1 -- For ghost cells on "right".
   end
   return Range.Range(lv, uv)
end

function Bc:_advance(tCurr, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local dir, edge = self._dir, self._edge
   local vdir      = self._vdir
   local global    = qOut:globalRange()

   if self._isFirst then
      local globalExt     = qOut:globalExtRange()
      local localExtRange = qOut:localExtRange()
      self._ghostRng      = localExtRange:intersect(
   	 self:getGhostRange(global, globalExt)) -- Range spanning ghost cells.
      -- Decompose ghost region into threads.
      self._ghostRangeDecomp = LinearDecomp.LinearDecompRange {
   	 range = self._ghostRng, numSplit = grid:numSharedProcs() }

      -- Get the unction onto the boundary grid
      if self._evaluateFn then
         self._ghostFld = DataStruct.Field {
           onGrid = self._boundaryGrid,
           numComponents = qOut:numComponents(),
           ghost = {1,1},
           metaData = qOut:getMetaData(),
         }
      end
   end

   if self._evaluateFn and (self._isFirst or self._feedback) then
      local myEvaluateFn = self._evaluateFn
      if self._feedback then
         local M0, M1, M2 = inFld[1], inFld[2], inFld[3]
         local indexer = M0:genIndexer()
         local ptrM0, ptrM1, ptrM2 = M0:get(1), M1:get(1), M2:get(1)
         local edgeM0, edgeM1, edgeM2 = 0.0, 0.0, 0.0
         for idx in self._ghostRangeDecomp:rowMajorIter(tId) do
            idx[1] = self._edge == "lower" and global:lower(dir) or global:upper(dir)
            M0:fill(indexer(idx), ptrM0)
            M1:fill(indexer(idx), ptrM1)
            M2:fill(indexer(idx), ptrM2)
            for k = 1, self.numConfBasis do
               edgeM0 = edgeM0 + ptrM0[k]*self._confBasisEdge[k]
               edgeM1 = edgeM1 + ptrM1[k]*self._confBasisEdge[k]
               edgeM2 = edgeM2 + ptrM2[k]*self._confBasisEdge[k]
            end
            myEvaluateFn = function(t, z)
               return self._evaluateFn(t, z, edgeM0, edgeM1, edgeM2)
            end
            break -- Hack, this need just a conf space loop
         end
      end
      self._projectEvaluateFn = ProjectOnBasis {
         onGrid = self._boundaryGrid,
         basis = self._basis,
         evaluate = myEvaluateFn,
      }
   end
   
   if self._evaluateFn and (self._isFirst or self._evolveFn) then
      for idx in self._ghostRangeDecomp:rowMajorIter(tId) do
         self._projectEvaluateFn:advance(tCurr, {}, {self._ghostFld})
         break
      end
   end

   local tId = self._grid:subGridSharedId() -- Local thread ID.
   
   -- Get the in and out pointers
   local ptrOut, ptrIn = qOut:get(1), qOut:get(1)
   local indexerOut, indexerIn = qOut:genIndexer(), qOut:genIndexer()
   if self._evaluateFn then
      ptrIn = self._ghostFld:get(1)
      indexerIn = self._ghostFld:genIndexer()
   end

   local idxIn = Lin.IntVec(self._grid:ndim()) -- Prealloc this
   for idxOut in self._ghostRangeDecomp:rowMajorIter(tId) do 
      qOut:fill(indexerOut(idxOut), ptrOut)

      -- Copy out index into in index
      idxOut:copyInto(idxIn)
      -- Reverse the in velocity index if needed
      if self._skinLoop == "flip" then
         idxIn[vdir] = global:upper(vdir) + 1 - idxIn[vdir]
      end
      if self._evaluateFn then
         idxIn[dir] = 1 -- the boundaryGrid has only 1 cell in dir
         self._ghostFld:fill(indexerIn(idxIn), ptrIn)
      else
         idxIn[dir] = edge == "lower" and global:lower(dir) or global:upper(dir)
         qOut:fill(indexerIn(idxIn), ptrIn)
      end

      for _, bc in ipairs(self._bcList) do
         -- Apply the 'bc' function. This can represent many boundary
         -- condition types ranging from a simple copy or a reflection
         -- with the sign flit to QM based electron emission model.
         bc(dir, tCurr, idxIn, ptrIn, ptrOut)
      end
   end

   self._isFirst = false
end

function Bc:storeBoundaryFlux(tCurr, rkIdx, qOut)
   if self._isFirst then
      -- make enough boundary field copies for 4 RK stages
      self._boundaryFluxFields = {
         DataStruct.Field {
           onGrid        = self._boundaryGrid,
           numComponents = qOut:numComponents(),
           ghost         = {1,1},
           metaData = qOut:getMetaData(),
         },
         DataStruct.Field {
           onGrid        = self._boundaryGrid,
           numComponents = qOut:numComponents(),
           ghost         = {1,1},
           metaData = qOut:getMetaData(),
         },
         DataStruct.Field {
           onGrid        = self._boundaryGrid,
           numComponents = qOut:numComponents(),
           ghost         = {1,1},
           metaData = qOut:getMetaData(),
         },
         DataStruct.Field {
           onGrid        = self._boundaryGrid,
           numComponents = qOut:numComponents(),
           ghost         = {1,1},
           metaData = qOut:getMetaData(),
         },
      }
      self._boundaryFluxRate = DataStruct.Field {
         onGrid        = self._boundaryGrid,
         numComponents = qOut:numComponents(),
         ghost         = {1,1},
         metaData = qOut:getMetaData(),
      }
      self._boundaryFluxFieldPrev = DataStruct.Field {
         onGrid        = self._boundaryGrid,
         numComponents = qOut:numComponents(),
         ghost         = {1,1},
         metaData = qOut:getMetaData(),
      }
      self._boundaryFluxFieldPrev:clear(0.0)
      self._boundaryPtr = {self._boundaryFluxFields[1]:get(1), 
                           self._boundaryFluxFields[2]:get(1),
                           self._boundaryFluxFields[3]:get(1),
                           self._boundaryFluxFields[4]:get(1)
                          }
      self._boundaryIdxr = self._boundaryFluxFields[1]:genIndexer()

      local global = qOut:globalRange()
      local globalExt = qOut:globalExtRange()
      local localExtRange = qOut:localExtRange()
      self._ghostRng = localExtRange:intersect(
   	 self:getGhostRange(global, globalExt)) -- Range spanning ghost cells.
      -- Decompose ghost region into threads.
      self._ghostRangeDecomp = LinearDecomp.LinearDecompRange {
   	 range = self._ghostRng, numSplit = self._grid:numSharedProcs() }
   end

   local ptrOut = qOut:get(1) -- Get pointers to (re)use inside inner loop.
   local indexer = qOut:genIndexer()

   local tId = self._grid:subGridSharedId() -- Local thread ID.
   
   local idxIn = Lin.IntVec(self._grid:ndim())
   for idxOut in self._ghostRangeDecomp:rowMajorIter(tId) do
      idxOut:copyInto(idxIn)

      qOut:fill(indexer(idxOut), ptrOut) 

      -- Before operating on ghosts, store ghost values for later flux diagnostics
      idxIn[self._dir] = 1
      self._boundaryFluxFields[rkIdx]:fill(self._boundaryIdxr(idxIn), self._boundaryPtr[rkIdx])
      for c = 1, qOut:numComponents() do
         self._boundaryPtr[rkIdx][c] = ptrOut[c]
      end
   end
end

function Bc:evalOnConfBoundary(inFld)
   if self._isFirst then
      self._confBoundaryField = DataStruct.Field {
           onGrid        = self._confBoundaryGrid,
           numComponents = inFld:numComponents(),
           ghost         = {1,1},
           metaData = inFld:getMetaData(),
         }
      self._confBoundaryFieldPtr = self._confBoundaryField:get(1)
      self._confBoundaryIdxr = self._confBoundaryField:genIndexer()
   end
   local global = inFld:globalRange()
   local dir, edge = self._dir, self._edge
   
   local localRange = self._confBoundaryField:localRange()
   local localRangeDecomp = LinearDecomp.LinearDecompRange {
      range = localRange, numSplit = self._grid:numSharedProcs() }

   local inFldPtr = inFld:get(1)
   local indexer = inFld:genIndexer()

   local tId = self._grid:subGridSharedId() -- Local thread ID.
   local idxIn = Lin.IntVec(self._grid:ndim())
   for idxOut in localRangeDecomp:rowMajorIter(tId) do
      idxOut:copyInto(idxIn)
      if edge == "lower" then
         idxIn[dir] = global:lower(dir)-inFld:lowerGhost() 
      else
         idxIn[dir] = global:upper(dir)+inFld:upperGhost()
      end
      
      inFld:fill(indexer(idxIn), inFldPtr) 
      self._confBoundaryField:fill(self._confBoundaryIdxr(idxOut), self._confBoundaryFieldPtr)
       
      for c = 1, inFld:numComponents() do
         self._confBoundaryFieldPtr[c] = inFldPtr[c]
      end
   end

   return self._confBoundaryField
end

function Bc:getDir()
   return self._dir
end

function Bc:getEdge()
   return self._edge
end

function Bc:label()
   return "Flux"..self._dirlabel..self._edge
end

function Bc:getBoundaryFluxFields()
   return self._boundaryFluxFields
end

function Bc:getBoundaryFluxRate()
   return self._boundaryFluxRate
end

function Bc:getBoundaryFluxFieldPrev()
   return self._boundaryFluxFieldPrev
end

function Bc:getBoundaryGrid()
   return self._boundaryGrid
end

function Bc:getConfBoundaryGrid()
   return self._confBoundaryGrid
end

function Bc:initBcDiagnostics(cdim)
   if cdim == self._grid:ndim() then
      self._confBoundaryGrid = self._boundaryGrid
   else 
      -- Create reduced boundary config-space grid with 1 cell in dimension of self._dir.
      local reducedLower, reducedUpper, reducedNumCells, reducedCuts = {}, {}, {}, {}
      for d = 1, cdim do
         if d==self._dir then
            table.insert(reducedLower, -self._grid:dx(d)/2)
            table.insert(reducedUpper, self._grid:dx(d)/2)
            table.insert(reducedNumCells, 1)
            table.insert(reducedCuts, 1)
         else
            table.insert(reducedLower, self._grid:lower(d))
            table.insert(reducedUpper, self._grid:upper(d))
            table.insert(reducedNumCells, self._grid:numCells(d))
            table.insert(reducedCuts, self._grid:cuts(d))
         end
      end
      local reducedDecomp = CartDecomp.CartProd {
         comm = self._splitComm,
         writeRank = self.writeRank,
         cuts = reducedCuts,
         useShared = self._grid:isShared(),
      }

      self._confBoundaryGrid = Grid.RectCart {
         lower = reducedLower,
         upper = reducedUpper,
         cells = reducedNumCells,
         decomposition = reducedDecomp,
      }
   end
end

return Bc
