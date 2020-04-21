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
local CartDecomp   = require "Lib.CartDecomp"
local DataStruct   = require "DataStruct"
local Grid         = require "Grid"
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto        = require "Lib.Proto"
local Range        = require "Lib.Range"
local UpdaterBase  = require "Updater.Base"
local DistFuncMomentCalc = require "Updater.DistFuncMomentCalc"
local CartFieldBinOp = require "Updater.CartFieldBinOp"
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"
local Mpi = require "Comm.Mpi"

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
   if self._skinLoop == "flip" or self._skinLoop == "integrate" then
      self._cdim = assert(
	 tbl.cdim,
	 "Updater.Bc: Must specify configuration space dimensions to apply with 'cdim'")
   end
   if self._skinLoop == "flip" then
      self._vdir = assert(
	 tbl.vdir,
	 "Updater.Bc: Must specify velocity direction to flip with 'vdir'")
   end
   if self._skinLoop == "integrate" then
      self._vdim = assert(
	 tbl.vdim,
	 "Updater.Bc: Must specify velocity space dimensions to apply with 'vdim'")
      self._numComps = assert(
	 tbl.numComps,
	 "Updater.Bc: Must specify the number of components of the field using 'numComps'")
   end

   self._ghostRangeDecomp = nil -- Will be constructed on first call to advance.

   self.idxS = Lin.IntVec(self._grid:ndim()) -- Prealloc this.

   self.hasExtFld = xsys.pickBool(tbl.hasExtFld, false)

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
      if self._dir == 1 then 
         dirRank = nodeRank % (cutsX*cutsY) % cutsX
      elseif self._dir == 2 then 
         dirRank = math.floor(nodeRank / cutsX) % cutsY
      elseif self._dir == 3 then
         dirRank = math.floor(nodeRank/cutsX/cutsY)
      end
      self._splitComm = Mpi.Comm_split(worldComm, dirRank, nodeRank)
      
      local reducedDecomp = CartDecomp.CartProd {
         comm = self._splitComm,
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
   local qIn  = inFld[1]

   local dir, edge = self._dir, self._edge
   local vdir      = self._vdir
   local global    = qOut:globalRange()

   if self._isFirst then
      local globalExt     = qOut:globalExtRange()
      local localExtRange = qOut:localExtRange()
      self._ghostRng      = localExtRange:intersect(
   	 self:getGhostRange(global, globalExt)) -- Range spanning ghost cells.
      if self._skinLoop == "integrate" then
   	 self._skin = qOut:localRange():selectLast(self._vdim)
      end
      -- Decompose ghost region into threads.
      self._ghostRangeDecomp = LinearDecomp.LinearDecompRange {
   	 range = self._ghostRng, numSplit = grid:numSharedProcs() }
   end

   local qG, qS = qOut:get(1), qOut:get(1) -- Get pointers to (re)use inside inner loop [G: Ghost, S: Skin].
   if self.hasExtFld then qS = qIn:get(1) end
   self.idxS = Lin.IntVec(grid:ndim()) -- Prealloc this.
   local indexer = qOut:genIndexer()

   local tId = self._grid:subGridSharedId() -- Local thread ID.
   for idxG in self._ghostRangeDecomp:rowMajorIter(tId) do -- Loop, applying BCs.
      idxG:copyInto(self.idxS)

      qOut:fill(indexer(idxG), qG) 

      -- If an in-field is specified the same indexes are used (gS
      -- points to the ghost layer of the in-field); otherwise, move
      -- the ghost index to point into the skin layer.
      if not self.hasExtFld then
	 self.idxS[dir] = edge == "lower" and global:lower(dir) or global:upper(dir)
	 if self._skinLoop == "flip" then
	    self.idxS[vdir] = global:upper(vdir) + 1 - self.idxS[vdir]
	 end
      end

      if self._skinLoop == "integrate" then 
   	 for c = 1, self._numComponents do qG[c] = 0 end

         for idx in self._skin:rowMajorIter() do
   	    for d = 1, self._vdim do self.idxS[self._cdim + d] = idx[d] end
   	    qOut:fill(indexer(self.idxS), qS)
            for _, bc in ipairs(self._bcList) do
               bc(dir, tCurr, self.idxS, qS, qG, self._bcList)
            end
         end
      else
	 if not self.hasExtFld then
	    qOut:fill(indexer(self.idxS), qS)
	 else
	    qIn:fill(indexer(self.idxS), qS)
	 end
         for _, bc in ipairs(self._bcList) do
            bc(dir, tCurr, self.idxS, qS, qG) -- TODO: PASS COORDINATES.
         end
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

   local qG = qOut:get(1) -- Get pointers to (re)use inside inner loop [G: Ghost].
   local indexer = qOut:genIndexer()

   local tId = self._grid:subGridSharedId() -- Local thread ID.
   for idxG in self._ghostRangeDecomp:rowMajorIter(tId) do -- Loop, applying BCs.
      idxG:copyInto(self.idxS)

      qOut:fill(indexer(idxG), qG) 

      -- Before operating on ghosts, store ghost values for later flux diagnostics
      self.idxS[self._dir] = 1
      self._boundaryFluxFields[rkIdx]:fill(self._boundaryIdxr(self.idxS), self._boundaryPtr[rkIdx])
      for c = 1, qOut:numComponents() do
         self._boundaryPtr[rkIdx][c] = qG[c]
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
   for idxOut in localRangeDecomp:rowMajorIter(tId) do
      idxOut:copyInto(self.idxS)
      if edge == "lower" then
         self.idxS[dir] = global:lower(dir)-inFld:lowerGhost() 
      else
         self.idxS[dir] = global:upper(dir)+inFld:upperGhost()
      end
      
      inFld:fill(indexer(self.idxS), inFldPtr) 
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
