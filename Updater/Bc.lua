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
local CartDecomp     = require "Lib.CartDecomp"
local CartFieldBinOp = require "Updater.CartFieldBinOp"
local DataStruct     = require "DataStruct"
local Grid           = require "Grid"
local Lin            = require "Lib.Linalg"
local LinearDecomp   = require "Lib.LinearDecomp"
local Mpi            = require "Comm.Mpi"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local Proto          = require "Lib.Proto"
local Range          = require "Lib.Range"
local UpdaterBase    = require "Updater.Base"
local MomDecl        = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"

-- Boundary condition updater.
local Bc = Proto(UpdaterBase)

local dirlabel = {"X", "Y", "Z"}

local function createFieldFromField(grid, fld, ghostCells)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = fld:numComponents(),
      ghost         = ghostCells,
      metaData      = fld:getMetaData(),
   }
   fld:clear(0.0)
   return fld
end

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
      self._cDim = assert(
	 tbl.cdim,
	 "Updater.Bc: Must specify configuration space dimensions to apply with 'cdim'")
      self._vdir = assert(
	 tbl.vdir,
	 "Updater.Bc: Must specify velocity direction to flip with 'vdir'")
   end

   self._ghostRangeDecomp = nil -- Will be constructed on first call to advance.

   self._idxIn  = Lin.IntVec(self._grid:ndim())
   self._idxOut = Lin.IntVec(self._grid:ndim())
   self._xcIn   = Lin.Vec(self._grid:ndim())
   self._xcOut  = Lin.Vec(self._grid:ndim())

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
      local commSet   = self._grid:commSet()
      local worldComm = commSet.comm
      local nodeComm  = commSet.nodeComm
      local nodeRank  = Mpi.Comm_rank(nodeComm)
      local dirRank   = nodeRank
      local cuts      = {}
      for d=1,3 do cuts[d] = self._grid:cuts(d) or 1 end
      local writeRank = -1
      if self._dir == 1 then 
         dirRank = nodeRank % (cuts[1]*cuts[2]) % cuts[1]
      elseif self._dir == 2 then 
         dirRank = math.floor(nodeRank/cuts[1]) % cuts[2]
      elseif self._dir == 3 then
         dirRank = math.floor(nodeRank/cuts[1]/cuts[2])
      end
      self._splitComm = Mpi.Comm_split(worldComm, dirRank, nodeRank)
      -- Set up which ranks to write from.
      if self._edge == "lower" and dirRank == 0 then 
         writeRank = nodeRank
      elseif self._edge == "upper" and dirRank == self._grid:cuts(self._dir)-1 then
         writeRank = nodeRank
      end
      self.writeRank = writeRank
      
      local reducedDecomp = CartDecomp.CartProd {
         comm      = self._splitComm,
         writeRank = writeRank,
         cuts      = reducedCuts,
         useShared = self._grid:isShared(),
      }

      self._boundaryGrid = Grid.RectCart {
         lower = reducedLower,
         upper = reducedUpper,
         cells = reducedNumCells,
         decomposition = reducedDecomp,
      }

      if tbl.inField then
         local inFld = tbl.inField
         -- make enough boundary field copies for 4 RK stages
         self._boundaryFluxFields = {createFieldFromField(self._boundaryGrid, inFld, {1,1}),
                                     createFieldFromField(self._boundaryGrid, inFld, {1,1}),
                                     createFieldFromField(self._boundaryGrid, inFld, {1,1}),
                                     createFieldFromField(self._boundaryGrid, inFld, {1,1}),}
         self._boundaryFluxRate = createFieldFromField(self._boundaryGrid, inFld, {1,1})
         self._boundaryFluxFieldPrev = createFieldFromField(self._boundaryGrid, inFld, {1,1})
         self._boundaryFluxFieldPrev:clear(0.0)
         self._boundaryPtr = {self._boundaryFluxFields[1]:get(1), 
                              self._boundaryFluxFields[2]:get(1),
                              self._boundaryFluxFields[3]:get(1),
                              self._boundaryFluxFields[4]:get(1),}
         self._boundaryIdxr = self._boundaryFluxFields[1]:genIndexer()

         local global = inFld:globalRange()
         local globalExt = inFld:globalExtRange()
         local localExtRange = inFld:localExtRange()
         self._ghostRng = self._ghostRng or localExtRange:intersect( self:getGhostRange(global, globalExt) ) -- Range spanning ghost cells.
         -- Decompose ghost region into threads.
         self._ghostRangeDecomp = self._ghostRangeDecomp or LinearDecomp.LinearDecompRange {
      	    range = self._ghostRng, numSplit = self._grid:numSharedProcs() }
      end
   end

   -- A function can be specied in the ghost cell layer to be used as
   -- a boundary condition. Additionally, a feedback function of fluid
   -- moments at the boundary can be set up.
   self._evaluateFn = tbl.evaluate
   if self._evaluateFn then
      self._basis    = assert(tbl.basis, "Bc.init: Evaluate is currently implemented only for DG; 'basis' must be specified.")
      self._evolveFn = xsys.pickBool(tbl.evolveFn, false)
      self._feedback = xsys.pickBool(tbl.feedback, false)
      if self._feedback then
         assert(tbl.cdim == 1, "Bc.init: Feedback boundary condition is implemented only for 1X simulations.")
         local confBasis   = assert(tbl.confBasis, "Bc.init: Must specify 'confBasis' when 'feedback' is true.")
         local cDim        = confBasis:ndim()
         self.numConfBasis = confBasis:numBasis()
         local node        = Lin.Vec(cDim)
         if self._edge == "lower" then node[1] = -1.0 else node[1] = 1.0 end
         self._confBasisEdge = Lin.Vec(self.numConfBasis)
         confBasis:evalBasis(node, self._confBasisEdge)
         -- Need moments, and we need the momentum of the particles directed
         -- towards the boundary only. We will create moment calculators here.
         -- Instead of calling DistFuncMomentCalc over the entire confSpace grid
         -- we would ideally just compute the moment in the skin cell.
         local pDim = self._grid:ndim()
         local vDim = pDim - cDim
         self.idxP, self.xcP, self.dxP = Lin.IntVec(pDim), Lin.Vec(pDim), Lin.Vec(pDim)
         self.mom0  = Lin.Vec(self.numConfBasis)
         self.pMom1 = Lin.Vec(self.numConfBasis*vDim)
         self.mom2  = Lin.Vec(self.numConfBasis)
         -- Compute the offsets used to shorten the velocity range. For now assume that the zero
         -- along any velocity dimension is located at a cell boundary and not inside of a cell.
         local partialMomReg    = self._edge == "lower" and "N" or "P"
         local partialMomDirP   = cDim + self._dir
         self.partialMomDirExts = {0,0}
         local partialMomDirCells = self._grid:numCells(partialMomDirP)
         for d = 1,pDim do self.idxP[d]=1 end   -- Could be any cell in other directions.
         for idx=1,partialMomDirCells do
            self.idxP[partialMomDirP] = idx
            self._grid:setIndex(self.idxP)
            self._grid:cellCenter(self.xcP)
            if (partialMomReg == "P") and (self.xcP[partialMomDirP] > 0.0) then
               self.partialMomDirExts[1] = -(idx-1)
               break
            elseif (partialMomReg == "N") and (self.xcP[partialMomDirP] > 0.0) then
               self.partialMomDirExts[2] = -(partialMomDirCells-(idx-1))
               break
            end
         end
         self._mom0Calc  = MomDecl.selectMomCalc("M0",  self._basis:id(), cDim, vDim, self._basis:polyOrder(), false)
         self._pMom1Calc = MomDecl.selectMomCalc("M1i", self._basis:id(), cDim, vDim, self._basis:polyOrder(), false)
         self._mom2Calc  = MomDecl.selectMomCalc("M2",  self._basis:id(), cDim, vDim, self._basis:polyOrder(), false)
      end
   end
end

function Bc:getGhostRange(global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if self._edge == "lower" then
      uv[self._dir] = global:lower(self._dir)-1   -- For ghost cells on "left".
   else
      lv[self._dir] = global:upper(self._dir)+1   -- For ghost cells on "right".
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
      self._ghostRng      = self._ghostRng or localExtRange:intersect(
   	 self:getGhostRange(global, globalExt)) -- Range spanning ghost cells.
      -- Decompose ghost region into threads.
      self._ghostRangeDecomp = self._ghostRangeDecomp or LinearDecomp.LinearDecompRange {
   	 range = self._ghostRng, numSplit = grid:numSharedProcs() }

      -- Get the function onto the boundary grid.
      if self._evaluateFn then
         self._ghostFld = createFieldFromField(self._boundaryGrid, qOut, {1,1})
      end
   end

   if self._evaluateFn and (self._isFirst or self._feedback) then
      local myEvaluateFn = self._evaluateFn
      if self._feedback then
         local distf = inFld[1]

         -- Mini version of DistFuncMomentCalc to compute moments in the skin cell.
         local pDim, cDim = self._grid:ndim(), 1
         local vDim = pDim - cDim
         local phaseRange = distf:localRange()
         --for dir = 1, cDim do
         --   phaseRange = phaseRange:extendDir(dir, distf:lowerGhost(), distf:upperGhost())
         --end
         local phaseIndexer = distf:genIndexer()
         local distfItr = distf:get(1)
         -- Construct ranges for nested loops.
         local confSkinRange = self._edge == "lower" and phaseRange:selectFirst(cDim):lowerSkin(dir,1) or phaseRange:selectFirst(cDim):upperSkin(dir,1)
         local confSkinRangeDecomp = LinearDecomp.LinearDecompRange {
            range = confSkinRange, numSplit = grid:numSharedProcs() }
         local velRange = phaseRange:selectLast(vDim)
         local tId = grid:subGridSharedId()    -- Local thread ID.
         velRange = velRange:extendDir(dir,self.partialMomDirExts[1],self.partialMomDirExts[2])
         local phaseIndexer = distf:genIndexer()
         local edgeM0, edgeM1, edgeM2 = 0.0, 0.0, 0.0
         -- Outer loop is threaded and over configuration space.
         for cIdx in confSkinRangeDecomp:rowMajorIter(tId) do
            cIdx:copyInto(self.idxP)
            for k = 1, self.numConfBasis do
               self.mom0[k], self.pMom1[k], self.mom2[k] = 0.0, 0.0, 0.0
            end
            -- Inner loop is over velocity space: no threading to avoid race conditions.
            for vIdx in velRange:rowMajorIter() do
               for d = 1, vDim do self.idxP[cDim+d] = vIdx[d] end

               grid:setIndex(self.idxP)
               grid:cellCenter(self.xcP)
               grid:getDx(self.dxP)

               distf:fill(phaseIndexer(self.idxP), distfItr)

               self._mom0Calc(self.xcP:data(), self.dxP:data(), distfItr:data(),  self.mom0:data())
               self._pMom1Calc(self.xcP:data(), self.dxP:data(), distfItr:data(), self.pMom1:data())
               self._mom2Calc(self.xcP:data(), self.dxP:data(), distfItr:data(),  self.mom2:data())
            end
            -- This currently only works for 1x. In higher dimensions one gets an expansion
            -- in describing the variation parallel to the boundary. How to handle that expansion
            -- requires a better implementation of the feedback BC.
            for k = 1, self.numConfBasis do
               edgeM0 = edgeM0 + self.mom0[k]*self._confBasisEdge[k]
               edgeM1 = edgeM1 + self.pMom1[k]*self._confBasisEdge[k]
               edgeM2 = edgeM2 + self.mom2[k]*self._confBasisEdge[k]
            end

            myEvaluateFn = function(t, z)
               return self._evaluateFn(t, z, edgeM0, edgeM1, edgeM2)
            end
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

   for idxOut in self._ghostRangeDecomp:rowMajorIter(tId) do 
      qOut:fill(indexerOut(idxOut), ptrOut)

      -- Copy out index into in index
      idxOut:copyInto(self._idxIn)
      -- Reverse the in velocity index if needed
      if self._skinLoop == "flip" then
         self._idxIn[vdir] = global:upper(vdir) + 1 - self._idxIn[vdir]
      end
      if self._evaluateFn then
         self._idxIn[dir] = 1 -- the boundaryGrid has only 1 cell in dir
         self._ghostFld:fill(indexerIn(self._idxIn), ptrIn)
      else
         self._idxIn[dir] = edge == "lower" and global:lower(dir) or global:upper(dir)
         qOut:fill(indexerIn(self._idxIn), ptrIn)
      end

      self._grid:setIndex(self._idxIn)
      self._grid:cellCenter(self._xcIn)
      self._grid:setIndex(idxOut)
      self._grid:cellCenter(self._xcOut)

      for _, bc in ipairs(self._bcList) do
         -- Apply the 'bc' function. This can represent many boundary
         -- condition types ranging from a simple copy or a reflection
         -- with the sign flit to QM based electron emission model.
         bc(dir, tCurr, self._idxIn, ptrIn, ptrOut, self._xcOut, self._xcIn)
      end
   end

   self._isFirst = false
end

function Bc:storeBoundaryFlux(tCurr, rkIdx, qOut)
   local ptrOut = qOut:get(1) -- Get pointers to (re)use inside inner loop.
   local indexer = qOut:genIndexer()

   local tId = self._grid:subGridSharedId() -- Local thread ID.
   for idxOut in self._ghostRangeDecomp:rowMajorIter(tId) do
      idxOut:copyInto(self._idxIn)

      qOut:fill(indexer(idxOut), ptrOut) 

      -- Before operating on ghosts, store ghost values for later flux diagnostics
      self._idxIn[self._dir] = 1
      self._boundaryFluxFields[rkIdx]:fill(self._boundaryIdxr(self._idxIn), self._boundaryPtr[rkIdx])
      for c = 1, qOut:numComponents() do
         self._boundaryPtr[rkIdx][c] = ptrOut[c]
      end
   end
end

function Bc:evalOnConfBoundary(inFld)
   if self._isFirst then
      self._confBoundaryField    = self._confBoundaryField or createFieldFromField(self._confBoundaryGrid, inFld, {1,1})
      self._confBoundaryFieldPtr = self._confBoundaryFieldPtr or self._confBoundaryField:get(1)
      self._confBoundaryIdxr     = self._confBoundaryIdxr or self._confBoundaryField:genIndexer()

      local confGlobal        = inFld:globalRange()
      local confGlobalExt     = inFld:globalExtRange()
      local confLocalExtRange = inFld:localExtRange()
      self._confGhostRng      = self._confGhostRng or confLocalExtRange:intersect(
   	 self:getGhostRange(confGlobal, confGlobalExt) ) -- Range spanning ghost cells.
      -- Decompose ghost region into threads.
      self._confGhostRangeDecomp = self._confGhostRangeDecomp or LinearDecomp.LinearDecompRange {
   	 range = self._confGhostRng, numSplit = self._grid:numSharedProcs() }
   end

   local inFldPtr = inFld:get(1)
   local indexer = inFld:genIndexer()

   local tId = self._grid:subGridSharedId() -- Local thread ID.
   for idxIn in self._confGhostRangeDecomp:rowMajorIter(tId) do
      idxIn:copyInto(self._idxOut)
      self._idxOut[self._dir] = 1
      
      inFld:fill(indexer(idxIn), inFldPtr) 
      self._confBoundaryField:fill(self._confBoundaryIdxr(self._idxOut), self._confBoundaryFieldPtr)
       
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

function Bc:initBcDiagnostics(cDim)
   if cDim == self._grid:ndim() then
      self._confBoundaryGrid = self._boundaryGrid
   else 
      -- Create reduced boundary config-space grid with 1 cell in dimension of self._dir.
      local reducedLower, reducedUpper, reducedNumCells, reducedCuts = {}, {}, {}, {}
      for d = 1, cDim do
         if d==self._dir then
            table.insert(reducedLower, -self._grid:dx(d)/2)
            table.insert(reducedUpper,  self._grid:dx(d)/2)
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
         comm      = self._splitComm,
         writeRank = self.writeRank,
         cuts      = reducedCuts,
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
