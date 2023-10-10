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
local DataStruct     = require "DataStruct"
local Grid           = require "Grid"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local Proto          = require "Lib.Proto"
local Range          = require "Lib.Range"
local UpdaterBase    = require "Updater.Base"
local MomDecl        = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"

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
   self._cDim = tbl.cdim
   if self._skinLoop == "flip" then
      assert(self._cDim, "Updater.Bc: Must specify configuration space dimensions to apply with 'cdim'")
      self._vdir = assert(tbl.vdir, "Updater.Bc: Must specify velocity direction to flip with 'vdir'")
   end

   local advArgs       = tbl.advanceArgs  -- Sample arguments for advance method.
   self._doDiagnostics = xsys.pickBool(tbl.doDiagnostics, false)   -- Indicate whether boundary diagnostics are needed.
   if self._doDiagnostics then
      assert(advArgs, "Updater.Bc: requesting boundary diagnostics requires passing sample arguments for the advance method through 'advanceArgs'.")
      assert(self._cDim, "Updater.Bc: Must specify configuration space dimensions to apply with 'cdim'")
   end

   -- A function can be specied in the ghost cell layer to be used as a boundary condition.
   -- Additionally, a feedback function of fluid moments at the boundary can be set up.
   self._evaluateFn = tbl.evaluate
   if self._evaluateFn then
      self._basis    = assert(tbl.basis, "Bc.init: Evaluate is currently implemented only for DG; 'basis' must be specified.")
      self._feedback = xsys.pickBool(tbl.feedback, false)
      self._evolveFn = xsys.pickBool(tbl.evolveFn, self._feedback)
      if self._feedback then
         assert(self._evolveFn, "Updater.Bc: if feedback=true 'evolveFn' must also be true.")
         assert(tbl.cdim == 1, "Bc.init: Feedback boundary condition is implemented only for 1X simulations.")
      end
   end

   self._idxIn  = Lin.IntVec(self._grid:ndim())
   self._idxOut = Lin.IntVec(self._grid:ndim())
   self._xcIn   = Lin.Vec(self._grid:ndim())
   self._xcOut  = Lin.Vec(self._grid:ndim())

   -- A function can be specied in the ghost cell layer to be used as
   -- a boundary condition. Additionally, a feedback function of fluid
   -- moments at the boundary can be set up.
   if self._evaluateFn then
      if self._feedback then
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

         self.advImpl = function(tCurr, inFld, outFld) Bc["_advanceFeedback"](self, tCurr, inFld, outFld) end
      else
         self.advImpl = function(tCurr, inFld, outFld) Bc["_advanceBCFunc"](self, tCurr, inFld, outFld) end
      end
   else
      self.advImpl = function(tCurr, inFld, outFld) Bc["_advanceBasic"](self, tCurr, inFld, outFld) end
   end

   if self._doDiagnostics or self._evaluateFn then 
       -- Create reduced boundary grid with 1 cell in dimension of self._dir.
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
         comm      = self._splitComm,  cuts = reducedCuts,
         writeRank = writeRank,
      }
      self._boundaryGrid = Grid.RectCart {
         lower = reducedLower,  cells = reducedNumCells,
         upper = reducedUpper,  decomposition = reducedDecomp,
      }

      if self._doDiagnostics then
         if self._cDim == self._grid:ndim() then
            self._confBoundaryGrid = self._boundaryGrid
         else 
            -- Create reduced boundary config-space grid with 1 cell in dimension of self._dir.
            local reducedLower, reducedUpper, reducedNumCells, reducedCuts = {}, {}, {}, {}
            for d = 1, self._cDim do
               if d==self._dir then
                  table.insert(reducedLower, -self._grid:dx(d)/2)
                  table.insert(reducedUpper,  self._grid:dx(d)/2)
                  table.insert(reducedNumCells, 1)
                  table.insert(reducedCuts, 1)
               else
                  table.insert(reducedLower,    self._grid:lower(d))
                  table.insert(reducedUpper,    self._grid:upper(d))
                  table.insert(reducedNumCells, self._grid:numCells(d))
                  table.insert(reducedCuts,     self._grid:cuts(d))
               end
            end
            local reducedDecomp = CartDecomp.CartProd {
               comm      = self._splitComm,  cuts      = reducedCuts,
               writeRank = self.writeRank,
            }
            self._confBoundaryGrid = Grid.RectCart {
               lower = reducedLower,  cells = reducedNumCells,
               upper = reducedUpper,  decomposition = reducedDecomp,
            }
         end
      end
   end

   -- Initialize tools constructed from fields (e.g. ranges).
   self.fldTools = advArgs and self:initFldTools(advArgs[1],advArgs[2]) or nil

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

function Bc:initFldTools(inFld, outFld)
   -- Pre-initialize tools (ranges, pointers, etc) depending on fields and used in the advance method.
   local tools = {}

   local distf = inFld[1]
   local qOut  = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local grid = self._grid

   local global     = qOut:globalRange()
   local globalExt  = qOut:globalExtRange()
   local localExt   = qOut:localExtRange()
   tools.ghostRange = localExt:intersect(self:getGhostRange(global, globalExt))   -- Range spanning ghost cells.

   -- Get the in and out indexers. 
   tools.indexerOut, tools.indexerIn = qOut:genIndexer(), qOut:genIndexer()

   self.flipIdx = self._skinLoop == "flip" 
      and function(idxIn) idxIn[self._vdir] = global:upper(self._vdir) + 1 - idxIn[self._vdir] end
      or function(idxIn) end

   if self._doDiagnostics then
      -- Make enough boundary field copies for 4 RK stages.
      self._boundaryFluxFields = {createFieldFromField(self._boundaryGrid, qOut, {1,1}),
                                  createFieldFromField(self._boundaryGrid, qOut, {1,1}),
                                  createFieldFromField(self._boundaryGrid, qOut, {1,1}),
                                  createFieldFromField(self._boundaryGrid, qOut, {1,1}),}
      self._boundaryFluxRate      = createFieldFromField(self._boundaryGrid, qOut, {1,1})
      self._boundaryFluxFieldPrev = createFieldFromField(self._boundaryGrid, qOut, {1,1})
      self._boundaryPtr = {self._boundaryFluxFields[1]:get(1), 
                           self._boundaryFluxFields[2]:get(1),
                           self._boundaryFluxFields[3]:get(1),
                           self._boundaryFluxFields[4]:get(1),}
      self._boundaryIdxr = self._boundaryFluxFields[1]:genIndexer()
   end

   if self._evaluateFn then
      tools.ghostFld        = createFieldFromField(self._boundaryGrid, qOut, {1,1})
      tools.ghostFldIndexer = tools.ghostFld:genIndexer()

      local myEvaluateFn = self._evaluateFn
      self._projectEvaluateFn = ProjectOnBasis {
         onGrid   = self._boundaryGrid,
         basis    = self._basis,
         evaluate = myEvaluateFn,
      }

      if self._feedback then
         local pDim, cDim = grid:ndim(), 1
         local vDim       = pDim - cDim
         local phaseRange, phaseIndexer = distf:localRange(), distf:genIndexer()
         -- Construct ranges for nested loops.
         tools.confSkinRange = self._edge == "lower"
            and phaseRange:selectFirst(cDim):lowerSkin(self._dir,1) or phaseRange:selectFirst(cDim):upperSkin(self._dir,1)
         tools.velRange = phaseRange:selectLast(vDim)
         tools.velRange = tools.velRange:extendDir(self._dir,self.partialMomDirExts[1],self.partialMomDirExts[2])
         tools.phaseIndexer = distf:genIndexer()
      else
         self._projectEvaluateFn:advance(tCurr, {}, {tools.ghostFld})
      end
   end
   
   return tools
end

function Bc:_advanceBasic(tCurr, inFld, outFld)
   -- Advance method for the basic BCs (copy, absorb, reflect, open, extern).
   local qOut = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local grid      = self._grid
   local dir, edge = self._dir, self._edge
   local global    = qOut:globalRange()

   -- Get the in and out pointers.
   local ptrOut, ptrIn = qOut:get(1), qOut:get(1)

   for idxOut in self.fldTools.ghostRange:rowMajorIter() do 
      qOut:fill(self.fldTools.indexerOut(idxOut), ptrOut)

      -- Copy out index into in index
      idxOut:copyInto(self._idxIn)
      self.flipIdx(self._idxIn)   -- Reverse the in velocity index if needed
      self._idxIn[dir] = edge=="lower" and global:lower(dir) or global:upper(dir)

      qOut:fill(self.fldTools.indexerIn(self._idxIn), ptrIn)

      self._grid:setIndex(self._idxIn);  self._grid:cellCenter(self._xcIn);
      self._grid:setIndex(idxOut);       self._grid:cellCenter(self._xcOut)

      for _, bc in ipairs(self._bcList) do
         -- Apply the 'bc' function. This can represent many boundary
         -- condition types ranging from a simple copy or a reflection
         -- with the sign flit to QM based electron emission model.
         bc(dir, tCurr, self._idxIn, ptrIn, ptrOut, self._xcOut, self._xcIn)
      end
   end
end

function Bc:_advanceBCFunc(tCurr, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local dir = self._dir

   if self._evolveFn then
      self._projectEvaluateFn:advance(tCurr, {}, {self.fldTools.ghostFld})
   end
   
   -- Get the in and out pointers
   local ptrOut, ptrIn = qOut:get(1), self.fldTools.ghostFld:get(1)

   for idxOut in self.fldTools.ghostRange:rowMajorIter() do 
      qOut:fill(self.fldTools.indexerOut(idxOut), ptrOut)

      -- Copy out index into in index
      idxOut:copyInto(self._idxIn)

      self._idxIn[dir] = 1 -- The boundaryGrid has only 1 cell in dir.
      self.fldTools.ghostFld:fill(self.fldTools.ghostFldIndexer(self._idxIn), ptrIn)

      grid:setIndex(self._idxIn)
      grid:cellCenter(self._xcIn)
      grid:setIndex(idxOut)
      grid:cellCenter(self._xcOut)

      for _, bc in ipairs(self._bcList) do
         -- Apply the 'bc' function. This can represent many boundary
         -- condition types ranging from a simple copy or a reflection
         -- with the sign flit to QM based electron emission model.
         bc(dir, tCurr, self._idxIn, ptrIn, ptrOut, self._xcOut, self._xcIn)
      end
   end
end

function Bc:_advanceFeedback(tCurr, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local dir = self._dir

   -- Mini version of DistFuncMomentCalc to compute moments in the skin cell.
   local distf    = inFld[1]
   local distfItr = distf:get(1)

   local cDim = 1
   local vDim = self._grid:ndim() - cDim
   local edgeM0, edgeM1, edgeM2 = 0.0, 0.0, 0.0
   -- Outer loop is threaded and over configuration space.
   for cIdx in self.fldTools.confSkinRange:rowMajorIter() do
      cIdx:copyInto(self.idxP)
      for k = 1, self.numConfBasis do
         self.mom0[k], self.pMom1[k], self.mom2[k] = 0.0, 0.0, 0.0
      end
      -- Inner loop is over velocity space: no threading to avoid race conditions.
      for vIdx in self.fldTools.velRange:rowMajorIter() do
         for d = 1, vDim do self.idxP[cDim+d] = vIdx[d] end

         grid:setIndex(self.idxP)
         grid:cellCenter(self.xcP)
         grid:getDx(self.dxP)

         distf:fill(self.fldTools.phaseIndexer(self.idxP), distfItr)

         self._mom0Calc( self.xcP:data(), self.dxP:data(), distfItr:data(),  self.mom0:data())
         self._pMom1Calc(self.xcP:data(), self.dxP:data(), distfItr:data(), self.pMom1:data())
         self._mom2Calc( self.xcP:data(), self.dxP:data(), distfItr:data(),  self.mom2:data())
      end
      -- This currently only works for 1x. In higher dimensions one gets an expansion
      -- in describing the variation parallel to the boundary. How to handle that expansion
      -- requires a better implementation of the feedback BC.
      for k = 1, self.numConfBasis do
         edgeM0 = edgeM0 + self.mom0[k]*self._confBasisEdge[k]
         edgeM1 = edgeM1 + self.pMom1[k]*self._confBasisEdge[k]
         edgeM2 = edgeM2 + self.mom2[k]*self._confBasisEdge[k]
      end

      myEvaluateFn = function(t, z) return self._evaluateFn(t, z, edgeM0, edgeM1, edgeM2) end
   end
   self._projectEvaluateFn:setFunc(myEvaluateFn)
   
   self._projectEvaluateFn:advance(tCurr, {}, {self.fldTools.ghostFld})

   -- Get the in and out pointers
   local ptrOut, ptrIn = qOut:get(1), self.fldTools.ghostFld:get(1)

   for idxOut in self.fldTools.ghostRange:rowMajorIter() do 
      qOut:fill(self.fldTools.indexerOut(idxOut), ptrOut)

      -- Copy out index into in index
      idxOut:copyInto(self._idxIn)
      self._idxIn[dir] = 1 -- the boundaryGrid has only 1 cell in dir
      self.fldTools.ghostFld:fill(self.fldTools.ghostFldIndexer(self._idxIn), ptrIn)

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
end

function Bc:_advance(tCurr, inFld, outFld)
   self.fldTools = self.fldTools or self:initFldTools(inFld,outFld)

   self.advImpl(tCurr, inFld, outFld)
end

function Bc:storeBoundaryFlux(tCurr, rkIdx, qOut)
   local ptrOut = qOut:get(1) -- Get pointers to (re)use inside inner loop.
   local indexer = qOut:genIndexer()

   for idxOut in self.fldTools.ghostRange:rowMajorIter() do
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
   if self._confGhostRange==nil then
      self._confBoundaryField    = createFieldFromField(self._confBoundaryGrid, inFld, {1,1})
      self._confBoundaryFieldPtr = self._confBoundaryField:get(1)
      self._confBoundaryIdxr     = self._confBoundaryField:genIndexer()

      local confGlobal        = inFld:globalRange()
      local confGlobalExt     = inFld:globalExtRange()
      local confLocalExtRange = inFld:localExtRange()
      self._confGhostRange    = confLocalExtRange:intersect(
         self:getGhostRange(confGlobal, confGlobalExt) ) -- Range spanning ghost cells.
   end

   local inFldPtr = inFld:get(1)
   local indexer = inFld:genIndexer()

   for idxIn in self._confGhostRange:rowMajorIter() do
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

function Bc:getDir() return self._dir end

function Bc:getEdge() return self._edge end

function Bc:label() return "Flux"..self._dirlabel..self._edge end

function Bc:getBoundaryFluxFields() return self._boundaryFluxFields end

function Bc:getBoundaryFluxRate() return self._boundaryFluxRate end

function Bc:getBoundaryFluxFieldPrev() return self._boundaryFluxFieldPrev end

function Bc:getBoundaryGrid() return self._boundaryGrid end

function Bc:getConfBoundaryGrid() return self._confBoundaryGrid end

return Bc
