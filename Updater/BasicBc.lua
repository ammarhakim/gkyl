-- Gkyl ------------------------------------------------------------------------
--
-- Apply a basic boundary condition, in which the function in the ghost cell
-- is only a function of the skin cell next to it.
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

-- Boundary condition updater.
local BasicBc = Proto(UpdaterBase)

local dirlabel = {"X", "Y", "Z"}

function BasicBc:init(tbl)
   BasicBc.super.init(self, tbl) -- Setup base object.

   self._grid = assert(tbl.onGrid, "Updater.BasicBc: Must specify grid to use with 'onGrid''")
   self._dir  = assert(tbl.dir, "Updater.BasicBc: Must specify direction to apply BCs with 'dir'")
   self._dirlabel = dirlabel[self._dir]

   self._edge = assert(
      tbl.edge, "Updater.BasicBc: Must specify edge to apply BCs with 'edge'. Must be one of 'upper', or 'lower'.")
   if self._edge ~= "lower" and self._edge ~= "upper" then
      error("Updater.BasicBc: 'edge' must be one of 'lower' or 'upper'. Was " .. self._edge .. " instead")
   end
   self._bcList = assert(
      tbl.boundaryConditions, "Updater.BasicBc: Must specify boundary conditions to apply with 'boundaryConditions'")

   self._skinLoop = tbl.skinLoop and tbl.skinLoop or "pointwise"
   self._cDim = tbl.cdim
   if self._skinLoop == "flip" then
      assert(self._cDim, "Updater.BasicBc: Must specify configuration space dimensions to apply with 'cdim'")
      self._vdir = assert(tbl.vdir, "Updater.BasicBc: Must specify velocity direction to flip with 'vdir'")
   end

   local advArgs = tbl.advanceArgs  -- Sample arguments for advance method.

   self._idxIn  = Lin.IntVec(self._grid:ndim())
   self._idxOut = Lin.IntVec(self._grid:ndim())
   self._xcIn   = Lin.Vec(self._grid:ndim())
   self._xcOut  = Lin.Vec(self._grid:ndim())

   -- Initialize tools constructed from fields (e.g. ranges).
   self.fldTools = advArgs and self:initFldTools(advArgs[1],advArgs[2]) or nil

end

function BasicBc:getGhostRange(global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if self._edge == "lower" then
      uv[self._dir] = global:lower(self._dir)-1   -- For ghost cells on "left".
   else
      lv[self._dir] = global:upper(self._dir)+1   -- For ghost cells on "right".
   end
   return Range.Range(lv, uv)
end

function BasicBc:initFldTools(inFld, outFld)
   -- Pre-initialize tools (ranges, pointers, etc) depending on fields and used in the advance method.
   local tools = {}

   local distf = inFld[1]
   local qOut  = assert(outFld[1], "BasicBc.advance: Must-specify an output field")

   local grid = self._grid

   local global     = qOut:globalRange()
   local globalExt  = qOut:globalExtRange()
   local localExt   = qOut:localExtRange()
   local ghostRange = localExt:intersect(self:getGhostRange(global, globalExt))   -- Range spanning ghost cells.
   -- Decompose ghost region into threads.
   tools.ghostRangeDecomp = LinearDecomp.LinearDecompRange {
      range = ghostRange, numSplit = grid:numSharedProcs() }

   -- Get the in and out indexers. 
   tools.indexerOut, tools.indexerIn = qOut:genIndexer(), qOut:genIndexer()

   self.flipIdx = self._skinLoop == "flip" 
      and function(idxIn) idxIn[self._vdir] = global:upper(self._vdir) + 1 - idxIn[self._vdir] end
      or function(idxIn) end

   return tools
end

function BasicBc:_advanceBasic(tCurr, inFld, outFld)
   -- Advance method for the basic BCs (copy, absorb, reflect, open, extern).
   local qOut = assert(outFld[1], "BasicBc.advance: Must-specify an output field")

   local grid      = self._grid
   local dir, edge = self._dir, self._edge
   local global    = qOut:globalRange()

   -- Get the in and out pointers.
   local ptrOut, ptrIn = qOut:get(1), qOut:get(1)

   local tId = grid:subGridSharedId() -- Local thread ID.
   for idxOut in self.fldTools.ghostRangeDecomp:rowMajorIter(tId) do 
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

function BasicBc:_advance(tCurr, inFld, outFld)
   self.fldTools = self.fldTools or self:initFldTools(inFld,outFld)

   self:_advanceBasic(tCurr, inFld, outFld)
end

function BasicBc:getDir() return self._dir end

function BasicBc:getEdge() return self._edge end

function BasicBc:label() return "Flux"..self._dirlabel..self._edge end

return BasicBc
