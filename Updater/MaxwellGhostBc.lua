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
local CartFieldBinOp= require "Updater.CartFieldBinOp"
local DistFuncMomentCalc  = require "Updater.DistFuncMomentCalc"

-- Boundary condition updater.
local MaxwellGhostBc = Proto(UpdaterBase)
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


function MaxwellGhostBc:init(tbl)
   MaxwellGhostBc.super.init(self, tbl) -- Setup base object.

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

   -- A function can be specied in the ghost cell layer to be used as a boundary condition.
   -- Additionally, a feedback function of fluid moments at the boundary can be set up.
   self._idxIn  = Lin.IntVec(self._grid:ndim())
   self._idxOut = Lin.IntVec(self._grid:ndim())
   self._xcIn   = Lin.Vec(self._grid:ndim())
   self._xcOut  = Lin.Vec(self._grid:ndim())

  self._boundaryGrid = tbl.boundaryGrid
  self._confBoundaryGrid = tbl.confBoundaryGrid

   -- Initialize tools constructed from fields (e.g. ranges).
   self.fldTools = advArgs and self:initFldTools(advArgs[1],advArgs[2]) or nil
end


function MaxwellGhostBc:getGhostRange(global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if self._edge == "lower" then
      uv[self._dir] = global:lower(self._dir)-1   -- For ghost cells on "left".
   else
      lv[self._dir] = global:upper(self._dir)+1   -- For ghost cells on "right".
   end
   return Range.Range(lv, uv)
end


function MaxwellGhostBc:initFldTools(inFld, outFld)
   -- Pre-initialize tools (ranges, pointers, etc) depending on fields and used in the advance method.
   local tools = {}

   local distf = inFld[1]
   local inGhostRange = inFld[2] -- Optional range on which we wish to apply BCs.
   local qOut  = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local grid = self._grid

   local global     = qOut:globalRange()
   local globalExt  = qOut:globalExtRange()
   local localExt   = qOut:localExtRange()
   local ghostRange = localExt:intersect(self:getGhostRange(global, globalExt))   -- Range spanning ghost cells.
   if inGhostRange then
      local ghostRangeAll = localExt:intersect(self:getGhostRange(global, globalExt)) 
      ghostRange = ghostRangeAll:intersect(inGhostRange)
   end
   -- Decompose ghost region into threads.
   tools.ghostRangeDecomp = LinearDecomp.LinearDecompRange {
      range = ghostRange, numSplit = grid:numSharedProcs() }

   -- Get the in and out indexers. 
   tools.indexerOut, tools.indexerIn = qOut:genIndexer(), qOut:genIndexer()

   self.flipIdx = self._skinLoop == "flip" 
      and function(idxIn) idxIn[self._vdir] = global:upper(self._vdir) + 1 - idxIn[self._vdir] end
      or function(idxIn) end
  
   --Use projMaxwell instead of evaluate functon
   tools.ghostFld        = createFieldFromField(self._boundaryGrid, qOut, {1,1})
   tools.ghostFldIndexer = tools.ghostFld:genIndexer()
   return tools
end

-- The advance method
function MaxwellGhostBc:advance(tCurr, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local dir = self._dir

   -- Get the in and out pointers
   local ptrOut, ptrIn = qOut:get(1), self.fldTools.ghostFld:get(1)

   local tId = grid:subGridSharedId() -- Local thread ID.
   for idxOut in self.fldTools.ghostRangeDecomp:rowMajorIter(tId) do
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


function MaxwellGhostBc:getDir() return self._dir end

function MaxwellGhostBc:getEdge() return self._edge end

function MaxwellGhostBc:label() return "Flux"..self._dirlabel..self._edge end

function MaxwellGhostBc:getBoundaryGrid() return self._boundaryGrid end

return MaxwellGhostBc
