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
  --First Section is standard from Bc
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
   --End First Section

   --Second Section passes new properties needed for projection
   self._basis = tbl.basis
   self.numDens = tbl.numDens
   self.uPar = tbl.uPar
   self.vtSq = tbl.vtSq
   self.maxwellianGhost = tbl.maxwellianGhost
   self.confGrid=tbl.confGrid
   self.confBasis=tbl.confBasis
   self.mass = tbl.mass
   self.projMaxwell=tbl.projMaxwell
   self.bmag=tbl.bmag
   self.distf = tbl.distf
   self._boundaryGrid=tbl.boundaryGrid
   self._confBoundaryGrid=tbl.confBoundaryGrid
   self.species = tbl.species
   self.bcFunc1=tbl.bcFunc1
   self.jacobGeo = tbl.jacobGeo
   self.numDensity = tbl.numDensity
   --Set Up weakMultiplier
   self.weakMultiplyConfPhase = CartFieldBinOp {
      onGrid     = self._boundaryGrid,
      weakBasis  = self._basis,
      fieldBasis = self.confBasis,
      operation  = "Multiply",
      onGhosts   = true,
   }
   --End Second section

   --Third section is standard
   local advArgs       = tbl.advanceArgs  -- Sample arguments for advance method.

   -- A function can be specied in the ghost cell layer to be used as a boundary condition.
   -- Additionally, a feedback function of fluid moments at the boundary can be set up.
   self._idxIn  = Lin.IntVec(self._grid:ndim())
   self._idxOut = Lin.IntVec(self._grid:ndim())
   self._xcIn   = Lin.Vec(self._grid:ndim())
   self._xcOut  = Lin.Vec(self._grid:ndim())

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
        comm      = self._splitComm,  cuts      = reducedCuts,
        writeRank = writeRank,        useShared = self._grid:isShared(),
     }
     self._boundaryGrid = Grid.RectCart {
        lower = reducedLower,  cells = reducedNumCells,
        upper = reducedUpper,  decomposition = reducedDecomp,
     }
  end

   -- Initialize tools constructed from fields (e.g. ranges).
   self.fldTools = advArgs and self:initFldTools(advArgs[1],advArgs[2]) or nil
  --End Third Section
   print("Initialized Bc")
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
   self.projMaxwell:advance(time,{self.numDens,self.uPar,self.vtSq,self.bmag},{tools.ghostFld})
   
   return tools
end

-- 2 Additional Methods needed for later b/c species doesn't have them and need fields to be on boundary grid
function MaxwellGhostBc:allocCartField()
   local f = DataStruct.Field {
        onGrid        = self._confBoundaryGrid,
        numComponents = self.confBasis:numBasis(),
        ghost         = {1, 1},
        metaData      = {polyOrder = self.confBasis:polyOrder(),
                         basisType = self.confBasis:id(),
                         charge = self.species.charge,
                         mass=self.species.mass},
   }
   f:clear(0.0)
   return f
end

function MaxwellGhostBc:allocMoment()
   return self:allocCartField()
end


-- Method Needed to scale density exactly
function MaxwellGhostBc:scaleDensity(distf)
   local M0e, M0 = self:allocMoment(), self:allocMoment()
   local M0mod   = self:allocMoment()

   self.numDensityCalc = DistFuncMomentCalc {
      onGrid     = self._boundaryGrid,   confBasis  = self.confBasis,
      phaseBasis = self._basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM0", -- GkM0 = < f >
   }
   self.numDensityCalc:advance(0.0, {distf}, {M0})
   local func = function (t, zn)
      return self.bcFunc1(t,zn,self.species)
   end
   local project = ProjectOnBasis {
      onGrid   = self._confBoundaryGrid,
      basis    = self.confBasis,
      evaluate = func,
      onGhosts = true,
   }
   project:advance(0.0, {}, {M0e})

   local weakDivision = CartFieldBinOp {
      --onGrid    = self.confGrid,
      onGrid    = self._confBoundaryGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
      onGhosts  = true,
   }

   -- Calculate M0mod = M0e / M0.
   weakDivision:advance(0.0, {M0, M0e}, {M0mod})
   -- Calculate distff = M0mod * distf.
   self.weakMultiplyConfPhase:advance(0.0, {M0mod, distf}, {distf})
end


-- The advance method
function MaxwellGhostBc:advance(tCurr, inFld, outFld)
   local grid = self._grid
   local qOut = assert(outFld[1], "Bc.advance: Must-specify an output field")

   local dir = self._dir

   --First Section to project into the ghost field
   self.projMaxwell:advance(time,{self.numDens,self.uPar,self.vtSq,self.bmag},{self.fldTools.ghostFld})
   --scale Density and multiply by conf space Jacobian
   self:scaleDensity(self.fldTools.ghostFld)
   self.weakMultiplyConfPhase:advance(0, {self.fldTools.ghostFld, self.jacobGeo}, {self.fldTools.ghostFld})
   --End first section
   
   --Second Section standard from Bc for function type advance
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
   --End Second section
end


function MaxwellGhostBc:getDir() return self._dir end

function MaxwellGhostBc:getEdge() return self._edge end

function MaxwellGhostBc:label() return "Flux"..self._dirlabel..self._edge end

function MaxwellGhostBc:getBoundaryGrid() return self._boundaryGrid end

return MaxwellGhostBc
