-- Gkyl ------------------------------------------------------------------------
--
-- Apply twist shift BCs for flux-tube simulations.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase      = require "App.BCs.BCsBase"
local DataStruct   = require "DataStruct"
local Updater      = require "Updater"
local Mpi          = require "Comm.Mpi"
local Proto        = require "Lib.Proto"
local Time         = require "Lib.Time"
local Range        = require "Lib.Range"
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local CartDecomp   = require "Lib.CartDecomp"
local Grid         = require "Grid"
local DiagsApp     = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags      = require "App.Diagnostics.GkDiagnostics"
local xsys         = require "xsys"
local lume         = require "Lib.lume"

local TwistShiftBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function TwistShiftBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function TwistShiftBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.yShiftFuncIn = assert(tbl.shiftFunction, "TwistShiftBC: must provide the function that computes the y-shift in 'shiftFunction'.")
   
   self.yShiftPolyOrder = tbl.shiftPolyOrder

   self.saveFlux = tbl.saveFlux or false
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
         self.saveFlux       = true
      end
   end
end

function TwistShiftBC:setName(nm) self.name = self.speciesName.."_"..nm end

function TwistShiftBC:createSolver(mySpecies, field, externalField)

--   print("arrived", self.bcEdge)
   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   if self.yShiftPolyOrder == nil then self.yShiftPolyOrder = self.basis:polyOrder() end

   if self.bcEdge=="lower" then
      self.yShiftFunc = function(t,xn) return -self.yShiftFuncIn(t,xn) end
   else
      self.yShiftFunc = self.yShiftFuncIn
   end

   self.bcSolver = Updater.TwistShiftBC {
      onGrid    = self.grid,       yShiftFunc      = self.yShiftFunc,
      basis     = self.basis,      yShiftPolyOrder = self.yShiftPolyOrder,
      confBasis = self.confBasis,  edge            = self.bcEdge,
   }

   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   self:createBoundaryGrid()

   -- Create a boundary grid that is global in Y: used to broadcast skin data.
   self:createBoundaryGridGlobalY()

   local distf = mySpecies["getDistF"] and mySpecies:getDistF() or mySpecies:getMoments()
   -- Define methods to allocate fields defined on boundary grid (used by e.g. diagnostics).
   self.allocCartField = function(self, grid, nComp, ghosts, metaData)
      local f = DataStruct.Field {
         onGrid        = grid,   ghost    = ghosts,
         numComponents = nComp,  metaData = metaData,
      }
      f:clear(0.0)
      return f
   end
   local allocDistf = function()
      return self:allocCartField(self.boundaryGrid, distf:numComponents(),
                                 {distf:lowerGhost(),distf:upperGhost()}, distf:getMetaData())
   end
   local allocDistfGlobalY = function()
      return self:allocCartField(self.boundaryGridGlobalY, distf:numComponents(),
                                 {distf:lowerGhost(),distf:upperGhost()}, distf:getMetaData())
   end
   self.boundaryField    = allocDistf()
   self.boundaryFieldPtr = self.boundaryField:get(1)
   self.boundaryFieldGlobalY = allocDistfGlobalY()

   -- Create the communicator used to send information along z.
   self:createGraphComm()

   -- Create MPI derived data type or receiving data.
   self:createSyncMPIdataTypes(distf)

   -- Create the range needed to loop over ghosts.
   local global, globalExt, localExtRange = distf:globalRange(), distf:globalExtRange(), distf:localExtRange()
   self.ghostRange = localExtRange:intersect(self:getGhostRange(global, globalExt))
   -- Decompose ghost region into threads.
   self.ghostRangeDecomp = LinearDecomp.LinearDecompRange{range=self.ghostRange, numSplit=self.grid:numSharedProcs()}

   self.boundaryIdxr = self.boundaryField:genIndexer()

   self.idxOut = Lin.IntVec(self.grid:ndim())

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then
      -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
      self:createConfBoundaryGrid()

      local numDensity = mySpecies:getNumDensity()
      self.allocMoment = function(self)
         return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                    {numDensity:lowerGhost(),numDensity:upperGhost()}, numDensity:getMetaData())
      end
      self.allocVectorMoment = function(self, dim)
         return self:allocCartField(self.confBoundaryGrid, dim*self.confBasis:numBasis(),
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

      self.tId = self.grid:subGridSharedId() -- Local thread ID.

      -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
      self.confBoundaryField    = self:allocMoment()
      self.confBoundaryFieldPtr = self.confBoundaryField:get(1)
      self.confBoundaryIdxr     = self.confBoundaryField:genIndexer()
      local confGlobal        = numDensity:globalRange()
      local confGlobalExt     = numDensity:globalExtRange()
      local confLocalExtRange = numDensity:localExtRange()
      self.confGhostRange = confLocalExtRange:intersect(self:getGhostRange(confGlobal, confGlobalExt)) -- Range spanning ghost cells.
      -- Decompose ghost region into threads.
      self.confGhostRangeDecomp = LinearDecomp.LinearDecompRange {range=self.confGhostRange, numSplit=self.grid:numSharedProcs()}

      -- Evaluate the magnetic field and jacobGeo in the boundary (needed by diagnostics).
      local bmag = externalField.geo.bmag 
      self.bmag = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                      {bmag:lowerGhost(),bmag:upperGhost()}, bmag:getMetaData())
      self.bmag:copy(self:evalOnConfBoundary(bmag))
      local bmagInvSq = externalField.geo.bmagInvSq
      self.bmagInvSq = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                          {bmagInvSq:lowerGhost(),bmagInvSq:upperGhost()}, bmagInvSq:getMetaData())
      self.bmagInvSq:copy(self:evalOnConfBoundary(bmagInvSq))
      local jacobGeo = externalField.geo.jacobGeo
      if jacobGeo then
         self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                             {jacobGeo:lowerGhost(),jacobGeo:upperGhost()}, jacobGeo:getMetaData())
         self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo))
      end
      local jacobGeoInv = externalField.geo.jacobGeoInv
      if jacobGeoInv then
         self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(),
                                                {jacobGeoInv:lowerGhost(),jacobGeoInv:upperGhost()}, jacobGeoInv:getMetaData())
         self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv))
      end

      -- Declare methods/functions needed for handling saved fluxes and needed by diagnostics.
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

      -- Number density calculator. Needed regardless of diagnostics (for recycling BCs).
      local mass = mySpecies.mass
      self.numDensityCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
         phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
         moment     = "GkM0", -- GkM0 = < f >
      }

      if not self.anyDiagnostics then
         self.calcBoundaryFluxRateFunc = function(dtIn) end
      else
         self.calcBoundaryFluxRateFunc = function(dtIn)
            -- Compute boundary flux rate ~ (fGhost_new - fGhost_old)/dt.
            self.boundaryFluxRate:combine( 1.0/dtIn, self.boundaryFluxFields[1],
                                          -1.0/dtIn, self.boundaryFluxFieldPrev)
            self.boundaryFluxFieldPrev:copy(self.boundaryFluxFields[1])
         end
         -- Set up weak multiplication and division operators (for diagnostics).
         self.confWeakMultiply = Updater.CartFieldBinOp {
            onGrid    = self.confBoundaryGrid,  operation = "Multiply",
            weakBasis = self.confBasis,         onGhosts  = true,
         }
         self.confWeakDivide = Updater.CartFieldBinOp {
            onGrid    = self.confBoundaryGrid,  operation = "Divide",
            weakBasis = self.confBasis,         onGhosts  = true,
         }
         -- Volume integral operator (for diagnostics).
         self.volIntegral = {
            scalar = Updater.CartFieldIntegratedQuantCalc {
               onGrid = self.confBoundaryGrid,  numComponents = 1,
               basis  = self.confBasis,         quantity      = "V",
            }
         }
         -- Moment calculators (for diagnostics).
         local mass = mySpecies.mass
         self.momDensityCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
            moment     = "GkM1", -- GkM1 = < v_parallel f >
         }
         self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
            moment     = "GkM2", -- GkM2 = < (v_parallel^2 + 2*mu*B/m) f >
         }
         self.M2parCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
            moment     = "GkM2par", -- GkM2par = < v_parallel^2 f >
         }
         self.M3parCalc = Updater.DistFuncMomentCalc {
            onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
            phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
            moment     = "GkM3par", -- GkM3par = < v_parallel^3 f >
         }
         if self.vdim > 1 then
            self.M2perpCalc = Updater.DistFuncMomentCalc {
               onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
               phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
               moment     = "GkM2perp", -- GkM2 = < (mu*B/m) f >
            }
            self.M3perpCalc = Updater.DistFuncMomentCalc {
               onGrid     = self.boundaryGrid,  confBasis = self.confBasis,
               phaseBasis = self.basis,         gkfacs    = {mass, self.bmag},
               moment     = "GkM3perp", -- GkM3perp = < vpar*(mu*B/m) f >
            }
         end
         self.divideByJacobGeo = self.jacobGeoInv
            and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeoInv}, {fldOut}) end
            or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
         self.multiplyByJacobGeo = self.jacobGeo
            and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeo}, {fldOut}) end
            or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
      end
   else
      self.storeBoundaryFluxFunc        = function(tCurr, rkIdx, qOut) end
      self.copyBoundaryFluxFieldFunc    = function(inIdx, outIdx) end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...) end
      self.calcBoundaryFluxRateFunc     = function(dtIn) end
   end
end

function TwistShiftBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function TwistShiftBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function TwistShiftBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function TwistShiftBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function TwistShiftBC:createDiagnostics(mySpecies, field)
   -- Create BC diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = GkDiags()}
      self.diagnostics:fullInit(mySpecies, field, self)
      -- Presently boundary diagnostics are boundary flux diagnostics. Append 'flux' to the diagnostic's
      -- name so files are named accordingly. Re-design this when non-flux diagnostics are implemented
      self.diagnostics.name = self.diagnostics.name..'_flux'
   end
   return self.diagnostics
end

function TwistShiftBC:createBoundaryGridGlobalY()
   -- Create a ghost boundary grid with only one cell in the direction of the BC
   -- and global in y (i.e. not decomposed along y).

   if self.grid:ndim() < 3 then return end

   local reducedLower, reducedUpper, reducedNumCells, reducedCuts = {}, {}, {}, {}
   for d = 1, self.grid:ndim() do
      if d==self.bcDir then
         table.insert(reducedLower, -self.grid:dx(d)/2.)
         table.insert(reducedUpper,  self.grid:dx(d)/2.)
         table.insert(reducedNumCells, 1)
         table.insert(reducedCuts, 1)
      else
         table.insert(reducedLower,    self.grid:lower(d))
         table.insert(reducedUpper,    self.grid:upper(d))
         table.insert(reducedNumCells, self.grid:numCells(d))
         table.insert(reducedCuts,     self.grid:cuts(d))
      end
   end
   reducedCuts[2] = 1 -- Do not decompose along y.
   local commSet  = self.grid:commSet()
   local worldComm, nodeComm = commSet.comm, commSet.nodeComm
   local nodeRank = Mpi.Comm_rank(nodeComm)
   local dirRank  = nodeRank
   local cuts     = {}
   for d=1,3 do cuts[d] = self.grid:cuts(d) or 1 end
   -- For self.bcDir == 3 only (assumes column-major ordering of MPI ranks) ...................
   dirRank = math.floor(nodeRank/(cuts[1]*cuts[2]))*cuts[2]+math.floor((nodeRank % (cuts[1]*cuts[2]))/cuts[1])
   -- ............................................
   self._splitCommGlobalY = Mpi.Comm_split(worldComm, dirRank, nodeRank)
   -- Set up which ranks to write from.
   local writeRank = -1
   if self.bcEdge == "lower" and dirRank == 0 then
      writeRank = nodeRank
   elseif self.bcEdge == "upper" and dirRank == self.grid:cuts(self.bcDir)-1 then
      writeRank = nodeRank
   end
   self.writeRankGlobalY = writeRank
   local reducedDecomp = CartDecomp.CartProd {
      comm      = self._splitCommGlobalY,  cuts      = reducedCuts,
      writeRank = writeRank,               useShared = self.grid:isShared(),
   }
   self.boundaryGridGlobalY = Grid.RectCart {
      lower = reducedLower,  cells         = reducedNumCells,
      upper = reducedUpper,  decomposition = reducedDecomp,
   }
end

function TwistShiftBC:createGraphComm()
   -- Create a graph communicator which connects all the upper ranks in z (sources)
   -- with all each lower rank in z (destination) if bcEdge=lower, and another graph
   -- connecting all the lower ranks in z (sources) with each upper rank in z (destination)
   -- for bcEdge=upper.

   if self.grid:ndim() < 3 then return end

   local xDir, yDir = 1, 2

   -- Create a list of all the lower/upper ranks in z, but at the same location in x.
   -- Some of this code follows DecomposedRange in Lib/CartDecomp.
   local decompRange = self.grid:decomposedRange()

   local cutsIdxr    = decompRange:cutsIndexer()
   local cutsInvIdxr = decompRange:cutsInvIndexer()
   local myId        = self.grid:subGridId() -- Grid ID on this processor.
   local myCutsIdx   = {}
   for d=1,decompRange:ndim() do myCutsIdx[d] = -9 end
   cutsInvIdxr(myId, myCutsIdx)

   -- Loop over "shortened" range which are basically
   -- sub-domains that lie on the lower domain boundary.
   local ones, upper = {}, {}
   for d=1,decompRange:ndim() do ones[d], upper[d] = 1, decompRange:cuts(d) end
   local cutsRange = Range.Range(ones, upper)
   local xyRecvRange, xySendRange
   if self.bcEdge == "lower" then
      xyRecvRange = cutsRange:shorten(self.bcDir)
      xySendRange = cutsRange:shortenFromBelow(self.bcDir)
   elseif self.bcEdge == "upper" then
      xyRecvRange = cutsRange:shortenFromBelow(self.bcDir)
      xySendRange = cutsRange:shorten(self.bcDir)
   end

   -- Sources only depend on x, so we can create a single source list for all y ranks.
   local xLoYRecvRange = xyRecvRange:shorten(yDir)
   -- y range with the same ID along x:
   local yRecvRange = xyRecvRange:extendDir(xDir,-myCutsIdx[xDir]+xyRecvRange:lower(xDir),
                                                  myCutsIdx[xDir]-xyRecvRange:upper(xDir))
   local srcNum, destNum, isSendRank
   if xyRecvRange:contains(myCutsIdx) then  -- This is a receiving rank.
      srcNum, destNum = yRecvRange:volume(), 0
      isSendRank = false
   elseif xySendRange:contains(myCutsIdx) then -- This is a sending rank.
      srcNum, destNum = 0, yRecvRange:volume()
      isSendRank = true
   end

   if isSendRank ~= nil then  -- Only sending/receiving ranks participate.

      -- List ranks along z boundaries starting with the lower-z boundary.
      local zSkinGroupRanks = Lin.IntVec(xyRecvRange:volume()+xySendRange:volume())
      local zLoRange = self.bcEdge == "lower" and xyRecvRange or xySendRange
      local zUpRange = self.bcEdge == "upper" and xyRecvRange or xySendRange
      local i = 0
      for idx in zLoRange:colMajorIter() do
         i = i+1;  zSkinGroupRanks[i] = cutsIdxr(idx)-1
      end
      for idx in zUpRange:colMajorIter() do
         i = i+1;  zSkinGroupRanks[i] = cutsIdxr(idx)-1
      end

      local nodeComm  = self.grid:commSet().comm
      local nodeGroup = Mpi.Comm_group(nodeComm)
      -- Create a group/comm of the ranks along the z-boundary.
      self.zSkinGroup = Mpi.Group_incl(nodeGroup, #zSkinGroupRanks, zSkinGroupRanks:data());
      self.zSkinComm  = Mpi.Comm_create_group(nodeComm, self.zSkinGroup, self.bcEdge == "lower" and 0 or 1);

      -- Determine the sources and destinations in the graph we'll use for the sync.
      local j = 0
      local srcIDs  = Lin.IntVec(srcNum)  -- Source rank IDs along y (for this x).
      local destIDs = Lin.IntVec(destNum) -- destination rank IDs along y (for each this x).
      for yIdx in yRecvRange:colMajorIter() do
         local idxRecv, idxSend = yIdx, yIdx:copy()
         idxSend[self.bcDir] = self.bcEdge == "lower" and idxRecv[self.bcDir]+decompRange:cuts(self.bcDir)-1
                                                       or idxRecv[self.bcDir]-decompRange:cuts(self.bcDir)+1
         j = j+1
         if isSendRank==true then
            destIDs[j] = cutsIdxr(idxRecv)-1
         elseif isSendRank==false then
            srcIDs[j]  = cutsIdxr(idxSend)-1
         end
      end

      -- Translate these ranks from the nodeComm to the zSkinComm.
      local destRanks, srcRanks
      if isSendRank==true then
         destRanks = Mpi.Group_translate_ranks(nodeGroup, destIDs, self.zSkinGroup)
         srcRanks  = Lin.IntVec(srcNum)
      elseif isSendRank==false then
         destRanks = Lin.IntVec(destNum)
         srcRanks  = Mpi.Group_translate_ranks(nodeGroup, srcIDs, self.zSkinGroup)
      end

      local reorder = 0
      local srcW, destW = Lin.IntVec(srcNum), Lin.IntVec(destNum)  -- Weights (not used).
      for i = 1, srcNum do srcW[i] = 1 end
      for i = 1, destNum do destW[i] = 1 end
      -- Create a group/comm with only the lower and upper in z ranks.
      self.graphComm = Mpi.Dist_graph_create_adjacent(self.zSkinComm, srcNum, srcRanks, srcW,
                                                      destNum, destRanks, destW, Mpi.INFO_NULL, reorder)
   else
      self.zSkinComm, self.graphComm = Mpi.COMM_NULL, Mpi.COMM_NULL
   end
end

function TwistShiftBC:createSyncMPIdataTypes(fIn)
   -- Create the MPI derived data (MPIDD) types to sync the field along z.

   if not Mpi.Is_comm_valid(self.graphComm) then return end

   -- Create the receiving MPIDD, which will place the z-skin (including x-ghost
   -- cells for correct data layout) from all ranks along y into a single buffer
   -- that is global in Y.

   local indexer       = self.boundaryFieldGlobalY:genIndexer()
   local numComponents = self.boundaryField:numComponents()
   local localExtRange = self.boundaryField:localExtRange()
   local layout        = self.boundaryField:layout()
   local elctCommType  = self.boundaryField:elemTypeMpi()

   -- Obtain the subDomain ID of the rank that has the same x-location but it's
   -- the first along y. We'll use a collective call to gather all the data along y
   -- and the MPIDD should use the data layout of the first-in-y data chunk.
   local decompRange = self.boundaryGrid:decomposedRange()

   local cutsIdxr, cutsInvIdxr = decompRange:cutsIndexer(), decompRange:cutsInvIndexer()
   local myId        = self.boundaryGrid:subGridId() -- Grid ID on this processor.
   local myCutsIdx   = {}
   for d=1,decompRange:ndim() do myCutsIdx[d] = -9 end
   cutsInvIdxr(myId, myCutsIdx)

   local ones, upper = {}, {}
   for d=1,decompRange:ndim() do ones[d], upper[d] = 1, decompRange:cuts(d) end
   local cutsRange = Range.Range(ones, upper)
   local xLoYRecvRange = cutsRange:shorten(2)
   local yLoIdx, yLoId
   for idx in xLoYRecvRange:colMajorIter() do
      if myCutsIdx[1]==idx[1] and myCutsIdx[self.bcDir]==idx[self.bcDir] then
         yLoIdx, yLoId = idx:copy(), cutsIdxr(idx)
         break
      end
   end

   local srcNum, destNum, _ = Mpi.Dist_graph_neighbors_count(self.graphComm)
   self.isRecvRank = srcNum > 0 and true or false

   self.recvPerMPIDataType, self.recvPerMPILoc = self.boundaryFieldGlobalY:elemTypeMpi(), 0
   local skelIds = decompRange:boundarySubDomainIds(self.bcDir)
--   local myRank        = self.grid:subGridId() -- Whole grid ID on this processor.
   for i = 1, #skelIds do
      local cId = skelIds[i].lower -- Should be = to skelIds[i].upper.

      -- Only create if we are on proper ranks.
      -- Note that if the node communicator has rank size of 1, then we can access all the
      -- memory needed for periodic boundary conditions and we do not need MPI Datatypes.
      if myId == cId and self.bcEdge == "lower" and self.isRecvRank then
         local rgnRecv = decompRange:subDomain(yLoId):lowerSkin(self.bcDir, self.boundaryField:lowerGhost()):extendDir(1,1,1)
         local idx     = rgnRecv:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.recvPerMPILoc      = (indexer(idx)-1)*numComponents
--         print(string.format("%s r:%d id=%d | rgnRecvLower=(%d,%d,%d) rgnRecvUpper=(%d,%d,%d) | idx=(%d,%d,%d) | loc=%d | yLoId=%d",self.bcEdge,myRank,myId,rgnRecv:lower(1),rgnRecv:lower(2),rgnRecv:lower(3),rgnRecv:upper(1),rgnRecv:upper(2),rgnRecv:upper(3),idx[1],idx[2],idx[3],self.recvPerMPILoc/numComponents,yLoId))
         self.recvPerMPIDataType = Mpi.createDataTypeFromRangeAndSubRange(
            rgnRecv, localExtRange, numComponents, layout, elctCommType)
      end
      if myId == cId and self.bcEdge == "upper" and self.isRecvRank then
         local rgnRecv = decompRange:subDomain(yLoId):upperSkin(self.bcDir, self.boundaryField:upperGhost()):extendDir(1,1,1)
         local idx     = rgnRecv:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.recvPerMPILoc      = (indexer(idx)-1)*numComponents
--         print(string.format("%s r:%d id=%d | rgnRecvLower=(%d,%d,%d) rgnRecvUpper=(%d,%d,%d) | idx=(%d,%d,%d) | loc=%d | yLoId=%d",self.bcEdge,myRank,myId,rgnRecv:lower(1),rgnRecv:lower(2),rgnRecv:lower(3),rgnRecv:upper(1),rgnRecv:upper(2),rgnRecv:upper(3),idx[1],idx[2],idx[3],self.recvPerMPILoc/numComponents,yLoId))
         self.recvPerMPIDataType = Mpi.createDataTypeFromRangeAndSubRange(
            rgnRecv, localExtRange, numComponents, layout, elctCommType)
      end
   end

   -- Create sending MPIDDs.
   -- These have to be slightly different than the CartField's own MPIDDs
   -- because in order to use Neighborhood_allgather we need to incorporate
   -- the left x-ghost layer, otherwise when one appends the next chunk of
   -- data along y it'll be in the wrong place.
   local myId            = self.grid:subGridId() -- Whole grid ID on this processor.
   local decomposedRange = self.grid:decomposedRange()
   local skelIds         = decomposedRange:boundarySubDomainIds(self.bcDir)
   local indexer                      = fIn:genIndexer()
   local numComponents, localExtRange = fIn:numComponents(), fIn:localExtRange()
   local layout, elctCommType         = fIn:layout(), fIn:elemTypeMpi()
   self.sendPerMPIDataType, self.sendPerMPILoc = fIn:elemTypeMpi(), 0 
   for i = 1, #skelIds do
      local loId, upId = skelIds[i].lower, skelIds[i].upper

      -- Only create if we are on proper ranks.
      -- Note that if the node communicator has rank size of 1, then we can access all the
      -- memory needed for periodic boundary conditions and we do not need MPI Datatypes.
      if myId == loId and self.bcEdge == "upper" and (not self.isRecvRank) then
         local rgnSend = decomposedRange:subDomain(loId):lowerSkin(self.bcDir, fIn:upperGhost()):extendDir(1,1,1)
         local idx = rgnSend:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.sendPerMPILoc      = (indexer(idx)-1)*numComponents
         self.sendPerMPIDataType = Mpi.createDataTypeFromRangeAndSubRange(
            rgnSend, localExtRange, numComponents, layout, elctCommType)
--         print(string.format("%s s:%d | sendLoc=%d",self.bcEdge, myId, self.sendPerMPILoc/numComponents))
      end
      if myId == upId and self.bcEdge == "lower" and (not self.isRecvRank) then
         local rgnSend = decomposedRange:subDomain(upId):upperSkin(self.bcDir, fIn:lowerGhost()):extendDir(1,1,1)
         local idx = rgnSend:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.sendPerMPILoc      = (indexer(idx)-1)*numComponents
         self.sendPerMPIDataType = Mpi.createDataTypeFromRangeAndSubRange(
            rgnSend, localExtRange, numComponents, layout, elctCommType)
--         print(string.format("%s s:%d | sendLoc=%d",self.bcEdge, myId, self.sendPerMPILoc/numComponents))
      end
   end
end

function TwistShiftBC:zSync(fIn, bufferOut)
   -- bcEdge=lower: each upper-z rank sends their upper-z skin cell data in the fIn field
   --               to all lower-z ranks who receive them into a buffer with 1-cell in z (+ghosts).
   -- bcEdge=upper: each lower-z rank sends their lower-z skin cell data in the fIn field
   --               to all lower-z ranks who receive them into a buffer with 1-cell in z (+ghosts).

   if not Mpi.Is_comm_valid(self.graphComm) then return end

   Mpi.Neighbor_allgather(fIn:dataPointer()+self.sendPerMPILoc, 1, self.sendPerMPIDataType, 
                          bufferOut:dataPointer()+self.recvPerMPILoc, 1, self.recvPerMPIDataType,
                          self.graphComm)
end

-- These are needed to recycle the GkDiagnostics with TwistShiftBC.
function TwistShiftBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                                self.boundaryFluxRate, self.boundaryFluxRate} end
function TwistShiftBC:getFlucF() return self.boundaryFluxRate end

function TwistShiftBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)

   local fIn = mySpecies:rkStepperFields()[outIdx] 

   self:zSync(fIn, self.boundaryFieldGlobalY)

--   -- Apply periodic BCs in z direction, but only on the first edge (lower or upper, whichever is first)
--   -- First check if sync has been called on this step by any BC.
--   local doSync = true
--   for _, bc in lume.orderedIter(mySpecies.nonPeriodicBCs) do 
--      if bc.synced then doSync = false end
--   end
--   if doSync then
--      -- Fetch skin cells from opposite z-edge (donor) and copy to ghost cells (target).
--      local fldGrid = fIn:grid()
--      local periodicDirs = fldGrid:getPeriodicDirs()
--      local modifiedPeriodicDirs = {3}
--      fldGrid:setPeriodicDirs(modifiedPeriodicDirs)
--      fIn:sync()
--      fldGrid:setPeriodicDirs(periodicDirs)
--      self.synced = true
--   else
--      -- Don't sync but reset synced flags for next step.
--      for _, bc in lume.orderedIter(mySpecies.nonPeriodicBCs) do bc.synced = false end
--   end
--
--   -- Copy field ghost cells to boundary field (buffer).
--   self:evalOnBoundary(fIn)
--
--   self.bcSolver:advance(tCurr, {self.boundaryField}, {fIn})
end

function TwistShiftBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return TwistShiftBC
