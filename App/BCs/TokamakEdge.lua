-- Gkyl ------------------------------------------------------------------------
--
-- A BC object for the gyrokinetic species that applies twist-shift BCs up to
-- a given radial location and sheath BCs after that. Inteded to simulate the
-- tokamak edge with open and closed field lines.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase      = require "App.BCs.BCsBase"
local DataStruct   = require "DataStruct"
local Updater      = require "Updater"
local Mpi          = require "Comm.Mpi"
local Proto        = require "Lib.Proto"
local Range        = require "Lib.Range"
local Grid         = require "Grid"
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local CartDecomp   = require "Lib.CartDecomp"
local Grid         = require "Grid"
local DiagsApp     = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags      = require "App.Diagnostics.GkDiagnostics"
local xsys         = require "xsys"
local lume         = require "Lib.lume"
local ffi          = require "ffi"
local GyrokineticModDecl = require "Eq.gkData.GyrokineticModDecl"

local sizeof = xsys.from(ffi, "sizeof")

local TokamakEdgeBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function TokamakEdgeBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function TokamakEdgeBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.xLCFS = assert(tbl.xLCFS, "TokamakEdgeBC: must specify x-location of the LCFS in 'xLCFS', and it must be at a cell boundary.")

   self.yShiftFuncIn = assert(tbl.shiftFunction, "TokamakEdgeBC: must provide the function that computes the y-shift in 'shiftFunction'.")

   self.yShiftPolyOrder = tbl.shiftPolyOrder

   self.bcKind = "sheath"
   self.phiWallFunc = tbl.phiWall
   if self.phiWallFunc then assert(type(self.phiWallFunc)=="function", "TokamakEdgeBC: phiWall must be a function (t, xn).") end

   self.evolve = false

   self.saveFlux = tbl.saveFlux or false
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
         self.saveFlux       = true
      end
   end
end

function TokamakEdgeBC:setName(nm) self.name = self.speciesName.."_"..nm end

function TokamakEdgeBC:bcReflect(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "flip".
   self.basis:flipSign(dir, fIn, fOut)
   local vparDir = self.cdim+1
   self.basis:flipSign(vparDir, fOut, fOut)
end
function TokamakEdgeBC:calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_, idx, f, fRefl)
   self.phi:fill(self.phiIdxr(idx), self.phiPtr)
   self.phiWallFld:fill(self.phiWallFldIdxr(idx), self.phiWallFldPtr)
   return self._calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, q_, m_,
                                     self.phiPtr:data(), self.phiWallFldPtr:data(), f:data(), fRefl:data())
end
function TokamakEdgeBC:bcSheath(dir, tm, idxIn, fIn, fOut)
   -- Note that GK reflection only valid in z-vpar.
   -- This is checked when bc is created.

   -- Need to figure out if we are on lower or upper domain edge
   local edgeVal
   local globalRange = self.grid:globalRange()
   if idxIn[dir] == globalRange:lower(dir) then
      -- This means we are at lower domain edge,
      -- so we need to evaluate basis functions at z=-1.
      edgeVal = -1
   else
      -- This means we are at upper domain edge
      -- so we need to evaluate basis functions at z=1.
      edgeVal = 1
   end
   -- Get vpar limits of cell.
   local vpardir = self.cdim+1
   local grid    = self.grid
   grid:setIndex(idxIn)
   local vL = grid:cellLowerInDir(vpardir)
   local vR = grid:cellUpperInDir(vpardir)
   local vlowerSq, vupperSq
   -- This makes it so that we only need to deal with absolute values of vpar.
   if math.abs(vR)>=math.abs(vL) then
      vlowerSq = vL*vL
      vupperSq = vR*vR
   else
      vlowerSq = vR*vR
      vupperSq = vL*vL
   end
   local w  = grid:cellCenterInDir(vpardir)
   local dv = grid:dx(vpardir)
   -- calculate reflected distribution function fhat
   -- note: reflected distribution can be
   -- 1) fhat=0 (no reflection, i.e. absorb),
   -- 2) fhat=f (full reflection)
   -- 3) fhat=c*f (partial reflection)
   self:calcSheathReflection(w, dv, vlowerSq, vupperSq, edgeVal, self.charge, self.mass, idxIn, fIn, self.fhatSheath)
   -- reflect fhat into skin cells
   self:bcReflect(dir, tm, nil, self.fhatSheath, fOut)
end

local function getGhostRange(edge, dir, global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if edge == "lower" then
      uv[dir] = global:lower(dir)-1   -- For ghost cells on "left".
   else
      lv[dir] = global:upper(dir)+1   -- For ghost cells on "right".
   end
   return Range.Range(lv, uv)
end
local function localExtRangeInDir(fIn, dir)
   local localrng = fIn:localRange()
   return localrng:extendDir(dir, fIn:lowerGhost(), fIn:upperGhost())
end
local function getSkinRange(edge, dir, global)
   local lv, uv = global:lowerAsVec(), global:upperAsVec()
   if edge == "lower" then
      uv[dir] = global:lower(dir)
   else
      lv[dir] = global:upper(dir)
   end
   return Range.Range(lv, uv)
end

function TokamakEdgeBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   -- Sheath BCs use phi and phiWall.
   self.setPhiWall = {advance = function(tCurr,inFlds,OutFlds) end}
   self.getPhi     = function(fieldIn, inIdx) return nil end

   assert(self.bcDir==self.cdim, "TokamakEdgeBC: sheath BC can only be used along the last/parallel configuration space dimension.")

   self.charge, self.mass = mySpecies.charge, mySpecies.mass

   self.fhatSheath = Lin.Vec(self.basis:numBasis())

   self.getPhi = function(fieldIn, inIdx) return fieldIn:rkStepperFields()[inIdx].phi end
   -- Pre-create pointer and indexer for phi potential.
   local phi = field:rkStepperFields()[1].phi
   self.phiPtr, self.phiIdxr = phi:get(1), phi:genIndexer()

   -- Create field and function for calculating wall potential according to user-provided function.
   self.phiWallFld = DataStruct.Field {
      onGrid   = field.grid,  numComponents    = field.basis:numBasis(),
      ghost    = {1,1},       syncPeriodicDirs = false,
      metaData = {polyOrder = field.basis:polyOrder(),
                  basisType = field.basis:id()},
   }
   self.phiWallFld:clear(0.)
   self.phiWallFldPtr, self.phiWallFldIdxr = self.phiWallFld:get(1), self.phiWallFld:genIndexer()
   if self.phiWallFunc then
      self.setPhiWall = Updater.EvalOnNodes {
         onGrid = field.grid,   evaluate = self.phiWallFunc,
         basis  = field.basis,  onGhosts = true,
      }
      self.setPhiWall:advance(0.0, {}, {self.phiWallFld})
   end
   self.phiWallFld:sync(false)

   self._calcSheathReflection = GyrokineticModDecl.selectSheathReflection(self.basis:id(), self.cdim,
                                                                          self.vdim, self.basis:polyOrder())

   -- Create ghost range for the sheath BC.
   local distf             = mySpecies["getDistF"] and mySpecies:getDistF() or mySpecies:getMoments()
   local global, globalExt = distf:globalRange(), distf:globalExtRange()
   local localExt          = distf:localExtRange()
   local ghostGlobal       = getGhostRange(self.bcEdge, self.bcDir, global, globalExt)
   -- Reduce this ghost range to the part of the x-domain with sheath BCs.
   -- Assume the split happens at a cell boundary and within the domain.
   assert(self.grid:lower(1)<self.xLCFS and self.xLCFS<self.grid:upper(1), "TokamakEdgeBC: 'xLCFS' coordinate must be within the x-domain.") 
   local needint = (self.xLCFS-self.grid:lower(1))/self.grid:dx(1)
   assert(math.floor(math.abs(needint-math.floor(needint))) < 1., "TokamakEdgeBC: 'xLCFS' must fall on a cell boundary along x.")
   -- Determine the index of the cell that abuts xLCFS from below.
   local coordLCFS = {self.xLCFS-1.e-7}
   self.idxLCFS    = {-9}
   local xGridIngr = self.grid:childGrid({1})
   local xGrid = Grid.RectCart {
      lower = xGridIngr.lower,  periodicDirs  = xGridIngr.periodicDirs,
      upper = xGridIngr.upper,  decomposition = xGridIngr.decomposition,
      cells = xGridIngr.cells,
   }
   xGrid:findCell(coordLCFS, self.idxLCFS) 
   local ghostGlobalSOL = ghostGlobal:shortenFromBelow(1, self.grid:numCells(1)-self.idxLCFS[1]+1)
   local localExtInDir  = localExtRangeInDir(distf, self.bcDir)
   local ghostLocalSOL  = localExtInDir:intersect(ghostGlobalSOL)

   -- BC updater for sheath BC.
   self.bcSolverSOL = Updater.Bc{
      onGrid   = self.grid,   edge               = self.bcEdge,
      cdim     = self.cdim,   boundaryConditions = {function(...) return self:bcSheath(...) end},
      dir      = self.bcDir,  vdir               = self.cdim+1,
      skinLoop = "flip",      confBasis          = self.confBasis,
      basis    = self.basis,
      advanceArgs = {{mySpecies:rkStepperFields()[1], ghostLocalSOL},
                     {mySpecies:rkStepperFields()[1]}},
   }

   -- Create the twist-shift BC solver for the closed field lines.
   if self.yShiftPolyOrder == nil then self.yShiftPolyOrder = self.basis:polyOrder() end

   if self.bcEdge=="lower" then
      self.yShiftFunc = function(t,xn) return -self.yShiftFuncIn(t,xn) end
   else
      self.yShiftFunc = self.yShiftFuncIn
   end

   local Ny = self.grid:numCells(2)
   if Ny == 1 then  -- Axisymmetric case. Assume periodicity in the core.
      self.zSync = function(fIn, bufferOut) self:zSyncCore(fIn) end
      self.bcSolverCore = { advance = function(...) end }
   else
      self.zSync = function(fIn, bufferOut) self:zSyncGlobalY(fIn, bufferOut) end
      self.bcSolverCore = Updater.TwistShiftBC {
         onGrid    = self.grid,       yShiftFunc      = self.yShiftFunc,
         basis     = self.basis,      yShiftPolyOrder = self.yShiftPolyOrder,
         confBasis = self.confBasis,  edge            = self.bcEdge,
      }
   end
   self.ghostGlobalCore = ghostGlobal:shorten(1, self.idxLCFS[1]+1)
   self.ghostLocalCore  = localExtInDir:intersect(self.ghostGlobalCore)

   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   self:createBoundaryGrid()

   -- Create a boundary grid that is global in Y: used to broadcast skin data.
   self:createBoundaryGridGlobalY()

   -- Need to define methods to allocate fields defined on boundary grid (used by diagnostics).
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

   -- Create communicator and MPI derived data types to send info along z.
   if Ny > 1 then
      self:createGraphComm()
      self:createSyncMPIdataTypes(distf)
   else
      self:createPeriodicGraphComm(distf)
      self:createPeriodicMPIdataTypes(distf)
   end

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

function TokamakEdgeBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function TokamakEdgeBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function TokamakEdgeBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function TokamakEdgeBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function TokamakEdgeBC:createDiagnostics(mySpecies, field)
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

function TokamakEdgeBC:createBoundaryGridGlobalY()
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

function TokamakEdgeBC:createGraphComm()
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
   local srcNum, destNum = 0, 0
   local isRecvRank, isSendRank
   if xyRecvRange:contains(myCutsIdx) then  -- This is a receiving rank.
      srcNum = yRecvRange:volume()
      isRecvRank = true
   end
   if xySendRange:contains(myCutsIdx) then -- This is a sending rank.
      destNum = yRecvRange:volume()
      isSendRank = true
   end

   if isRecvRank or isSendRank then  -- Only sending/receiving ranks participate.

      -- List ranks along z boundaries starting with the lower-z boundary.
      local zSkinRanksNum = decompRange:cuts(self.bcDir) > 1 and xyRecvRange:volume()+xySendRange:volume()
                                                              or xyRecvRange:volume()
      local zSkinGroupRanks = Lin.IntVec(zSkinRanksNum)
      local zLoRange = self.bcEdge == "lower" and xyRecvRange or xySendRange
      local zUpRange = self.bcEdge == "upper" and xyRecvRange or xySendRange
      local i = 0
      for idx in zLoRange:colMajorIter() do
         i = i+1;  zSkinGroupRanks[i] = cutsIdxr(idx)-1
      end
      if decompRange:cuts(self.bcDir) > 1 then
         for idx in zUpRange:colMajorIter() do
            i = i+1;  zSkinGroupRanks[i] = cutsIdxr(idx)-1
         end
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
         if isSendRank then
            destIDs[j] = cutsIdxr(idxRecv)-1
         end
         if isRecvRank then
            srcIDs[j]  = cutsIdxr(idxSend)-1
         end
      end

      -- Translate these ranks from the nodeComm to the zSkinComm.
      local destRanks, srcRanks
      if isSendRank then
         destRanks = Mpi.Group_translate_ranks(nodeGroup, destIDs, self.zSkinGroup)
      end
      if isRecvRank then
         srcRanks  = Mpi.Group_translate_ranks(nodeGroup, srcIDs, self.zSkinGroup)
      end
      if srcRanks==nil  then srcRanks  = Lin.IntVec(srcNum)  end -- Needed even if srcNum=0.
      if destRanks==nil then destRanks = Lin.IntVec(destNum) end -- Needed even if destNum=0.

      local reorder = 0
      -- Create a group/comm with only the lower and upper in z ranks.
      self.graphComm = Mpi.Dist_graph_create_adjacent(self.zSkinComm, srcNum, srcRanks, Mpi.UNWEIGHTED,
                                                      destNum, destRanks, Mpi.UNWEIGHTED, Mpi.INFO_NULL, reorder)
   else
      self.zSkinComm, self.graphComm = Mpi.COMM_NULL, Mpi.COMM_NULL
   end
end

function TokamakEdgeBC:createZskinComm()
   -- Create a group/communicator of z-boundary ranks only.
   local decompRange = self.grid:decomposedRange()

   local ones, upper = {}, {}
   for d=1,decompRange:ndim() do ones[d], upper[d] = 1, decompRange:cuts(d) end
   local cutsRange = Range.Range(ones, upper)

   local myId          = self.grid:subGridId() -- Grid ID on this processor.
   local zSkinRanksNum = cutsRange:shorten(self.bcDir):volume()
   if decompRange:cuts(self.bcDir) > 1 then zSkinRanksNum = zSkinRanksNum*2 end
   local zSkinGroupRanks = Lin.IntVec(zSkinRanksNum)
   local skelIds = decompRange:boundarySubDomainIds(self.bcDir)
   local c, isOnZboundary = 0, false
   for i = 1, #skelIds do
      local loId = skelIds[i].lower
      c = c+1;  zSkinGroupRanks[c] = loId-1
      if myId == loId then isOnZboundary = true end
   end
   if decompRange:cuts(self.bcDir) > 1 then
      for i = 1, #skelIds do
         local upId = skelIds[i].upper
         c = c+1;  zSkinGroupRanks[c] = upId-1
         if myId == upId then isOnZboundary = true end
      end
   end
   if isOnZboundary then
      local str = ""
      for i=1,#zSkinGroupRanks do str=str..tostring(zSkinGroupRanks[i]).."," end
      local nodeComm  = self.grid:commSet().comm
      local nodeGroup = Mpi.Comm_group(nodeComm)
      self.zSkinGroup = Mpi.Group_incl(nodeGroup, #zSkinGroupRanks, zSkinGroupRanks:data());
      self.zSkinComm  = Mpi.Comm_create_group(nodeComm, self.zSkinGroup, self.bcEdge == "lower" and 0 or 1);
   else
      self.zSkinComm = Mpi.COMM_NULL
   end
end

function TokamakEdgeBC:createSyncMPIdataTypes(fIn)
   -- Create the MPI derived data (MPIDD) types to sync the field along z.
   -- The idea is to use MPI_Neighbor_alltoallw with MPIDDs. The sending MPIDD
   -- will be the same for all ranks and will be based on the localRange of the
   -- simulation, as they are all just sending their local data. The receiving
   -- MPIDD has similar dimensions in x-y, but has different offsets because the
   -- data will be placed in a buffer with 1 cell in z.
   -- Futhermore, we'll use the same receiving type for all data received along y
   -- but we'll use displacements to put the data in the correct place in the
   -- global-along-y buffer.

   if not Mpi.Is_comm_valid(self.graphComm) then return end

   -- Create the receiving MPIDD, which will place the z-skin (including x-ghost
   -- cells for correct data layout) from all ranks along y into a single buffer
   -- that is global in Y.

   local indexer       = self.boundaryFieldGlobalY:genIndexer()
   local numComponents = self.boundaryFieldGlobalY:numComponents()
   local localExtRange = self.boundaryFieldGlobalY:localExtRange()
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
   isRecvRank = srcNum > 0 and true or false
   isSendRank = destNum > 0 and true or false

   self.recvDataType, self.recvLoc = Mpi.MPI_Datatype_vec(srcNum), 0
   for i = 1, srcNum do self.recvDataType[i-1] = elctCommType end

   -- The receive displacements ensure data is put in the place along y.
   self.recvDispl = Mpi.MPI_Aint_vec(srcNum)

   local skelIds = decompRange:boundarySubDomainIds(self.bcDir)
   local idxStart = Lin.IntVec(self.boundaryFieldGlobalY:ndim())
   for d = 1,self.boundaryFieldGlobalY:ndim() do idxStart[d] = -9 end
   for i = 1, #skelIds do
      local cId = skelIds[i].lower -- Should be = to skelIds[i].upper.

      -- Only create if we are on proper ranks.
      -- Note that if the node communicator has rank size of 1, then we can access all the
      -- memory needed for periodic boundary conditions and we do not need MPI Datatypes.
      if myId == cId and self.bcEdge == "lower" and isRecvRank then
         local rgnRecv = decompRange:subDomain(yLoId):lowerSkin(self.bcDir, self.boundaryField:lowerGhost()):extendDir(1,1,1)
         local idx     = rgnRecv:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.recvLoc  = (indexer(idx)-1)*numComponents
         for i = 1, srcNum do
            self.recvDataType[i-1] = Mpi.createDataTypeFromRangeAndSubRange(
               rgnRecv, localExtRange, numComponents, layout, elctCommType)[0]
            idx:copyInto(idxStart)
            idxStart[2] = idx[2]+(i-1)*decompRange:subDomain(yLoId):shape(2)
            self.recvDispl[i-1] = (indexer(idxStart)-1)*numComponents*sizeof(elctCommType)
         end
      end
      if myId == cId and self.bcEdge == "upper" and isRecvRank then
         local rgnRecv = decompRange:subDomain(yLoId):upperSkin(self.bcDir, self.boundaryField:upperGhost()):extendDir(1,1,1)
         local idx     = rgnRecv:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.recvLoc  = (indexer(idx)-1)*numComponents
         local str = ""
         for i = 1, srcNum do
            self.recvDataType[i-1] = Mpi.createDataTypeFromRangeAndSubRange(
               rgnRecv, localExtRange, numComponents, layout, elctCommType)[0]
            idx:copyInto(idxStart)
            idxStart[2] = idx[2]+(i-1)*decompRange:subDomain(yLoId):shape(2)
            self.recvDispl[i-1] = (indexer(idxStart)-1)*numComponents*sizeof(elctCommType)
            str = str .. tostring((indexer(idxStart)-1)) .. ","
         end
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

   self.sendDataType, self.sendLoc = Mpi.MPI_Datatype_vec(destNum), 0
   for i = 1, destNum do self.sendDataType[i-1] = elctCommType end

   for i = 1, #skelIds do
      local loId, upId = skelIds[i].lower, skelIds[i].upper

      -- Only create if we are on proper ranks.
      -- Note that if the node communicator has rank size of 1, then we can access all the
      -- memory needed for periodic boundary conditions and we do not need MPI Datatypes.
      if myId == loId and self.bcEdge == "upper" and isSendRank then
         local rgnSend = decomposedRange:subDomain(loId):lowerSkin(self.bcDir, fIn:upperGhost()):extendDir(1,1,1)
         local idx = rgnSend:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.sendLoc      = (indexer(idx)-1)*numComponents
         for i = 1, destNum do
            self.sendDataType[i-1] = Mpi.createDataTypeFromRangeAndSubRange(
               rgnSend, localExtRange, numComponents, layout, elctCommType)[0]
         end
      end
      if myId == upId and self.bcEdge == "lower" and isSendRank then
         local rgnSend = decomposedRange:subDomain(upId):upperSkin(self.bcDir, fIn:lowerGhost()):extendDir(1,1,1)
         local idx = rgnSend:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.sendLoc      = (indexer(idx)-1)*numComponents
         for i = 1, destNum do
            self.sendDataType[i-1] = Mpi.createDataTypeFromRangeAndSubRange(
               rgnSend, localExtRange, numComponents, layout, elctCommType)[0]
         end
      end
   end

   -- The send displacements are zero because the sendLoc will
   -- place the data pointer in the desired location.
   self.sendDispl = Mpi.MPI_Aint_vec(destNum)
   for i = 1, destNum do self.sendDispl[i-1] = 0 end

   -- MPI alltoall requires an array with the number of elements sent
   -- to/received from each rank. Set them to 1 as we use MPIDDs.
   self.recvCount, self.sendCount = Lin.IntVec(srcNum), Lin.IntVec(destNum)
   for i = 1, srcNum  do self.recvCount[i] = 1 end
   for i = 1, destNum do self.sendCount[i] = 1 end
end

function TokamakEdgeBC:createPeriodicGraphComm(fIn)
   -- Create a graph communicator which connects all the upper ranks in z (sources)
   -- with their corresponding lower rank in z (destination) if bcEdge=lower, or
   -- the lower ranks in z (sources) with their corresponding upper rank in z (destination)
   -- for bcEdge=upper, assuming periodicity in z.

   if self.grid:ndim() < 3 then return end

   local xDir, yDir = 1, 2

   -- Create a group/communicator of z-boundary ranks only.
   self:createZskinComm()

   -- Create graph communicator. Only ranks w/ core (twist-shift) BCs will send/receive.
   if Mpi.Is_comm_valid(self.zSkinComm) then

      -- Identify the corresponding rank at the other z boundary.
      local decompRange, myId = self.grid:decomposedRange(), self.grid:subGridId()
      local skelIds           = decompRange:boundarySubDomainIds(self.bcDir)
      local zOppositeRank
      for i = 1, #skelIds do
         local loId, upId = skelIds[i].lower, skelIds[i].upper
         if myId == loId then zOppositeRank = upId end
         if myId == upId then zOppositeRank = loId end
      end
      
      local global, globalExt = fIn:globalRange(), fIn:globalExtRange()
      local localExtInDir     = localExtRangeInDir(fIn, self.bcDir)
      local ghostGlobalOppositeZ = getGhostRange(self.bcEdge=="lower" and "upper" or "lower",
                                                 self.bcDir, global, globalExt)
      local ghostGlobalCoreOppositeZ = ghostGlobalOppositeZ:shorten(xDir, self.idxLCFS[xDir]+1)
      local ghostLocalCoreOppositeZ  = localExtInDir:intersect(ghostGlobalCoreOppositeZ)

      self.hasCoreBC          = self.ghostLocalCore:volume() > 0
      self.hasCoreBCOppositeZ = ghostLocalCoreOppositeZ:volume() > 0 -- =true if this is a core rank on the opposite edge.
        
      local srcNum, destNum
      local srcIDs   -- Source rank IDs along y (for this x).
      local destIDs  -- destination rank IDs along y (for each this x).
      if self.hasCoreBC or self.hasCoreBCOppositeZ then
         srcNum  = self.hasCoreBC and 1 or 0
         destNum = self.hasCoreBC and 0 or 1
         srcIDs  = Lin.IntVec(srcNum)  -- Source rank IDs along y (for this x).
         destIDs = Lin.IntVec(destNum) -- destination rank IDs along y (for each this x).

         if self.hasCoreBC then srcIDs[1] = zOppositeRank-1 end
         if self.hasCoreBCOppositeZ then destIDs[1] = zOppositeRank-1 end
      else
         -- These ranks do not need to enforce periodicity in z, they are SOL-only ranks.
         srcNum, destNum = 0, 0
         srcIDs  = Lin.IntVec(srcNum)  -- Source rank IDs along y (for this x).
         destIDs = Lin.IntVec(destNum) -- destination rank IDs along y (for each this x).
      end

      local nodeGroup = Mpi.Comm_group(self.grid:commSet().comm)

      -- Translate these ranks from the nodeComm to the zSkinComm.
      local destRanks, srcRanks
      if destNum > 0 then
         destRanks = Mpi.Group_translate_ranks(nodeGroup, destIDs, self.zSkinGroup)
      end
      if srcNum > 0 then
         srcRanks  = Mpi.Group_translate_ranks(nodeGroup, srcIDs, self.zSkinGroup)
      end
      if srcRanks==nil  then srcRanks  = Lin.IntVec(srcNum)  end -- Needed even if srcNum=0.
      if destRanks==nil then destRanks = Lin.IntVec(destNum) end -- Needed even if destNum=0.

      local reorder = 0
      -- Create a group/comm with only the lower and upper in z ranks.
      self.graphComm = Mpi.Dist_graph_create_adjacent(self.zSkinComm, srcNum, srcRanks, Mpi.UNWEIGHTED,
                                                      destNum, destRanks, Mpi.UNWEIGHTED, Mpi.INFO_NULL, reorder)
   else
      self.graphComm = Mpi.COMM_NULL
   end
end

function TokamakEdgeBC:createPeriodicMPIdataTypes(fIn)
   -- Create the MPI derived data (MPIDD) types to sync the field along z
   -- for the (axisymmetric) case in which the core is periodic.
   -- These are similar to the CartField MPIDDs but restricted in x according
   -- to the location of the LCFS.

   if not Mpi.Is_comm_valid(self.graphComm) then return end

   local srcNum, destNum, _ = Mpi.Dist_graph_neighbors_count(self.graphComm)
   local isRecvRank = srcNum > 0 and true or false
   local isSendRank = destNum > 0 and true or false

   local myId        = self.grid:subGridId() -- Whole grid ID on this processor.
   local decompRange = self.grid:decomposedRange()
   local skelIds     = decompRange:boundarySubDomainIds(self.bcDir)
   local indexer                      = fIn:genIndexer()
   local numComponents, localExtRange = fIn:numComponents(), fIn:localExtRange()
   local layout, elctCommType         = fIn:layout(), fIn:elemTypeMpi()

   -- Create sending MPIDDs. This is modeled after code in CartField.
   self.sendDataType, self.sendLoc = Mpi.MPI_Datatype_vec(destNum), 0
   for i = 1, destNum do self.sendDataType[i-1] = elctCommType end

   -- The send displacements are zero because the sendLoc will
   -- place the data pointer in the desired location.
   self.sendDispl = Mpi.MPI_Aint_vec(destNum)
   for i = 1, destNum do self.sendDispl[i-1] = 0 end

   local skinGlobalOppositeZ = getSkinRange(self.bcEdge=="lower" and "upper" or "lower", self.bcDir, fIn:globalRange())
   local skinGlobalCoreOppositeZ = skinGlobalOppositeZ:shorten(1, self.idxLCFS[1])

   for i = 1, #skelIds do
      local loId, upId = skelIds[i].lower, skelIds[i].upper
      if myId == loId and self.bcEdge == "upper" then
         local rgnSend = decompRange:subDomain(loId):lowerSkin(self.bcDir, fIn:upperGhost())
         if self.hasCoreBCOppositeZ then rgnSend = skinGlobalCoreOppositeZ:intersect(rgnSend) end 
         local idx = rgnSend:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.sendLoc = (indexer(idx)-1)*numComponents
         for i = 1, destNum do 
            self.sendDataType[i-1] = Mpi.createDataTypeFromRangeAndSubRange(
               rgnSend, localExtRange, numComponents, layout, elctCommType)[0]
         end
      end
      if myId == upId and self.bcEdge == "lower" then
         local rgnSend = decompRange:subDomain(upId):upperSkin(self.bcDir, fIn:lowerGhost())
         if self.hasCoreBCOppositeZ then rgnSend = skinGlobalCoreOppositeZ:intersect(rgnSend) end 
         local idx = rgnSend:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.sendLoc = (indexer(idx)-1)*numComponents
         for i = 1, destNum do 
            self.sendDataType[i-1] = Mpi.createDataTypeFromRangeAndSubRange(
               rgnSend, localExtRange, numComponents, layout, elctCommType)[0]
         end
      end
   end

   -- Create receving MPIDDs.
   self.recvDataType, self.recvLoc = Mpi.MPI_Datatype_vec(srcNum), 0
   for i = 1, srcNum do self.recvDataType[i-1] = elctCommType end

   -- The receive displacements ensure data is put in the place along y.
   self.recvDispl = Mpi.MPI_Aint_vec(srcNum)

   for i = 1, #skelIds do
      local loId, upId = skelIds[i].lower, skelIds[i].upper
      if myId == loId and self.bcEdge == "lower" then
         local rgnRecv = decompRange:subDomain(loId):lowerGhost(self.bcDir, fIn:lowerGhost())
         if self.hasCoreBC then rgnRecv = self.ghostGlobalCore:intersect(rgnRecv) end
         local idx = rgnRecv:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.recvLoc = (indexer(idx)-1)*numComponents
         for i = 1, srcNum do
            self.recvDispl[i-1] = self.recvLoc*sizeof(elctCommType)
            self.recvDataType[i-1] = Mpi.createDataTypeFromRangeAndSubRange(
               rgnRecv, localExtRange, numComponents, layout, elctCommType)[0]
         end
      end
      if myId == upId and self.bcEdge == "upper" then
         local rgnRecv = decompRange:subDomain(upId):upperGhost(self.bcDir, fIn:upperGhost())
         if self.hasCoreBC then rgnRecv = self.ghostGlobalCore:intersect(rgnRecv) end
         local idx = rgnRecv:lowerAsVec()
         -- Set idx to starting point of region you want to recv.
         self.recvLoc = (indexer(idx)-1)*numComponents
         for i = 1, srcNum do
            self.recvDispl[i-1] = self.recvLoc*sizeof(elctCommType)
            self.recvDataType[i-1] = Mpi.createDataTypeFromRangeAndSubRange(
               rgnRecv, localExtRange, numComponents, layout, elctCommType)[0]
         end
      end
   end

   -- MPI alltoall requires an array with the number of elements sent
   -- to/received from each rank. Set them to 1 as we use MPIDDs.
   self.recvCount, self.sendCount = Lin.IntVec(srcNum), Lin.IntVec(destNum)
   for i = 1, srcNum  do self.recvCount[i] = 1 end
   for i = 1, destNum do self.sendCount[i] = 1 end

end

function TokamakEdgeBC:zSyncGlobalY(fIn, bufferOut)
   -- bcEdge=lower: each upper-z rank sends their upper-z skin cell data in the fIn field
   --               to all lower-z ranks who receive them into a buffer with 1-cell in z (+ghosts).
   -- bcEdge=upper: each lower-z rank sends their lower-z skin cell data in the fIn field
   --               to all lower-z ranks who receive them into a buffer with 1-cell in z (+ghosts).

   if not Mpi.Is_comm_valid(self.graphComm) then return end

   Mpi.Neighbor_alltoallw(fIn:dataPointer()+self.sendLoc, self.sendCount:data(), self.sendDispl, self.sendDataType,
                          bufferOut:dataPointer(), self.recvCount:data(), self.recvDispl, self.recvDataType,
                          self.graphComm)
end

function TokamakEdgeBC:zSyncCore(fIn)
   -- Sync along z assuming periodicity in the core only.

   if not Mpi.Is_comm_valid(self.graphComm) then return end

   Mpi.Neighbor_alltoallw(fIn:dataPointer()+self.sendLoc, self.sendCount:data(), self.sendDispl, self.sendDataType,
                          fIn:dataPointer(), self.recvCount:data(), self.recvDispl, self.recvDataType,
                          self.graphComm)
end

-- These are needed to recycle the GkDiagnostics with TokamakEdgeBC.
function TokamakEdgeBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                                 self.boundaryFluxRate, self.boundaryFluxRate} end
function TokamakEdgeBC:getFlucF() return self.boundaryFluxRate end

function TokamakEdgeBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)

   local fIn = mySpecies:rkStepperFields()[outIdx]

   -- .......... Twist-shift BCs in the core (closed-flux region) ............... --
   self.zSync(fIn, self.boundaryFieldGlobalY)
   self.bcSolverCore:advance(tCurr, {self.boundaryFieldGlobalY, self.ghostRangeCore}, {fIn})

   -- .......... Sheath BCs in the SOL (open field-line region) ............... --
   self.setPhiWall:advance(tCurr, {}, {self.phiWallFld}) -- Compute wall potential if needed (i.e. sheath BC).
   self.phi = self.getPhi(field, inIdx)              -- If needed get the current plasma potential (for sheath BC).
   self.bcSolverSOL:advance(tCurr, {fIn}, {fIn})
end

function TokamakEdgeBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return TokamakEdgeBC