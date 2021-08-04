-- Gkyl ------------------------------------------------------------------------
--
-- Basic boundary condition for a gyrofluid species, i.e. those that can be
-- applied with Updater/Bc.lua.
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
local GyrofluidDiags = require "App.Diagnostics.GyrofluidDiagnostics"
local Constants      = require "Lib.Constants"

local GyrofluidBasicBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function GyrofluidBasicBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GyrofluidBasicBC:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.bcKind      = assert(tbl.kind, "GyrofluidBasicBC: must specify the type of BC in 'kind'.")
   self.diagnostics = tbl.diagnostics or {}
   self.saveFlux    = tbl.saveFlux or false

   self.saveFlux = tbl.saveFlux or false
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
         self.saveFlux       = true
      end
   end
end

function GyrofluidBasicBC:setName(nm) self.name = self.speciesName.."_"..nm end
function GyrofluidBasicBC:setConfBasis(basis) self.basis = basis end
function GyrofluidBasicBC:setConfGrid(grid) self.grid = grid end

function GyrofluidBasicBC:bcCopy(dir, tm, idxIn, fIn, fOut)
   for i = 1, self.nMoments*self.basis:numBasis() do fOut[i] = fIn[i] end
end

function GyrofluidBasicBC:bcAbsorb(dir, tm, idxIn, fIn, fOut)
   -- The idea is that by setting the plasma quantities to zero in the
   -- ghost cell nothing is transported into the domain, and whatever is transported
   -- out is lost. We can't set them to exactly zero or else the sound speed
   -- and drift velocity would diverge, so we set them to something small.
   local numB = self.basis:numBasis()
   for i = 1, numB do fOut[0*numB+i] = 1.e-10*fIn[0*numB+i] end   -- Mass density.
   for i = 1, numB do fOut[1*numB+i] = 0. end                     -- Momentum density.
   for i = 1, numB do fOut[2*numB+i] = 1.e-10*fIn[2*numB+i] end   -- Energy density.
   for i = 1, numB do fOut[3*numB+i] = 1.e-10*fIn[3*numB+i] end   -- Perpendicular pressure (divided by B).
end

function GyrofluidBasicBC:bcSheath(dir, tm, idxIn, fIn, fOut)
   local function evAtBoundary(ptrIn, cOff)
      -- Given the pointer to the function in the skin cell, evaluate it
      -- at the lower/upper boundary.
      if self.basis:polyOrder()==1 then
         return self.bcEdge=='lower'
            and (ptrIn[cOff+1]-math.sqrt(3)*ptrIn[cOff+2])/math.sqrt(2)
            or  (ptrIn[cOff+1]+math.sqrt(3)*ptrIn[cOff+2])/math.sqrt(2)
      end
   end

   local function zerothCoeff(valIn)
      -- Assuming a cell-wise constant function, take the value in the cell
      -- and return the 0th DG coefficient. This would normally be done in a kernel.
      if self.basis:polyOrder()==1 then
         return valIn*math.sqrt(2.)
      elseif self.basis:polyOrder()==2 then
         return valIn*2.
      end
   end

   local numB = self.basis:numBasis()

   local mJacM0_b = evAtBoundary(fIn,0)   -- m*n*Jac at the boundary.

   fOut[0*numB+1] = zerothCoeff(mJacM0_b)   -- Mass density (times Jacobian).
   for i = 2, numB do fOut[0*numB+i] = 0. end

   -- Potential, magnetic field and Jacobian at the boundary.
   self.phi:fill(self.indexer(idxIn), self.phiPtr)
   self.bmag:fill(self.indexer(idxIn), self.bmagPtr)
   self.jacob:fill(self.indexer(idxIn), self.jacobPtr)
   local phi_b   = evAtBoundary(self.phiPtr,0)
   local bmag_b  = evAtBoundary(self.bmagPtr,0)
   local jacob_b = evAtBoundary(self.jacobPtr,0)

   -- Temperature and sound speed at the boundary.
   self.primMomSelf:fill(self.indexer(idxIn), self.primMomSelfPtr)
   self.cSound:fill(self.indexer(idxIn), self.cSoundPtr)
   local Tpar_b  = evAtBoundary(self.primMomSelfPtr,1*numB)
   local Tperp_b = evAtBoundary(self.primMomSelfPtr,2*numB)
   local cs_b    = evAtBoundary(self.cSoundPtr,0)

   local upar_b
   if self.charge > 0. then
      -- upar_i = max(c_s, upar).
      upar_b = math.max(cs_b, math.abs(evAtBoundary(self.primMomSelfPtr,0)))
   elseif self.charge < 0. then
      -- upar_e = upar_i*exp(Lambda - max(0, e*phi/Te)).
      local Te_b = (Tpar_b+2.*Tperp_b)/3.
      local eV   = Constants.ELEMENTARY_CHARGE
      self.ionPrimMomSelf:fill(self.indexer(idxIn), self.ionPrimMomSelfPtr)
      local upari_b = math.max(cs_b, math.abs(evAtBoundary(self.ionPrimMomSelfPtr,0)))

      upar_b = upari_b*math.exp( self.Lambda - math.max(0., eV*phi_b/Te_b) )
   end
   upar_b = self.bcEdge=='lower' and -upar_b or upar_b
   fOut[1*numB+1] = zerothCoeff(mJacM0_b*upar_b)   -- Momentum density (times Jacobian).
   for i = 2, numB do fOut[1*numB+i] = 0. end
   
   local mJacM1_b = evAtBoundary(fOut,1*numB)   -- New (ghost cell) m*n*upar*Jac at the boundary.
   
   -- Perpendicular and parallel pressures (times Jacobian) at the boundary.
   local pPerpJac_b = (mJacM0_b/self.mass)*Tperp_b
   local pParJac_b  = (mJacM0_b/self.mass)*Tpar_b

   -- Kinetic energy density (times Jacobian) at the boundary.
   fOut[2*numB+1] = zerothCoeff(0.5*(pParJac_b + upar_b*mJacM1_b) + pPerpJac_b)
   for i = 2, numB do fOut[2*numB+i] = 0. end

   fOut[3*numB+1] = zerothCoeff(pPerpJac_b/bmag_b)    -- Jacobian*pperp/B.
   for i = 2, numB do fOut[3*numB+i] = 0. end
end

function GyrofluidBasicBC:createSolver(mySpecies, field, externalField)

   self.nMoments = mySpecies.nMoments

   -- Bohm sheath BC uses phi.
   self.getPhi = function(fieldIn, inIdx) return nil end 

   local bcFunc, skinType
   if self.bcKind == "copy" then
      bcFunc   = function(...) return self:bcCopy(...) end
      skinType = "pointwise"
   elseif self.bcKind == "absorb" then
      bcFunc   = function(...) return self:bcAbsorb(...) end
      skinType = "pointwise"
   elseif self.bcKind == "sheath" then

      self.charge, self.mass = mySpecies.charge, mySpecies.mass
      self.Lambda = math.log(mySpecies.ionMass/(2.*math.pi*self.mass)) -- Only used by electrons. 
      
      self.getPhi = function(fieldIn, inIdx) return fieldIn:rkStepperFields()[inIdx].phi end
      local phi   = field:rkStepperFields()[1].phi
      self.bmag   = externalField.geo.bmag
      self.jacob  = externalField.geo.jacobGeo
      if self.jacob==nil then  -- In order to support simulations without jacobGeo.
         local evOnNodes = Updater.EvalOnNodes {
            onGrid = self.grid,   evaluate = function(t, xn) return 1. end,
            basis  = self.basis,  onGhosts = true,
         }
         self.jacob = DataStruct.Field {
            onGrid        = self.grid,
            numComponents = self.basis:numBasis(),
            ghost         = {1, 1},
            metaData      = {polyOrder = self.basis:polyOrder(),
                             basisType = self.basis:id(),},
         }
         evOnNodes:advance(0., {}, {self.jacob})
      end
      -- Pre-create pointers and indexers.
      self.phiPtr, self.phiIdxr = phi:get(1), phi:genIndexer()
      self.indexer = phi:genIndexer()
      self.bmagPtr, self.jacobPtr = self.bmag:get(1), self.jacob:get(1)

      self.primMomSelf    = mySpecies.primMomSelf   -- Primitive moments: upar, Tpar, Tperp.
      self.primMomSelfPtr = self.primMomSelf:get(1)

      self.ionPrimMomSelf    = mySpecies.ionPrimMomSelf
      self.ionPrimMomSelfPtr = self.ionPrimMomSelf:get(1)

      self.cSound    = mySpecies.cSound   -- Sound speed.
      self.cSoundPtr = self.cSound:get(1)

      bcFunc   = function(...) return self:bcSheath(...) end
      skinType = "pointwise"
   else
      assert(false, "GyrofluidBasicBC: BC kind not recognized.")
   end

   self.bcSolver = Updater.Bc {
      onGrid             = self.grid,
      cdim               = self.grid:ndim(),
      dir                = self.bcDir,
      edge               = self.bcEdge,
      boundaryConditions = {bcFunc},
      skinLoop           = skinType,
   }

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then
      -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
      self:createBoundaryGrid()

      local moms = mySpecies:getMoments()
      -- Need to define methods to allocate fields defined on boundary grid (used by diagnostics).
      self.allocCartField = function(self, grid, nComp, ghosts, metaData)
         local f = DataStruct.Field {
            onGrid        = grid,   ghost    = ghosts,
            numComponents = nComp,  metaData = metaData,
         }
         f:clear(0.0)
         return f
      end
      self.allocMoment = function(self)
         return self:allocCartField(self.boundaryGrid, self.basis:numBasis(), 
                                    {moms:lowerGhost(),moms:upperGhost()}, moms:getMetaData())
      end
      self.allocVectorMoment = function(self, dim)
         return self:allocCartField(self.boundaryGrid, dim*self.basis:numBasis(), 
                                    {moms:lowerGhost(),moms:upperGhost()}, moms:getMetaData())
      end

      -- Allocate fields needed.
      self.boundaryFluxFields = {}  -- Fluxes through the boundary, into ghost region, from each RK stage.
      self.boundaryPtr        = {}
      self.momInIdxr          = moms:genIndexer()
      for i = 1, #mySpecies:rkStepperFields() do
         self.boundaryFluxFields[i] = self:allocMoment()
         self.boundaryPtr[i]        = self.boundaryFluxFields[i]:get(1)
      end
      self.boundaryFluxRate      = self:allocMoment()
      self.boundaryFluxFieldPrev = self:allocMoment()
      self.boundaryIdxr          = self.boundaryFluxFields[1]:genIndexer()

      self.idxOut = Lin.IntVec(self.grid:ndim())

      -- Create the range needed to loop over ghosts.
      local global, globalExt, localExtRange = moms:globalRange(), moms:globalExtRange(), moms:localExtRange()
      self.ghostRng = localExtRange:intersect(self:getGhostRange(global, globalExt))
      -- Decompose ghost region into threads.
      self.ghostRangeDecomp = LinearDecomp.LinearDecompRange{range=self.ghostRng, numSplit=self.grid:numSharedProcs()}
      self.tId              = self.grid:subGridSharedId() -- Local thread ID.

      self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut)
         local ptrOut = qOut:get(1)
         for idx in self.ghostRangeDecomp:rowMajorIter(self.tId) do
            idx:copyInto(self.idxOut)
            qOut:fill(self.momInIdxr(idx), ptrOut)

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
      if not self.anyDiagnostics then
         self.calcBoundaryFluxRateFunc = function(dtIn) end
      else
         self.calcBoundaryFluxRateFunc = function(dtIn)
            -- Compute boundary flux rate ~ (fGhost_new - fGhost_old)/dt.
            self.boundaryFluxRate:combine( 1.0/dtIn, self.boundaryFluxFields[1],
                                          -1.0/dtIn, self.boundaryFluxFieldPrev)
            self.boundaryFluxFieldPrev:copy(self.boundaryFluxFields[1])
         end
      end
   else
      self.storeBoundaryFluxFunc        = function(tCurr, rkIdx, qOut) end
      self.copyBoundaryFluxFieldFunc    = function(inIdx, outIdx) end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...) end
      self.calcBoundaryFluxRateFunc     = function(dtIn) end
   end
end

function GyrofluidBasicBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end
function GyrofluidBasicBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function GyrofluidBasicBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function GyrofluidBasicBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function GyrofluidBasicBC:createDiagnostics(mySpecies, field)
   -- Create BC diagnostics.
   self.diagnostics = nil
   if self.tbl.diagnostics then
      self.diagnostics = DiagsApp{implementation = GyrofluidDiags()}
      self.diagnostics:fullInit(mySpecies, field, self)
      -- Presently boundary diagnostics are boundary flux diagnostics. Append 'flux' to the diagnostic's
      -- name so files are named accordingly. Re-design this when non-flux diagnostics are implemented
      self.diagnostics.name = self.diagnostics.name..'_flux'
   end
   return self.diagnostics
end

function GyrofluidBasicBC:getNoJacMoments() return self.boundaryFluxRate end  -- Used by diagnostics.

function GyrofluidBasicBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   self.phi = self.getPhi(field, inIdx)   -- If needed get the current plasma potential (for sheath BC).

   local fIn = mySpecies:rkStepperFields()[outIdx]
   self.bcSolver:advance(tCurr, {fIn}, {fIn})
end

function GyrofluidBasicBC:getBoundaryFluxFields()
   return self.boundaryFluxFields
end

-- ................... Classes meant as aliases to simplify input files ...................... --
local GyrofluidAbsorbBC = Proto(GyrofluidBasicBC)
function GyrofluidAbsorbBC:fullInit(mySpecies)
   self.tbl.kind = "absorb"
   GyrofluidAbsorbBC.super.fullInit(self, mySpecies)
end

local GyrofluidCopyBC = Proto(GyrofluidBasicBC)
function GyrofluidCopyBC:fullInit(mySpecies)
   self.tbl.kind = "copy"
   GyrofluidCopyBC.super.fullInit(self, mySpecies)
end

local GyrofluidSheathBC = Proto(GyrofluidBasicBC)
function GyrofluidSheathBC:fullInit(mySpecies)
   self.tbl.kind = "sheath"
   GyrofluidSheathBC.super.fullInit(self, mySpecies)
end
-- ................... End of GyrofluidBasicBC alias classes .................... --

return {GyrofluidBasic  = GyrofluidBasicBC,
        GyrofluidAbsorb = GyrofluidAbsorbBC,
        GyrofluidCopy   = GyrofluidCopyBC,
        GyrofluidSheath = GyrofluidSheathBC,}
