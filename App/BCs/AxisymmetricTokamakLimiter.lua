-- Gkyl ------------------------------------------------------------------------
--
-- A BC object for the gyrokinetic species that applies periodic BCs up to
-- a given radial location and sheath BCs after that. Inteded to simulate the
-- limited tokamak edge with open and closed field lines.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase    = require "App.BCs.BCsBase"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"
local Proto      = require "Lib.Proto"
local Grid       = require "Grid"
local BCtools    = require "App.BCs.GkBCtools"
local DiagsApp   = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags    = require "App.Diagnostics.GkDiagnostics"
local xsys       = require "xsys"
local lume       = require "Lib.lume"
local ffi        = require "ffi"

local AxiTokLimBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function AxiTokLimBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function AxiTokLimBC:fullInit(mySpecies)
   local tbl = self.tbl -- Previously stored table.

   self.xLCFS = assert(tbl.xLCFS, "TokamakEdgeBC: must specify x-location of the LCFS in 'xLCFS', and it must be at a cell boundary.")

   self.phiWallIn = tbl.phiWall
   if self.phiWallIn then assert(type(self.phiWallIn)=="table" or type(self.phiWallIn)=="function", "GkSheathBC: phiWall must be a function or a table with 'spatialDependence', and optionally 'timeDependence', functions.") end

   self.evolve = xsys.pickBool(tbl.evolve, self.feedback or false)

   self.saveFlux = tbl.saveFlux or false
   self.anyDiagnostics = false
   if tbl.diagnostics then
      if #tbl.diagnostics>0 then
         self.anyDiagnostics = true
         self.saveFlux       = true
      end
   end
end

function AxiTokLimBC:setName(nm) self.name = self.speciesName.."_"..nm end

function AxiTokLimBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()

   local Ny = self.grid:numCells(2)
   assert(self.cdim==3 and Ny==1, "AxiTokLimBC: this BC is currently for 3D simulations with 1 cell along y.") 

   -- Create the boundary grid and other boundary tools.
   BCtools.createBoundaryTools(mySpecies, field, externalField, self)

   local charge, mass = mySpecies.charge, mySpecies.mass

   -- Create field and function for calculating wall potential according to user-provided function.
   self.phiWallFld = self:allocCartField(field.grid, field.basis:numBasis(), {1,1},
                                         {polyOrder=field.basis:polyOrder(), basisType=field.basis:id()},false)
   self.phiWallFld:clear(0.)
   if self.phiWallIn then
      local projPhiWall = Updater.EvalOnNodes {
         onGrid = field.grid,   evaluate = function(t,xn) return 1. end,
         basis  = field.basis,  onGhosts = true,
      }
      if type(self.phiWallIn) == 'function' then
         self.setPhiWall = projPhiWall
         self.setPhiWall:setFunc(self.phiWallIn)
         self.setPhiWall:advance(0.0, {}, {self.phiWallFld})
      elseif type(self.phiWallIn) == 'table' then
         projPhiWall:setFunc(function(t,xn) return self.phiWallIn["spatialDependence"](xn) end)
         if self.phiWallIn["timeDependence"] then
            self.phiWallFldInit = self:allocCartField(field.grid, field.basis:numBasis(), {1,1},
                                                      {polyOrder=field.basis:polyOrder(), basisType=field.basis:id()},false)
            projPhiWall:advance(self.phiWallIn["timeDependence"](0.), {}, {self.phiWallFldInit})
            self.setPhiWall = {advance = function(tCurr,inFlds,OutFlds)
               self.phiWallFld:combine(self.phiWallIn["timeDependence"](tCurr), self.phiWallFldInit)
            end}
         else
            -- No time dependence
            projPhiWall:advance(0., {}, {self.phiWallFld})
            self.setPhiWall = {advance = function(tCurr,inFlds,OutFlds) end}
         end
      end
   else
      self.setPhiWall = {advance = function(tCurr,inFlds,OutFlds) end}
   end
   self.phiWallFld:sync(false)

   -- Part of global ghost range this rank owns.
   local distf = mySpecies:getDistF()
   self.myGlobalGhostRange = self.bcEdge=="lower" and distf:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                   or distf:localGlobalGhostRangeIntersectUpper()[self.bcDir]
   self.myGlobalSkinRange = self.bcEdge=="lower" and distf:localGlobalSkinRangeIntersectLower()[self.bcDir]
                                                  or distf:localGlobalSkinRangeIntersectUpper()[self.bcDir]
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
   local myGlobalGhostRangeSOL = self.myGlobalGhostRange:shortenFromBelow(1, self.grid:numCells(1)-self.idxLCFS[1]+1)
   self.myGlobalGhostRangeSOL  = distf:localExtRange():subRange(myGlobalGhostRangeSOL:lowerAsVec(),myGlobalGhostRangeSOL:upperAsVec())
   local myGlobalSkinRangeSOL  = self.myGlobalSkinRange:shortenFromBelow(1, self.grid:numCells(1)-self.idxLCFS[1]+1)
   self.myGlobalSkinRangeSOL   = distf:localExtRange():subRange(myGlobalSkinRangeSOL:lowerAsVec(),myGlobalSkinRangeSOL:upperAsVec())

   self.bcSolverSOL = Updater.GkSheathBc{
      onGrid  = self.grid,   edge    = self.bcEdge,  
      cdim    = self.cdim,   basis   = self.basis,
      dir     = self.bcDir,  phiWall = self.phiWallFld,
      charge  = charge,      mass    = mass,
      skinRange = self.myGlobalSkinRangeSOL,  ghostRange = self.myGlobalGhostRangeSOL,
   }
   self.bcSolverAdvanceSOL = function(tm, inFlds, outFlds)
      self.bcSolverSOL:advance(tm, {inFlds[1]}, outFlds)
   end

   -- The saveFlux option is used for boundary diagnostics, or BCs that require
   -- the fluxes through a boundary (e.g. neutral recycling).
   if self.saveFlux then

      -- Create boundary tools for saving fluxes.
      BCtools.createFluxTools(mySpecies, field, externalField, self)

      if not self.anyDiagnostics then
         self.calcBoundaryFluxRateFunc = function(dtIn) end
      else
         -- Create boundary tools for diagnostics.
         BCtools.createDiagnosticTools(mySpecies, field, externalField, self)
      end
   else
      self.storeBoundaryFluxFunc        = function(tCurr, rkIdx, qOut) end
      self.copyBoundaryFluxFieldFunc    = function(inIdx, outIdx) end
      self.combineBoundaryFluxFieldFunc = function(outIdx, a, aIdx, ...) end
      self.calcBoundaryFluxRateFunc     = function(dtIn) end
   end
end

function AxiTokLimBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function AxiTokLimBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function AxiTokLimBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function AxiTokLimBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function AxiTokLimBC:createDiagnostics(mySpecies, field)
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

-- These are needed to recycle the GkDiagnostics with AxiTokLimBC.
function AxiTokLimBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                               self.boundaryFluxRate, self.boundaryFluxRate} end
function AxiTokLimBC:getFlucF() return self.boundaryFluxRate end

function AxiTokLimBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)

   local fIn = mySpecies:rkStepperFields()[outIdx] 

   -- .......... Periodic BCs in the core (closed-flux region) ............... --
   -- Apply periodic BCs in z direction, but only by the first edge (lower or upper, whichever is called first)
   -- First check if sync has been called on this step by any BC.
   local doSync = true
   for _, bc in lume.orderedIter(mySpecies.nonPeriodicBCs) do 
      if bc.synced then doSync = false end
   end
   if doSync then
      -- Fetch skin cells from opposite z-edge (donor) and copy to ghost cells (target).
      local fldGrid = fIn:grid()
      local periodicDirs = fldGrid:getPeriodicDirs()
      local modifiedPeriodicDirs = {3}
      fldGrid:setPeriodicDirs(modifiedPeriodicDirs)
      fIn:sync()
      fldGrid:setPeriodicDirs(periodicDirs)
      self.synced = true
   else
      -- Don't sync but reset synced flags for next step.
      for _, bc in lume.orderedIter(mySpecies.nonPeriodicBCs) do bc.synced = false end
   end

   -- .......... Sheath BCs in the SOL (open field-line region) ............... --
   self.setPhiWall:advance(tCurr, {}, {self.phiWallFld})  -- Compute wall potential if needed.
   local phi = field:rkStepperFields()[inIdx].phi         -- Get the current plasma potential.
   self.bcSolverAdvanceSOL(tCurr, {phi}, {fIn})
end

function AxiTokLimBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return AxiTokLimBC
