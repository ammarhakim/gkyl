-- Gkyl ------------------------------------------------------------------------
--
-- Maxwellian Boundary Condition for a Gyrokinetic species
-- Applied with Updater/MaxwellGhostBc
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase    = require "App.BCs.BCsBase"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"
local Mpi        = require "Comm.Mpi"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"
local Range      = require "Lib.Range"
local Lin        = require "Lib.Linalg"
local Grid       = require "Grid"
local DiagsApp   = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags    = require "App.Diagnostics.GkDiagnostics"
local xsys       = require "xsys"

local MaxwellianBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function MaxwellianBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function MaxwellianBC:fullInit(mySpecies)
   self.maxwellianKind = self.tbl.maxwellianKind
   self.densityGhost   = self.tbl.densityGhost
   self.uParGhost      = self.tbl.uParGhost
   self.tempGhost      = self.tbl.tempGhost
   self.initFunc       = self.tbl.initFunc
   self.fromFile       = self.tbl.fromFile
end

function MaxwellianBC:setName(nm) self.name = self.speciesName.."_"..nm end

function MaxwellianBC:createSolver(mySpecies, field, externalField)

   self.basis, self.grid = mySpecies.basis, mySpecies.grid
   self.ndim, self.cdim, self.vdim = self.grid:ndim(), self.confGrid:ndim(), self.grid:ndim()-self.confGrid:ndim()



   self.ghostFld = self:allocCartField(self.boundaryGrid, self.basis:numBasis(), {1,1}, distf:getMetaData())
   self.bcSolver = Updater.MaxwellGhostBc{
      edge               = self.bcEdge,
      dir      = self.bcDir,
      boundaryGrid = self.boundaryGrid,
      confBoundaryGrid = self.confBoundaryGrid,
      ghostFld         = self.ghostFld,
      --localGhostRangeWOcorners = self.localGhostRangeWOcorners,
      myGlobalGhostRange = self.myGlobalGhostRange,
   }
   self.bcSolverAdvance = function(tm, inFlds, outFlds)
      self.bcSolver:advance(tm, {inFlds[1]}, outFlds)
   end
   self.phaseFieldIo = AdiosCartFieldIo {
      elemType   = mySpecies:getDistF():elemType(),
      method     = "MPI",
      writeGhost = false,
      metaData   = {polyOrder = self.basis:polyOrder(),
                    basisType = self.basis:id(),
                    charge    = self.charge,
                    mass      = self.mass,},
   }
   if not self.fromFile then
      if self.maxwellianKind == 'local' then
         self.projMaxwell = Updater.MaxwellianOnBasis{
            onGrid     = self.boundaryGrid,     confBasis = self.confBasis,
            phaseBasis = self.basis ,           mass      = self.mass,
            confGrid   = self.confBoundaryGrid, onGhosts  = true,
         }
         local confProject = Updater.ProjectOnBasis {
            onGrid   = self.confBoundaryGrid,
            basis    = self.confBasis,
            evaluate = function(t, xn) return 0. end,   -- Set below.
            onGhosts = false
         }
         local numDens = self:allocMoment()
         local uPar    = self:allocMoment()
         local vtSq    = self:allocMoment()

         confProject:setFunc(function(t, xn) return 1. end)
         confProject:advance(time, {}, {numDens})
         confProject:setFunc(self.uParGhost)
         confProject:advance(time, {}, {uPar})
         confProject:setFunc(self.tempGhost)
         confProject:advance(time, {}, {vtSq})
         vtSq:scale(1./self.mass)

         self.projMaxwell:advance(time,{numDens,uPar,vtSq,self.bmag},{self.ghostFld})
         local M0e, M0 = self:allocMoment(), self:allocMoment()
         local M0mod   = self:allocMoment()
         self.numDensityCalc:advance(0.0, {self.ghostFld}, {M0})
         confProject:setFunc(self.densityGhost)
         confProject:advance(0.0, {}, {M0e})
         self.confWeakDivide:advance(0.0, {M0, M0e}, {M0mod})
         self.phaseWeakMultiply:advance(0.0, {M0mod,self.ghostFld}, {self.ghostFld})

      elseif self.maxwellianKind == 'canonical' then
         if mySpecies.jacobPhaseFunc and self.vdim > 1 then
            local initFuncWithoutJacobian = self.initFunc
            self.initFunc = function (t, xn)
               local xconf = {}
               for d = 1, self.cdim do xconf[d] = xn[d] end
               local J = mySpecies.jacobPhaseFunc(t,xconf)
               local f = initFuncWithoutJacobian(t,xn)
               return J*f
            end
         end
         local project = Updater.ProjectOnBasis {
            onGrid   = self.boundaryGrid,
            basis    = self.basis,
            evaluate = self.initFunc,
            onGhosts = false
         }
         project:advance(time, {},{self.ghostFld})
      end
   end
end

function MaxwellianBC:createCouplingSolver(species,field,externalField)
   if self.bcKind == 'maxwellianGhost' then
      if self.fromFile then
         local tm, fr = self.phaseFieldIo:read(self.ghostFld, self.fromFile)
      else
         if self.maxwellianKind == 'canonical' then
            local ionName = nil
            for nm, s in lume.orderedIter(species) do
               if 0.0 < s.charge then ionName = nm end
            end

            local numDens = self:allocMoment()
            if self.species.charge < 0.0 then
               local bcName = self.name:gsub(self.speciesName .. "_","")
               local numDensScaleTo = species[ionName]['nonPeriodicBCs'][bcName]:allocMoment()
               self.numDensityCalc:advance(0.0, {self.ghostFld}, {numDens})
               species[ionName]['nonPeriodicBCs'][bcName].numDensityCalc:advance(0.0,
                  {species[ionName]['nonPeriodicBCs'][bcName].ghostFld}, {numDensScaleTo})

               local M0mod = self:allocMoment()
               self.confWeakDivide:advance(0.0,{numDens, numDensScaleTo} , {M0mod})
               self.phaseWeakMultiply:advance(0.0, {M0mod,self.ghostFld}, {self.ghostFld})
            end
            self.numDensityCalc:advance(0.0, {self.ghostFld}, {numDens})
            numDens:write(string.format("%s_M0_%d.bp", self.name, 0), 0.0, 0, false)
         end
         -- Multiply by conf space Jacobian Now
         self.phaseWeakMultiply:advance(0, {self.ghostFld, self.jacobGeo}, {self.ghostFld})
         self.phaseFieldIo:write(self.ghostFld, string.format("%s_Maxwellian_%d.bp", self.name, 0), 0., 0, false)
      end
   end
end

function MaxwellianBC:createBoundaryTools(mySpecies,field,externalField)
   -- Create reduced boundary grid with 1 cell in dimension of self.bcDir.
   local distf, numDensity = mySpecies:getDistF(), mySpecies:getNumDensity()
   --Make the ghost range without corners
   --local global, globalExt, localExtRange = distf:globalRange(), distf:globalExtRange(), distf:localExtRange()
   --local globalGhostRange = self:getGhostRange(global, globalExt)
   --local lv, uv = globalGhostRange:lowerAsVec(), globalGhostRange:upperAsVec()
   --for d = 1,self.grid:ndim() do
   --   if d ~= self.bcDir then
   --      uv[d] = uv[d]-distf:upperGhost()
   --      lv[d] = lv[d]+distf:lowerGhost()
   --   end
   --end
   --local ghostRangeWOcorners = Range.Range(lv, uv)
   --self.localGhostRangeWOcorners = localExtRange:intersect(ghostRangeWOcorners)

   local globalGhostRange = self.bcEdge=="lower" and distf:localGhostRangeLower()[self.bcDir]
                                                  or distf:localGhostRangeUpper()[self.bcDir]
   self:createBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())
   -- Create reduced boundary config-space grid with 1 cell in dimension of self.bcDir.
   self:createConfBoundaryGrid(globalGhostRange, self.bcEdge=="lower" and distf:lowerGhostVec() or distf:upperGhostVec())

   -- Need to define methods to allocate fields defined on boundary grid (used by diagnostics).
   self.allocCartField = function(self, grid, nComp, ghosts, metaData)
      local f = DataStruct.Field {
         onGrid        = grid,   ghost    = ghosts,
         numComponents = nComp,  metaData = metaData,
      }
      f:clear(0.0)
      return f
   end

   self.allocDistf = function()
      return self:allocCartField(self.boundaryGrid, self.basis:numBasis(), {0,0}, distf:getMetaData())
   end

   self.allocMoment = function(self)
      return self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
   end

   self.allocVectorMoment = function(self, dim)
      return self:allocCartField(self.confBoundaryGrid, dim*self.confBasis:numBasis(), {0,0}, numDensity:getMetaData())
   end
   self.allocIntMoment = function(self, comp)
      local metaData = {charge = self.charge,  mass = self.mass,}
      local ncomp = comp or 1
      local f = DataStruct.DynVector{numComponents = ncomp,     writeRank = self.confBoundaryGrid:commSet().writeRank,
                                     metaData      = metaData,  comm      = self.confBoundaryGrid:commSet().comm,}
      return f
   end

   -- Allocate fields needed.
   self.boundaryFluxFields = {}  -- Fluxes through the boundary, into ghost region, from each RK stage.
   self.distfInIdxr        = distf:genIndexer()
   for i = 1, #mySpecies:rkStepperFields() do
      self.boundaryFluxFields[i] = self.allocDistf()
   end
   self.boundaryFluxRate      = self.allocDistf()
   self.boundaryFluxFieldPrev = self.allocDistf()

   -- Part of global ghost range this rank owns.
   self.myGlobalGhostRange = self.bcEdge=="lower" and distf:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                   or distf:localGlobalGhostRangeIntersectUpper()[self.bcDir]

   -- The following are needed to evaluate a conf-space CartField on the confBoundaryGrid.
   self.confBoundaryField = self:allocMoment()
   -- Range spanning ghost cells.
   self.myGlobalConfGhostRange = self.bcEdge=="lower" and numDensity:localGlobalGhostRangeIntersectLower()[self.bcDir]
                                                       or numDensity:localGlobalGhostRangeIntersectUpper()[self.bcDir]

   -- Evaluate the magnetic field and jacobGeo in the boundary (needed by diagnostics).
   local bmag = externalField.geo.bmag
   self.bmag = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, bmag:getMetaData())
   self.bmag:copy(self:evalOnConfBoundary(bmag))
   local bmagInvSq = externalField.geo.bmagInvSq
   self.bmagInvSq = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, bmagInvSq:getMetaData())
   self.bmagInvSq:copy(self:evalOnConfBoundary(bmagInvSq))
   local jacobGeo = externalField.geo.jacobGeo
   if jacobGeo then
      self.jacobGeo = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeo:getMetaData())
      self.jacobGeo:copy(self:evalOnConfBoundary(jacobGeo))
   end
   local jacobGeoInv = externalField.geo.jacobGeoInv
   if jacobGeoInv then
      self.jacobGeoInv = self:allocCartField(self.confBoundaryGrid, self.confBasis:numBasis(), {0,0}, jacobGeoInv:getMetaData())
      self.jacobGeoInv:copy(self:evalOnConfBoundary(jacobGeoInv))
   end

   -- Declare methods/functions needed for handling saved fluxes and needed by diagnostics.
   self.storeBoundaryFluxFunc = function(tCurr, rkIdx, qOut)
      self.boundaryFluxFields[rkIdx]:copyRangeToRange(qOut, self.boundaryFluxFields[rkIdx]:localRange(), self.myGlobalGhostRange)
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
         weakBasis = self.confBasis,  operation = "Multiply",
         onGhosts  = true,
      }
      self.confWeakDivide = Updater.CartFieldBinOp {
         weakBasis = self.confBasis,  operation = "Divide",
         onRange   = self.confBoundaryField:localRange(),  onGhosts = false,
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
end

function MaxwellianBC:storeBoundaryFlux(tCurr, rkIdx, qOut)
   self.storeBoundaryFluxFunc(tCurr, rkIdx, qOut)
end

function MaxwellianBC:copyBoundaryFluxField(inIdx, outIdx)
   self.copyBoundaryFluxFieldFunc(inIdx, outIdx)
end
function MaxwellianBC:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   self.combineBoundaryFluxFieldFunc(outIdx, a, aIdx, ...)
end
function MaxwellianBC:computeBoundaryFluxRate(dtIn)
   self.calcBoundaryFluxRateFunc(dtIn)
end

function MaxwellianBC:createDiagnostics(mySpecies, field)
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

-- These are needed to recycle the GkDiagnostics with MaxwellianBC.
function MaxwellianBC:rkStepperFields() return {self.boundaryFluxRate, self.boundaryFluxRate,
                                             self.boundaryFluxRate, self.boundaryFluxRate} end
function MaxwellianBC:getFlucF() return self.boundaryFluxRate end

function MaxwellianBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   local fIn = mySpecies:rkStepperFields()[outIdx]
   self.bcSolverAdvance(tCurr, {self.ghostFld}, {fIn})
end

function MaxwellianBC:getBoundaryFluxFields() return self.boundaryFluxFields end

return MaxwellianBC
