-- Gkyl ------------------------------------------------------------------------
--
-- VlasovOnCartGrid support code: Species object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis = require "Basis"
local BoundaryCondition = require "Updater.BoundaryCondition"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local LinearTrigger = require "Lib.LinearTrigger"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local Updater = require "Updater"
local date = require "Lib.date"
local xsys = require "xsys"

-- function to create basis functions
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   end   
end

-- function to check if moment name is correct
local function isMomentNameGood(nm)
   if nm == "M0" or nm == "M1i" or nm == "M2ij" or nm == "M2" or nm == "M3i" then
      return true
   end
   return false
end

-- Class to store species-specific info
local Species = Proto()

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function Species:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function Species:fullInit(vlasovTbl)
   local tbl = self.tbl -- previously store table
   
   self.name = "name"
   self.cfl =  0.1
   self.charge, self.mass = tbl.charge, tbl.mass
   self.lower, self.upper = tbl.lower, tbl.upper
   self.cells = tbl.cells
   self.vdim = #self.cells -- velocity dimensions
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default, evolve species
   self.confBasis = nil -- Will be set later

   assert(#self.lower == self.vdim, "'lower' must have " .. self.vdim .. " entries")
   assert(#self.upper == self.vdim, "'upper' must have " .. self.vdim .. " entries")

   self.decompCuts = {}
   -- parallel decomposition stuff
   if tbl.decompCuts then
      assert(self.vdim == #tbl.decompCuts, "decompCuts should have exactly " .. self.vdim .. " entries")
      self.decompCuts = tbl.decompCuts
   else
      -- if not specified, use 1 processor
      for d = 1, self.vdim do self.decompCuts[d] = 1 end
   end

   -- create triggers to write distribution functions and moments
   if tbl.nDistFuncFrame then
      self.distIoTrigger = LinearTrigger(0, vlasovTbl.tEnd, tbl.nDistFuncFrame)
   else
      self.distIoTrigger = LinearTrigger(0, vlasovTbl.tEnd, vlasovTbl.nFrame)
   end
   
   if tbl.nDiagnosticFrame then
      self.diagIoTrigger = LinearTrigger(0, vlasovTbl.tEnd, tbl.nDiagnosticFrame)
   else
      self.diagIoTrigger = LinearTrigger(0, vlasovTbl.tEnd, vlasovTbl.nFrame)
   end

   self.distIoFrame = 0 -- frame number for distrubution fucntion
   self.diagIoFrame = 0 -- frame number for diagnostcia

   -- read in which diagnostic moments to compute on output
   self.diagnosticMoments = { }
   if tbl.diagnosticMoments then
      for i, nm in ipairs(tbl.diagnosticMoments) do
	 if isMomentNameGood(nm) then
	    self.diagnosticMoments[i] = nm
	 else
	    assert(false, string.format("Moment %s not valid", nm))
	 end
      end
   end

   -- for storing integrated moments
   self.integratedMoments = DataStruct.DynVector { numComponents = 5 }

   -- store initial condition function
   self.initFunc = tbl.init
end

-- methods for species object
function Species:ndim()
   return self.vdim
end
function Species:setName(nm)
   self.name = nm
end
function Species:setCfl(cfl)
   self.cfl = cfl
end
function Species:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end
function Species:setConfBasis(basis)
   self.confBasis = basis
end
function Species:setConfGrid(cgrid)
   self.confGrid = cgrid
end

function Species:createGrid(cLo, cUp, cCells, cDecompCuts, cPeriodicDirs)
   self.cdim = #cCells
   self.ndim = self.cdim+self.vdim

   -- create decomposition
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, cDecompCuts[d]) end
   for d = 1, self.vdim do table.insert(decompCuts, self.decompCuts[d]) end
   self.decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      shared = false,
   }

   -- create computational domain
   local lower, upper, cells = {}, {}, {}
   for d = 1, self.cdim do
      table.insert(lower, cLo[d])
      table.insert(upper, cUp[d])
      table.insert(cells, cCells[d])
   end
   for d = 1, self.vdim do
      table.insert(lower, self.lower[d])
      table.insert(upper, self.upper[d])
      table.insert(cells, self.cells[d])
   end
   self.grid = Grid.RectCart {
      lower = lower,
      upper = upper,
      cells = cells,
      periodicDirs = cPeriodicDirs,
      decomposition = self.decomp,
   }
end

function Species:createBasis(nm, polyOrder)
   self.basis = createBasis(nm, self.ndim, polyOrder)
end

function Species:alloc(nFields)
   -- allocate fields needed in RK update
   self.distf = {}
   for i = 1, nFields do
      self.distf[i] = DataStruct.Field {
	 onGrid = self.grid,
	 numComponents = self.basis:numBasis(),
	 ghost = {1, 1}
      }
   end
   -- create Adios object for field I/O
   self.distIo = AdiosCartFieldIo {
      elemType = self.distf[1]:elemType(),
      method = self.ioMethod,
   }

   -- allocate field to store coupling moments (for use in
   -- charge calculations) 
   self.momDensity = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.confBasis:numBasis()*self.vdim,
      ghost = {1, 1}
   }
end

function Species:allocMomCouplingFields()
   -- only need currents for coupling to fields (returning a table
   -- with single entry, i.e. space to store currents)
   return {
      DataStruct.Field {
	 onGrid = self.confGrid,
	 numComponents = self.confBasis:numBasis()*self.vdim,
	 ghost = {1, 1}
      }
   }
end

function Species:createSolver(hasE, hasB)
   -- create updater to advance solution by one time-step
   self.vlasovSlvr = Updater.VlasovDisCont {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      charge = self.charge,
      mass = self.mass,
      cfl = self.cfl,
      hasElectricField = hasE,
      hasMagneticField = hasB,
   }
   -- create updater to compute momentum density (for use in current
   -- calculations)
   self.momDensityCalc = Updater.DistFuncMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
      moment = "M1i",
   }

   -- create updater to compute integrated moments
   self.intMomentCalc = Updater.DistFuncIntegratedMomentCalc {
      onGrid = self.grid,
      phaseBasis = self.basis,
      confBasis = self.confBasis,
   }
end

function Species:createDiagnostics()
   local numComp = {}
   numComp["M0"] = 1
   numComp["M1i"] = self.vdim
   numComp["M2ij"] = self.vdim*(self.vdim+1)/2
   numComp["M2"] = 1
   numComp["M3i"] = self.vdim
   
   self.diagnosticMomentFields = { }
   self.diagnosticMomentUpdaters = { } 
   -- allocate space to store moments and create moment updater
   for i, mom in ipairs(self.diagnosticMoments) do
      self.diagnosticMomentFields[i] = DataStruct.Field {
	 onGrid = self.confGrid,
	 numComponents = self.confBasis:numBasis()*numComp[mom],
	 ghost = {1, 1}
      }
      
      self.diagnosticMomentUpdaters[i] = Updater.DistFuncMomentCalc {
	 onGrid = self.grid,
	 phaseBasis = self.basis,
	 confBasis = self.confBasis,
	 moment = mom,
      }
   end
end

function Species:calcDiagnosticMoments()
   local numMoms = #self.diagnosticMoments
   for i = 1, numMoms do
      self.diagnosticMomentUpdaters[i]:advance(
	 0.0, 0.0, {self.distf[1]}, {self.diagnosticMomentFields[i]})
   end
end

function Species:initDist()
   local project = Updater.ProjectOnBasis {
      onGrid = self.grid,
      basis = self.basis,
      evaluate = self.initFunc
   }
   project:advance(0.0, 0.0, {}, {self.distf[1]})
   self:applyBc(0.0, 0.0, self.distf[1])
end

function Species:write(tm)
   if self.evolve then
      -- compute integrated diagnostics
      self.intMomentCalc:advance(tm, 0.0, { self.distf[1] }, { self.integratedMoments })
      
      -- only write stuff if triggered
      if self.distIoTrigger(tm) then
	 self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, self.distIoFrame), tm)
	 self.distIoFrame = self.distIoFrame+1
      end

      if self.diagIoTrigger(tm) then
	 -- compute moments and write them out
	 self:calcDiagnosticMoments()
	 for i, mom in ipairs(self.diagnosticMoments) do
	    -- should one use AdiosIo object for this?
	    self.diagnosticMomentFields[i]:write(
	       string.format("%s_%s_%d.bp", self.name, mom, self.diagIoFrame), tm)
	 end
	 self.integratedMoments:write(
	    string.format("%s_intMom_%d.bp", self.name, self.diagIoFrame), tm)
	 self.diagIoFrame = self.diagIoFrame+1
      end
   else
      -- if not evolving species, don't write anything except initial conditions
      if self.distIoFrame == 0 then
	 self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, 0), tm)

	 -- compute moments and write them out
	 self:calcDiagnosticMoments()
	 for i, mom in ipairs(self.diagnosticMoments) do
	    -- should one use AdiosIo object for this?
	    self.diagnosticMomentFields[i]:write(string.format("%s_%s_%d.bp", self.name, mom, 0), tm)
	 end
      end
      self.distIoFrame = self.distIoFrame+1
   end
end

function Species:rkStepperFields()
   return self.distf
end

   
function Species:forwardEuler(tCurr, dt, fIn, emIn, fOut)
   if self.evolve then
      return self.vlasovSlvr:advance(tCurr, dt, {fIn, emIn}, {fOut})
   else
      fOut:copy(fIn) -- just copy stuff over
      return true, GKYL_MAX_DOUBLE
   end
end

function Species:incrementCouplingMoments(tCurr, dt, fIn, momOut)
   -- compute momentum and increment as current into supplied field
   self.momDensityCalc:advance(tCurr, dt, {fIn}, { self.momDensity })
   momOut[1]:accumulate(self.charge, self.momDensity)
end

function Species:applyBc(tCurr, dt, fIn)
   fIn:sync()
end

function Species:totalSolverTime()
   return self.vlasovSlvr.totalTime
end

function Species:streamTime()
   return self.vlasovSlvr:streamTime()
end

function Species:forceTime()
   return self.vlasovSlvr:forceTime()
end

function Species:incrementTime()
   return self.vlasovSlvr:incrementTime()
end

function Species:momCalcTime()
   local tm = self.momDensityCalc.totalTime
   for i, mom in ipairs(self.diagnosticMoments) do
      tm = tm + self.diagnosticMomentUpdaters[i].totalTime
   end
   return tm
end

function Species:intMomCalcTime()
   return self.intMomentCalc.totalTime
end

return Species
