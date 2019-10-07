-- Gkyl ------------------------------------------------------------------------
--
-- App support code: FluidSpecies object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo      = require "Io.AdiosCartFieldIo"
local Basis                 = require "Basis"
local Collisions            = require "App.Collisions"
local DataStruct            = require "DataStruct"
local DecompRegionCalc      = require "Lib.CartDecomp"
local Grid                  = require "Grid"
local LinearTrigger         = require "Lib.LinearTrigger"
local Mpi                   = require "Comm.Mpi"
local Proto                 = require "Lib.Proto"
local Projection            = require "App.Projection"
local ProjectionBase        = require "App.Projection.ProjectionBase"
local SpeciesBase           = require "App.Species.SpeciesBase"
local Time                  = require "Lib.Time"
local Updater               = require "Updater"
local ffi                   = require "ffi"
local xsys                  = require "xsys"
local ConstDiffusionModDecl = require "Eq.constDiffusionData.ConstDiffusionModDecl"

-- Function to create basis functions.
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   end
end

-- Base class for kinetic species.
local FluidSpecies = Proto(SpeciesBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function FluidSpecies:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function FluidSpecies:fullInit(appTbl)
   local tbl = self.tbl -- Previously store table.

   self.cfl      =  0.1
   self.charge   = tbl.charge and tbl.charge or 1.0
   self.mass     = tbl.mass and tbl.mass or 1.0
   self.ioMethod = "MPI"

   self.evolve              = xsys.pickBool(tbl.evolve, true) -- By default, evolve species.
   -- By default, do not write species if it is not evolved.
   self.forceWrite          = xsys.pickBool(tbl.forceWrite, false)
   self.evolveCollisionless = xsys.pickBool(tbl.evolveCollisionless,
                                            self.evolve)
   self.evolveCollisions    = xsys.pickBool(tbl.evolveCollisions, self.evolve)
   self.evolveSources       = xsys.pickBool(tbl.evolveSources, self.evolve)

   self.confBasis = nil -- Will be set later

   -- Create triggers to write diagnostics.
   if tbl.nDiagnosticFrame then
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nDiagnosticFrame)
   else
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.diagIoFrame = 0 -- Frame number for diagnostics.

   -- For storing integrated moments.
   self.integratedMoments = nil -- Allocated in alloc() method.

   -- Store initial condition function.
   self.initFunc = tbl.init

   -- Default to a single moment.
   self.nMoments = 1
   self.nGhost   = 1 -- Default is 1 ghost-cell in each direction.

   self.hasNonPeriodicBc        = false -- To indicate if we have non-periodic BCs.
   self.bcx, self.bcy, self.bcz = { }, { }, { }

   -- Read in boundary conditions.
   -- Check to see if bc type is good is now done in createBc.
   if tbl.bcx then
      self.bcx[1], self.bcx[2] = tbl.bcx[1], tbl.bcx[2]
      self.hasNonPeriodicBc    = true
   end
   if tbl.bcy then
      self.bcy[1], self.bcy[2] = tbl.bcy[1], tbl.bcy[2]
      self.hasNonPeriodicBc    = true
   end
   if tbl.bcz then
      self.bcz[1], self.bcz[2] = tbl.bcz[1], tbl.bcz[2]
      self.hasNonPeriodicBc    = true
   end
   
   self.ssBc = {}
   if tbl.ssBc then
      self.ssBc[1] = tbl.ssBc[1]
   end

   self.boundaryConditions   = { } -- List of Bcs to apply.
   self.ssBoundaryConditions = { } -- List of stair-stepped Bcs to apply.
   self.zeroFluxDirections   = {}

   self.bcTime = 0.0 -- Timer for BCs.

   -- Collisions: currently used for a diffusion term.
   self.collisions = {}
   for nm, val in pairs(tbl) do
      if Collisions.CollisionsBase.is(val) then
         self.collisions[nm] = val
         self.collisions[nm]:setName(nm)
         val:setSpeciesName(self.name)
         val:fullInit(tbl)    -- Initialize collisions (diffusion).
      end
   end

   -- Initialization.
   self.projections = {}
   for nm, val in pairs(tbl) do
      if ProjectionBase.is(val) then
         self.projections[nm] = val
      end
   end
   if tbl.sourceTimeDependence then
      self.sourceTimeDependence = tbl.sourceTimeDependence
   else
      self.sourceTimeDependence = function (t) return 1.0 end
   end
   -- It is possible to use the keyword 'initSource' to specify a
   -- function directly without using a Projection object.
   if type(tbl.init) == "function" then
      self.projections["init"] = Projection.FluidProjection.FunctionProjection {
         func = function (t, zn)
            return tbl.init(t, zn, self)
         end,
         isInit = true,
      }
   elseif type(tbl.init) == "string" then
      -- Specify the suffix of the file with the initial condition (including the extension).
      -- The prefix is assumed to be the name of the input file.
      self.projections["init"] = Projection.FluidProjection.ReadInput {
         inputFile = tbl.init,
         isInit    = true,
      }
   end
   if type(tbl.source) == "function" then
      self.projections["initSource"] = Projection.FluidProjection.FunctionProjection {
         func = function (t, zn)
            return tbl.source(t, zn, self)
         end,
         isInit   = false,
         isSource = true,
      }
   elseif type(tbl.source) == "string" then
      -- Specify the suffix of the file with the source (including the extension).
      -- The prefix is assumed to be the name of the input file.
      self.projections["initSource"] = Projection.FluidProjection.ReadInput {
         inputFile = tbl.source,
         isInit    = false,
         isSource  = true,
      }
   end

   self.useShared         = xsys.pickBool(appTbl.useShared, false)
   self.positivity        = xsys.pickBool(tbl.applyPositivity, false)
   self.positivityDiffuse = xsys.pickBool(tbl.positivityDiffuse, self.positivity)
   self.positivityRescale = xsys.pickBool(tbl.positivityRescale, false)
   self.deltaF            = xsys.pickBool(appTbl.deltaF, false)

   self.tCurr = 0.0
end

function FluidSpecies:getCharge() return self.charge end
function FluidSpecies:getMass() return self.mass end
function FluidSpecies:getEvolve() return self.evolve end

function FluidSpecies:getNdim()
   return self.ndim
end
function FluidSpecies:vdim()
   return 0
end
function FluidSpecies:setName(nm)
   self.name = nm
end
function FluidSpecies:setCfl(cfl)
   self.cfl = cfl
   for _, c in pairs(self.collisions) do
      c:setCfl(cfl)
   end
end
function FluidSpecies:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end
function FluidSpecies:setConfBasis(basis)
   self.confBasis = basis
   for _, c in pairs(self.collisions) do
      c:setConfBasis(basis)
   end
end
function FluidSpecies:setConfGrid(cgrid)
   self.confGrid = cgrid
   for _, c in pairs(self.collisions) do
      c:setConfGrid(cgrid)
   end
end

function FluidSpecies:createGrid(cLo, cUp, cCells, cDecompCuts, cPeriodicDirs)
   self.cdim = #cCells
   self.ndim = self.cdim

   -- Create decomposition.
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, cDecompCuts[d]) end
   self.decomp = DecompRegionCalc.CartProd {
      cuts      = decompCuts,
      useShared = self.useShared,
   }

   -- Create computational domain.
   local lower, upper, cells = {}, {}, {}
   for d = 1, self.cdim do
      table.insert(lower, cLo[d])
      table.insert(upper, cUp[d])
      table.insert(cells, cCells[d])
   end
   self.grid = Grid.RectCart {
      lower         = lower,
      upper         = upper,
      cells         = cells,
      periodicDirs  = cPeriodicDirs,
      decomposition = self.decomp,
   }
end

function FluidSpecies:createBasis(nm, polyOrder)
   self.basis = createBasis(nm, self.ndim, polyOrder)
end

function FluidSpecies:allocMoment()
   local m = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {self.nGhost, self.nGhost}
   }
   return m
end
function FluidSpecies:allocVectorMoment(dim)
   local m = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis()*dim,
      ghost         = {self.nGhost, self.nGhost}
   }
   return m
end

function FluidSpecies:allocMomCouplingFields()
   return {self:allocVectorMoment(self.nMoments)}
end

function FluidSpecies:bcAbsorbFunc(dir, tm, idxIn, fIn, fOut, fBC)
   -- Note that for bcAbsorb there is no operation on fIn,
   -- so skinLoop (which determines indexing of fIn) does not matter.
   for i = 1, self.nMoments*self.basis:numBasis() do
      fOut[i] = 0.0
   end
end

function FluidSpecies:bcCopyFunc(dir, tm, idxIn, fIn, fOut, fBC)
   for i = 1, self.nMoments*self.basis:numBasis() do
      fOut[i] = fIn[i]
   end
end

function FluidSpecies:bcDirichletFunc(dir, tm, idxIn, fIn, fOut, fBC)
   -- Impose f=fBC at the boundary.
--   fOut[1] = 2^(3/2)*fBC-fIn[1]
--   fOut[2] = fIn[2]
   if (idxIn == 1) then
      self.constDiffDirichletBCs[dir][1](self.grid:dx(dir),fIn:data(), fBC, fOut:data())
   else
      self.constDiffDirichletBCs[dir][2](self.grid:dx(dir),fIn:data(), fBC, fOut:data())
   end
end

function FluidSpecies:bcNeumannFunc(dir, tm, idxIn, fIn, fOut, fpBC)
   -- Impose f'=fpBC at the boundary.
--   fOut[1] = (2^(9/2)*self.grid:dx(dir)*fpBC+10*math.sqrt(3)*fIn[2]+21*fIn[1])/21
--   fOut[2] = (2^(5/2)*math.sqrt(3)*self.grid:dx(dir)*fpBC-3*fIn[2])/21
   if (idxIn == 1) then
      self.constDiffNeumannBCs[dir][1](self.grid:dx(dir),fIn:data(), fpBC, fOut:data())
   else
      self.constDiffNeumannBCs[dir][2](self.grid:dx(dir),fIn:data(), fpBC, fOut:data())
   end
end

-- Function to construct a BC updater.
function FluidSpecies:makeBcUpdater(dir, edge, bcList, skinLoop,
                                    hasExtFld)

   -- If BC is Dirichlet or Neumann select appropriate kernels.
   if (bcList[3] == 5) then
      local nm, ndim, p = self.basis:id(), self.basis:ndim(), self.basis:polyOrder()
      self.constDiffDirichletBCs = ConstDiffusionModDecl.selectBCs(nm, ndim, p, "Dirichlet")
   elseif (bcList[3] == 6) then
      local nm, ndim, p = self.basis:id(), self.basis:ndim(), self.basis:polyOrder()
      self.constDiffNeumannBCs = ConstDiffusionModDecl.selectBCs(nm, ndim, p, "Neumann")
   end

   return Updater.Bc {
      onGrid             = self.grid,
      boundaryConditions = bcList,
      dir                = dir,
      edge               = edge,
      skinLoop           = skinLoop,
      cdim               = self.cdim,
      vdim               = self.vdim,
      hasExtFld          = hasExtFld,
   }
end

-- Function to construct a stair-stepped BC updater.
function FluidSpecies:makeSsBcUpdater(dir, inOut, bcList)
   return Updater.StairSteppedBc {
      onGrid             = self.grid,
      inOut              = inOut,
      boundaryConditions = bcList,
      dir                = dir,
   }
end

function FluidSpecies:createBCs()
   -- Functions to make life easier while reading in BCs to apply.
   -- Note: appendBoundaryConditions defined in sub-classes.
   local function handleBc(dir, bc)
      if bc[1] then
	 self:appendBoundaryConditions(dir, 'lower', bc[1])
      end
      if bc[2] then
	 self:appendBoundaryConditions(dir, 'upper', bc[2])
      end
   end

   -- Add various BCs to list of BCs to apply.
   handleBc(1, self.bcx)
   handleBc(2, self.bcy)
   handleBc(3, self.bcz)
end

function FluidSpecies:createSolver(funcField)
   -- Create solvers for collisions (diffusion).
   for _, c in pairs(self.collisions) do
      c:createSolver(funcField)
   end

   if self.positivity then
      self.posChecker = Updater.PositivityCheck {
         onGrid = self.grid,
         basis  = self.basis,
      }

      self.posRescaler = Updater.PositivityRescale {
         onGrid = self.grid,
         basis  = self.basis,
      }
   end
end

function FluidSpecies:alloc(nRkDup)
   -- Allocate fields needed in RK update.
   self.moments = {}
   for i = 1, nRkDup do
      self.moments[i] = self:allocVectorMoment(self.nMoments)
      self.moments[i]:clear(0.0)
   end
   -- Create Adios object for field I/O.
   self.momIo = AdiosCartFieldIo {
      elemType = self.moments[1]:elemType(),
      method   = self.ioMethod,
      metaData = {
         polyOrder = self.basis:polyOrder(),
         basisType = self.basis:id()
      },
   }
   self.couplingMoments   = self:allocVectorMoment(self.nMoments)
   self.integratedMoments = DataStruct.DynVector { numComponents = self.nMoments }

   if self.positivity then
      self.fPos = self:allocVectorMoment(self.nMoments)
   end

   -- Array with one component per cell to store cflRate in each cell.
   self.cflRateByCell = DataStruct.Field {
      onGrid        = self.grid,
      numComponents = 1,
      ghost         = {1, 1},
   }
   self.cflRateByCell:clear(0.0)
   self.cflRatePtr  = self.cflRateByCell:get(1)
   self.cflRateIdxr = self.cflRateByCell:genIndexer()
   self.dt          = ffi.new("double[2]")
   self.dtGlobal    = ffi.new("double[2]")

   self:createBCs()
end

function FluidSpecies:initDist()

   local initCnt = 0
   for _, pr in pairs(self.projections) do
      pr:fullInit(self)
      pr:run(0.0, self.moments[2])
      -- This barrier is needed as when using MPI-SHM some
      -- processes will get to accumulate before projection is finished.
      Mpi.Barrier(self.grid:commSet().sharedComm)
      if pr.isInit then
         self.moments[1]:accumulate(1.0, self.moments[2])
         initCnt = initCnt + 1
      end
      if pr.isSource then
         if not self.mSource then
            self.mSource = self:allocVectorMoment(self.nMoments)
         end
         self.mSource:accumulate(1.0, self.moments[2])
         if self.positivityRescale then
            self.posRescaler:advance(0.0, {self.mSource}, {self.mSource})
         end
      end
   end
   assert(initCnt > 0,
          string.format("FluidSpecies: Species '%s' not initialized!", self.name))
   self.moments[2]:clear(0.0)

   if self.positivityRescale or self.positivityDiffuse then
      self.posRescaler:advance(0.0, {self.moments[1]}, {self.moments[1]}, false)
   end
end

function FluidSpecies:rkStepperFields()
   return self.moments
end

-- For RK timestepping.
function FluidSpecies:copyRk(outIdx, aIdx)
   self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx])
end
-- For RK timestepping. 
function FluidSpecies:combineRk(outIdx, a, aIdx, ...)
   local args  = {...} -- Package up rest of args as table.
   local nFlds = #args/2
   self:rkStepperFields()[outIdx]:combine(a, self:rkStepperFields()[aIdx])
   for i = 1, nFlds do -- Accumulate rest of the fields.
      self:rkStepperFields()[outIdx]:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]])
   end
end

function FluidSpecies:suggestDt()
   return GKYL_MAX_DOUBLE
end

function FluidSpecies:clearCFL()
end

function FluidSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   self.tCurr = tCurr
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   if self.evolveCollisionless then
      self.solver:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      local em = emIn[1]:rkStepperFields()[inIdx]
      if self.positivityRescale then
         self.posRescaler:advance(tCurr, {fIn}, {self.fPos}, false)
         self.solver:advance(tCurr, {self.fPos, em}, {fRhsOut})
      else
         self.solver:advance(tCurr, {fIn, em}, {fRhsOut})
      end
   else
      fRhsOut:clear(0.0) -- No RHS.
   end

   -- Perform the collision (diffusion) update.
   if self.evolveCollisions then
      for _, c in pairs(self.collisions) do
         c.diffusionSlvr:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
         c:advance(tCurr, fIn, species, fRhsOut)
         -- The full 'species' list is needed for the cross-species
         -- collisions.
      end
   end

   if self.mSource and self.evolveSources then
      -- Add source term to the RHS.
      -- Barrier over shared communicator before accumulate.
      Mpi.Barrier(self.grid:commSet().sharedComm)
      fRhsOut:accumulate(self.sourceTimeDependence(tCurr), self.mSource)
   end
end

function FluidSpecies:checkPositivity(tCurr, idx)
  local status = true
  if self.positivity then
     status = self.posChecker:advance(tCurr, {self:rkStepperFields()[idx]}, {})
  end
  return status
end

function FluidSpecies:applyBcIdx(tCurr, idx, isFirstRk)
  if self.positivityDiffuse then
     self.posRescaler:advance(tCurr, {self:rkStepperFields()[idx]}, {self:rkStepperFields()[idx]}, true, isFirstRk)
  end
  for dir = 1, self.ndim do
     self:applyBc(tCurr, self:rkStepperFields()[idx], dir)
  end
  if self.positivity then
     self:checkPositivity(tCurr, idx)
  end
end

function FluidSpecies:applyBc(tCurr, fIn, dir)
   local tmStart = Time.clock()
   if self.evolve then
      if self.hasNonPeriodicBc then
         for _, bc in ipairs(self.boundaryConditions) do
            if (not dir) or dir == bc:getDir() then
               bc:advance(tCurr, {}, {fIn})
            end
         end
      end
      for _, bc in ipairs(self.ssBoundaryConditions) do
         if (not dir) or dir == bc:getDir() then
             bc:advance(tCurr, {}, {fIn})
          end
       end
      fIn:sync()
   end
   self.bcTime = self.bcTime + Time.clock()-tmStart
end

function FluidSpecies:createDiagnostics()
   -- Create updater to compute volume-integrated moments.
   self.intMom2Calc = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = self.grid,
      basis         = self.basis,
      numComponents = self.nMoments,
      quantity      = "V"
   }
end

function FluidSpecies:write(tm, force)
   if self.evolve or self.forceWrite then
      -- Compute integrated diagnostics.
      self.intMom2Calc:advance(tm, { self.moments[1] }, { self.integratedMoments })
      
      -- Only write stuff if triggered.
      if self.diagIoTrigger(tm) or force then
	 self.momIo:write(
	    self.moments[1], string.format("%s_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame)
         self.integratedMoments:write(
            string.format("%s_intMom_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame)

         if self.positivityDiffuse then
            self.posRescaler:write(tm, self.diagIoFrame, self.name)
         end

         if tm == 0.0 and self.mSource then
            self.momIo:write(self.mSource, string.format("%s_mSource_0.bp", self.name), tm, self.diagIoFrame)
         end
 
         self.diagIoFrame = self.diagIoFrame+1
      end
   else
      -- If not evolving species, don't write anything except initial conditions.
      if self.diagIoFrame == 0 then
         self.momIo:write(self.moments[1], string.format("%s_%d.bp", self.name, 0), tm, 0)
      end
      self.diagIoFrame = self.diagIoFrame+1
   end
end

function FluidSpecies:writeRestart(tm)
   self.momIo:write(
      self.moments[1], string.format("%s_restart.bp", self.name), tm, self.diagIoFrame)
   self.integratedMoments:write(
      string.format("%s_intMom_restart.bp", self.name), tm, self.diagIoFrame, false)
end

function FluidSpecies:readRestart()
   local tm, fr = self.momIo:read(self.moments[1], string.format("%s_restart.bp", self.name))
   self.diagIoFrame = fr -- Reset internal frame counter.
   self.integratedMoments:read(string.format("%s_intMom_restart.bp", self.name))   
   
   for dir = 1, self.ndim do
      self:applyBc(tm, self.moments[1], dir)
   end
   self.moments[1]:sync() -- Must get all ghost-cell data correct.

   -- Iterate triggers.
   self.diagIoTrigger(tm)

   return tm
end

-- Timers.
function FluidSpecies:totalSolverTime()
   if self.solver then
      return self.solver.totalTime
   end
   return 0
end
function FluidSpecies:totalBcTime()
   return self.bcTime
end
function FluidSpecies:momCalcTime()
   return 0
end
function FluidSpecies:intMomCalcTime()
   return self.intMom2Calc.totalTime
end

return FluidSpecies
