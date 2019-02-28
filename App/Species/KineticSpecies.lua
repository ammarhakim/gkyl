-- Gkyl ------------------------------------------------------------------------
--
-- App support code: KineticSpecies object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
-- System imports
local xsys = require "xsys"

-- Gkeyll imports
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis = require "Basis"
local Collisions = require "App.Collisions"
local DataStruct = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local GaussQuadRules = require "Lib.GaussQuadRules"
local Grid = require "Grid"
local LinearDecomp = require "Lib.LinearDecomp"
local LinearTrigger = require "Lib.LinearTrigger"
local Mpi = require "Comm.Mpi"
local Projection = require "App.Projection"
local ProjectionBase = require "App.Projection.ProjectionBase"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local SpeciesBase = require "App.Species.SpeciesBase"
local Time = require "Lib.Time"
local Updater = require "Updater"
local ffi = require "ffi"

-- function to create basis functions
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "tensor" then
      return Basis.CartModalTensor { ndim = ndim, polyOrder = polyOrder }
   end
end

-- base class for kinetic species
local KineticSpecies = Proto(SpeciesBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function KineticSpecies:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization
function KineticSpecies:fullInit(appTbl)
   local tbl = self.tbl -- previously store table

   self.cfl =  0.1
   self.charge = tbl.charge and tbl.charge or 1.0
   self.mass = tbl.mass and tbl.mass or 1.0
   self.n0 = tbl.n0 or n0
   self.lower, self.upper = tbl.lower, tbl.upper
   self.cells = tbl.cells
   self.vdim = #self.cells -- velocity dimensions
   self.ioMethod = "MPI"
   self.evolve = xsys.pickBool(tbl.evolve, true) -- by default, evolve species
   self.evolveCollisionless = xsys.pickBool(tbl.evolveCollisionless,
					    self.evolve) 
   self.evolveCollisions = xsys.pickBool(tbl.evolveCollisions, self.evolve) 
   self.evolveSources = xsys.pickBool(tbl.evolveSources, self.evolve) 
   self.confBasis = nil -- Will be set later

   assert(#self.lower == self.vdim, "'lower' must have " .. self.vdim .. " entries")
   assert(#self.upper == self.vdim, "'upper' must have " .. self.vdim .. " entries")
   self.coordinateMap = tbl.coordinateMap

   self.useShared = xsys.pickBool(appTbl.useShared, false)
   
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
      self.distIoTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nDistFuncFrame)
   else
      self.distIoTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   if tbl.nDiagnosticFrame then
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nDiagnosticFrame)
   else
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   -- create trigger for how frequently to compute integrated moments
   self.calcIntQuantFlag = false
   if appTbl.calcIntQuantEvery then
      self.calcIntQuantTrigger = LinearTrigger(0, appTbl.tEnd,  math.floor(1/appTbl.calcIntQuantEvery))
   else
      self.calcIntQuantFlag = true
   end

   self.distIoFrame = 0 -- frame number for distribution function
   self.diagIoFrame = 0 -- frame number for diagnostics

   self.writeGhost = xsys.pickBool(appTbl.writeGhost, false)

   -- write perturbed moments by subtracting background before moment calc.. false by default
   self.perturbedMoments = false
   -- read in which diagnostic moments to compute on output
   self.diagnosticMoments = { }
   if tbl.diagnosticMoments then
      for i, nm in pairs(tbl.diagnosticMoments) do
         if i == "perturbed" and nm == true then 
            self.perturbedMoments = true
         elseif type(i) == "number" then
	    self.diagnosticMoments[i] = nm
         end
      end
   end

   -- read in which integrated diagnostic moments to compute on output
   self.diagnosticIntegratedMoments = { }
   if tbl.diagnosticIntegratedMoments then
      for i, nm in ipairs(tbl.diagnosticIntegratedMoments) do
         self.diagnosticIntegratedMoments[i] = nm
      end
   end

   -- get a random seed for random initial conditions
   self.randomseed = tbl.randomseed

   -- Initialization
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
   -- it is possible to use the keyword 'init', 'initBackground', and
   -- 'initSource' to specify a function directly without using a
   -- Projection object
   if type(tbl.init) == "function" then
      self.projections["init"] = Projection.KineticProjection.FunctionProjection {
	 func = function (t, zn)
	    return tbl.init(t, zn, self)
	 end,
	 isInit = true,
      }
   end
   if type(tbl.initBackground) == "function" then
      self.projections["initBackground"] = Projection.KineticProjection.FunctionProjection {
	 func = function (t, zn)
	    return tbl.initBackground(t, zn, self)
	 end,
	 isInit = false,
	 isBackground = true,
      }
   end
   if type(tbl.source) == "function" then
      self.projections["initSource"] = Projection.KineticProjection.FunctionProjection {
	 func = function (t, zn)
	    return tbl.source(t, zn, self)
	 end,
	 isInit = false,
	 isSource = true,
      }
   end
   -- >> LEGACY CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   if type(tbl.init) == "table" and tbl.init[1] == "maxwellian" then 
      self.projections["init"] = Projection.GkProjection.MaxwellianProjection {
	 density = tbl.init.density,
	 driftSpeed = tbl.init.driftSpeed,
	 temperature = tbl.init.temperature,
	 exactScaleM0 = true,
	 exactLagFixM012 = false,
	 isInit = true,
      }
   end 
   if type(tbl.initBackground) == "table" and tbl.initBackground[1] == "maxwellian" then 
      self.projections["initBackground"] = Projection.GkProjection.MaxwellianProjection {
	 density = tbl.initBackground.density,
	 driftSpeed = tbl.initBackground.driftSpeed,
	 temperature = tbl.initBackground.temperature,
	 exactScaleM0 = true,
	 exactLagFixM012 = false,
	 isInit = false,
	 isBackground = true,
      }
   end 
   if type(tbl.source) == "table" and tbl.source[1] == "maxwellian" then 
      self.projections["initSource"] = Projection.GkProjection.MaxwellianProjection {
	 density = tbl.source.density,
	 driftSpeed = tbl.source.driftSpeed,
	 temperature = tbl.source.temperature,
	 exactScaleM0 = true,
	 exactLagFixM012 = false,
	 isInit = false,
	 isSource = true,
      }
   end 
   -- <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   self.deltaF = xsys.pickBool(appTbl.deltaF, false)
   self.fluctuationBCs = xsys.pickBool(tbl.fluctuationBCs, false)
   if self.deltaF then self.fluctuationBCs = true end

   self.zeroFluxDirections = {}

   self.hasNonPeriodicBc = false -- to indicate if we have non-periodic BCs
   self.bcx, self.bcy, self.bcz = { }, { }, { }

   -- read in boundary conditions
   -- check to see if bc type is good is now done in createBc
   if tbl.bcx then
      self.bcx[1], self.bcx[2] = tbl.bcx[1], tbl.bcx[2]
      if self.bcx[1] == nil or self.bcx[2] == nil then assert(false, "KineticSpecies: unsupported BC type") end
      self.hasNonPeriodicBc = true
   end
   if tbl.bcy then
      self.bcy[1], self.bcy[2] = tbl.bcy[1], tbl.bcy[2]
      if self.bcy[1] == nil or self.bcy[2] == nil then assert(false, "KineticSpecies: unsupported BC type") end
      self.hasNonPeriodicBc = true
   end
   if tbl.bcz then
      self.bcz[1], self.bcz[2] = tbl.bcz[1], tbl.bcz[2]
      if self.bcz[1] == nil or self.bcz[2] == nil then assert(false, "KineticSpecies: unsupported BC type") end
      self.hasNonPeriodicBc = true
   end

   self.boundaryConditions = { } -- list of Bcs to apply

   self.bcTime = 0.0 -- timer for BCs
   self.integratedMomentsTime = 0.0 -- timer for integrated moments

   -- Collisions/Sources
   self.collisions = {}
   for nm, val in pairs(tbl) do
      if Collisions.CollisionsBase.is(val) then
	 self.collisions[nm] = val
	 self.collisions[nm]:setName(nm)
	 val:setSpeciesName(self.name)
	 val:fullInit(tbl) -- initialize collisions
      end
   end

   self.positivity = xsys.pickBool(tbl.applyPositivity, false)
   self.positivityDiffuse = xsys.pickBool(tbl.positivityDiffuse, self.positivity)
   self.positivityRescale = xsys.pickBool(tbl.positivityRescale, false)
   
   -- for GK only: flag for gyroaveraging
   self.gyavg = xsys.pickBool(tbl.gyroaverage, false)

   -- gravity table for running Vlasov simulations with constant gravity
   self.constGravity = tbl.constGravity

   self.tCurr = 0.0
end

function KineticSpecies:getCharge() return self.charge end
function KineticSpecies:getMass() return self.mass end

function KineticSpecies:getNdim()
   return self.ndim
end
function KineticSpecies:getVdim()
   return self.vdim
end
function KineticSpecies:setName(nm) -- Needs to be called before fullInit().
   self.name = nm
end
function KineticSpecies:setCfl(cfl)
   self.cfl = cfl
   for _, c in pairs(self.collisions) do
      c:setCfl(cfl)
   end   
end
function KineticSpecies:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end
function KineticSpecies:setConfBasis(basis)
   self.confBasis = basis
   for _, c in pairs(self.collisions) do
      c:setConfBasis(basis)
   end
end
function KineticSpecies:setConfGrid(grid)
   self.confGrid = grid
   for _, c in pairs(self.collisions) do
      c:setConfGrid(grid)
   end
end

function KineticSpecies:createGrid(cLo, cUp, cCells, cDecompCuts,
				   cPeriodicDirs, cMap)
   self.cdim = #cCells
   self.ndim = self.cdim+self.vdim

   -- create decomposition
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, cDecompCuts[d]) end
   for d = 1, self.vdim do table.insert(decompCuts, self.decompCuts[d]) end
   self.decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      useShared = self.useShared,
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

   local GridConstructor = Grid.RectCart
   local coordinateMap = {} -- table of functions
   -- construct comp -> phys mappings if they exist
   if self.coordinateMap or cMap then
      if cMap and self.coordinateMap then
         for d = 1, self.cdim do
            table.insert(coordinateMap, cMap[d])
         end
         for d = 1, self.vdim do
            table.insert(coordinateMap, self.coordinateMap[d])
         end
      elseif cMap then
         for d = 1, self.cdim do
            table.insert(coordinateMap, cMap[d])
         end
         for d = 1, self.vdim do
            table.insert(coordinateMap, function (z) return z end)
         end
      else
         for d = 1, self.cdim do
            table.insert(coordinateMap, function (z) return z end)
         end
         for d = 1, self.vdim do
            table.insert(coordinateMap, self.coordinateMap[d])
         end
      end
      GridConstructor = Grid.NonUniformRectCart
   end
   self.grid = GridConstructor {
      lower = lower,
      upper = upper,
      cells = cells,
      periodicDirs = cPeriodicDirs,
      decomposition = self.decomp,
      mappings = coordinateMap,
   }

   for _, c in pairs(self.collisions) do
      c:setPhaseGrid(self.grid)
   end
end

function KineticSpecies:createBasis(nm, polyOrder)
   self.basis = createBasis(nm, self.ndim, polyOrder)
   for _, c in pairs(self.collisions) do
      c:setPhaseBasis(self.basis)
   end
end

function KineticSpecies:allocDistf()
   local f = DataStruct.Field {
	onGrid = self.grid,
	numComponents = self.basis:numBasis(),
	ghost = {1, 1},
   }
   f:clear(0.0)
   return f
end
function KineticSpecies:allocMoment()
   local m = DataStruct.Field {
	onGrid = self.confGrid,
	numComponents = self.confBasis:numBasis(),
	ghost = {1, 1},
   }
   m:clear(0.0)
   return m
end
function KineticSpecies:allocVectorMoment(dim)
   local m = DataStruct.Field {
	onGrid = self.confGrid,
	numComponents = self.confBasis:numBasis()*dim,
	ghost = {1, 1},
   }
   m:clear(0.0)
   return m
end

-- various functions to apply BCs of different types
function KineticSpecies:bcAbsorbFunc(dir, tm, idxIn, fIn, fOut)
   -- note that for bcAbsorb there is no operation on fIn,
   -- so skinLoop (which determines indexing of fIn) does not matter 
   for i = 1, self.basis:numBasis() do
      fOut[i] = 0.0
   end
end
function KineticSpecies:bcOpenFunc(dir, tm, idxIn, fIn, fOut)
   -- requires skinLoop = "pointwise"
   self.basis:flipSign(dir, fIn, fOut)
end
function KineticSpecies:bcCopyFunc(dir, tm, idxIn, fIn, fOut)
   -- requires skinLoop = "pointwise"
   for i = 1, self.basis:numBasis() do
      fOut[i] = fIn[i]
   end
end

-- function to construct a BC updater
function KineticSpecies:makeBcUpdater(dir, vdir, edge, bcList, skinLoop,
				      hasExtFld)
   return Updater.Bc {
      onGrid = self.grid,
      boundaryConditions = bcList,
      dir = dir,
      vdir = vdir,
      edge = edge,
      skinLoop = skinLoop,
      cdim = self.cdim,
      vdim = self.vdim,
      hasExtFld = hasExtFld,
   }
end

function KineticSpecies:createBCs()
   -- functions to make life easier while reading in BCs to apply
   -- note: appendBoundaryConditions defined in sub-classes
   local function handleBc(dir, bc)
      if bc[1] then
	 self:appendBoundaryConditions(dir, "lower", bc[1])
      end
      if bc[2] then
	 self:appendBoundaryConditions(dir, "upper", bc[2])
      end
   end

   -- add various BCs to list of BCs to apply
   handleBc(1, self.bcx)
   handleBc(2, self.bcy)
   handleBc(3, self.bcz)
end

function KineticSpecies:createSolver(funcField)
   -- create solvers for collisions
   for _, c in pairs(self.collisions) do
      c:createSolver(funcField)
   end
end

function KineticSpecies:alloc(nRkDup)
   -- allocate fields needed in RK update
   self.distf = {}
   for i = 1, nRkDup do
      self.distf[i] = self:allocDistf()
   end
   -- create Adios object for field I/O
   self.distIo = AdiosCartFieldIo {
      elemType = self.distf[1]:elemType(),
      method = self.ioMethod,
      writeGhost = self.writeGhost
   }

   if self.positivity then
      self.fPos = self:allocDistf()
   end

   -- array with one component per cell to store cflRate in each cell
   self.cflRateByCell = DataStruct.Field {
	onGrid = self.grid,
	numComponents = 1,
	ghost = {1, 1},
   }
   self.cflRateByCell:clear(0.0)
   self.cflRatePtr = self.cflRateByCell:get(1)
   self.cflRateIdxr = self.cflRateByCell:genIndexer()
   self.dt = ffi.new("double[2]")
   self.dtGlobal = ffi.new("double[2]")

   self:createBCs()
end

-- note: do not call applyBc here. it is called later in initialization sequence.
function KineticSpecies:initDist()
   if self.randomseed then 
      math.randomseed(self.randomseed) 
   else
      math.randomseed(47*Mpi.Comm_rank(Mpi.COMM_WORLD)+os.time())
   end

   local syncPeriodicDirs = true
   if self.fluctuationBCs then syncPeriodicDirs = false end

   local initCnt, backgroundCnt = 0, 0
   for _, pr in pairs(self.projections) do
      pr:fullInit(self)
      pr:run(0.0, self.distf[2])
      if pr.isInit then
	 self.distf[1]:accumulate(1.0, self.distf[2])
	 initCnt = initCnt + 1
      end
      if pr.isBackground then
	 if not self.f0 then 
	    self.f0 = self:allocDistf()
	 end
	 self.f0:accumulate(1.0, self.distf[2])
	 self.f0:sync(syncPeriodicDirs)
	 backgroundCnt = backgroundCnt + 1
      end
      if pr.isSource then
	 if not self.fSource then 
	    self.fSource = self:allocDistf()
	 end
	 self.fSource:accumulate(1.0, self.distf[2])
         if self.positivityRescale then
           self.posRescaler:advance(0.0, {self.fSource}, {self.fSource})
         end
      end
      if pr.isReservoir then
	 if not self.fReservoir then 
	    self.fReservoir = self:allocDistf()
	 end
	 self.fReservoir:accumulate(1.0, self.distf[2])
      end
   end
   assert(initCnt > 0,
	  string.format("KineticSpecies: Species '%s' not initialized!", self.name))
   if self.f0 and backgroundCnt == 0 then 
      self.f0:copy(self.distf[1])
   end

   if self.fluctuationBCs then 
      assert(backgroundCnt > 0, "KineticSpecies: must specify an initial background distribution with 'initBackground' in order to use fluctuation-only BCs") 
   end

   self.distf[2]:clear(0.0)

   -- calculate initial density averaged over simulation domain
   --self.n0 = nil
   --local dens0 = self:allocMoment()
   --self.numDensityCalc:advance(0, {self.distf[1]}, {dens0})
   --local data
   --local dynVec = DataStruct.DynVector { numComponents = 1 }
   ---- integrate 
   --local calcInt = Updater.CartFieldIntegratedQuantCalc {
   --   onGrid = self.confGrid,
   --   basis = self.confBasis,
   --   numComponents = 1,
   --   quantity = "V"
   --}
   --calcInt:advance(0.0, {dens0}, {dynVec})
   --_, data = dynVec:lastData()
   --self.n0 = data[1]/self.confGrid:gridVolume()
   --print("Average density is " .. self.n0)
end

function KineticSpecies:rkStepperFields()
   return self.distf
end

-- for RK timestepping 
function KineticSpecies:copyRk(outIdx, aIdx)
   self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx])
end
-- for RK timestepping 
function KineticSpecies:combineRk(outIdx, a, aIdx, ...)
   local args = {...} -- package up rest of args as table
   local nFlds = #args/2
   self:rkStepperFields()[outIdx]:combine(a, self:rkStepperFields()[aIdx])
   for i = 1, nFlds do -- accumulate rest of the fields
      self:rkStepperFields()[outIdx]:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]])
   end

   -- Barrier after accumulate since applyBc has different loop structure
   Mpi.Barrier(self.grid:commSet().sharedComm)
 
   if a<=self.dtGlobal[0] then -- this should be sufficient to determine if this combine is a forwardEuler step
      -- only applyBc on forwardEuler combine
      self:applyBc(nil, self:rkStepperFields()[outIdx])
      -- only positivity diffuse on forwardEuler combine
      if self.positivityDiffuse then
         self.posRescaler:advance(self.tCurr, {self:rkStepperFields()[outIdx]}, {self:rkStepperFields()[outIdx]})
      end
   end
end

function KineticSpecies:suggestDt()
   -- loop over local region 
   local grid = self.grid
   self.dt[0] = GKYL_MAX_DOUBLE

   local tId = grid:subGridSharedId() -- local thread ID
   local localRange = self.cflRateByCell:localRange()
   local localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = localRange, numSplit = grid:numSharedProcs() }

   for idx in localRangeDecomp:rowMajorIter(tId) do
      -- calculate local min dt from local cflRates on each proc
      self.cflRateByCell:fill(self.cflRateIdxr(idx), self.cflRatePtr)
      self.dt[0] = math.min(self.dt[0], self.cfl/self.cflRatePtr:data()[0])
   end

   -- all reduce to get global min dt
   Mpi.Allreduce(self.dt, self.dtGlobal, 1, Mpi.DOUBLE, Mpi.MIN, grid:commSet().comm)

   return math.min(self.dtGlobal[0], GKYL_MAX_DOUBLE)
end

function KineticSpecies:clearCFL()
   -- clear cflRateByCell for next cfl calculation
   self.cflRateByCell:clear(0.0)
end

function KineticSpecies:applyBc(tCurr, fIn)
   -- fIn is total distribution function
   local tmStart = Time.clock()

   if self.evolve then 
      local syncPeriodicDirsTrue = true

      if self.fluctuationBCs then
        -- if fluctuation-only BCs, subtract off background before applying BCs
        fIn:accumulate(-1.0, self.f0)

        -- possibly needed Barrier for fluctuation-only BCs with shared memory on
        -- there is a barrier before applyBc is entered, but need a barrier before sync()
        -- needs to be tested, but I think this is needed -- Jimmy Juno 02/28/19
        Mpi.Barrier(self.grid:commSet().sharedComm)
      end

      -- apply non-periodic BCs (to only fluctuations if fluctuation BCs)
      if self.hasNonPeriodicBc then
         for _, bc in ipairs(self.boundaryConditions) do
            bc:advance(tCurr, {self.fReservoir}, {fIn})
         end
      end

      -- apply periodic BCs (to only fluctuations if fluctuation BCs)
      fIn:sync(syncPeriodicDirsTrue)

      if self.fluctuationBCs then
        -- put back together total distribution
        fIn:accumulate(1.0, self.f0)

        -- possibly needed Barrier for fluctuation-only BCs with shared memory on
        -- there is a barrier at the end of sync(), but need a barrier between accumulate and sync()
        -- needs to be tested, but I think this is needed -- Jimmy Juno 02/28/19
        Mpi.Barrier(self.grid:commSet().sharedComm)

        -- update ghosts in total distribution, without enforcing periodicity
        fIn:sync(not syncPeriodicDirsTrue)
      end
   end

   self.bcTime = self.bcTime + (Time.clock()-tmStart)
end

function KineticSpecies:createDiagnostics()
end

function KineticSpecies:calcDiagnosticMoments()
   local f = self.distf[1]
   if self.f0 and self.perturbedMoments then f:accumulate(-1, self.f0) end
   for i, mom in pairs(self.diagnosticMoments) do
      self.diagnosticMomentUpdaters[mom]:advance(
	 0.0, {f}, {self.diagnosticMomentFields[mom]})
   end
   if self.f0 and self.perturbedMoments then f:accumulate(1, self.f0) end
end

function KineticSpecies:calcDiagnosticWeakMoments()
   for i, mom in pairs(self.diagnosticWeakMoments) do
      self.weakDivision:advance(0.0, self.weakMomentOpFields[mom], {self.diagnosticMomentFields[mom]})
      if self.weakMomentScaleFac[mom] then self.diagnosticMomentFields[mom]:scale(self.weakMomentScaleFac[mom]) end
   end
end

-- species-specific
function KineticSpecies:calcDiagnosticAuxMoments()
end

-- species-specific
function KineticSpecies:calcDiagnosticIntegratedMoments()
end

function KineticSpecies:calcAndWriteDiagnosticMoments(tm)
    self:calcDiagnosticMoments()
    for i, mom in ipairs(self.diagnosticMoments) do
       -- should one use AdiosIo object for this?
       self.diagnosticMomentFields[mom]:write(
          string.format("%s_%s_%d.bp", self.name, mom, self.diagIoFrame), tm, self.diagIoFrame, self.writeGhost)
    end

    if self.diagnosticWeakMoments then 
       self:calcDiagnosticWeakMoments()
       for i, mom in ipairs(self.diagnosticWeakMoments) do
          -- should one use AdiosIo object for this?
          self.diagnosticMomentFields[mom]:write(
             string.format("%s_%s_%d.bp", self.name, mom, self.diagIoFrame), tm, self.diagIoFrame, self.writeGhost)
       end
    end

    if self.diagnosticAuxMoments then 
    self:calcDiagnosticAuxMoments()
       for i, mom in ipairs(self.diagnosticAuxMoments) do
          -- should one use AdiosIo object for this?
          self.diagnosticMomentFields[mom]:write(
             string.format("%s_%s_%d.bp", self.name, mom, self.diagIoFrame), tm, self.diagIoFrame, self.writeGhost)
       end
    end

    -- write integrated moments
    for i, mom in ipairs(self.diagnosticIntegratedMoments) do
       self.diagnosticIntegratedMomentFields[mom]:write(
          string.format("%s_%s_%d.bp", self.name, mom, self.diagIoFrame), tm, self.diagIoFrame)
    end
end

function KineticSpecies:isEvolving()
   return self.evolve
end

function KineticSpecies:write(tm, force)
   if self.evolve then
      local tmStart = Time.clock()
      -- compute integrated diagnostics
      if self.calcIntQuantFlag == false then
         if self.calcIntQuantTrigger(tm) then
            self:calcDiagnosticIntegratedMoments(tm)
         end
      else
         self:calcDiagnosticIntegratedMoments(tm)
      end
      -- time computation of integrated moments
      self.integratedMomentsTime = self.integratedMomentsTime + Time.clock() - tmStart

      -- only write stuff if triggered
      if self.distIoTrigger(tm) or force then
	 self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, self.distIoFrame), tm, self.distIoFrame)
         if self.f0 then
            if tm == 0.0 then
	       self.distIo:write(self.f0, string.format("%s_f0_%d.bp", self.name, self.distIoFrame), tm, self.distIoFrame)
	    end
            self.distf[1]:accumulate(-1, self.f0)
            self.distIo:write(self.distf[1], string.format("%s_f1_%d.bp", self.name, self.distIoFrame), tm, self.distIoFrame)
            self.distf[1]:accumulate(1, self.f0)
         end
         if tm == 0.0 and self.fSource then
            self.distIo:write(self.fSource, string.format("%s_fSource_0.bp", self.name), tm, self.distIoFrame)
         end
	 self.distIoFrame = self.distIoFrame+1
      end


      if self.diagIoTrigger(tm) or force then
         -- compute moments and write them out
         self:calcAndWriteDiagnosticMoments(tm)

         if self.evolveCollisions then
            for _, c in pairs(self.collisions) do
               c:write(tm, self.diagIoFrame)
            end
         end

         self.diagIoFrame = self.diagIoFrame+1
      end
   else
      -- if not evolving species, don't write anything except initial conditions
      if self.distIoFrame == 0 then
	 self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, 0), tm, 0)

	 -- compute moments and write them out
	 self:calcAndWriteDiagnosticMoments(tm)
      end
      self.distIoFrame = self.distIoFrame+1
   end

end

function KineticSpecies:writeRestart(tm)
   -- (the final "false" prevents writing of ghost cells)
   self.distIo:write(self.distf[1], string.format("%s_restart.bp", self.name), tm, self.distIoFrame, false)
   for i, mom in ipairs(self.diagnosticMoments) do
      self.diagnosticMomentFields[mom]:write(
	 string.format("%s_%s_restart.bp", self.name, mom), tm, self.diagIoFrame, false)
   end   

   -- (the final "false" prevents flushing of data after write)
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
      self.diagnosticIntegratedMomentFields[mom]:write(
         string.format("%s_%s_restart.bp", self.name, mom), tm, self.diagIoFrame, false)
   end
end

function KineticSpecies:readRestart()
   local tm, fr = self.distIo:read(self.distf[1], string.format("%s_restart.bp", self.name))

   self:applyBc(tm, self.distf[1]) -- apply BCs and set ghost-cell data
   
   self.distIoFrame = fr -- reset internal frame counter
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
      local _, dfr = self.diagnosticIntegratedMomentFields[mom]:read(
         string.format("%s_%s_restart.bp", self.name, mom))
      self.diagIoFrame = dfr -- reset internal diagnostic IO frame counter
   end

   for i, mom in ipairs(self.diagnosticMoments) do
      local _, dfr = self.diagnosticMomentFields[mom]:read(
         string.format("%s_%s_restart.bp", self.name, mom))
      self.diagIoFrame = dfr -- reset internal diagnostic IO frame counter
   end

   -- iterate triggers
   self.distIoTrigger(tm)
   self.diagIoTrigger(tm)
   
   return tm
end

-- timers
function KineticSpecies:totalSolverTime()
   return self.solver.totalTime
end
function KineticSpecies:totalBcTime()
   return self.bcTime
end
function KineticSpecies:intMomCalcTime()
   return self.integratedMomentsTime
end

return KineticSpecies
