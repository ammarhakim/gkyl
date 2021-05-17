-- Gkyl ------------------------------------------------------------------------
--
-- App support code: KineticSpecies object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
-- System imports.
local xsys = require "xsys"

-- Gkeyll imports.
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis            = require "Basis"
local Collisions       = require "App.Collisions"
local DataStruct       = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid             = require "Grid"
local LinearTrigger    = require "Lib.LinearTrigger"
local Mpi              = require "Comm.Mpi"
local Projection       = require "App.Projection"
local ProjectionBase   = require "App.Projection.ProjectionBase"
local Proto            = require "Lib.Proto"
local SpeciesBase      = require "App.Species.SpeciesBase"
local SourceBase       = require "App.Sources.SourceBase"
local Time             = require "Lib.Time"
local Updater          = require "Updater"
local ffi              = require "ffi"
local lume             = require "Lib.lume"

-- Function to create basis functions.
local function createBasis(nm, ndim, polyOrder)
   if nm == "serendipity" then
      return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "maximal-order" then
      return Basis.CartModalMaxOrder { ndim = ndim, polyOrder = polyOrder }
   elseif nm == "tensor" then
      return Basis.CartModalTensor { ndim = ndim, polyOrder = polyOrder }
   end
end

-- Base class for kinetic species.
local KineticSpecies = Proto(SpeciesBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function KineticSpecies:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function KineticSpecies:fullInit(appTbl)
   local tbl = self.tbl -- Previously store table.

   self.charge = tbl.charge or 1.0
   self.mass   = tbl.mass or 1.0
   self.n0     = tbl.n0 or n0
   self.lower, self.upper = tbl.lower, tbl.upper
   self.cells  = tbl.cells
   self.vdim   = #self.cells -- Velocity dimensions.

   self.evolve              = xsys.pickBool(tbl.evolve, true) -- By default, evolve species.
   self.evolveCollisionless = xsys.pickBool(tbl.evolveCollisionless, self.evolve) 
   self.evolveCollisions    = xsys.pickBool(tbl.evolveCollisions, self.evolve) 

   assert(#self.lower == self.vdim, "'lower' must have " .. self.vdim .. " entries")
   assert(#self.upper == self.vdim, "'upper' must have " .. self.vdim .. " entries")
   self.coordinateMap = tbl.coordinateMap

   self.useShared = xsys.pickBool(appTbl.useShared, false)
   
   self.decompCuts = {}
   -- WE DO NOT ALLOW DECOMPOSITION IN VELOCITY SPACE
   for d = 1, self.vdim do self.decompCuts[d] = 1 end

   local nFrame = tbl.nDiagnosticFrame and tbl.nDiagnosticFrame or appTbl.nFrame
   -- Create triggers to write distribution functions and moments.
   if tbl.nDistFuncFrame then
      self.distIoTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nDistFuncFrame)
   else
      self.distIoTrigger = LinearTrigger(0, appTbl.tEnd, nFrame)
   end
   self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, nFrame)

   -- Create trigger for how frequently to compute integrated moments.
   -- Do not compute the integrated diagnostics less frequently than we output data.
   if appTbl.calcIntQuantEvery then
      self.calcIntQuantTrigger = LinearTrigger(0, appTbl.tEnd,  math.max(nFrame,math.floor(1/appTbl.calcIntQuantEvery)))
   else
      self.calcIntQuantTrigger = function(t) return true end
   end

   -- Write ghost cells on boundaries of global domain (for BCs).
   self.writeGhost = xsys.pickBool(appTbl.writeGhost, false)

   -- Write perturbed moments by subtracting background before moment calc.. false by default.
   self.perturbedMoments = false
   -- Read in which diagnostic moments to compute on output.
   self.requestedDiagnosticMoments = tbl.diagnosticMoments or {}
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

   -- Read in which integrated diagnostic moments to compute on output.
   self.diagnosticIntegratedMoments = { }
   if tbl.diagnosticIntegratedMoments then
      for i, nm in ipairs(tbl.diagnosticIntegratedMoments) do
         self.diagnosticIntegratedMoments[i] = nm
      end
   end

   -- Read in which boundary diagnostic moments to compute on output.
   self.boundaryFluxDiagnostics       = false
   self.diagnosticBoundaryFluxMoments = { }
   self.requestedDiagnosticBoundaryFluxMoments = tbl.diagnosticBoundaryFluxMoments or {}
   if tbl.diagnosticBoundaryFluxMoments then
      self.boundaryFluxDiagnostics = true
      for i, nm in pairs(tbl.diagnosticBoundaryFluxMoments) do
         self.diagnosticBoundaryFluxMoments[i] = nm
      end
   end

   -- Read in which boundary diagnostic moments to compute on output.
   self.boundaryFluxDiagnostics = false
   self.diagnosticIntegratedBoundaryFluxMoments = { }
   if tbl.diagnosticIntegratedBoundaryFluxMoments then
      self.boundaryFluxDiagnostics = true
      for i, nm in pairs(tbl.diagnosticIntegratedBoundaryFluxMoments) do
         self.diagnosticIntegratedBoundaryFluxMoments[i] = nm
      end
   end

   -- Get a random seed for random initial conditions.
   self.randomseed = tbl.randomseed

   -- Initialize table containing sources (if any).
   self.sources = {} 
   for nm, val in pairs(tbl) do
      if SourceBase.is(val) or string.find(nm,"source") then
         if ProjectionBase.is(val) then val = self:projToSource(val) end
	 self.sources[nm] = val
	 self.sources[nm]:setName(nm)
	 val:setSpeciesName(self.name)
	 val:fullInit(tbl) -- Initialize sources
      end
   end
   lume.setOrder(self.sources)  -- Save order in metatable to loop in the same order (w/ orderedIter, better for I/O).

   self.projections = {}
   for nm, val in pairs(tbl) do
      if ProjectionBase.is(val) and not string.find(nm,"source") then
         self.projections[nm] = val
      end
   end
   -- It is possible to use the keywords 'init' and 'background'
   -- to specify a function directly without using a Projection object.
   if type(tbl.init) == "function" then
      self.projections["init"] = Projection.KineticProjection.FunctionProjection {
	 func = function(t, zn) return tbl.init(t, zn, self) end,
      }
   end
   if type(tbl.background) == "function" then
      self.projections["background"] = Projection.KineticProjection.FunctionProjection {
	 func = function(t, zn) return tbl.background(t, zn, self) end,
      }
   end
   lume.setOrder(self.projections)  -- Save order in metatable to loop in the same order (w/ orderedIter, better for I/O).

   self.deltaF         = xsys.pickBool(appTbl.deltaF, false)
   self.fluctuationBCs = xsys.pickBool(tbl.fluctuationBCs, false)
   if self.deltaF then self.fluctuationBCs = true end

   self.zeroFluxDirections = {}

   self.hasNonPeriodicBc   = false -- To indicate if we have non-periodic BCs.
   self.boundaryConditions = { }   -- list of Bcs to apply
   self.bcx, self.bcy, self.bcz = { }, { }, { }
   -- Functional BCs
   self.evolveFnBC = xsys.pickBool(tbl.evolveFnBC, true)
   self.feedbackBC = xsys.pickBool(tbl.feedbackBC, false)

   -- Read in boundary conditions.
   -- Check to see if bc type is good is now done in createBc.
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

   -- Collisions.
   self.collisions = {}
   for nm, val in pairs(tbl) do
      if Collisions.CollisionsBase.is(val) then
	 self.collisions[nm] = val
	 self.collisions[nm]:setName(nm)
	 val:setSpeciesName(self.name)
	 val:fullInit(tbl) -- Initialize collisions
      end
   end

   self.positivity        = xsys.pickBool(tbl.applyPositivity, false)
   self.positivityDiffuse = xsys.pickBool(tbl.positivityDiffuse, self.positivity)
   self.positivityRescale = xsys.pickBool(tbl.positivityRescale, false)
   
   -- for GK only: flag for gyroaveraging.
   self.gyavg = xsys.pickBool(tbl.gyroaverage, false)

   self.ioMethod           = "MPI"
   self.distIoFrame        = 0 -- Frame number for distribution function.
   self.diagIoFrame        = 0 -- Frame number for diagnostics.
   self.dynVecRestartFrame = 0 -- Frame number of restarts (for DynVectors only).
   self.cfl    =  0.1
   self.nGhost = 1   -- Default is 1 ghost-cell in each direction.

   self.tCurr = 0.0

   self.integratedMomentsTime = 0.0 -- Timer for integrated moments.
   self.bcTime = 0.0   -- Timer for BCs.

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
   for _, c in pairs(self.collisions) do c:setCfl(cfl) end   
end
function KineticSpecies:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end
function KineticSpecies:setConfBasis(basis)
   self.confBasis = basis
   for _, c in pairs(self.collisions) do c:setConfBasis(basis) end
   for _, src in pairs(self.sources) do src:setConfBasis(basis) end
end
function KineticSpecies:setConfGrid(grid)
   self.confGrid = grid
   for _, c in pairs(self.collisions) do c:setConfGrid(grid) end
   for _, src in pairs(self.sources) do src:setConfGrid(grid) end
end

function KineticSpecies:createGrid(confGridIn)
   local confGrid = assert(confGridIn or self.confGrid, "KineticSpecies:createGrid ... must pass in confGrid or call setConfGrid prior to createGrid") 
 
   self.cdim = confGrid:ndim()
   self.ndim = self.cdim+self.vdim

   -- Create decomposition.
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, confGrid:cuts(d)) end
   for d = 1, self.vdim do table.insert(decompCuts, self.decompCuts[d]) end
   self.decomp = DecompRegionCalc.CartProd {
      cuts      = decompCuts,
      useShared = self.useShared,
   }

   -- Create computational domain.
   local lower, upper, cells = {}, {}, {}
   for d = 1, self.cdim do
      table.insert(lower, confGrid:lower(d))
      table.insert(upper, confGrid:upper(d))
      table.insert(cells, confGrid:numCells(d))
   end
   for d = 1, self.vdim do
      table.insert(lower, self.lower[d])
      table.insert(upper, self.upper[d])
      table.insert(cells, self.cells[d])
   end

   local GridConstructor = Grid.RectCart
   local coordinateMap = {} -- Table of functions
   -- Construct comp -> phys mappings if they exist
   if self.coordinateMap or confGrid:getMappings() then
      if confGrid:getMappings() and self.coordinateMap then
         for d = 1, self.cdim do
            lower[d], upper[d] = confGrid:logicalLower(d), confGrid:logicalUpper(d)
            table.insert(coordinateMap, confGrid:getMappings(d))
         end
         for d = 1, self.vdim do
            table.insert(coordinateMap, self.coordinateMap[d])
         end
      elseif confGrid:getMappings() then
         for d = 1, self.cdim do
            lower[d], upper[d] = confGrid:logicalLower(d), confGrid:logicalUpper(d)
            table.insert(coordinateMap, confGrid:getMappings(d))
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
      lower         = lower,
      upper         = upper,
      cells         = cells,
      periodicDirs  = confGrid:getPeriodicDirs(),
      decomposition = self.decomp,
      mappings      = coordinateMap,
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

   -- Output of grid file is placed here because as the file name is associated
   -- with a species, we wish to save the basisID and polyOrder in it. But these
   -- can only be extracted from self.basis after this is created.
   if self.grid:getMappings() then
      local metaData = {polyOrder = self.basis:polyOrder(),
                        basisType = self.basis:id(),
                        charge    = self.charge,
                        mass      = self.mass,
                        grid      = GKYL_OUT_PREFIX .. "_" .. self.name .. "_grid.bp"}
      self.grid:write(self.name .. "_grid.bp", 0.0, metaData)
   end
end

-- Field allocation in the species objects should be performed with one
-- of the following four functions instead of calling DataStruct directly.
function KineticSpecies:allocCartField(grid,nComp,ghosts,metaData)
   local f = DataStruct.Field {
      onGrid        = grid,
      numComponents = nComp,
      ghost         = ghosts,
      metaData      = metaData,
   }
   f:clear(0.0)
   return f
end
function KineticSpecies:allocDistf()
   local metaData = {polyOrder = self.basis:polyOrder(),
                     basisType = self.basis:id(),
                     charge    = self.charge,
                     mass      = self.mass,
                     grid      = GKYL_OUT_PREFIX .. "_" .. self.name .. "_grid.bp"}
   return self:allocCartField(self.grid,self.basis:numBasis(),{self.nGhost,self.nGhost},metaData)
end
function KineticSpecies:allocMoment()
   local metaData = {polyOrder = self.basis:polyOrder(),
                     basisType = self.basis:id(),
                     charge    = self.charge,
                     mass      = self.mass,}
   return self:allocCartField(self.confGrid,self.confBasis:numBasis(),{self.nGhost,self.nGhost},metaData)
end
function KineticSpecies:allocVectorMoment(dim)
   local metaData = {polyOrder = self.basis:polyOrder(),
                     basisType = self.basis:id(),
                     charge    = self.charge,
                     mass      = self.mass,}
   return self:allocCartField(self.confGrid,dim*self.confBasis:numBasis(),{self.nGhost,self.nGhost},metaData)
end

-- Various functions to apply BCs of different types.
function KineticSpecies:bcAbsorbFunc(dir, tm, idxIn, fIn, fOut)
   -- Note that for bcAbsorb there is no operation on fIn,
   -- so skinLoop (which determines indexing of fIn) does not matter 
   for i = 1, self.basis:numBasis() do fOut[i] = 0.0 end
end
function KineticSpecies:bcOpenFunc(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "pointwise".
   self.basis:flipSign(dir, fIn, fOut)
end
function KineticSpecies:bcCopyFunc(dir, tm, idxIn, fIn, fOut)
   -- Requires skinLoop = "pointwise".
   for i = 1, self.basis:numBasis() do fOut[i] = fIn[i] end
end

-- Function to construct a BC updater.
function KineticSpecies:makeBcUpdater(dir, vdir, edge, bcList, skinLoop,
				      evaluateFn)
   return Updater.Bc {
      onGrid             = self.grid,
      boundaryConditions = bcList,
      dir                = dir,
      vdir               = vdir,
      edge               = edge,
      skinLoop           = skinLoop,
      cdim               = self.cdim,
      vdim               = self.vdim,
      basis              = self.basis,
      evaluate           = evaluateFn,
      evolveFn           = self.evolveFnBC,
      feedback           = self.feedbackBC,
      confBasis          = self.confBasis,
      confGrid           = self.confGrid,
   }
end

function KineticSpecies:createBCs()
   -- Functions to make life easier while reading in BCs to apply.
   -- Note: appendBoundaryConditions defined in sub-classes.
   local function handleBc(dir, bc)
      if bc[1] then
	 self:appendBoundaryConditions(dir, "lower", bc[1])
      end
      if bc[2] then
	 self:appendBoundaryConditions(dir, "upper", bc[2])
      end
   end

   -- Add various BCs to list of BCs to apply.
   handleBc(1, self.bcx)
   handleBc(2, self.bcy)
   handleBc(3, self.bcz)

   -- Calculate external boundary condition if applicable
   if self.tbl.computeExternalBC then self:initExternalBC() end
end

function KineticSpecies:createSolver(externalField)
   -- Create solvers for collisions.
   for _, c in pairs(self.collisions) do c:createSolver(externalField) end
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
   if self.positivityRescale then
      self.prePosM0     = self:allocMoment()
      self.postPosM0    = self:allocMoment()
      self.delPosM0     = self:allocMoment()
      self.intDelPosM0  = DataStruct.DynVector{numComponents = 1}
      self.calcIntPosM0 = Updater.CartFieldIntegratedQuantCalc {
	    onGrid        = self.confGrid,
	    basis         = self.confBasis,
	    numComponents = 1,
	    quantity      = "V",
	    timeIntegrate = timeIntegrate,
      }
   end
end

function KineticSpecies:alloc(nRkDup)
   -- Allocate fields needed in RK update.
   self.distf = {}
   for i = 1, nRkDup do
      self.distf[i] = self:allocDistf()
      self.distf[i]:clear(0.0)
   end

   if self.positivity then
      self.fDelPos = {}
      for i = 1, nRkDup do
         self.fDelPos[i] = self:allocDistf()
         self.fDelPos[i]:clear(0.0)
      end
   end

   -- Create Adios object for field I/O.
   self.distIo = AdiosCartFieldIo {
      elemType   = self.distf[1]:elemType(),
      method     = self.ioMethod,
      writeGhost = self.writeGhost,
      metaData   = {polyOrder = self.basis:polyOrder(),
                    basisType = self.basis:id(),
                    charge    = self.charge,
                    mass      = self.mass,    
                    grid      = GKYL_OUT_PREFIX .. "_" .. self.name .. "_grid.bp"},
   }

   if self.positivity then self.fPos = self:allocDistf() end

   self.fPrev = self:allocDistf()
   self.fPrev:clear(0.0)

   -- Array with one component per cell to store cflRate in each cell.
   self.cflRateByCell = self:allocCartField(self.grid, 1, {1,1})
   self.cflRatePtr    = self.cflRateByCell:get(1)
   self.cflRateIdxr   = self.cflRateByCell:genIndexer()
   self.dtGlobal      = ffi.new("double[2]")

   self:createBCs()

   -- Create a table of flags to indicate whether moments have been computed.
   -- At first we consider 6 flags: coupling moments (M0, M1i, M2)
   -- boundary corrections (m1Correction, m2Correction), star moments
   -- (m0Star, m1Star, m2Star), self primitive moments (uSelf, vtSqSelf),
   -- cross primitive moments (uCross, vtSqCross), and spatially varying
   -- cross-species collisionality (varNuXCross).
   self.momentFlags = {}
   for iF = 1,4 do self.momentFlags[iF] = false end
   -- The fifth and sixth entries need a table to store 
   -- a flag for each pair of species colliding.
   self.momentFlags[5] = {}  -- Corresponds to uCross and vtSqCross.
   self.momentFlags[6] = {}  -- Corresponds to varNuXCross.
end

-- Note: do not call applyBc here. It is called later in initialization sequence.
function KineticSpecies:initDist(extField)
   if self.randomseed then 
      math.randomseed(self.randomseed) 
   else
      math.randomseed(47*Mpi.Comm_rank(Mpi.COMM_WORLD)+os.time())
   end

   local syncPeriodicDirs = true
   if self.fluctuationBCs then syncPeriodicDirs = false end

   local initCnt, backgroundCnt = 0, 0
   for nm, pr in lume.orderedIter(self.projections) do
      pr:fullInit(self)
      pr:advance(0.0, {extField}, {self.distf[2]})
      -- This barrier is needed as when using MPI-SHM some
      -- processes will get to accumulate before projection is finished.
      Mpi.Barrier(self.grid:commSet().sharedComm)
      if string.find(nm,"init") then
	 self.distf[1]:accumulate(1.0, self.distf[2])
	 initCnt = initCnt + 1
         if pr.scaleWithSourcePower then self.scaleInitWithSourcePower = true end
      end
      if string.find(nm,"background") then
	 if not self.f0 then self.f0 = self:allocDistf() end
	 self.f0:accumulate(1.0, self.distf[2])
	 self.f0:sync(syncPeriodicDirs)
	 backgroundCnt = backgroundCnt + 1
      end
      -- if pr.isReservoir then
      --    if not self.fReservoir then 
      --       self.fReservoir = self:allocDistf()
      --    end
      --    self.fReservoir:accumulate(1.0, self.distf[2])
      -- end
   end
   
   -- Set up profile function for species sources.
   for nm, src in lume.orderedIter(self.sources) do src:createSolver(self, extField) end

   assert(initCnt>0, string.format("KineticSpecies: Species '%s' not initialized!", self.name))
   if self.f0 and backgroundCnt == 0 then self.f0:copy(self.distf[1]) end

   if self.fluctuationBCs then 
      assert(backgroundCnt > 0, "KineticSpecies: must specify an initial background distribution with 'background' in order to use fluctuation-only BCs") 
   end

   self.distf[2]:clear(0.0)
   
   self:setActiveRKidx(1)

   -- Calculate initial density averaged over simulation domain.
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

function KineticSpecies:setActiveRKidx(rkIdx)
   self.activeRKidx = rkIdx
end

function KineticSpecies:rkStepperFields() return self.distf end

function KineticSpecies:getDistF(rkIdx)
   if rkIdx == nil then
      return self:rkStepperFields()[self.activeRKidx]
   else
      return self:rkStepperFields()[rkIdx]
   end
end

function KineticSpecies:copyRk(outIdx, aIdx)
   self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx])

   if self.hasNonPeriodicBc and self.boundaryFluxDiagnostics then
      for _, bc in ipairs(self.boundaryConditions) do
         bc:getBoundaryFluxFields()[outIdx]:copy(bc:getBoundaryFluxFields()[aIdx])
      end
   end

   if self.positivity then
      self.fDelPos[outIdx]:copy(self.fDelPos[aIdx])
   end
end

function KineticSpecies:combineRk(outIdx, a, aIdx, ...)
   local args = {...} -- Package up rest of args as table.
   local nFlds = #args/2
   self:rkStepperFields()[outIdx]:combine(a, self:rkStepperFields()[aIdx])
   for i = 1, nFlds do -- Accumulate rest of the fields.
      self:rkStepperFields()[outIdx]:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]])
   end

   if self.hasNonPeriodicBc and self.boundaryFluxDiagnostics then
      for _, bc in ipairs(self.boundaryConditions) do
         bc:getBoundaryFluxFields()[outIdx]:combine(a, bc:getBoundaryFluxFields()[aIdx])
         for i = 1, nFlds do -- Accumulate rest of the fields.
            bc:getBoundaryFluxFields()[outIdx]:accumulate(args[2*i-1], bc:getBoundaryFluxFields()[args[2*i]])
         end
      end
   end

   if self.positivity then
      self.fDelPos[outIdx]:combine(a, self.fDelPos[aIdx])
      for i = 1, nFlds do -- Accumulate rest of the fields.
         self.fDelPos[outIdx]:accumulate(args[2*i-1], self.fDelPos[args[2*i]])
      end
   end
end

function KineticSpecies:suggestDt()
   if not self.evolve then return GKYL_MAX_DOUBLE end

   local dtSuggested = math.min(self.cfl/self.cflRateByCell:reduce('max')[1], GKYL_MAX_DOUBLE)

   -- If dtSuggested == GKYL_MAX_DOUBLE, it is likely because of NaNs. 
   -- If so, return 0 so that no timestep is taken, and we will abort the simulation.
   if dtSuggested == GKYL_MAX_DOUBLE then dtSuggested = 0.0 end

   return dtSuggested
end

function KineticSpecies:setDtGlobal(dtGlobal) self.dtGlobal[0] = dtGlobal end

function KineticSpecies:clearCFL()
   -- Clear cflRateByCell for next cfl calculation.
   self.cflRateByCell:clear(0.0)
end

function KineticSpecies:clearMomentFlags(species)
   -- Clear the momentFlags table to indicate that moments (and other
   -- quantities) need to be computed again.
   for iF = 1,4 do self.momentFlags[iF] = false end
   for sN, _ in lume.orderedIter(species) do
      if sN ~= self.name then
         self.momentFlags[5][sN] = false
         self.momentFlags[6][sN] = false
      end
   end
end

function KineticSpecies:checkPositivity(tCurr, idx)
  local status = true
  if self.positivity then
     --status = self.posChecker:advance(tCurr, {self:rkStepperFields()[idx]}, {})
  end
  return status
end

function KineticSpecies:applyBcIdx(tCurr, idx, isFirstRk)
  if self.positivityDiffuse then
     self.fDelPos[idx]:combine(-1.0, self:rkStepperFields()[idx])
     self.posRescaler:advance(tCurr, {self:rkStepperFields()[idx]}, {self:rkStepperFields()[idx]}, true, isFirstRk)
     self.fDelPos[idx]:accumulate(1.0, self:rkStepperFields()[idx])
  end
  self:applyBc(tCurr, self:rkStepperFields()[idx])
  if self.positivity then
     self:checkPositivity(tCurr, idx)
  end
  if self.positivityRescale then
     self.numDensityCalc:advance(tCurr, {self:rkStepperFields()[idx]}, {self.prePosM0})
     
     self.posRescaler:advance(tCurr, {self:rkStepperFields()[idx]}, {self:rkStepperFields()[idx]}, false)
     
     self.numDensityCalc:advance(tCurr, {self:rkStepperFields()[idx]}, {self.postPosM0})
     self.delPosM0:combine(1.0, self.postPosM0, -1.0, self.prePosM0)
     self.calcIntPosM0:advance(tCurr, {self.delPosM0}, {self.intDelPosM0})
  end
end

function KineticSpecies:applyBc(tCurr, fIn)
   -- fIn is total distribution function.
   local tmStart = Time.clock()

   if self.evolve then 
      local syncPeriodicDirsTrue = true

      if self.fluctuationBCs then
         -- If fluctuation-only BCs, subtract off background before applying BCs.
         fIn:accumulate(-1.0, self.f0)
      end

      -- Apply non-periodic BCs (to only fluctuations if fluctuation BCs).
      if self.hasNonPeriodicBc then
         if self.feedbackBC then
            for _, bc in ipairs(self.boundaryConditions) do
               bc:advance(tCurr, {fIn}, {fIn})
            end
         else
            for _, bc in ipairs(self.boundaryConditions) do
               bc:advance(tCurr, {}, {fIn})
            end
         end
      end

      -- Apply periodic BCs (to only fluctuations if fluctuation BCs)
      fIn:sync(syncPeriodicDirsTrue)

      if self.fluctuationBCs then
         -- Put back together total distribution
         fIn:accumulate(1.0, self.f0)
 
         -- Update ghosts in total distribution, without enforcing periodicity.
         fIn:sync(not syncPeriodicDirsTrue)
      end
   end

   self.bcTime = self.bcTime + (Time.clock()-tmStart)
end

function KineticSpecies:createDiagnostics()
   -- Set up weak multiplication and division operators.
   self.weakMultiplication = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Multiply",
      onGhosts  = true,
   }
   self.weakDivision = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
      onGhosts  = true,
   }
end

function KineticSpecies:calcDiagnosticMoments(tm)
   local fIn   = self:rkStepperFields()[1]
   local tCurr = tm
   if self.f0 and self.perturbedMoments then 
      fIn:accumulate(-1, self.f0)
      tCurr = -tm-1   -- Setting tCurr = -tm-1 forces the updater to recompute moments on this timestep.
   end
   for i, mom in pairs(self.diagnosticMoments) do
      self.diagnosticMomentUpdaters[mom]:advance(tCurr, {fIn}, {self.diagnosticMomentFields[mom]})
      -- Remove geometric jacobian factor.
      if self.jacobGeoInv then
         self.weakMultiplication:advance(0.0, {self.diagnosticMomentFields[mom], self.jacobGeoInv},
                                              {self.diagnosticMomentFields[mom]})
      end
   end
   if self.f0 and self.perturbedMoments then fIn:accumulate(1, self.f0) end
end

-- Some species-specific parts, but this function still gets called.
function KineticSpecies:calcDiagnosticWeakMoments(tm, weakMoments, bc)
   local fIn
   local tCurr = tm
   local label = ""
   if bc then
      label = bc:label() 
      fIn   = bc:getBoundaryFluxRate()
   else
      fIn = self:rkStepperFields()[1]
      if self.f0 and self.perturbedMoments then
         fIn:accumulate(-1, self.f0)
         tCurr = -tm-1   -- Setting tCurr = -tm-1 forces the updater to recompute moments on this timestep.
      end
   end

   for i, mom in ipairs(weakMoments) do
      self.diagnosticMomentUpdaters[mom..label].advance(self, tCurr, {fIn}, {self.diagnosticMomentFields[mom..label]})
   end

   if bc==nil and self.f0 and self.perturbedMoments then fIn:accumulate(1, self.f0) end
end

function KineticSpecies:calcDiagnosticBoundaryFluxMoments(tm)
   for _, bc in ipairs(self.boundaryConditions) do
      for i, mom in ipairs(self.diagnosticBoundaryFluxMoments) do
         self.diagnosticMomentUpdaters[mom..bc:label()]:advance(
            tm, {bc:getBoundaryFluxRate()}, {self.diagnosticMomentFields[mom..bc:label()]})
      end 
   end
end

-- Species-specific.
function KineticSpecies:calcDiagnosticIntegratedMoments(tm) end

function KineticSpecies:calcAndWriteDiagnosticMoments(tm)
    self:calcDiagnosticMoments(tm)
    if self.diagnosticWeakMoments then 
       self:calcDiagnosticWeakMoments(tm, self.diagnosticWeakMoments)
    end
    if self.diagnosticBoundaryFluxMoments then
       self:calcDiagnosticBoundaryFluxMoments(tm)
    end
    if self.diagnosticWeakBoundaryFluxMoments then
       for _, bc in ipairs(self.boundaryConditions) do
          self:calcDiagnosticWeakMoments(tm, self.diagnosticWeakBoundaryFluxMoments, bc)
       end
    end

    for i, mom in ipairs(self.requestedDiagnosticMoments) do
       local fldNm = mom
       if self.diagnosticMomentFields[fldNm]==nil then
          -- Cross-species diagnostics have another species name appended to them, so the name
          -- in requestedDiagnosticMoments does not match the name in diagnosticMomentFields.
          for nm, _ in pairs(self.diagnosticMomentFields) do
             if string.find(nm, "uCross") or string.find(nm, "vtSqCross") or
                string.find(nm, "GkUparCross") or string.find(nm, "GkVtSqCross") then fldNm=nm end
          end
       end
       self.diagnosticMomentFields[fldNm]:write(
          string.format("%s_%s_%d.bp", self.name, fldNm, self.diagIoFrame), tm, self.diagIoFrame, self.writeGhost)
    end

    for i, mom in ipairs(self.requestedDiagnosticBoundaryFluxMoments) do
       for _, bc in ipairs(self.boundaryConditions) do
          self.diagnosticMomentFields[mom..bc:label()]:write(
             string.format("%s_%s_%d.bp", self.name, mom..bc:label(), self.diagIoFrame), tm, self.diagIoFrame, self.writeGhost)
       end
    end

    -- Write integrated moments.
    for i, mom in ipairs(self.diagnosticIntegratedMoments) do
       -- These moments are handled in src:writeDiagnosticIntegratedMoments 
       if not (mom == "intSrcM0" or mom == "intSrcM1" or mom == "intSrcM2" or mom == "intSrcKE") then
          self.diagnosticIntegratedMomentFields[mom]:write(
             string.format("%s_%s.bp", self.name, mom), tm, self.diagIoFrame)
       end
    end

    -- Write source integrated diagnostics.
    for nm, src in lume.orderedIter(self.sources) do
       src:writeDiagnosticIntegratedMoments(tm, self.diagIoFrame)
    end

    for i, mom in ipairs(self.diagnosticIntegratedBoundaryFluxMoments) do
       for _, bc in ipairs(self.boundaryConditions) do
          self.diagnosticIntegratedMomentFields[mom..bc:label()]:write(
             string.format("%s_%s.bp", self.name, mom..bc:label()), tm, self.diagIoFrame)
       end
    end

    -- Write ionization diagnostics
    if self.calcReactRate then
       local sourceIz = self.collisions[self.collNmIoniz]:getIonizSrc()
       self.fMaxwellIz:write(string.format("%s_fMaxwell_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame, self.writeSkin)
       self.vtSqIz:write(string.format("%s_vtSqIz_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame, self.writeSkin)
       self.voronovReactRate:write(string.format("%s_coefIz_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame, self.writeSkin)
       sourceIz:write(string.format("%s_sourceIz_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame, self.writeSkin)
       -- include dynvector for zeroth vector of ionization source
       tmStart = Time.clock()
       self.intSrcIzM0:write(
          string.format("%s_intSrcIzM0.bp", self.name), tm, self.diagIoFrame)
       self.integratedMomentsTime = self.integratedMomentsTime + Time.clock() - tmStart       
    end

    if self.calcIntSrcIz then
       tmStart = Time.clock()
       local sourceIz = self.collisions[self.collNmIoniz]:getIonizSrc()
       sourceIz:write(string.format("%s_sourceIz_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame, self.writeSkin)
       self.intSrcIzM0:write(
          string.format("%s_intSrcIzM0.bp", self.name), tm, self.diagIoFrame)
       self.integratedMomentsTime = self.integratedMomentsTime + Time.clock() - tmStart    
    end
       
    -- Write CX diagnostics
    if self.calcCXSrc then
       self.vSigmaCX:write(string.format("%s_vSigmaCX_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame, self.writeSkin)
       self.collisions[self.collNmCX].sourceCX:write(string.format("%s_sourceCX_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame, self.writeSkin)
    end

    -- Write recycling diagnostics
    if self.hasRecycleBcs then
        for _, bc in ipairs(self.boundaryConditions) do
	   label = bc:label()
	   if self.cdim == 1 or (self.cdim == 3 and string.match(label,"Z")) then
	      wlabel = (label):gsub("Flux","")
	      self.recycleCoef[label]:write(string.format("%s%s_%d.bp", 'recycleCoef', wlabel, self.diagIoFrame), tm, self.diagIoFrame, false)
	      self.recycleDistF[label]:write(string.format("%s_%s%s_%d.bp", self.name, 'recycleDistF', wlabel, self.diagIoFrame), tm, self.diagIoFrame, false)
	   end
	end
    end

    -- Vlasov positivity diagnostics
    if self.positivityRescale then
       self.intDelPosM0:write( string.format("%s_intPosDelM0.bp", self.name), tm, self.diagIoFrame)
    end
end

function KineticSpecies:isEvolving()
   return self.evolve
end

function KineticSpecies:write(tm, force)
   if self.evolve then
      if self.hasNonPeriodicBc and self.boundaryFluxDiagnostics and tm > 0 then
         for _, bc in ipairs(self.boundaryConditions) do
            -- compute boundary flux rate ~ (fGhost_new - fGhost_old)/dt
            bc:getBoundaryFluxRate():combine(1.0/self.dtGlobal[0], bc:getBoundaryFluxFields()[1], -1.0/self.dtGlobal[0], bc:getBoundaryFluxFieldPrev())
            bc:getBoundaryFluxFieldPrev():copy(bc:getBoundaryFluxFields()[1])
         end
      end

      local tmStart = Time.clock()
      -- Compute integrated diagnostics.
      if self.calcIntQuantTrigger(tm) then
         self:calcDiagnosticIntegratedMoments(tm)
      end
      self.integratedMomentsTime = self.integratedMomentsTime + Time.clock() - tmStart

      -- Only write stuff if triggered.
      if self.distIoTrigger(tm) or force then
         self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, self.distIoFrame), tm, self.distIoFrame)
         if self.f0 then
            if tm == 0.0 then
	       self.f0:write(string.format("%s_f0_%d.bp", self.name, self.distIoFrame), tm, self.distIoFrame, true)
	    end
            self.distf[1]:accumulate(-1, self.f0)
            self.distIo:write(self.distf[1], string.format("%s_f1_%d.bp", self.name, self.distIoFrame), tm, self.distIoFrame)
            self.distf[1]:accumulate(1, self.f0)
         end
         for _, src in lume.orderedIter(self.sources) do src:write(tm, self.distIoFrame) end
	 self.distIoFrame = self.distIoFrame+1
      end


      if self.diagIoTrigger(tm) or force then
         -- Compute moments and write them out.
         self:calcAndWriteDiagnosticMoments(tm)

         if self.evolveCollisions then
            for _, c in pairs(self.collisions) do
               c:write(tm, self.diagIoFrame)
            end
         end

         if self.positivityDiffuse then
            self.posRescaler:write(tm, self.diagIoFrame, self.name)
         end

         self.diagIoFrame = self.diagIoFrame+1
      end
   else
      -- If not evolving species, don't write anything except initial conditions.
      if self.distIoFrame == 0 then

         local tmStart = Time.clock()
         self:calcDiagnosticIntegratedMoments(tm)   -- Compute integrated diagnostics.
         self.integratedMomentsTime = self.integratedMomentsTime + Time.clock() - tmStart

	 self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, 0), tm, 0)

	 -- Compute moments and write them out.
	 self:calcAndWriteDiagnosticMoments(tm)
      end
      self.distIoFrame = self.distIoFrame+1
   end

end

function KineticSpecies:writeRestart(tm)
   -- (The final "true/false" determines writing of ghost cells).
   local writeGhost = false
   if self.hasSheathBCs or self.fluctuationBCs then writeGhost = true end

   self.distIo:write(self.distf[1], string.format("%s_restart.bp", self.name), tm, self.distIoFrame, writeGhost)

   for i, mom in pairs(self.diagnosticMoments) do
      self.diagnosticMomentFields[mom]:write(
	 string.format("%s_%s_restart.bp", self.name, mom), tm, self.diagIoFrame, false)
   end   

   -- Write restart files for integrated moments. Note: these are only needed for the rare case that the
   -- restart write frequency is higher than the normal write frequency from nFrame.
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
      -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
      self.diagnosticIntegratedMomentFields[mom]:write(
         string.format("%s_%s_restart.bp", self.name, mom), tm, self.dynVecRestartFrame, false, false)
   end

   if self.calcReactRate then
      self.intSrcIzM0:write(
	 string.format("%s_intSrcIzM0_restart.bp", self.name), tm, self.dynVecRestartFrame, false, false)
   end
   if self.calcIntSrcIz then
      self.intSrcIzM0:write(
	 string.format("%s_intSrcIzM0_restart.bp", self.name), tm, self.dynVecRestartFrame, false, false)
   end

   self.dynVecRestartFrame = self.dynVecRestartFrame + 1
end

function KineticSpecies:readRestart()
   local readGhost = false
   if self.hasSheathBCs or self.fluctuationBCs then readGhost = true end

   local tm, fr = self.distIo:read(self.distf[1], string.format("%s_restart.bp", self.name), readGhost)
   self.distIoFrame = fr -- Reset internal frame counter.

   -- Set ghost cells.
   self.distf[1]:sync()

   -- Apply BCs (unless skin cells have been read because of special BCs).
   if not self.hasSheathBCs and not self.fluctuationBCs then 
      self:applyBc(tm, self.distf[1]) 
   end 
   
   for i, mom in ipairs(self.diagnosticIntegratedMoments) do
      self.diagnosticIntegratedMomentFields[mom]:read(
         string.format("%s_%s_restart.bp", self.name, mom))
   end
   
   if self.calcReactRate then
      self.intSrcIzM0:read(
	 string.format("%s_intSrcIzM0_restart.bp", self.name))
   end
   if self.calcIntSrcIz then
      self.intSrcIzM0:read(
	 string.format("%s_intSrcIzM0_restart.bp", self.name))
   end

   for i, mom in pairs(self.diagnosticMoments) do
      local _, dfr = self.diagnosticMomentFields[mom]:read(
         string.format("%s_%s_restart.bp", self.name, mom))
      self.diagIoFrame = dfr -- Reset internal diagnostic IO frame counter.
   end

   -- Iterate triggers.
   self.distIoTrigger(tm)
   self.diagIoTrigger(tm)
   
   return tm
end

-- Timers.
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
