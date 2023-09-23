-- Gkyl ------------------------------------------------------------------------
--
-- Vlasov species
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto            = require "Lib.Proto"
local SpeciesBase      = require "App.Species.SpeciesBase"
local ProjectionBase   = require "App.Projection.ProjectionBase"
local CollisionsBase   = require "App.Collisions.CollisionsBase"
local SourceBase       = require "App.Sources.SourceBase"
local BCsBase          = require "App.BCs.BCsBase"
local Projection       = require "App.Projection"
local Source           = require "App.Sources.VmSource"
local DiagsApp         = require "App.Diagnostics.SpeciesDiagnostics"
local VlasovDiags      = require "App.Diagnostics.VlasovDiagnostics"
local BasicBC          = require ("App.BCs.VlasovBasic").VlasovBasic
local Basis            = require "Basis"
local DataStruct       = require "DataStruct"
local Grid             = require "Grid"
local DecompRegionCalc = require "Lib.CartDecomp"
local LinearTrigger    = require "Lib.LinearTrigger"
local Mpi              = require "Comm.Mpi"
local Time             = require "Lib.Time"
local Updater          = require "Updater"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local xsys             = require "xsys"
local lume             = require "Lib.lume"
local ffi              = require "ffi"

local VlasovSpecies = Proto(SpeciesBase)

local SP_BC_EXTERN  = 4
local SP_BC_RECYCLE = 7
VlasovSpecies.bcExternal = SP_BC_EXTERN     -- Load external BC file.
VlasovSpecies.bcRecycle  = SP_BC_RECYCLE

-- ............. Backwards compatible treatment of BCs .....................--
-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB    = 1
local SP_BC_REFLECT   = 3
local SP_BC_COPY      = 5
-- AHH: This was 2 but seems that is unstable. So using plain copy.
local SP_BC_OPEN      = SP_BC_COPY
local SP_BC_ZEROFLUX  = 6
VlasovSpecies.bcCopy      = SP_BC_COPY       -- Copy stuff.
VlasovSpecies.bcAbsorb    = SP_BC_ABSORB     -- Absorb all particles.
VlasovSpecies.bcOpen      = SP_BC_OPEN       -- Zero gradient.
VlasovSpecies.bcReflect   = SP_BC_REFLECT    -- Specular reflection.
VlasovSpecies.bcZeroFlux  = SP_BC_ZEROFLUX

function VlasovSpecies:makeBcApp(bcIn, dir, edge)
   local bcOut
   if type(bcIn) == "function" then
      bcOut = BasicBC{kind="function", bcFunction=bcIn, diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_COPY then
      print("VlasovSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="copy", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_ABSORB then
      print("VlasovSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="absorb", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_OPEN then
      print("VlasovSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      -- AHH: open seems unstable. So using plain copy.
      bcOut = BasicBC{kind="copy", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_REFLECT then
      print("VlasovSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="reflect", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_ZEROFLUX or bcIn.tbl.kind=="zeroFlux" then
      bcOut = "zeroFlux"
      table.insert(self.zeroFluxDirections, dir)
   end
   return bcOut
end

-- ............. End of backwards compatibility for BCs .....................--

-- This ctor stores what is passed to it and defers most of the
-- construction to the fullInit() method below.
-- We also place here the things we want every species to know about
-- every other species (when we parallelize over species).
function VlasovSpecies:init(tbl) self.tbl = tbl end

-- Function to create basis functions.
local function createBasis(nm, cdim, vdim, polyOrder)
   local ndim = cdim+vdim
   if nm == "serendipity" then
      if polyOrder == 1 then
         return Basis.CartModalHybrid { cdim = cdim, vdim = vdim }
      else
         return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
      end
   elseif nm == "tensor" then
      return Basis.CartModalTensor { ndim = ndim, polyOrder = polyOrder }
   end
end

-- Function to create velocity space basis functions.
local function createVelBasis(nm, vdim, polyOrder)
   if nm == "serendipity" then
      if polyOrder == 1 then
         return Basis.CartModalSerendipity { ndim = vdim, polyOrder = 2 }
      else
         return Basis.CartModalSerendipity { ndim = vdim, polyOrder = polyOrder }
      end
   elseif nm == "tensor" then
      return Basis.CartModalTensor { ndim = vdim, polyOrder = polyOrder }
   end
end

function VlasovSpecies:createBasis(nm, polyOrder)
   self.basis = createBasis(nm, self.cdim, self.vdim, polyOrder)
   self.velBasis = createVelBasis(nm, self.vdim, polyOrder)
   for _, c in lume.orderedIter(self.collisions) do c:setPhaseBasis(self.basis) end

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

function VlasovSpecies:alloc(nRkDup)
   -- Allocate distribution function and other quantities.

   -- Allocate fields needed in RK update.
   self.distf = {}
   for i = 1, nRkDup do
      self.distf[i] = self:allocDistf()
      self.distf[i]:clear(0.0)
   end
   self:setActiveRKidx(1)

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

   -- Array with one component per cell to store cflRate in each cell.
   self.cflRateByCell = self:allocCartField(self.grid, 1, {1,1})
   self.dtGlobal      = ffi.new("double[2]")
   self.dtGlobal[0], self.dtGlobal[1] = 1.0, 1.0   -- Temporary value (so diagnostics at t=0 aren't inf).

   -- Allocate fields to store coupling moments (for use in coupling
   -- to field and collisions).
   self.numDensity  = self:allocMoment()
   self.momDensity  = self:allocVectorMoment(self.vdim)
   self.fiveMoments = self:allocVectorMoment(self.vdim+2)

   self.ptclEnergyAux = self:allocMoment()

   -- Special relativistic only arrays
   -- Allocate fields to store velocity space arrays:
   -- p/gamma (relativistic velocity)
   -- gamma (particle Lorentz boost factor = sqrt(1 + p^2))
   -- gamma_inv = 1/gamma
   self.p_over_gamma  = self:allocVectorVelMoment(self.vdim)
   self.gamma  = self:allocVelMoment()
   self.gamma_inv  = self:allocVelMoment()

   -- Allocate additional auxiliary fields needed for certain relativistic moments
   -- V_drift (bulk velocity)
   -- GammaV2 (bulk velocity Lorentz boost factor squared = 1/(1 - V_drift^2))
   -- GammaV_inv = 1/GammaV = sqrt(1 - V_drift^2)
   self.V_drift = self:allocVectorMoment(self.vdim)
   self.GammaV2 = self:allocMoment()
   self.GammaV_inv = self:allocMoment()

   -- Allocate field for external forces if any.
   if self.hasExtForce then 
      self.vExtForce = self:allocVectorMoment(self.vdim)
      self.vExtFptr, self.vExtFidxr = self.vExtForce:get(1), self.vExtForce:genIndexer()
   end
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VlasovSpecies:fullInit(appTbl)
   local tbl = self.tbl -- Previously store table.

   self.charge = tbl.charge or 1.0
   self.mass   = tbl.mass or 1.0
   self.lower, self.upper = tbl.lower, tbl.upper
   self.cells = tbl.cells
   self.vdim  = #self.cells -- Velocity dimensions.

   self.evolve              = xsys.pickBool(tbl.evolve, true) -- By default, evolve species.
   self.evolveCollisionless = xsys.pickBool(tbl.evolveCollisionless, self.evolve)

   assert(#self.lower == self.vdim, "'lower' must have " .. self.vdim .. " entries")
   assert(#self.upper == self.vdim, "'upper' must have " .. self.vdim .. " entries")
   self.coordinateMap = tbl.coordinateMap

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

   self.diagnostics = {}  -- Table in which we'll place diagnostic objects.

   -- Write ghost cells on boundaries of global domain (for BCs).
   self.writeGhost = xsys.pickBool(appTbl.writeGhost, false)

   -- Option to group diagnostics (i.e. it writes one file for all grid diags, and one file
   -- for all integrated diags) in all diagnostic Apps, rather than writing one file for each.
   self.groupDiags = appTbl.groupDiagnostics or false

   -- Get a random seed for random initial conditions.
   self.randomseed = tbl.randomseed

   -- Initialize table containing sources (if any).
   self.sources = {}
   for nm, val in pairs(tbl) do
      if SourceBase.is(val) or string.find(nm,"source") then
         if ProjectionBase.is(val) then val = self:projToSource(val) end
         self.sources[nm] = val
         val:setSpeciesName(self.name)
         val:setName(nm)   -- Do :setName after :setSpeciesName for sources.
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
   lume.setOrder(self.projections)  -- Save order in metatable to loop in the same order (w/ orderedIter, better for I/O).

   self.zeroFluxDirections = {}

   -- Read in boundary conditions.
   self.bcInDir = {{ }, { }, { }}   -- List of BCs to apply.
   if tbl.bcx then
      if tbl.bcx[1] == nil or tbl.bcx[2] == nil then assert(false, "VlasovSpecies: unsupported BC type") end
      self.bcInDir[1] = {tbl.bcx[1], tbl.bcx[2]}
   end
   if tbl.bcy then
      if tbl.bcy[1] == nil or tbl.bcy[2] == nil then assert(false, "VlasovSpecies: unsupported BC type") end
      self.bcInDir[2] = {tbl.bcy[1], tbl.bcy[2]}
   end
   if tbl.bcz then
      if tbl.bcz[1] == nil or tbl.bcz[2] == nil then assert(false, "VlasovSpecies: unsupported BC type") end
      self.bcInDir[3] = {tbl.bcz[1], tbl.bcz[2]}
   end
   -- Initialize boundary conditions.
   self.nonPeriodicBCs = {}
   local dirLabel  = {'X','Y','Z'}
   local edgeLabel = {'lower','upper'}
   for d, bcsTbl in ipairs(self.bcInDir) do
      for e, bcOb in ipairs(bcsTbl) do
         local goodBC = false
         local val    = bcOb
         if not BCsBase.is(val) then val = self:makeBcApp(bcOb, d, e) end
         if BCsBase.is(val) then
            local nm = 'bc'..dirLabel[d]..edgeLabel[e]
            self.nonPeriodicBCs[nm] = val
            val:setSpeciesName(self.name)
            val:setName(nm)   -- Do :setName after :setSpeciesName for BCs.
            val:setDir(d)
            val:setEdge(edgeLabel[e])
            val:fullInit(tbl)
            goodBC = true
         elseif val=="zeroFlux" then
            goodBC = true
         end
         assert(goodBC, "VlasovSpecies: bc not recognized.")
      end
   end
   lume.setOrder(self.nonPeriodicBCs)  -- Save order in metatable to loop in the same order (w/ orderedIter, better for I/O).

   -- Collisions.
   self.collisions = {}
   for nm, val in pairs(tbl) do
      if CollisionsBase.is(val) then
         self.collisions[nm] = val
         val:setSpeciesName(self.name)
         val:setName(nm) -- Do :setName after :setSpeciesName for collisions.
         val:fullInit(tbl) -- Initialize collisions
      end
   end
   lume.setOrder(self.collisions)

   -- If there is an external force, get the force function.
   self.hasExtForce = false
   if tbl.vlasovExtForceFunc then
      self.vlasovExtForceFunc = tbl.vlasovExtForceFunc
      self.hasExtForce = true
   end

   self.ioMethod           = "MPI"
   self.distIoFrame        = 0 -- Frame number for distribution function.
   self.diagIoFrame        = 0 -- Frame number for diagnostics.
   self.dynVecRestartFrame = 0 -- Frame number of restarts (for DynVectors only).
   self.cfl    =  0.1
   self.nGhost = 1   -- Default is 1 ghost-cell in each direction.

   self.tCurr = 0.0

   self.timers = {
      mom = 0.,   momcross = 0.,   advance = 0.,
      advancecross = 0.,   collisions = 0.,   collisionless = 0.,
      boundflux = 0.,  sources = 0.,   bc = 0.,
   }
end

function VlasovSpecies:setCfl(cfl)
   self.cfl = cfl
   for _, c in lume.orderedIter(self.collisions) do c:setCfl(cfl) end
end

function VlasovSpecies:setIoMethod(ioMethod) self.ioMethod = ioMethod end

function VlasovSpecies:setConfBasis(basis)
   self.confBasis = basis
   for _, c in lume.orderedIter(self.collisions) do c:setConfBasis(basis) end
   for _, src in lume.orderedIter(self.sources) do src:setConfBasis(basis) end
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:setConfBasis(basis) end
end
function VlasovSpecies:setConfGrid(grid)
   self.confGrid = grid
   for _, c in lume.orderedIter(self.collisions) do c:setConfGrid(grid) end
   for _, src in lume.orderedIter(self.sources) do src:setConfGrid(grid) end
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:setConfGrid(grid) end
end

function VlasovSpecies:createGrid(confGridIn)
   local confGrid = assert(confGridIn or self.confGrid, "VlasovSpecies:createGrid ... must pass in confGrid or call setConfGrid prior to createGrid")

   self.cdim = confGrid:ndim()
   self.ndim = self.cdim+self.vdim

   -- Create decomposition.
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, confGrid:cuts(d)) end
   for d = 1, self.vdim do table.insert(decompCuts, 1) end -- No decomposition in v-space.
   self.decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts, comm = confGrid:commSet().comm,
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
      lower     = lower,  periodicDirs  = confGrid:getPeriodicDirs(),
      upper     = upper,  decomposition = self.decomp,
      cells     = cells,  mappings      = coordinateMap,
      messenger = confGrid:getMessenger(),
   }

   for _, c in lume.orderedIter(self.collisions) do c:setPhaseGrid(self.grid) end

   -- Construct velocity space grid from phase space grid
   local dimsV = {}
   for d = 1, self.vdim do
      table.insert(dimsV, self.cdim+d)
   end
   -- Get the ingredients of the velocity space grid from the phase space grid
   local velGridIngr = self.grid:childGrid(dimsV)
   self.velGrid = GridConstructor {
      lower = velGridIngr.lower,  periodicDirs  = velGridIngr.periodicDirs,
      upper = velGridIngr.upper,  decomposition = velGridIngr.decomposition,
      cells = velGridIngr.cells,      
   }
end

-- Field allocation in the species objects should be performed with one
-- of the following functions instead of calling DataStruct directly.
function VlasovSpecies:allocCartField(grid,nComp,ghosts,metaData)
   local f = DataStruct.Field {
      onGrid        = grid,   ghost    = ghosts,
      numComponents = nComp,  metaData = metaData,
   }
   f:clear(0.0)
   return f
end
function VlasovSpecies:allocDistf()
   local metaData = {polyOrder = self.basis:polyOrder(),
                     basisType = self.basis:id(),
                     charge    = self.charge,
                     mass      = self.mass,
                     grid      = GKYL_OUT_PREFIX .. "_" .. self.name .. "_grid.bp"}
   return self:allocCartField(self.grid,self.basis:numBasis(),{self.nGhost,self.nGhost},metaData)
end
function VlasovSpecies:allocMoment()
   local metaData = {polyOrder = self.confBasis:polyOrder(),
                     basisType = self.confBasis:id(),
                     charge    = self.charge,
                     mass      = self.mass,}
   return self:allocCartField(self.confGrid,self.confBasis:numBasis(),{self.nGhost,self.nGhost},metaData)
end
function VlasovSpecies:allocVectorMoment(dim)
   local metaData = {polyOrder = self.confBasis:polyOrder(),
                     basisType = self.confBasis:id(),
                     charge    = self.charge,
                     mass      = self.mass,}
   return self:allocCartField(self.confGrid,dim*self.confBasis:numBasis(),{self.nGhost,self.nGhost},metaData)
end
-- Velocity space arrays (no ghost cells in velocity space arrays)
function VlasovSpecies:allocVelMoment()
   local metaData = {polyOrder = self.velBasis:polyOrder(),
                     basisType = self.velBasis:id(),
                     charge    = self.charge,
                     mass      = self.mass,}
   return self:allocCartField(self.velGrid,self.velBasis:numBasis(),{0,0},metaData)
end
function VlasovSpecies:allocVectorVelMoment(dim)
   local metaData = {polyOrder = self.velBasis:polyOrder(),
                     basisType = self.velBasis:id(),
                     charge    = self.charge,
                     mass      = self.mass,}
   return self:allocCartField(self.velGrid,dim*self.velBasis:numBasis(),{0,0},metaData)
end
function VlasovSpecies:allocIntMoment(comp)
   local metaData = {charge = self.charge,
                     mass   = self.mass,}
   local ncomp = comp or 1
   local f = DataStruct.DynVector {numComponents = ncomp,     writeRank = 0,
                                   metaData      = metaData,  comm      = self.confGrid:commSet().comm,}
   return f
end

function VlasovSpecies:createSolver(field, externalField)
   -- Set up weak multiplication and division operators.
   self.confWeakMultiply = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "Multiply",
      onGhosts  = true,
   }
   local dummyConfField = self:allocMoment()
   self.confWeakDivide = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "Divide",
      onRange   = dummyConfField:localRange(),  onGhosts = false,
   }
   self.confWeakDotProduct = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,  operation = "DotProduct",
      onGhosts  = true,
   }
   self.confPhaseWeakMultiply = Updater.CartFieldBinOp {
      weakBasis  = self.basis,  operation = "Multiply",
      fieldBasis = self.confBasis,
   }

   if self.evolve then
      self.applyBcFunc = function(tCurr, field, externalField, inIdx, outIdx)
         return VlasovSpecies["applyBcEvolve"](self, tCurr, field, externalField, inIdx, outIdx)
      end
   else
      self.applyBcFunc = function(tCurr, field, externalField, inIdx, outIdx)
         return VlasovSpecies["applyBcDontEvolve"](self, tCurr, field, externalField, inIdx, outIdx)
      end
   end

   -- Create solvers for collisions.
   for _, c in lume.orderedIter(self.collisions) do c:createSolver(self, externalField) end

   -- Create BC solvers.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:createSolver(self, field, externalField) end

   self.qbym = self.charge/self.mass

   local plasmaE, plasmaB = field:hasEB()
   local extHasE, extHasB = externalField:hasEB()

   local hasE = plasmaE or extHasE
   local hasB = plasmaB or extHasB

   -- External forces are accumulated to the electric field part of
   -- totalEmField
   if self.hasExtForce then
      hasE, hasB = true, true
   end

   if hasB then
      self.totalEmField = self:allocVectorMoment(8)     -- 8 components of EM field.
   else
      --self.totalEmField = self:allocVectorMoment(3)     -- Electric field only.
      self.totalEmField = self:allocMoment()  -- Phi only (Vlasov-Poisson)
   end
   if self.hasExtForce then
      self.totEmFptr, self.totEmFidxr = self.totalEmField:get(1), self.totalEmField:genIndexer()
   end

   self.computePlasmaB = true and plasmaB or extHasB

   if self.model_id == "GKYL_MODEL_SR" then 
      -- Initialize velocity-space arrays for relativistic Vlasov
      -- Only need to do this once, so we don't need to store the updater
      local initSRVarsCalc = Updater.CalcSRVars {
         velGrid = self.velGrid, 
         confBasis = self.confBasis, 
         velBasis = self.velBasis, 
         phaseBasis = self.basis, 
         operation = "init", 
      }
      initSRVarsCalc:advance(0.0, {}, {self.p_over_gamma, self.gamma, self.gamma_inv})
      -- Create table of pointers to fields needed in update
      self.fldPtrs = {self.totalEmField, self.p_over_gamma}
      self.momPtrs = {self.p_over_gamma, self.gamma, self.gamma_inv}
   else
      -- Create table of pointers to fields needed in update
      self.fldPtrs = {self.totalEmField}
      -- No auxiliary fields for moments for non-relativistic Vlasov
      self.momPtrs = {}
   end

   -- Create updater to advance solution by one time-step.
   if self.evolveCollisionless then
      self.solver = Updater.VlasovDG {
         onGrid    = self.grid,                      confBasis = self.confBasis,                phaseBasis = self.basis, 
         confRange = self.totalEmField:localRange(), velRange = self.p_over_gamma:localRange(), phaseRange = self.distf[1]:localRange(), 
         model_id  = self.model_id,                  field_id = self.field_id,                  fldPtrs    = self.fldPtrs, 
      }
      self.collisionlessAdvance = function(tCurr, inFlds, outFlds)
         self.solver:advance(tCurr, inFlds, outFlds)
      end

      -- Boundary flux updater.
      self.boundaryFluxSlvr = Updater.BoundaryFluxCalc {
         onGrid = self.grid,  equation    = self.solver:getEquation(), 
         cdim   = self.cdim,  equation_id = "vlasov",
      }
      self.collisionlessBoundaryAdvance = function(tCurr, inFlds, outFlds)
         local tmStart = Time.clock()
         self.boundaryFluxSlvr:advance(tCurr, inFlds, outFlds)
         self.timers.boundflux = self.timers.boundflux + Time.clock() - tmStart
      end
   else
      self.solver = {totalTime = 0.}
      self.collisionlessAdvance = function(tCurr, inFlds, outFlds) end
      self.collisionlessBoundaryAdvance = function(tCurr, inFlds, outFlds) end
   end

   -- Create updaters to compute various moments.
   -- M0 and M1i are common to both non-relativistic and relativistic Vlasov, so use common DistFuncMomentDG 
   self.numDensityCalc = Updater.DistFuncMomentDG {
      onGrid    = self.grid,                      confBasis = self.confBasis,                phaseBasis = self.basis, 
      confRange = self.totalEmField:localRange(), velRange = self.p_over_gamma:localRange(), 
      model_id  = self.model_id,                  isIntegrated = false,                      momPtrs    = self.momPtrs, 
      moment    = "M0", 
   }
   self.momDensityCalc = Updater.DistFuncMomentDG {
      onGrid    = self.grid,                      confBasis = self.confBasis,                phaseBasis = self.basis, 
      confRange = self.totalEmField:localRange(), velRange = self.p_over_gamma:localRange(), 
      model_id  = self.model_id,                  isIntegrated = false,                      momPtrs    = self.momPtrs,  
      moment    = "M1i", 
   }
   if self.model_id == "GKYL_MODEL_SR" then
      self.ptclEnergyCalc = Updater.DistFuncMomentDG {
         onGrid    = self.grid,                      confBasis = self.confBasis,                phaseBasis = self.basis, 
         confRange = self.totalEmField:localRange(), velRange = self.p_over_gamma:localRange(), 
         model_id  = self.model_id,                  isIntegrated = false,                      momPtrs    = self.momPtrs,   
         moment    = "Energy", 
      } 
   else
      self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis  = self.confBasis,
         phaseBasis = self.basis,  moment     = "M2",
      }
      self.M2ijCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis = self.confBasis,
         phaseBasis = self.basis,  moment    = "M2ij",
      }
      self.M3iCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis = self.confBasis,
         phaseBasis = self.basis,  moment    = "M3i",
      }
      self.calcMaxwell = Updater.MaxwellianOnBasis {
         onGrid     = self.grid,   confGrid  = self.confGrid,
         phaseBasis = self.basis,  confBasis = self.confBasis,
      }
   end
   if self.needFiveMoments then
      -- Create updater to compute M0, M1i, M2 moments.
      self.fiveMomentsCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis = self.confBasis,
         phaseBasis = self.basis,  moment    = "FiveMoments",
      }
      self.calcSelfCouplingMom = self.computePlasmaB
         and function(tCurr, fIn)
            -- Compute M0, M1i and M2.
            self.fiveMomentsCalc:advance(tCurr, {fIn}, {self.fiveMoments})
            -- Copy momentum density to its own field.
            self.momDensity:combineOffset(1., self.fiveMoments, 1*self.confBasis:numBasis())
         end
         or function(tCurr, fIn)
            -- Compute M0, M1i and M2.
            self.fiveMomentsCalc:advance(tCurr, {fIn}, {self.fiveMoments})
            -- Copy number density to its own field.
            self.numDensity:combineOffset(1., self.fiveMoments, 0)
         end
   else
      self.calcSelfCouplingMom = self.computePlasmaB
         and function(tCurr, fIn)
            self.momDensityCalc:advance(tCurr, {fIn}, { self.momDensity })
         end
         or function(tCurr, fIn)
            self.numDensityCalc:advance(tCurr, {fIn}, { self.numDensity })
         end
   end

   if self.hasExtForce then
      self.evalVlasovExtForce = Updater.ProjectOnBasis {
         onGrid = self.confGrid,   evaluate = self.vlasovExtForceFunc,
         basis  = self.confBasis,  onGhosts = false
      }

      self.accumulateExtForce = function(tCurr, totalEmField)
         local vExtForce  = self.vExtForce
         local vItr, eItr = self.vExtFptr, self.totEmFptr
         self.evalVlasovExtForce:advance(tCurr, {}, {vExtForce})

         -- Analogous to the current, the external force only gets accumulated onto the electric field.
         for idx in totalEmField:localRangeIter() do
            vExtForce:fill(self.vExtFidxr(idx), vItr)
            totalEmField:fill(self.totEmFidxr(idx), eItr)
            for i = 1, vExtForce:numComponents() do eItr[i] = eItr[i]+vItr[i] end
         end
      end
   else
      self.accumulateExtForce = function(tCurr, totalEmField) end
   end

   -- Create an updater for volume integrals. Used by diagnostics.
   self.volIntegral = {
      scalar = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.confGrid,   numComponents = 1,
         basis  = self.confBasis,  quantity      = "V",
      },
      vector = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.confGrid,   numComponents = self.vdim,
         basis  = self.confBasis,  quantity      = "V",
      },
   }

   -- Create species source solvers.
   for _, src in lume.orderedIter(self.sources) do src:createSolver(self, externalField) end
end

function VlasovSpecies:initDist(extField)
   -- Note: do not call applyBc here. It is called later in initialization sequence.

   if self.randomseed then
      math.randomseed(self.randomseed)
   else
      math.randomseed(47*Mpi.Comm_rank(Mpi.COMM_WORLD)+os.time())
   end

   local syncPeriodicDirs = true

   local initCnt = 0
   for nm, pr in lume.orderedIter(self.projections) do
      pr:fullInit(self)
      pr:advance(0.0, {extField}, {self.distf[2]})
      if string.find(nm,"init") then
         self.distf[1]:accumulate(1.0, self.distf[2])
         initCnt = initCnt + 1
      end
   end

   assert(initCnt>0, string.format("VlasovSpecies: Species '%s' not initialized!", self.name))

   self.distf[2]:clear(0.0)
end

function VlasovSpecies:createCouplingSolver(population, field, externalField)
   -- After all species have called their createSolver methods, we here create the objects
   -- needed for cross-species solves (e.g. cross-species collisions).

   local species = population:getSpecies()

   -- Create cross collision solvers.
   for _, c in lume.orderedIter(self.collisions) do c:createCouplingSolver(population, field, externalField) end

   -- Create source solvers.
   for _, src in lume.orderedIter(self.sources) do src:createCouplingSolver(population, field, externalField) end

   -- Create BC solvers.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:createCouplingSolver(species, field, externalField) end
end

function VlasovSpecies:initCrossSpeciesCoupling(population)
   -- This method establishes the interaction between different
   -- species that is not mediated by the field (solver), like
   -- collisions.

   local species = population:getSpecies()

   -- Determine if M0, M1i and M2 are needed.
   self.needFiveMoments = false

   -- Create a double nested table of colliding species.
   -- In this table we will encode information about that collition such as:
   --   * does the collision take place?
   --   * Operator modeling the collision.
   -- Other features of a collision may be added in the future.
   self.collPairs = {}
   for sN, _ in lume.orderedIter(species) do
      self.collPairs[sN] = {}
      for sO, _ in lume.orderedIter(species) do
         self.collPairs[sN][sO] = {}
         -- Need next below because species[].collisions is created as an empty table. 
         if species[sN].collisions and next(species[sN].collisions) then 
            for collNm, _ in pairs(species[sN].collisions) do
               -- This species collides with someone.
               self.collPairs[sN][sO].on = lume.any(species[sN].collisions[collNm].collidingSpecies,
                                                    function(e) return e==sO end)
               if self.collPairs[sN][sO].on then
                  self.collPairs[sN][sO].kind = species[sN].collisions[collNm].collKind
                  self.needFiveMoments = true  -- MF 2022/09/16: currently all collision models need M0, M1, M2.
               end
            end
         else
            -- This species does not collide with anyone.
            self.collPairs[sN][sO].on = false
         end
      end
   end

   local isThisSpeciesMine = population:isSpeciesMine(self.name)
   -- Allocate fiveMoments if we collide with other species not in this rank.
   if self.needFiveMoments and (not isThisSpeciesMine) then
      self.fiveMoments = self:allocVectorMoment(self.vdim+2)
   end

   -- Create list of ranks we need to send/recv local fiveMoments to/from.
   self.fiveMomentsXfer = {}
   self.fiveMomentsXfer.destRank, self.fiveMomentsXfer.srcRank  = {}, {}
   self.fiveMomentsXfer.sendReqStat, self.fiveMomentsXfer.recvReqStat = nil, nil
   for sO, info in pairs(self.collPairs[self.name]) do
      local sOrank = population:getSpeciesOwner(sO)
      local selfRank = population:getSpeciesOwner(self.name)
      if sO~=self.name and info.on then
         if isThisSpeciesMine then
            -- Only species owned by this rank send fiveMoments to other ranks.
            if #self.fiveMomentsXfer.destRank == 0 and (not population:isSpeciesMine(sO)) then
               table.insert(self.fiveMomentsXfer.destRank, sOrank)
               self.fiveMomentsXfer.sendReqStat = Mpi.RequestStatus()
            end
         else
            -- Only species not owned by this rank receive fiveMoments from other ranks.
            if #self.fiveMomentsXfer.srcRank == 0 and (not population:isSpeciesMine(self.name)) then
               table.insert(self.fiveMomentsXfer.srcRank, selfRank)
               self.fiveMomentsXfer.recvReqStat = Mpi.RequestStatus()
            end
         end
      end
   end

   -- Initialize the BC cross-coupling interactions.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:initCrossSpeciesCoupling(species) end

   -- Initialize the source cross-coupling interactions.
   for _, src in lume.orderedIter(self.sources) do src:initCrossSpeciesCoupling(population) end
end

function VlasovSpecies:advance(tCurr, population, emIn, inIdx, outIdx)
   local tmStart = Time.clock()

   self:setActiveRKidx(inIdx)
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]
   fRhsOut:clear(0.0)

   -- Accumulate functional Maxwell fields (if needed).
   local emField         = emIn[1]:rkStepperFields()[inIdx]
   local emExternalField = emIn[2]:rkStepperFields()[1]
   local totalEmField    = self.totalEmField
   totalEmField:clear(0.0)

   if emField then totalEmField:accumulate(self.qbym, emField) end
   if emExternalField then totalEmField:accumulate(self.qbym, emExternalField) end

   -- If external force present (gravity, body force, etc.) accumulate it to electric field.
   self.accumulateExtForce(tCurr, totalEmField)

   self.collisionlessAdvance(tCurr, {fIn}, {fRhsOut, self.cflRateByCell})
   self.timers.collisionless = self.solver.totalTime

   -- Perform the collision update.
   self.timers.collisions = 0.
   for _, c in lume.orderedIter(self.collisions) do
      c:advance(tCurr, fIn, population, {fRhsOut, self.cflRateByCell})   -- 'population' needed for cross-species collisions.
      self.timers.collisions = self.timers.collisions + c:getTimer('advance')
   end

   self.collisionlessBoundaryAdvance(tCurr, {fIn}, {fRhsOut})

   self.timers.sources = 0.
   for _, src in lume.orderedIter(self.sources) do
      src:advance(tCurr, fIn, population:getSpecies(), fRhsOut)
      self.timers.sources = self.timers.sources + src:getTimer('advance')
   end
   
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do
      bc:storeBoundaryFlux(tCurr, outIdx, fRhsOut)   -- Save boundary fluxes.
   end

   self.timers.advance = self.timers.advance + Time.clock() - tmStart
end

function VlasovSpecies:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx)
   -- Perform some operations after the updates have been computed, but before
   -- the combine RK (in PlasmaOnCartGrid) is called.
   local tmStart = Time.clock()

   local species = population:getSpecies()

   -- Wait to finish sending fiveMoments if needed.
   population:speciesXferField_waitSend(self.fiveMomentsXfer)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx)
   end

   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:advanceCrossSpeciesCoupling(tCurr, species, outIdx) end

   for _, src in lume.orderedIter(self.sources) do src:advanceCrossSpeciesCoupling(tCurr, species, outIdx) end

   self.timers.advancecross = self.timers.advancecross + Time.clock() - tmStart
end

function VlasovSpecies:calcCouplingMoments(tCurr, rkIdx, species)
   -- Compute moments needed in coupling to fields and collisions.
   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]

   -- Compute M0, M1i and/or M2 depending on what fields and collisions need.
   self.calcSelfCouplingMom(tCurr, fIn)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:calcCouplingMoments(tCurr, rkIdx, species)
   end

   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:calcCouplingMoments(tCurr, rkIdx, species) end

   self.timers.mom = self.timers.mom + Time.clock() - tmStart
end

function VlasovSpecies:calcCrossCouplingMoments(tCurr, rkIdx, population)
   -- Perform cross-species calculation related to coupling moments that require the
   -- self-species coupling moments.
   local tmStart = Time.clock()

   -- Begin sending/receiving fiveMoments if needed.
   population:speciesXferField_begin(self.fiveMomentsXfer, self.fiveMoments, 22)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:calcCrossCouplingMoments(tCurr, rkIdx, population)
   end

   self.timers.momcross = self.timers.momcross + Time.clock() - tmStart
end

function VlasovSpecies:copyRk(outIdx, aIdx)
   self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx])

   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:copyBoundaryFluxField(aIdx, outIdx) end

   if self.positivity then
      self.fDelPos[outIdx]:copy(self.fDelPos[aIdx])
   end
end

function VlasovSpecies:combineRk(outIdx, a, aIdx, ...)
   local args = {...} -- Package up rest of args as table.
   local nFlds = #args/2
   self:rkStepperFields()[outIdx]:combine(a, self:rkStepperFields()[aIdx])
   for i = 1, nFlds do -- Accumulate rest of the fields.
      self:rkStepperFields()[outIdx]:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]])
   end

   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do
      bc:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   end

   if self.positivity then
      self.fDelPos[outIdx]:combine(a, self.fDelPos[aIdx])
      for i = 1, nFlds do -- Accumulate rest of the fields.
         self.fDelPos[outIdx]:accumulate(args[2*i-1], self.fDelPos[args[2*i]])
      end
   end
end

function VlasovSpecies:applyBcIdx(tCurr, field, externalField, inIdx, outIdx, isFirstRk)
   self:applyBc(tCurr, field, externalField, inIdx, outIdx)
end

function VlasovSpecies:applyBcDontEvolve(tCurr, field, externalField, inIdx, outIdx) end
function VlasovSpecies:applyBcEvolve(tCurr, field, externalField, inIdx, outIdx)
   -- fIn is total distribution function.
   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[outIdx]

   local syncPeriodicDirsTrue = true

   -- Apply non-periodic BCs (to only fluctuations if fluctuation BCs).
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:advance(tCurr, self, field, externalField, inIdx, outIdx) end

   -- Apply periodic BCs (to only fluctuations if fluctuation BCs)
   fIn:sync(syncPeriodicDirsTrue)

   self.timers.bc = self.timers.bc + Time.clock() - tmStart
end
function VlasovSpecies:applyBc(tCurr, field, externalField, inIdx, outIdx)
   self.applyBcFunc(tCurr, field, externalField, inIdx, outIdx)
end

function VlasovSpecies:createDiagnostics(field)
   -- Many diagnostics require dividing or multiplying by the Jacobian
   -- (if present). Predefine the functions that do that.
   self.divideByJacobGeo = self.jacobGeoInv
      and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeoInv}, {fldOut}) end
      or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
   self.multiplyByJacobGeo = self.jacobGeo
      and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeo}, {fldOut}) end
      or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end

   -- Create this species' diagnostics.
   if self.tbl.diagnostics then
      self.diagnostics[self.name] = DiagsApp{implementation = VlasovDiags()}
      self.diagnostics[self.name]:fullInit(self, field, self)
   end

   for srcNm, src in lume.orderedIter(self.sources) do
      self.diagnostics[self.name..srcNm] = src:createDiagnostics(self, field)
   end

   for bcNm, bc in lume.orderedIter(self.nonPeriodicBCs) do
      self.diagnostics[self.name..bcNm] = bc:createDiagnostics(self, field)
   end

   for collNm, coll in lume.orderedIter(self.collisions) do
      self.diagnostics[self.name..collNm] = coll:createDiagnostics(self, field)
   end
   lume.setOrder(self.diagnostics)
end

function VlasovSpecies:clearCFL()
   -- Clear cflRateByCell for next cfl calculation.
   self.cflRateByCell:clear(0.0)
end

function VlasovSpecies:suggestDt()
   if not self.evolve then return GKYL_MAX_DOUBLE end

   local dtSuggested = math.min(self.cfl/self.cflRateByCell:reduce('max')[1], GKYL_MAX_DOUBLE)

   -- If dtSuggested == GKYL_MAX_DOUBLE, it is likely because of NaNs.
   -- If so, return 0 so that no timestep is taken, and we will abort the simulation.
   if dtSuggested == GKYL_MAX_DOUBLE then dtSuggested = 0.0 end

   return dtSuggested
end

function VlasovSpecies:setDtGlobal(dtGlobal) self.dtGlobal[0] = dtGlobal end

function VlasovSpecies:rkStepperFields() return self.distf end

function VlasovSpecies:getDistF(rkIdx)
   if rkIdx == nil then
      return self:rkStepperFields()[self.activeRKidx]
   else
      return self:rkStepperFields()[rkIdx]
   end
end

function VlasovSpecies:fluidMoments() return self.fiveMoments end

function VlasovSpecies:getNumDensity(rkIdx)
   -- If no rkIdx specified, assume numDensity has already been calculated.
   if rkIdx == nil then return self.numDensity end 

   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]
   self.numDensityCalc:advance(nil, {fIn}, { self.numDensity })

   self.timers.mom = self.timers.mom + Time.clock() - tmStart
   return self.numDensity
end

function VlasovSpecies:getMomDensity(rkIdx)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then return self.momDensity end 

   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]
   self.momDensityCalc:advance(nil, {fIn}, { self.momDensity })

   self.timers.mom = self.timers.mom + Time.clock() - tmStart
   return self.momDensity
end

-- Please test this for higher than 1x1v... (MF: JJ?).
function VlasovSpecies:Maxwellian(xn, n0, vdnIn, T0)
   local vdn = vdnIn or {0, 0, 0}
   local vt2 = T0/self.mass
   local v2 = 0.0
   for d = self.cdim+1, self.cdim+self.vdim do
      v2 = v2 + (xn[d] - vdnIn[d-self.cdim])^2
   end
   return n0 / math.sqrt(2*math.pi*vt2)^self.vdim * math.exp(-v2/(2*vt2))
end

function VlasovSpecies:projToSource(proj)
   -- For backwards compatibility: in the case where the source is specified
   -- as a projection object in the input file, this function turns that
   -- projection object into a Source object.
   local tbl = proj.tbl
   local pow = tbl.power
   return Source { profile = proj, power = pow }
end

function VlasovSpecies:write(tm, field, force)
   if self.evolve then

      for _, dOb in lume.orderedIter(self.diagnostics) do
         dOb:resetState(tm)   -- Reset booleans indicating if diagnostic has been computed.
      end

      for _, bc in lume.orderedIter(self.nonPeriodicBCs) do
         bc:computeBoundaryFluxRate(self.dtGlobal[0])
      end

      -- Compute integrated diagnostics.
      if self.calcIntQuantTrigger(tm) then
         for _, dOb in lume.orderedIter(self.diagnostics) do
            dOb:calcIntegratedDiagnostics(tm, self, field)
         end
      end

      -- Only write stuff if triggered.
      if self.distIoTrigger(tm) or force then
         local fIn = self:rkStepperFields()[1]

         self.distIo:write(fIn, string.format("%s_%d.bp", self.name, self.distIoFrame), tm, self.distIoFrame)

         for _, src in lume.orderedIter(self.sources) do src:write(tm, self.distIoFrame) end

         self.distIoFrame = self.distIoFrame+1
      end

      if self.diagIoTrigger(tm) or force then
         -- Compute moments and write them out.

         -- Compute diagnostics defined on a grid.
         for _, dOb in lume.orderedIter(self.diagnostics) do
            dOb:calcGridDiagnostics(tm, self, field)
         end

         for _, dOb in lume.orderedIter(self.diagnostics) do   -- Write grid and integrated diagnostics.
            dOb:write(tm, self.diagIoFrame)
         end

         for _, c in lume.orderedIter(self.collisions) do
            c:write(tm, self.diagIoFrame)  -- MF: Preferably this method will go away. Use diagnostics instead.
         end

         self.diagIoFrame = self.diagIoFrame+1
      end
   else
      -- If not evolving species, don't write anything except initial conditions.
      if self.distIoFrame == 0 then

         -- Compute integrated diagnostics.
         for _, dOb in lume.orderedIter(self.diagnostics) do
            dOb:calcIntegratedDiagnostics(tm, self, field)
         end

         self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, 0), tm, 0)

         -- Compute diagnostics defined on a grid.
         for _, dOb in lume.orderedIter(self.diagnostics) do
            dOb:calcGridDiagnostics(tm, self, field)
         end

         for _, dOb in lume.orderedIter(self.diagnostics) do   -- Write grid and integrated diagnostics.
            dOb:write(tm, self.diagIoFrame)
         end
      end
      self.distIoFrame = self.distIoFrame+1
   end

end

function VlasovSpecies:writeRestart(tm)
   -- (The final "true/false" in calls to :write determines writing of ghost cells).
   local writeGhost = false

   self.distIo:write(self.distf[1], string.format("%s_restart.bp", self.name), tm, self.distIoFrame, writeGhost)

   for _, dOb in lume.orderedIter(self.diagnostics) do   -- Write restart diagnostics.
      dOb:writeRestart(tm, self.diagIoFrame, self.dynVecRestartFrame)
   end

   self.dynVecRestartFrame = self.dynVecRestartFrame + 1
end

function VlasovSpecies:readRestart(field, externalField)
   local readGhost = false

   local tm, fr = self.distIo:read(self.distf[1], string.format("%s_restart.bp", self.name), readGhost)
   self.distIoFrame = fr -- Reset internal frame counter.

   -- Set ghost cells.
   self.distf[1]:sync()

   -- Apply BCs (unless skin cells have been read because of special BCs).
   self:applyBc(tm, field, externalField, 1, 1)

   local diagIoFrame_new
   for _, dOb in lume.orderedIter(self.diagnostics) do   -- Read grid and integrated diagnostics.
      local _, dfr = dOb:readRestart()
      diagIoFrame_new = diagIoFrame_new or dfr
      if dfr then
         assert(diagIoFrame_new==dfr, string.format("VlasovSpecies:readRestart expected diagnostics from previous run to have the same last frame. Instead got %d and %d", diagIoFrame_new, dfr))
      end
   end
   -- The 'or self.distIoFrame' below is for sims without diagnostics, or when the first
   -- run didn't request diagnostics, but the latter (when the restart requests diagnostics
   -- while the first one didn't) requires commenting out the loop above (a hack, beware).
   self.diagIoFrame = diagIoFrame_new or self.distIoFrame

   -- Iterate triggers.
   self.distIoTrigger(tm)
   self.diagIoTrigger(tm)

   return tm
end

function VlasovSpecies:clearTimers()
   for nm, _ in pairs(self.timers) do self.timers[nm] = 0. end
   self.solver.totalTime = 0.
   for _, c in lume.orderedIter(self.collisions) do c:clearTimers() end
   for _, src in lume.orderedIter(self.sources) do src:clearTimers() end
end

-- ................... Classes meant as aliases to simplify input files ...................... --
-- Default: Non-Relativistic Vlasov-Maxwell (Cartesian geometry)
local VlasovMaxwellSpecies = Proto(VlasovSpecies)
function VlasovMaxwellSpecies:fullInit(mySpecies)
   self.model_id = "GKYL_MODEL_DEFAULT"
   self.field_id = "GKYL_FIELD_E_B"
   VlasovMaxwellSpecies.super.fullInit(self, mySpecies)
end

-- Vlasov-Poisson, only phi (Cartesian geometry)
local VlasovPoissonSpecies = Proto(VlasovSpecies)
function VlasovPoissonSpecies:fullInit(mySpecies)
   self.model_id = "GKYL_MODEL_DEFAULT"
   self.field_id = "GKYL_FIELD_PHI"
   VlasovPoissonSpecies.super.fullInit(self, mySpecies)
end

-- Vlasov-Poisson, phi + constant A (constant background magnetic field, Cartesian geometry)
local VlasovPoissonASpecies = Proto(VlasovSpecies)
function VlasovPoissonASpecies:fullInit(mySpecies)
   self.model_id = "GKYL_MODEL_DEFAULT"
   self.field_id = "GKYL_FIELD_PHI_A"
   VlasovPoissonASpecies.super.fullInit(self, mySpecies)
end

-- Neutral Non-Relativistic Vlasov (Cartesian geometry)
local VlasovNeutralSpecies = Proto(VlasovSpecies)
function VlasovNeutralSpecies:fullInit(mySpecies)
   self.model_id = "GKYL_MODEL_DEFAULT"
   self.field_id = "GKYL_FIELD_NULL"
   VlasovNeutralSpecies.super.fullInit(self, mySpecies)
end

-- Neutral Vlasov (General geometry)
local VlasovGenGeoNeutralSpecies = Proto(VlasovSpecies)
function VlasovGenGeoNeutralSpecies:fullInit(mySpecies)
   self.model_id = "GKYL_MODEL_GEN_GEO"
   self.field_id = "GKYL_FIELD_NULL"
   VlasovGenGeoNeutralSpecies.super.fullInit(self, mySpecies)
end

-- Special Relativistic Vlasov-Maxwell (Cartesian geometry)
local VlasovSRMaxwellSpecies = Proto(VlasovSpecies)
function VlasovSRMaxwellSpecies:fullInit(mySpecies)
   self.model_id = "GKYL_MODEL_SR"
   self.field_id = "GKYL_FIELD_E_B"
   VlasovSRMaxwellSpecies.super.fullInit(self, mySpecies)
end

-- Neutral Special Relativistic Vlasov (Cartesian geometry)
local VlasovSRNeutralSpecies = Proto(VlasovSpecies)
function VlasovSRNeutralSpecies:fullInit(mySpecies)
   self.model_id = "GKYL_MODEL_SR"
   self.field_id = "GKYL_FIELD_NULL"
   VlasovSRNeutralSpecies.super.fullInit(self, mySpecies)
end
-- ................... End of VlasovSpecies alias classes .................... --

return {VlasovMaxwell       = VlasovMaxwellSpecies,
        VlasovPoisson       = VlasovPoissonSpecies,
        VlasovPoissonA      = VlasovPoissonASpecies,
        VlasovNeutral       = VlasovNeutralSpecies,
        VlasovGenGeoNeutral = VlasovGenGeoNeutralSpecies, 
        VlasovSRMaxwell     = VlasovSRMaxwellSpecies,
        VlasovSRNeutral     = VlasovSRNeutralSpecies}
