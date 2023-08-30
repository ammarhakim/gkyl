-- Gkyl ------------------------------------------------------------------------
--
-- Gyrokinetic species object.
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
local Source           = require "App.Sources.GkSource"
local DiagsApp         = require "App.Diagnostics.SpeciesDiagnostics"
local GkDiags          = require "App.Diagnostics.GkDiagnostics"
local BasicBC          = require("App.BCs.GkBasic").GkBasic
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
local Constants        = require "Lib.Constants"

local GkSpecies = Proto(SpeciesBase)

-- ............. Backwards compatible treatment of BCs .....................--
-- Add constants to object indicate various supported boundary conditions.
local SP_BC_ABSORB   = 1
local SP_BC_OPEN     = 2
local SP_BC_REFLECT  = 3
local SP_BC_SHEATH   = 4
local SP_BC_ZEROFLUX = 5
local SP_BC_COPY     = 6
GkSpecies.bcAbsorb   = SP_BC_ABSORB      -- Absorb all particles.
GkSpecies.bcOpen     = SP_BC_OPEN        -- Zero gradient.
GkSpecies.bcReflect  = SP_BC_REFLECT     -- Specular reflection.
GkSpecies.bcSheath   = SP_BC_SHEATH      -- Sheath.
GkSpecies.bcZeroFlux = SP_BC_ZEROFLUX    -- Zero flux.
GkSpecies.bcCopy     = SP_BC_COPY        -- Copy stuff.

function GkSpecies:makeBcApp(bcIn, dir, edge)
   local bcOut
   if type(bcIn) == "function" then
      bcOut = BasicBC{kind="function", bcFunction=bcIn, diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_COPY then
      print("GkSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="copy", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_ABSORB then
      print("GkSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="absorb", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_OPEN then
      print("GkSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      -- AHH: open seems unstable. So using plain copy.
      bcOut = BasicBC{kind="copy", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_REFLECT then
      print("GkSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="reflect", diagnostics={}, saveFlux=false}
   elseif bcIn == SP_BC_SHEATH then
      print("GkSpecies: warning... old way of specifyin BCs will be deprecated. Use BC apps instead.")
      bcOut = BasicBC{kind="sheath", diagnostics={}, saveFlux=false}
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
function GkSpecies:init(tbl) self.tbl = tbl end

-- Function to create basis functions.
local function createBasis(nm, cdim, vdim, polyOrder)
   local ndim = cdim+vdim
   if nm == "serendipity" then
      if polyOrder == 1 then
         return Basis.CartModalGkHybrid { cdim = cdim, vdim = vdim }
      else
         return Basis.CartModalSerendipity { ndim = ndim, polyOrder = polyOrder }
      end
   elseif nm == "tensor" then
      return Basis.CartModalTensor { ndim = ndim, polyOrder = polyOrder }
   end
end

function GkSpecies:createBasis(nm, polyOrder)
   self.basis = createBasis(nm, self.cdim, self.vdim, polyOrder)
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

function GkSpecies:alloc(nRkDup)
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

   self.flucF = (self.fluctuationBCs or self.perturbedDiagnostics) and self:allocDistf() or nil   -- Fluctuation.

   -- Array with one component per cell to store cflRate in each cell.
   self.cflRateByCell = self:allocCartField(self.grid, 1, {1,1})
   self.dtGlobal      = ffi.new("double[2]")
   self.dtGlobal[0], self.dtGlobal[1] = 1.0, 1.0   -- Temporary value (so diagnostics at t=0 aren't inf).

   -- Allocate fields to store coupling moments (for use in coupling
   -- to field and collisions).
   self.numDensity    = self:allocMoment()
   self.numDensityAux = self:allocMoment()
   self.momDensity    = self:allocMoment()
   self.momDensityAux = self:allocMoment()
   self.ptclEnergy    = self:allocMoment()
   self.ptclEnergyAux = self:allocMoment()
   self.threeMoments  = self:allocVectorMoment(3)
			
   self.vDegFreedom = self.vdim == 1 and 1.0 or 3.0

   self.first = true
end

function GkSpecies:fullInit(appTbl)
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

   self.perturbedDiagnostics = false
   if tbl.diagnostics then
      if lume.any(tbl.diagnostics, function(e) return e=="perturbed" end) then
         lume.remove(tbl.diagnostics,"perturbed")
         self.perturbedDiagnostics = true
      end
   end

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

   -- Read in boundary conditions.
   self.bcInDir = {{ }, { }, { }}   -- List of BCs to apply.
   if tbl.bcx then
      if tbl.bcx[1] == nil or tbl.bcx[2] == nil then assert(false, "GkSpecies: unsupported BC type") end
      self.bcInDir[1] = {tbl.bcx[1], tbl.bcx[2]}
   end
   if tbl.bcy then
      if tbl.bcy[1] == nil or tbl.bcy[2] == nil then assert(false, "GkSpecies: unsupported BC type") end
      self.bcInDir[2] = {tbl.bcy[1], tbl.bcy[2]}
   end
   if tbl.bcz then
      if tbl.bcz[1] == nil or tbl.bcz[2] == nil then assert(false, "GkSpecies: unsupported BC type") end
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
         assert(goodBC, "GkSpecies: bc not recognized.")
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

function GkSpecies:setCfl(cfl)
   self.cfl = cfl
   for _, c in lume.orderedIter(self.collisions) do c:setCfl(cfl) end
end

function GkSpecies:setIoMethod(ioMethod) self.ioMethod = ioMethod end

function GkSpecies:setConfBasis(basis)
   self.confBasis = basis
   for _, c in lume.orderedIter(self.collisions) do c:setConfBasis(basis) end
   for _, src in pairs(self.sources) do src:setConfBasis(basis) end
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:setConfBasis(basis) end
end
function GkSpecies:setConfGrid(grid)
   self.confGrid = grid
   for _, c in lume.orderedIter(self.collisions) do c:setConfGrid(grid) end
   for _, src in pairs(self.sources) do src:setConfGrid(grid) end
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:setConfGrid(grid) end
end

function GkSpecies:createGrid(confGridIn)
   local confGrid = assert(confGridIn or self.confGrid, "GkSpecies:createGrid ... must pass in confGrid or call setConfGrid prior to createGrid")

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
end

-- Field allocation in the species objects should be performed with one
-- of the following four functions instead of calling DataStruct directly.
function GkSpecies:allocCartField(grid,nComp,ghosts,metaData)
   local f = DataStruct.Field {
      onGrid        = grid,   ghost    = ghosts,
      numComponents = nComp,  metaData = metaData,
   }
   f:clear(0.0)
   return f
end
function GkSpecies:allocDistf()
   local metaData = {polyOrder = self.basis:polyOrder(),
                     basisType = self.basis:id(),
                     charge    = self.charge,
                     mass      = self.mass,
                     grid      = GKYL_OUT_PREFIX .. "_" .. self.name .. "_grid.bp"}
   return self:allocCartField(self.grid,self.basis:numBasis(),{self.nGhost,self.nGhost},metaData)
end
function GkSpecies:allocMoment()
   local metaData = {polyOrder = self.confBasis:polyOrder(),
                     basisType = self.confBasis:id(),
                     charge    = self.charge,
                     mass      = self.mass,}
   return self:allocCartField(self.confGrid,self.confBasis:numBasis(),{self.nGhost,self.nGhost},metaData)
end
function GkSpecies:allocVectorMoment(dim)
   local metaData = {polyOrder = self.confBasis:polyOrder(),
                     basisType = self.confBasis:id(),
                     charge    = self.charge,
                     mass      = self.mass,}
   return self:allocCartField(self.confGrid,dim*self.confBasis:numBasis(),{self.nGhost,self.nGhost},metaData)
end
function GkSpecies:allocIntMoment(comp)
   local metaData = {charge = self.charge,
                     mass   = self.mass,}
   local ncomp = comp or 1
   local f = DataStruct.DynVector {numComponents = ncomp,     writeRank = 0,
                                   metaData      = metaData,  comm      = self.confGrid:commSet().comm,}
   return f
end

function GkSpecies:createSolver(field, externalField)
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

   -- Functions to compute fluctuations given the current moments and background,
   -- and the full-F moments given the fluctuations and background.
   self.getF_or_deltaF = self.deltaF and function(fIn)
      self.flucF:combine(1.0, fIn, -1.0, self.fBackground)
      return self.flucF
   end or function(fIn) return fIn end
   self.minusBackgroundF = self.fluctuationBCs
                          and function(fIn) fIn:accumulate(-1.0, self.fBackground) end
                          or function(fIn) end
   self.calcFullF = self.fluctuationBCs and function(fIn, syncFullFperiodicDirs)
      fIn:accumulate(1.0, self.fBackground)
      fIn:sync(syncFullFperiodicDirs)
   end or function(fIn, syncFullFperiodicDirs) end
   self.calcDeltaF = self.perturbedDiagnostics
                    and function(fIn) self.flucF:combine(1.0, fIn, -1.0, self.fBackground) end
                    or function(fIn) end

   if self.fluctuationBCs or self.perturbedDiagnostics then
      self.writeFluctuation = self.perturbedDiagnostics
         and function(tm, fr, fIn)
            self.distIo:write(self.flucF, string.format("%s_fluctuation_%d.bp", self.name, self.diagIoFrame), tm, fr)
         end
         or function(tm, fr, fIn)
            self.calcDeltaF(fIn)
            self.distIo:write(self.flucF, string.format("%s_fluctuation_%d.bp", self.name, self.diagIoFrame), tm, fr)
         end
   else
      self.writeFluctuation = function(tm, fr, fIn) end
   end

   if self.evolve then
      self.applyBcFunc = function(tCurr, field, externalField, inIdx, outIdx)
         return GkSpecies["applyBcEvolve"](self, tCurr, field, externalField, inIdx, outIdx)
      end
   else
      self.applyBcFunc = function(tCurr, field, externalField, inIdx, outIdx)
         return GkSpecies["applyBcDontEvolve"](self, tCurr, field, externalField, inIdx, outIdx)
      end
   end

   local hasE, hasB       = field:hasEB()
   local extHasE, extHasB = externalField:hasEB()

   local hasPhi  = hasE or extHasE
   local hasApar = hasB or extHasB

   -- Set up Jacobian.
   if externalField then
      self.bmagFunc = externalField.bmagFunc
      -- If vdim>1, get the phase-space Jacobian (=bmag) from geo.
      self.jacobPhaseFunc = self.bmagFunc
      self.jacobGeoFunc   = externalField.jacobGeoFunc

      self.bmag        = assert(externalField.geo.bmag, "nil bmag")
      self.bmagInv     = externalField.geo.bmagInv
      self.bmagInvSq   = externalField.geo.bmagInvSq
      self.jacobGeo    = externalField.geo.jacobGeo
      self.jacobGeoInv = externalField.geo.jacobGeoInv

      self.jacobGeoDbmagSq = self:allocMoment()
      self.confWeakMultiply:advance(0., {self.jacobGeo, self.bmagInvSq}, {self.jacobGeoDbmagSq})

      -- Compute the magnetic field in the center of the domain (e.g. for the Poisson equation).
      local xMid = {}
      for d = 1,self.cdim do xMid[d]=self.confGrid:mid(d) end
      self.bmagMid = self.bmagFunc(0.0, xMid)

   end

   -- Minimum vtSq supported by the grid (for p=1 only for now):
   local TparMin  = (self.mass/6.)*self.grid:dx(self.cdim+1)
   local TperpMin = self.vdim==1 and TparMin or (self.bmagMid/3.)*self.grid:dx(self.cdim+2)
   self.vtSqMinSupported = (TparMin + 2.*TperpMin)/(3.*self.mass)

   -- Create solvers for collisions.
   for _, c in lume.orderedIter(self.collisions) do c:createSolver(self, externalField) end

   -- Create BC solvers.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:createSolver(self, field, externalField) end

   self.hasSheathBCs = false
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do
      self.hasSheathBCs = self.hasSheathBCs or (bc.bcKind=="sheath" and true or false)
   end

   -- Create updater to advance solution by one time-step.
   if self.evolveCollisionless then
      self.solver = Updater.GyrokineticDG {
         onGrid     = self.grid,       confRange = self.bmag:localRange(), 
         confBasis  = self.confBasis,  charge    = self.charge,
         phaseBasis = self.basis,      mass      = self.mass,
      }
      self.collisionlessAdvance = function(tCurr, inFlds, outFlds)
         self.solver:advance(tCurr, inFlds, outFlds)
      end
      
--      if hasApar then
--         -- Set up solver that adds on volume term involving dApar/dt and the entire vpar surface term.
--         self.equationStep2 = GyrokineticEq.GkEqStep2 {
--            onGrid     = self.grid,       charge     = self.charge,
--            phaseBasis = self.basis,      mass       = self.mass,
--            confBasis  = self.confBasis,
--         }
--
--         -- Note that the surface update for this term only involves the vpar direction.
--         self.solverStep2 = Updater.HyperDisCont {
--            onGrid   = self.grid,           zeroFluxDirections = self.zeroFluxDirections,
--            basis    = self.basis,          updateDirections   = {self.cdim+1},
--            cfl      = self.cfl,            clearOut           = false,   -- Continue accumulating into output field.
--            equation = self.equationStep2,  globalUpwind       = false,   -- Don't reduce max speed.
--         }
--      end

      -- Boundary flux updater.
      self.boundaryFluxSlvr = Updater.BoundaryFluxCalc {
         onGrid = self.grid,  equation    = self.solver:getEquation(),
         cdim   = self.cdim,  equation_id = "gyrokinetic",
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
   self.numDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM0", -- GkM0 = < f >
   }
   self.momDensityCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM1", -- GkM1 = < v_parallel f > 
   }
   self.ptclEnergyCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM2", -- GkM2 = < (v_parallel^2 + 2*mu*B/m) f >
   }
   self.M2parCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM2par", -- GkM2par = < v_parallel^2 f >
   }
   self.M3parCalc = Updater.DistFuncMomentCalc {
      onGrid     = self.grid,   confBasis  = self.confBasis,
      phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
      moment     = "GkM3par", -- GkM3par = < v_parallel^3 f >
   }
   if self.vdim > 1 then
      self.M2perpCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis  = self.confBasis,
         phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
         moment     = "GkM2perp", -- GkM2 = < (mu*B/m) f >
      }
      self.M3perpCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis  = self.confBasis,
         phaseBasis = self.basis,  gkfacs     = {self.mass, self.bmag},
         moment     = "GkM3perp", -- GkM3perp = < vpar*(mu*B/m) f >
      }
   end
   self.calcMaxwell = Updater.MaxwellianOnBasis {
      onGrid     = self.grid,   confGrid  = self.confGrid,
      phaseBasis = self.basis,  confBasis = self.confBasis,
      mass       = self.mass,
   }
   if self.needThreeMoments then
      -- Create updater to compute M0, M1, M2 moments.
      self.threeMomentsCalc = Updater.DistFuncMomentCalc {
         onGrid     = self.grid,   confBasis = self.confBasis,
         phaseBasis = self.basis,  moment    = "GkThreeMoments",
         gkfacs     = {self.mass, self.bmag},
      }
      self.calcSelfCouplingMom = function(tCurr, fIn)
         -- Compute M0, M1i and M2.
         self.threeMomentsCalc:advance(tCurr, {fIn}, {self.threeMoments})
         -- Copy number density to its own field.
         self.numDensity:combineOffset(1., self.threeMoments, 0)
      end
   else
      self.calcSelfCouplingMom = function(tCurr, fIn)
         self.numDensityCalc:advance(tCurr, {fIn}, { self.numDensity })
      end
   end

   -- Create an updater for volume integrals. Used by diagnostics.
   -- Placed in a table with key 'scalar' to keep consistency with VlasovSpecies (makes diagnostics simpler).
   self.volIntegral = {
      scalar = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.confGrid,   numComponents = 1,
         basis  = self.confBasis,  quantity      = "V",
      }
   }

   -- Select the function that returns the mass density factor for the polarization (Poisson equation).
   -- Allow user to specify polarization density factor (n in polarization density).
   -- If the polarization is linearized, it should be specified in the input file (as a number, file, or function).
   self.polWeight = self:allocMoment() -- Polarization weight mass*jacobGeo*n0/B^2 (for Poisson equation).
   if field.linearizedPolarization then
      local evOnNodes = Updater.EvalOnNodes {
         onGrid = self.confGrid,   evaluate = function(t, xn) return 1. end,
         basis  = self.confBasis,  onGhosts = false, --true,
      }
      if self.tbl.polarizationDensityFactor == nil then
         print("*** App.Species.GkSpecies: WARNING... not specifying 'polarizationDensityFactor' and relying on n0 in the input file will be deprecated. Please change your input file to specify 'polarizationDensityFactor' (the density factor in the Poisson equation). ***")
         local den0 = self.tbl.n0 or n0
         evOnNodes:setFunc(function(t,xn) return den0*self.mass/(self.bmagMid^2) end)
         evOnNodes:advance(0., {}, {self.polWeight})
      else
         local polDenFacIn = self.tbl.polarizationDensityFactor
   
         if type(polDenFacIn) == "number" then
            evOnNodes:setFunc(function(t,xn) return polDenFacIn*self.mass/(self.bmagMid^2) end)
            evOnNodes:advance(0., {}, {self.polWeight})
         elseif type(polDenFacIn) == "string" or type(polDenFacIn) == "function" then
            if type(polDenFacIn) == "string" then
               self.distIo:read(self.polWeight, polDenFacIn) --, true)
            else
               evOnNodes:setFunc(polDenFacIn)
               evOnNodes:advance(0., {}, {self.polWeight})
            end
   
            -- Apply open BCs (although BCs here should not matter/be used).
            local function makeOpenBcUpdater(dir, edge)
               local bcOpen = function(dir, tm, idxIn, fIn, fOut)
                  self.confBasis:flipSign(dir, fIn, fOut)   -- Requires skinLoop = "pointwise".
               end
               return Updater.Bc {
                  onGrid   = self.confGrid,  edge = edge,
                  dir      = dir,            boundaryConditions = {bcOpen},
                  skinLoop = "pointwise",
               }
            end
            local openBCupdaters = {}
            for dir = 1, self.cdim do
               if not lume.any(self.confGrid:getPeriodicDirs(), function(t) return t==dir end) then
                  openBCupdaters["lower"] = makeOpenBcUpdater(dir, "lower")
                  openBCupdaters["upper"] = makeOpenBcUpdater(dir, "upper")
               end
            end
            for _, bc in pairs(openBCupdaters) do bc:advance(0., {}, {self.polWeight}) end
            self.polWeight:sync(true)
   
            self.distIo:write(self.polWeight, string.format("%s_polarizationDensityFactor_%d.bp", self.name, self.diagIoFrame), 0., self.diagIoFrame, false) --true)
   
            self.polWeight:scale(self.mass)
            self.confWeakMultiply:advance(0., {self.bmagInvSq, self.polWeight}, {self.polWeight})
         end
      end

      -- Include a factor of jacobGeo (included naturally when linearizedPolarization=false).
      self.confWeakMultiply:advance(0., {self.jacobGeo, self.polWeight}, {self.polWeight})

      self.getPolWeight = function() return self.polWeight end
   else
      self.getPolWeight = function()
         self.polWeight:combine(self.mass, self.numDensity)
         self.confWeakMultiply:advance(0., {self.jacobGeoDbmagSq, self.polWeight}, {self.polWeight})
         return self.polWeight
      end
   end

   -- Create species source solvers.
   for _, src in lume.orderedIter(self.sources) do src:createSolver(self, externalField) end
end

function GkSpecies:createCouplingSolver(population, field, externalField)
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

-- Note: do not call applyBc here. It is called later in initialization sequence.
function GkSpecies:initDist(extField)
   if self.randomseed then
      math.randomseed(self.randomseed)
   else
      math.randomseed(47*Mpi.Comm_rank(Mpi.COMM_WORLD)+os.time())
   end

   local syncPeriodicDirs = true
   if self.fluctuationBCs then syncPeriodicDirs = false end

   local initCnt, backgroundCnt = 0, 0
   local scaleInitWithSourcePower = false
   for nm, pr in lume.orderedIter(self.projections) do
      pr:fullInit(self)
      pr:advance(0.0, {extField}, {self.distf[2]})
      if string.find(nm,"init") then
         self.distf[1]:accumulate(1.0, self.distf[2])
         initCnt = initCnt + 1
         if pr.scaleWithSourcePower then scaleInitWithSourcePower = true end
      end
      if string.find(nm,"background") then
         if not self.fBackground then self.fBackground = self:allocDistf() end
         self.fBackground:accumulate(1.0, self.distf[2])
         self.fBackground:sync(syncPeriodicDirs)
         backgroundCnt = backgroundCnt + 1
      end
   end

   if scaleInitWithSourcePower then
      -- MF 2021/05/27: This assumes there's only one source object per species in the input file.
      for _, src in lume.orderedIter(self.sources) do self.distf[1]:scale(src.powerScalingFac) end
   end

   assert(initCnt>0, string.format("GkSpecies: Species '%s' not initialized!", self.name))
   if self.fBackground then
      if backgroundCnt == 0 then self.fBackground:copy(self.distf[1]) end
      self.fBackground:write(string.format("%s_background_%d.bp", self.name, self.diagIoFrame), 0., self.diagIoFrame, true)
   end

   if self.fluctuationBCs then
      assert(backgroundCnt > 0, "GkSpecies: must specify an initial background distribution with 'background' in order to use fluctuation-only BCs")
   end

   self.distf[2]:clear(0.0)
end

function GkSpecies:setActiveRKidx(rkIdx) self.activeRKidx = rkIdx end

function GkSpecies:initCrossSpeciesCoupling(population)
   -- This method establishes the interaction between different
   -- species that is not mediated by the field (solver), like
   -- collisions.

   local species = population:getSpecies()

   -- Determine if M0, M1i and M2 are needed.
   self.needThreeMoments = false

   -- Create a double nested table of colliding species.
   -- In this table we will encode information about that collition such as:
   --   * does the collision take place?
   --   * Operator modeling the collision.
   -- Other features of a collision may be added in the future
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
                  self.needThreeMoments = true  -- MF 2022/09/16: currently all collision models need M0, M1, M2.
               end
            end
         else
            -- This species does not collide with anyone.
            self.collPairs[sN][sO].on = false
         end
      end
   end

   local isThisSpeciesMine = population:isSpeciesMine(self.name)
   -- Allocate threeMoments if we collide with other species not in this rank.
   if self.needThreeMoments and (not isThisSpeciesMine) then
      self.threeMoments = self:allocVectorMoment(3)
   end

   -- Create list of ranks we need to send/recv local threeMoments to/from.
   self.threeMomentsXfer = {}
   self.threeMomentsXfer.destRank, self.threeMomentsXfer.srcRank  = {}, {}
   self.threeMomentsXfer.sendReqStat, self.threeMomentsXfer.recvReqStat = nil, nil
   for sO, info in pairs(self.collPairs[self.name]) do
      local sOrank = population:getSpeciesOwner(sO)
      local selfRank = population:getSpeciesOwner(self.name)
      if sO~=self.name and info.on then
         if isThisSpeciesMine then
            -- Only species owned by this rank send threeMoments to other ranks.
            if #self.threeMomentsXfer.destRank == 0 and (not population:isSpeciesMine(sO)) then
               table.insert(self.threeMomentsXfer.destRank, sOrank)
               self.threeMomentsXfer.sendReqStat = Mpi.RequestStatus()
            end
         else
            -- Only species not owned by this rank receive threeMoments from other ranks.
            if #self.threeMomentsXfer.srcRank == 0 and (not population:isSpeciesMine(self.name)) then
               table.insert(self.threeMomentsXfer.srcRank, selfRank)
               self.threeMomentsXfer.recvReqStat = Mpi.RequestStatus()
            end
         end
       end
    end

   -- Initialize the BC cross-coupling interactions.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:initCrossSpeciesCoupling(species) end

end

function GkSpecies:advance(tCurr, population, emIn, inIdx, outIdx)
   local tmStart = Time.clock()

   self:setActiveRKidx(inIdx)
   self.tCurr = tCurr
   local fIn     = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em     = emIn[1]:rkStepperFields()[inIdx] -- Dynamic fields (e.g. phi, Apar)
   local emFunc = emIn[2]:rkStepperFields()[1]     -- Geometry/external field.

   local dApardtProv = emIn[1].dApardtProv

   -- Clear RHS, because HyperDisCont set up with clearOut = false.
   fRhsOut:clear(0.0)

   -- Do collisions first so that collisions contribution to cflRate is included in GK positivity.
   self.timers.collisions = 0.
   for _, c in lume.orderedIter(self.collisions) do
      c:advance(tCurr, fIn, population, {fRhsOut, self.cflRateByCell})
      self.timers.collisions = self.timers.collisions + c:getTimer('advance')
   end

   self.collisionlessAdvance(tCurr, {fIn, em, emFunc, dApardtProv}, {fRhsOut, self.cflRateByCell})
   self.timers.collisionless = self.solver.totalTime

   self.collisionlessBoundaryAdvance(tCurr, {fIn}, {fRhsOut})

   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do
      bc:storeBoundaryFlux(tCurr, outIdx, fRhsOut)   -- Save boundary fluxes.
   end
   emIn[1]:useBoundaryFlux(tCurr, outIdx)  -- Some field objects need to use the boundary fluxes right away.

   self.timers.sources = 0.
   for _, src in lume.orderedIter(self.sources) do
      src:advance(tCurr, fIn, population:getSpecies(), fRhsOut)
      self.timers.sources = self.timers.sources + src:getTimer('advance')
   end

   self.timers.advance = self.timers.advance + Time.clock() - tmStart
end

function GkSpecies:advanceStep2(tCurr, population, emIn, inIdx, outIdx)
   local tmStart = Time.clock()
   local fIn = self:rkStepperFields()[inIdx]
   local fRhsOut = self:rkStepperFields()[outIdx]

   local em = emIn[1]:rkStepperFields()[inIdx]
   local dApardtProv = emIn[1].dApardtProv
   local emFunc = emIn[2]:rkStepperFields()[1]

   if self.evolveCollisionless then
      self.solverStep2:setDtAndCflRate(self.dtGlobal[0], self.cflRateByCell)
      self.solverStep2:advance(tCurr, {fIn, em, emFunc, dApardtProv}, {fRhsOut})
   end
   self.timers.advance = self.timers.advance + Time.clock() - tmStart
end

function GkSpecies:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx)
   -- Perform some operations after the updates have been computed, but before
   -- the combine RK (in PlasmaOnCartGrid) is called.
   local tmStart = Time.clock()

   local species = population:getSpecies()

   -- Wait to finish sending threeMoments if needed.
   population:speciesXferField_waitSend(self.threeMomentsXfer)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:advanceCrossSpeciesCoupling(tCurr, population, emIn, inIdx, outIdx)
   end

   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:advanceCrossSpeciesCoupling(tCurr, species, outIdx) end
   self.timers.advancecross = self.timers.advancecross + Time.clock() - tmStart
end

function GkSpecies:createDiagnostics(field)
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
      self.diagnostics[self.name] = DiagsApp{implementation = GkDiags()}
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

function GkSpecies:calcCouplingMoments(tCurr, rkIdx, species)
   -- Compute moments needed in coupling to fields and collisions.
   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[rkIdx]

   fIn = self.getF_or_deltaF(fIn)  -- Return full-F, or compute and return fluctuations.

   -- Compute M0, M1i and/or M2 depending on what fields and collisions need.
   self.calcSelfCouplingMom(tCurr, fIn)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:calcCouplingMoments(tCurr, rkIdx, species)
   end

   self.timers.mom = self.timers.mom + Time.clock() - tmStart
end

function GkSpecies:calcCrossCouplingMoments(tCurr, rkIdx, population)
   -- Perform cross-species calculation related to coupling moments that require the
   -- self-species coupling moments.
   local tmStart = Time.clock()

   -- Begin sending/receiving threeMoments if needed.
   population:speciesXferField_begin(self.threeMomentsXfer, self.threeMoments, 22)

   for _, coll in lume.orderedIter(self.collisions) do
      coll:calcCrossCouplingMoments(tCurr, rkIdx, population)
   end

   self.timers.momcross = self.timers.momcross + Time.clock() - tmStart
end

function GkSpecies:setActiveRKidx(rkIdx) self.activeRKidx = rkIdx end

function GkSpecies:rkStepperFields() return self.distf end

function GkSpecies:getDistF(rkIdx)
   if rkIdx == nil then
      return self:rkStepperFields()[self.activeRKidx]
   else
      return self:rkStepperFields()[rkIdx]
   end
end

function GkSpecies:getFlucF() return self.flucF end

function GkSpecies:copyRk(outIdx, aIdx)
   self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx])

   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:copyBoundaryFluxField(aIdx, outIdx) end
end

function GkSpecies:combineRk(outIdx, a, aIdx, ...)
   local args = {...} -- Package up rest of args as table.
   local nFlds = #args/2
   self:rkStepperFields()[outIdx]:combine(a, self:rkStepperFields()[aIdx])
   for i = 1, nFlds do -- Accumulate rest of the fields.
      self:rkStepperFields()[outIdx]:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]])
   end

   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do
      bc:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   end
end

function GkSpecies:suggestDt()
   if not self.evolve then return GKYL_MAX_DOUBLE end

   local dtSuggested = math.min(self.cfl/self.cflRateByCell:reduce('max')[1], GKYL_MAX_DOUBLE)

   -- If dtSuggested == GKYL_MAX_DOUBLE, it is likely because of NaNs.
   -- If so, return 0 so that no timestep is taken, and we will abort the simulation.
   if dtSuggested == GKYL_MAX_DOUBLE then dtSuggested = 0.0 end

   return dtSuggested
end

function GkSpecies:setDtGlobal(dtGlobal) self.dtGlobal[0] = dtGlobal end

function GkSpecies:clearCFL()
   -- Clear cflRateByCell for next cfl calculation.
   self.cflRateByCell:clear(0.0)
end

function GkSpecies:applyBcIdx(tCurr, field, externalField, inIdx, outIdx, isFirstRk)
   self:applyBc(tCurr, field, externalField, inIdx, outIdx)
end

function GkSpecies:applyBcDontEvolve(tCurr, field, externalField, inIdx, outIdx) end
function GkSpecies:applyBcEvolve(tCurr, field, externalField, inIdx, outIdx)
   -- fIn is total distribution function.
   local tmStart = Time.clock()

   local fIn = self:rkStepperFields()[outIdx]

   local syncPeriodicDirsTrue = true

   if self.fluctuationBCs then
      -- If fluctuation-only BCs, subtract off background before applying BCs.
      fIn:accumulate(-1.0, self.fBackground)
   end

   -- Apply non-periodic BCs (to only fluctuations if fluctuation BCs).
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:advance(tCurr, self, field, externalField, inIdx, outIdx) end

   -- Apply periodic BCs (to only fluctuations if fluctuation BCs)
   fIn:sync(syncPeriodicDirsTrue)

   if self.fluctuationBCs then
      -- Put back together total distribution
      fIn:accumulate(1.0, self.fBackground)

      -- Update ghosts in total distribution, without enforcing periodicity.
      fIn:sync(not syncPeriodicDirsTrue)
   end

   self.timers.bc = self.timers.bc + Time.clock() - tmStart
end
function GkSpecies:applyBc(tCurr, field, externalField, inIdx, outIdx)
   self.applyBcFunc(tCurr, field, externalField, inIdx, outIdx)
end

function GkSpecies:fluidMoments() return self.threeMoments end

function GkSpecies:vtSqMin() return self.vtSqMinSupported end

function GkSpecies:getNumDensity(rkIdx)
   -- If no rkIdx specified, assume numDensity has already been calculated.
   if rkIdx == nil then return self.numDensity end 
   local fIn = self:rkStepperFields()[rkIdx]

   local tmStart = Time.clock()

   fIn = self.getF_or_deltaF(fIn)
   self.numDensityCalc:advance(nil, {fIn}, { self.numDensityAux })

   self.timers.mom = self.timers.mom + Time.clock() - tmStart
   return self.numDensityAux
end

function GkSpecies:getMomDensity(rkIdx)
   -- If no rkIdx specified, assume momDensity has already been calculated.
   if rkIdx == nil then return self.momDensity end 
   local fIn = self:rkStepperFields()[rkIdx]
 
   local tmStart = Time.clock()

   fIn = self.getF_or_deltaF(fIn)
   self.momDensityCalc:advance(nil, {fIn}, { self.momDensityAux })

   self.timers.mom = self.timers.mom + Time.clock() - tmStart
   return self.momDensityAux
end

function GkSpecies:getPolarizationWeight() return self.getPolWeight() end

function GkSpecies:Maxwellian(xn, n0, vdIn, T0)
   local vd   = vdIn or 0.0
   local vt2  = T0/self.mass
   local vpar = xn[self.cdim+1]
   local v2   = (vpar-vd)^2
   if self.vdim > 1 then 
      local mu = xn[self.cdim+2]
      v2 = v2 + 2*math.abs(mu)*self.bmagFunc(0,xn)/self.mass
      return n0*(2*math.pi*vt2)^(-3/2)*math.exp(-v2/(2*vt2))
   else
      return n0*(2*math.pi*vt2)^(-1/2)*math.exp(-v2/(2*vt2))
   end
end

function GkSpecies:projToSource(proj)
   local tbl = proj.tbl
   local pow = tbl.power
   return Source { profile = proj, power = pow }
end

function GkSpecies:write(tm, field, force)
   if self.evolve then

      -- Compute delta-F (if perturbed diagnostics are requested) and put it in self.flucF.
      self.calcDeltaF(self:rkStepperFields()[1])

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
         self.writeFluctuation(tm, self.diagIoFrame, fIn)

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

function GkSpecies:writeRestart(tm)
   -- (The final "true/false" in calls to :write determines writing of ghost cells).
   local writeGhost = false
   if self.hasSheathBCs or self.fluctuationBCs then writeGhost = true end

   self.distIo:write(self.distf[1], string.format("%s_restart.bp", self.name), tm, self.distIoFrame, writeGhost)

   for _, dOb in lume.orderedIter(self.diagnostics) do   -- Write restart diagnostics.
      dOb:writeRestart(tm, self.diagIoFrame, self.dynVecRestartFrame)
   end

   self.dynVecRestartFrame = self.dynVecRestartFrame + 1
end

function GkSpecies:readRestart(field, externalField)
   local readGhost = false
   if self.hasSheathBCs or self.fluctuationBCs then readGhost = true end

   local tm, fr = self.distIo:read(self.distf[1], string.format("%s_restart.bp", self.name), readGhost)
   self.distIoFrame = fr -- Reset internal frame counter.

   -- Set ghost cells.
   self.distf[1]:sync()

   -- Apply BCs (unless skin cells have been read because of special BCs).
   if not self.hasSheathBCs and not self.fluctuationBCs then
      self:applyBc(tm, field, externalField, 1, 1)
   end

   local diagIoFrame_new
   for _, dOb in lume.orderedIter(self.diagnostics) do   -- Read grid and integrated diagnostics.
      local _, dfr = dOb:readRestart()
      diagIoFrame_new = diagIoFrame_new or dfr
      if dfr then
         assert(diagIoFrame_new==dfr, string.format("GkSpecies:readRestart expected diagnostics from previous run to have the same last frame. Instead got %d and %d", diagIoFrame_new, dfr))
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

function GkSpecies:clearTimers() 
   for nm, _ in pairs(self.timers) do self.timers[nm] = 0. end
   self.solver.totalTime = 0.
   for _, c in lume.orderedIter(self.collisions) do c:clearTimers() end
   for _, src in lume.orderedIter(self.sources) do src:clearTimers() end
end

return GkSpecies
