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
local BCs              = require "App.BCs"
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

   -- Determine if user wants diagnostics of the fluctuations.
   -- ~~~~~ Backwards compatibility with old diagnostic specification. To be removed in the future. ~~~~~~ --
   if tbl.diagnosticMoments then
      print("App.Species.KineticSpecies: warning... 'diagnosticMoments' will be deprecated. use 'diagnostics' instead.")
      tbl.diagnostics = tbl.diagnostics or {}
      for nm, v in pairs(tbl.diagnosticMoments) do
         if nm == "perturbed" and v == true then 
            table.insert(tbl.diagnostics, "perturbed")
         elseif type(nm) == "number" then
	    table.insert(tbl.diagnostics, v)
         end
      end
   end
   if tbl.diagnosticIntegratedMoments then
      print("App.Species.KineticSpecies: warning... 'diagnosticIntegratedMoments' will be deprecated. use 'diagnostics' instead.")
      tbl.diagnostics = tbl.diagnostics or {}
      for _, v in ipairs(tbl.diagnosticIntegratedMoments) do table.insert(tbl.diagnostics, v) end
   end
   -- ~~~~~~~~~~~~ End of diagnostics backwards compatibility code. ~~~~~~~~~~~~~~~~~~~~~~~ --

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
      if tbl.bcx[1] == nil or tbl.bcx[2] == nil then assert(false, "KineticSpecies: unsupported BC type") end
      self.bcInDir[1] = {tbl.bcx[1], tbl.bcx[2]}
   end
   if tbl.bcy then
      if tbl.bcy[1] == nil or tbl.bcy[2] == nil then assert(false, "KineticSpecies: unsupported BC type") end
      self.bcInDir[2] = {tbl.bcy[1], tbl.bcy[2]}
   end
   if tbl.bcz then
      if tbl.bcz[1] == nil or tbl.bcz[2] == nil then assert(false, "KineticSpecies: unsupported BC type") end
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
         if not BCs.BCsBase.is(val) then val = self:makeBcApp(bcOb, d, e) end
         if BCs.BCsBase.is(val) then
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
   self.nonconPositivity  = xsys.pickBool(tbl.nonconPositivity, false)
   
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
function KineticSpecies:getNdim() return self.ndim end
function KineticSpecies:getVdim() return self.vdim end
function KineticSpecies:setName(nm) self.name = nm end -- Needs to be called before fullInit().

function KineticSpecies:setCfl(cfl)
   self.cfl = cfl
   for _, c in pairs(self.collisions) do c:setCfl(cfl) end   
end

function KineticSpecies:setIoMethod(ioMethod) self.ioMethod = ioMethod end

function KineticSpecies:setConfBasis(basis)
   self.confBasis = basis
   for _, c in pairs(self.collisions) do c:setConfBasis(basis) end
   for _, src in pairs(self.sources) do src:setConfBasis(basis) end
   for _, bc in pairs(self.nonPeriodicBCs) do bc:setConfBasis(basis) end
end
function KineticSpecies:setConfGrid(grid)
   self.confGrid = grid
   for _, c in pairs(self.collisions) do c:setConfGrid(grid) end
   for _, src in pairs(self.sources) do src:setConfGrid(grid) end
   for _, bc in pairs(self.nonPeriodicBCs) do bc:setConfGrid(grid) end
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

function KineticSpecies:createSolver(field, externalField)
   -- Set up weak multiplication and division operators.
   self.confWeakMultiply = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,   operation = "Multiply",
      weakBasis = self.confBasis,  onGhosts  = true,
   }
   self.confWeakDivide = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,   operation = "Divide",
      weakBasis = self.confBasis,  onGhosts  = true,
   }
   self.confWeakDotProduct = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,   operation = "DotProduct",
      weakBasis = self.confBasis,  onGhosts  = true,
   }
   self.confPhaseWeakMultiply = Updater.CartFieldBinOp {
      onGrid    = self.grid,   fieldBasis = self.confBasis,
      weakBasis = self.basis,  operation  = "Multiply",
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
         return KineticSpecies["applyBcEvolve"](self, tCurr, field, externalField, inIdx, outIdx)
      end
   else
      self.applyBcFunc = function(tCurr, field, externalField, inIdx, outIdx)
         return KineticSpecies["applyBcDontEvolve"](self, tCurr, field, externalField, inIdx, outIdx)
      end
   end

   -- Create solvers for collisions.
   for _, c in pairs(self.collisions) do c:createSolver(externalField) end

   -- Create BC solvers.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:createSolver(self, field, externalField) end

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
   if self.nonconPositivity then
      self.prePosM0     = self:allocMoment()
      self.postPosM0    = self:allocMoment()
      self.delPosM0     = self:allocMoment()
      self.intDelPosM0  = DataStruct.DynVector{numComponents = 1}
      self.intPosM0     = DataStruct.DynVector{numComponents = 1}
      self.calcIntPosM0 = Updater.CartFieldIntegratedQuantCalc {
	    onGrid        = self.confGrid,
	    basis         = self.confBasis,
	    numComponents = 1,
	    quantity      = "V",
      }
      self.nonconPos = Updater.PositivityRescale {
         onGrid = self.grid,
         basis  = self.basis,
         nonconservative = true,
      }
   end
end

function KineticSpecies:createCouplingSolver(species, field, externalField)
   -- After all species have called their createSolver methods, we here create the objects
   -- needed for cross-species solves (e.g. cross-species collisions).

   -- Create BC solvers.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:createCouplingSolver(species, field, externalField) end
end

function KineticSpecies:alloc(nRkDup)
   -- Allocate fields needed in RK update.
   self.distf = {}
   for i = 1, nRkDup do
      self.distf[i] = self:allocDistf()
      self.distf[i]:clear(0.0)
   end
   self:setActiveRKidx(1)

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

   self.flucF = (self.fluctuationBCs or self.perturbedDiagnostics) and self:allocDistf() or nil   -- Fluctuation.

   if self.positivity then self.fPos = self:allocDistf() end

   self.fPrev = self:allocDistf()
   self.fPrev:clear(0.0)

   -- Array with one component per cell to store cflRate in each cell.
   self.cflRateByCell = self:allocCartField(self.grid, 1, {1,1})
   self.cflRatePtr    = self.cflRateByCell:get(1)
   self.cflRateIdxr   = self.cflRateByCell:genIndexer()
   self.dtGlobal      = ffi.new("double[2]")
   self.dtGlobal[0], self.dtGlobal[1] = 1.0, 1.0   -- Temporary value (so diagnostics at t=0 aren't inf).

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
   local scaleInitWithSourcePower = false
   for nm, pr in lume.orderedIter(self.projections) do
      pr:fullInit(self)
      pr:advance(0.0, {extField}, {self.distf[2]})
      -- This barrier is needed as when using MPI-SHM some
      -- processes will get to accumulate before projection is finished.
      Mpi.Barrier(self.grid:commSet().sharedComm)
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
      -- if pr.isReservoir then
      --    if not self.fReservoir then 
      --       self.fReservoir = self:allocDistf()
      --    end
      --    self.fReservoir:accumulate(1.0, self.distf[2])
      -- end
   end

   if scaleInitWithSourcePower then 
      -- MF 2021/05/27: This assumes there's only one source object per species in the input file.
      for _, src in lume.orderedIter(self.sources) do self.distf[1]:scale(src.powerScalingFac) end
   end
   
   assert(initCnt>0, string.format("KineticSpecies: Species '%s' not initialized!", self.name))
   if self.fBackground then
      if backgroundCnt == 0 then self.fBackground:copy(self.distf[1]) end
      self.fBackground:write(string.format("%s_background_%d.bp", self.name, self.diagIoFrame), 0., self.diagIoFrame, true)
   end

   if self.fluctuationBCs then 
      assert(backgroundCnt > 0, "KineticSpecies: must specify an initial background distribution with 'background' in order to use fluctuation-only BCs") 
   end

   self.distf[2]:clear(0.0)
   
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

function KineticSpecies:setActiveRKidx(rkIdx) self.activeRKidx = rkIdx end

function KineticSpecies:rkStepperFields() return self.distf end

function KineticSpecies:getDistF(rkIdx)
   if rkIdx == nil then
      return self:rkStepperFields()[self.activeRKidx]
   else
      return self:rkStepperFields()[rkIdx]
   end
end

function KineticSpecies:getFlucF() return self.flucF end

function KineticSpecies:copyRk(outIdx, aIdx)
   self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx])

   for _, bc in pairs(self.nonPeriodicBCs) do bc:copyBoundaryFluxField(aIdx, outIdx) end

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

   for _, bc in pairs(self.nonPeriodicBCs) do
      bc:combineBoundaryFluxField(outIdx, a, aIdx, ...)
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
     status = self.posChecker:advance(tCurr, {self:rkStepperFields()[idx]}, {})
  end
  return status
end

function KineticSpecies:applyBcIdx(tCurr, field, externalField, inIdx, outIdx, isFirstRk)
  if self.positivityDiffuse then
     self.fDelPos[outIdx]:combine(-1.0, self:rkStepperFields()[outIdx])
     self.posRescaler:advance(tCurr, {self:rkStepperFields()[outIdx]}, {self:rkStepperFields()[outIdx]}, true, isFirstRk)
     self.fDelPos[outIdx]:accumulate(1.0, self:rkStepperFields()[outIdx])
  end
  self:applyBc(tCurr, field, externalField, inIdx, outIdx)
  if self.positivity then
     self:checkPositivity(tCurr, outIdx)
  end
  if self.nonconPositivity then
     self.numDensityCalc:advance(tCurr, {self:rkStepperFields()[outIdx]}, {self.prePosM0})
     
     self.nonconPos:advance(tCurr, {self:rkStepperFields()[outIdx]}, {self:rkStepperFields()[outIdx]}, false)
     
     self.numDensityCalc:advance(tCurr, {self:rkStepperFields()[outIdx]}, {self.postPosM0})
     self.delPosM0:combine(1.0, self.postPosM0, -1.0, self.prePosM0)

     local tm, lv = self.intPosM0:lastData()
     self.calcIntPosM0:advance(tCurr, {self.delPosM0}, {self.intDelPosM0})
     local tmDel, lvDel = self.intDelPosM0:lastData()
     self.intPosM0:appendData(tmDel, {lv[1]+lvDel[1]})
  end
end

function KineticSpecies:applyBcDontEvolve(tCurr, field, externalField, inIdx, outIdx) end
function KineticSpecies:applyBcEvolve(tCurr, field, externalField, inIdx, outIdx)
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

   self.bcTime = self.bcTime + (Time.clock()-tmStart)
end
function KineticSpecies:applyBc(tCurr, field, externalField, inIdx, outIdx)
   self.applyBcFunc(tCurr, field, externalField, inIdx, outIdx)
end

function KineticSpecies:createDiagnostics(field)
   -- Many diagnostics require dividing or multiplying by the Jacobian
   -- (if present). Predefine the functions that do that.
   self.divideByJacobGeo = self.jacobGeoInv
      and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeoInv}, {fldOut}) end
      or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
   self.multiplyByJacobGeo = self.jacobGeo
      and function(tm, fldIn, fldOut) self.confWeakMultiply:advance(tm, {fldIn, self.jacobGeo}, {fldOut}) end
      or function(tm, fldIn, fldOut) fldOut:copy(fldIn) end
end

function KineticSpecies:calcAndWriteDiagnosticMoments(tm)
    -- IMPORTANT: do not use this method anymore. It should disappear. The stuff below will be moved elsewhere (MF).

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

    -- Vlasov positivity diagnostics
    if self.nonconPositivity then
       self.intPosM0:write( string.format("%s_intPosM0.bp", self.name), tm, self.diagIoFrame)
    end
end

function KineticSpecies:isEvolving() return self.evolve end

function KineticSpecies:write(tm, force)
   if self.evolve then

      -- Compute delta-F (if perturbed diagnostics are requested) and put it in self.flucF.
      self.calcDeltaF(self:rkStepperFields()[1])

      for _, dOb in lume.orderedIter(self.diagnostics) do
         dOb:resetState(tm)   -- Reset booleans indicating if diagnostic has been computed.
      end

      for _, bc in pairs(self.nonPeriodicBCs) do
         bc:computeBoundaryFluxRate(self.dtGlobal[0])
      end

      local tmStart = Time.clock()
      -- Compute integrated diagnostics.
      if self.calcIntQuantTrigger(tm) then
         self:calcDiagnosticIntegratedMoments(tm)   -- MF: to be removed. Only here for some neutral diagnostics (to be moved).
         for _, dOb in lume.orderedIter(self.diagnostics) do
            dOb:calcIntegratedDiagnostics(tm, self)   -- Compute integrated diagnostics (this species' and other objects').
         end
      end
      self.integratedMomentsTime = self.integratedMomentsTime + Time.clock() - tmStart

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
         self:calcAndWriteDiagnosticMoments(tm)   -- MF: to be removed. Only here for some neutral diagnostics (to be moved).

         for _, dOb in lume.orderedIter(self.diagnostics) do
            dOb:calcGridDiagnostics(tm, self)   -- Compute grid diagnostics (this species' and other objects').
         end

         for _, dOb in lume.orderedIter(self.diagnostics) do   -- Write grid and integrated diagnostics.
            dOb:write(tm, self.diagIoFrame)
         end

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
         self:calcDiagnosticIntegratedMoments(tm)   -- MF: to be removed. Only here for some neutral diagnostics (to be moved).
         for _, dOb in lume.orderedIter(self.diagnostics) do
            dOb:calcIntegratedDiagnostics(tm, self)   -- Compute integrated diagnostics (this species' and other objects').
         end
         self.integratedMomentsTime = self.integratedMomentsTime + Time.clock() - tmStart

	 self.distIo:write(self.distf[1], string.format("%s_%d.bp", self.name, 0), tm, 0)

	 self:calcAndWriteDiagnosticMoments(tm)   -- MF: to be removed. Only here for some neutral diagnostics (to be moved).

         for _, dOb in lume.orderedIter(self.diagnostics) do
            dOb:calcGridDiagnostics(tm, self)   -- Compute grid diagnostics (this species' and other objects').
         end

         for _, dOb in lume.orderedIter(self.diagnostics) do   -- Write grid and integrated diagnostics.
            dOb:write(tm, self.diagIoFrame)
         end
      end
      self.distIoFrame = self.distIoFrame+1
   end

end

function KineticSpecies:writeRestart(tm)
   -- (The final "true/false" in calls to :write determines writing of ghost cells).
   local writeGhost = false
   if self.hasSheathBCs or self.fluctuationBCs then writeGhost = true end

   self.distIo:write(self.distf[1], string.format("%s_restart.bp", self.name), tm, self.distIoFrame, writeGhost)

   for _, dOb in lume.orderedIter(self.diagnostics) do   -- Write restart diagnostics.
      dOb:writeRestart(tm, self.diagIoFrame, self.dynVecRestartFrame)
   end

   -- The following two should be moved elsehwere (MF).
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

function KineticSpecies:readRestart(field, externalField)
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
      assert(diagIoFrame_new==dfr, "KineticSpecies:readRestart expected diagnostics from previous run to have the same last frame.") 
   end
   -- The 'or self.distIoFrame' below is for sims without diagnostics, or when the first
   -- run didn't request diagnostics, but the latter (when the restart requests diagnostics
   -- while the first one didn't) requires commenting out the loop above (a hack, beware).
   self.diagIoFrame = diagIoFrame_new or self.distIoFrame
   
   -- The following two should be moved elsehwere (MF).
   if self.calcReactRate then
      self.intSrcIzM0:read(string.format("%s_intSrcIzM0_restart.bp", self.name))
   end
   if self.calcIntSrcIz then
      self.intSrcIzM0:read(string.format("%s_intSrcIzM0_restart.bp", self.name))
   end

   -- Iterate triggers.
   self.distIoTrigger(tm)
   self.diagIoTrigger(tm)
   
   return tm
end

-- Timers.
function KineticSpecies:totalSolverTime() return self.solver.totalTime end
function KineticSpecies:totalBcTime() return self.bcTime end
function KineticSpecies:intMomCalcTime() return self.integratedMomentsTime end

return KineticSpecies
