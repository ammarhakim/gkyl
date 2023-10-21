-- Gkyl ------------------------------------------------------------------------
--
-- App support code: FluidSpecies object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Collisions       = require "App.Collisions"
local DataStruct       = require "DataStruct"
local LinearTrigger    = require "Lib.LinearTrigger"
local Mpi              = require "Comm.Mpi"
local Proto            = require "Lib.Proto"
local Projection       = require "App.Projection"
local ProjectionBase   = require "App.Projection.ProjectionBase"
local SourceBase       = require "App.Sources.SourceBase"
local SpeciesBase      = require "App.Species.SpeciesBase"
local BCs              = require "App.BCs"
local DiagsApp         = require "App.Diagnostics.SpeciesDiagnostics"
local FluidDiags       = require "App.Diagnostics.FluidDiagnostics"
local Time             = require "Lib.Time"
local Updater          = require "Updater"
local ffi              = require "ffi"
local xsys             = require "xsys"
local lume             = require "Lib.lume"

-- Base class for kinetic species.
local FluidSpecies = Proto(SpeciesBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function FluidSpecies:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function FluidSpecies:fullInit(appTbl)
   local tbl = self.tbl -- Previously stored table.

   self.charge = tbl.charge or 1.0
   self.mass   = tbl.mass or 1.0

   self.evolve              = xsys.pickBool(tbl.evolve, true) -- Default: evolve species.
   self.evolveCollisionless = xsys.pickBool(tbl.evolveCollisionless, self.evolve)

   local nFrame = tbl.nDiagnosticFrame and tbl.nDiagnosticFrame or appTbl.nFrame
   -- Create triggers to write diagnostics.
   if tbl.nDiagnosticFrame then
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nDiagnosticFrame)
   else
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, nFrame)
   end

   -- Create trigger for how frequently to compute integrated moments.
   -- Do not compute the integrated diagnostics less frequently than we output data.
   if appTbl.calcIntQuantEvery then
      self.calcIntQuantTrigger = LinearTrigger(0, appTbl.tEnd,  math.max(nFrame,math.floor(1/appTbl.calcIntQuantEvery)))
   else
      self.calcIntQuantTrigger = function(t) return true end
   end

   -- Write ghost cells on boundaries of global domain (for BCs).
   self.writeGhost = xsys.pickBool(appTbl.writeGhost, false)

   -- Option to group diagnostics (i.e. it writes one file for all grid diags, and one file
   -- for all integrated diags) in all diagnostic Apps, rather than writing one file for each.
   self.groupDiags = appTbl.groupDiagnostics or false

   -- Get a random seed for random initial conditions.
   self.randomseed = tbl.randomseed

   self.diagnostics = {}  -- Table in which we'll place diagnostic objects.
   -- Determine if user wants diagnostics of the fluctuations.
   self.perturbedDiagnostics = false
   if tbl.diagnostics then
      if lume.any(tbl.diagnostics, function(e) return e=="perturbed" end) then
         lume.remove(tbl.diagnostics,"perturbed")
         self.perturbedDiagnostics = true
      end
   end

   -- Initialize table containing sources (if any).
   self.sources = {}
   for nm, val in pairs(tbl) do
      if SourceBase.is(val) or string.find(nm,"source") then
         self.sources[nm] = val
         val:setSpeciesName(self.name)
         val:setName(nm)   -- Do :setName after :setSpeciesName for sources.
         val:fullInit(tbl) -- Initialize sources
      end
   end
   lume.setOrder(self.sources)  -- Save order in metatable to loop in the same order (w/ orderedIter, better for I/O).

   self.projections = {}
   for nm, val in pairs(tbl) do
      if ProjectionBase.is(val) then self.projections[nm] = val end
   end
   -- The keys 'init' and 'background' can be used to speciy initial conditions
   -- or a background (for perturbed=true sims). These can be functions that would
   -- be projected with the Projection app, or a string which would specify a file to read.
   -- When providing a file/string: specify the suffix of the file (including the extension),
   -- and the prefix is assumed to be the name of the input file.
   if type(tbl.init) == "function" then
      self.projections["init"] = Projection.FluidProjection.FunctionProjection {
         func = function (t, zn) return tbl.init(t, zn, self) end,
      }
   elseif type(tbl.init) == "string" then
      self.projections["init"] = Projection.FluidProjection.ReadInput {
         inputFile = tbl.init,
      }
   end
   if type(tbl.background) == "function" then
      self.projections["init"] = Projection.FluidProjection.FunctionProjection {
         func = function (t, zn) return tbl.background(t, zn, self) end,
      }
   elseif type(tbl.background) == "string" then
      self.projections["init"] = Projection.FluidProjection.ReadInput {
         inputFile = tbl.background,
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
      if tbl.bcx[1] == nil or tbl.bcx[2] == nil then assert(false, "FluidSpecies: unsupported BC type") end
      self.bcInDir[1] = {tbl.bcx[1], tbl.bcx[2]}
   end
   if tbl.bcy then
      if tbl.bcy[1] == nil or tbl.bcy[2] == nil then assert(false, "FluidSpecies: unsupported BC type") end
      self.bcInDir[2] = {tbl.bcy[1], tbl.bcy[2]}
   end
   if tbl.bcz then
      if tbl.bcz[1] == nil or tbl.bcz[2] == nil then assert(false, "FluidSpecies: unsupported BC type") end
      self.bcInDir[3] = {tbl.bcz[1], tbl.bcz[2]}
   end
   if tbl.bcvx then
      if tbl.bcvx[1] == nil or tbl.bcvx[2] == nil then assert(false, "FluidSpecies: unsupported BC type") end
      self.bcInDir[4] = {tbl.bcvx[1], tbl.bcvx[2]}
   end
   if tbl.bcvy then
      if tbl.bcvy[1] == nil or tbl.bcvy[2] == nil then assert(false, "FluidSpecies: unsupported BC type") end
      self.bcInDir[5] = {tbl.bcvy[1], tbl.bcvy[2]}
   end
   if tbl.bcvz then
      if tbl.bcvz[1] == nil or tbl.bcvz[2] == nil then assert(false, "FluidSpecies: unsupported BC type") end
      self.bcInDir[6] = {tbl.bcvz[1], tbl.bcvz[2]}
   end
   -- Initialize boundary conditions.
   self.nonPeriodicBCs = {}
   local dirLabel  = {'X','Y','Z'}
   local edgeLabel = {'lower','upper'}
   for d, bcsTbl in ipairs(self.bcInDir) do
      for e, bcOb in ipairs(bcsTbl) do
         local goodBC = false
         local val    = bcOb
         if not BCs.BCsBase.is(val) then val = self:makeBcApp(bcOb) end
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
         assert(goodBC, "FluidSpecies: bc not recognized.")
      end
   end
   lume.setOrder(self.nonPeriodicBCs)  -- Save order in metatable to loop in the same order (w/ orderedIter, better for I/O).
   
   -- Collisions: currently used for a diffusion term.
   self.collisions = {}
   for nm, val in pairs(tbl) do
      if Collisions.CollisionsBase.is(val) then
         self.collisions[nm] = val
         self.collisions[nm]:setName(nm)
         val:setSpeciesName(self.name)
         val:fullInit(tbl)    -- Initialize collisions (e.g. diffusion).
      end
   end

   self.positivity        = xsys.pickBool(tbl.applyPositivity, false)
   self.positivityDiffuse = xsys.pickBool(tbl.positivityDiffuse, self.positivity)
   self.positivityRescale = xsys.pickBool(tbl.positivityRescale, false)

   self.diagIoFrame        = 0 -- Frame number for diagnostics.
   self.dynVecRestartFrame = 0 -- Frame number of restarts (for DynVectors only).
   self.cfl      = 0.1
   if (not self.nMoments) then self.nMoments = 1 end -- Default to a single moment.
   if (not self.nGhost) then self.nGhost   = 1 end  -- Default is 1 ghost-cell in each direction.
   
   self.integratedMomentsTime = 0.0 -- Timer for integrated moments.
   self.bcTime = 0.0   -- Timer for BCs.
end

function FluidSpecies:getCharge() return self.charge end
function FluidSpecies:getMass() return self.mass end
function FluidSpecies:getNdim() return self.ndim end
function FluidSpecies:setName(nm) self.name = nm end

function FluidSpecies:setCfl(cfl)
   self.cfl = cfl
   for _, c in pairs(self.collisions) do c:setCfl(cfl) end
   for _, src in pairs(self.sources) do src:setCfl(cfl) end
end

function FluidSpecies:setConfBasis(basis)
   self.basis = basis
   for _, c in pairs(self.collisions) do c:setConfBasis(basis) end
   for _, src in pairs(self.sources) do src:setConfBasis(basis) end
   for _, bc in pairs(self.nonPeriodicBCs) do bc:setConfBasis(basis) end
end
function FluidSpecies:setConfGrid(grid, myadios)
   self.grid    = grid
   self.myadios = myadios
   self.ndim    = self.grid:ndim()
   for _, c in pairs(self.collisions) do c:setConfGrid(grid) end
   for _, src in pairs(self.sources) do src:setConfGrid(grid) end
   for _, bc in pairs(self.nonPeriodicBCs) do bc:setConfGrid(grid) end
end

function FluidSpecies:createGrid(grid)
   self.grid = grid
   -- Output of grid file is placed here so we can associate a grid file
   -- with this species alone (e.g. what if one species is evolved on a
   -- coarser/finer mesh than another one).
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
function FluidSpecies:allocCartField(grid, nComp, ghosts, metaData)
   local f = DataStruct.Field {
      onGrid        = grid,   ghost    = ghosts,
      numComponents = nComp,  metaData = metaData,
   }
   f:clear(0.0)
   return f
end
function FluidSpecies:allocMoment()
   local metaData = {polyOrder = self.basis:polyOrder(),
                     basisType = self.basis:id(),
                     charge    = self.charge,
                     mass      = self.mass,
                     grid      = GKYL_OUT_PREFIX .. "_" .. self.name .. "_grid.bp"}
   return self:allocCartField(self.grid, self.basis:numBasis(), {self.nGhost,self.nGhost}, metaData)
end
function FluidSpecies:allocVectorMoment(dim)
   local metaData = {polyOrder = self.basis:polyOrder(),
                     basisType = self.basis:id(),
                     charge    = self.charge,
                     mass      = self.mass,
                     grid      = GKYL_OUT_PREFIX .. "_" .. self.name .. "_grid.bp"}
   return self:allocCartField(self.grid, dim*self.basis:numBasis(), {self.nGhost,self.nGhost}, metaData)
end

function FluidSpecies:createSolver(field, externalField)
   if externalField then
      -- Set up Jacobian.
      self.jacobFun = externalField.jacobGeoFunc
      if externalField.geo then
         self.jacob    = externalField.geo.jacobGeo
         self.jacobInv = externalField.geo.jacobGeoInv
      end
   end

   -- Operators needed for time-dependent calculation and diagnostics.
   if self.ndim <= 3 and self.jacob then
      self.weakMultiply = Updater.CartFieldBinOp {
         operation = "Multiply",  weakBasis = self.basis,
         onGhosts  = true,
      }
      self.weakDivide = Updater.CartFieldBinOp {
         operation = "Divide",  weakBasis = self.basis,
         onRange   = self.jacob:localExtRange(),  onGhosts = true,
      }
   end
   self.volIntegral = {
      scalar = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.grid,   numComponents = 1,
         basis  = self.basis,  quantity      = "V"
      },
      vector = Updater.CartFieldIntegratedQuantCalc {
         onGrid = self.grid,   numComponents = self.nMoments,
         basis  = self.basis,  quantity      = "V"
      }
   }

   -- Create solvers for collisions (e.g. diffusion).
   for _, c in pairs(self.collisions) do c:createSolver(self, externalField) end

   -- Create solvers for sources.
   for _, src in lume.orderedIter(self.sources) do src:createSolver(self, externalField) end

   -- Create BC solvers.
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:createSolver(self, field, externalField) end

   -- Create positivity functions.
   if self.positivity then
      self.posChecker = Updater.PositivityCheck {
         onGrid = self.grid,
         basis  = self.basis,
      }
      self.checkPositivity = function(tCurr, idx)
         return self.posChecker:advance(tCurr, {self:rkStepperFields()[idx]}, {})
      end
   else
      self.checkPositivity = function(tCurr, idx) end
   end
   if self.positivityRescale or self.positivityDiffuse then
      self.posRescaler = Updater.PositivityRescale {
         onGrid = self.grid,
         basis  = self.basis,
      }
      if self.positivityDiffuse then
         self.posRescalerDiffAdv = function(tCurr, rkIdx, computeDiagnostics, zeroOut)
            self.posRescaler:advance(tCurr, {self:rkStepperFields()[rkIdx]}, {self:rkStepperFields()[rkIdx]},
                                     computeDiagnostics, zeroOut)
         end
         self.posRescalerDiffWrite = function(tm, fr) self.posRescaler:write(tm, fr, self.name) end
      else
         self.posRescalerDiffAdv   = function(tCurr, rkIdx, computeDiagnostics, zeroOut) end
         self.posRescalerDiffWrite = function(tm, fr) end
      end
   else
      self.posRescaler          = {advance=function(tCurr, inFlds, outFlds, computeDiagnostics, zeroOut) end, totalTime=0.}
      self.posRescalerDiffAdv   = function(tCurr, rkIdx, computeDiagnostics, zeroOut) end
      self.posRescalerDiffWrite = function(tm, fr) end
   end

   -- Functions to compute fluctuations given the current moments and background,
   -- and the full-F moments given the fluctuations and background.
   self.getMom_or_deltaMom = self.deltaF and function(momIn)
      self.flucMom:combine(1.0, momIn, -1.0, self.momBackground)
      return self.flucMom
   end or function(momIn) return momIn end
   self.minusBackgroundMom = self.fluctuationBCs
      and function(momIn) momIn:accumulate(-1.0, self.momBackground) end
      or function(momIn) end
   self.calcFullMom = self.fluctuationBCs
      and function(momIn, syncFullFperiodicDirs)
         momIn:accumulate(1.0, self.momBackground)
         momIn:sync(syncFullFperiodicDirs)
      end or function(momIn, syncFullFperiodicDirs) end 
   self.calcDeltaMom = self.perturbedDiagnostics 
      and function(momIn) self.flucMom:combine(1.0, momIn, -1.0, self.momBackground) end
      or function(momIn) end

   if self.fluctuationBCs or self.perturbedDiagnostics then
      self.writeFluctuation = self.perturbedDiagnostics 
         and function(tm, fr, momIn)
            self.momIo:write(self.flucMom, string.format("%s_fluctuation_%d.bp", self.name, self.diagIoFrame), tm, fr)
         end
         or function(tm, fr, momIn)
            self.calcDeltaMom(momIn)
            self.momIo:write(self.flucMom, string.format("%s_fluctuation_%d.bp", self.name, self.diagIoFrame), tm, fr)
         end
   else
      self.writeFluctuation = function(tm, fr, momIn) end
   end

   if self.evolve then
      self.suggestDtFunc = function() return FluidSpecies["suggestDtEvolve"](self) end
      self.applyBcFunc   = function(tCurr, field, externalField, inIdx, outIdx)
         return FluidSpecies["applyBcEvolve"](self, tCurr, field, externalField, inIdx, outIdx)
      end
      self.calcCouplingMomentsFunc = function(tCurr, rkIdx, species)
         return self:calcCouplingMomentsEvolve(tCurr, rkIdx, species)
      end
   else
      self.suggestDtFunc = function() return FluidSpecies["suggestDtDontEvolve"](self) end
      self.applyBcFunc   = function(tCurr, field, externalField, inIdx, outIdx)
         return FluidSpecies["applyBcDontEvolve"](self, tCurr, field, externalField, inIdx, outIdx)
      end
      self.calcCouplingMomentsFunc = function(tCurr, rkIdx, species)
         return self:calcCouplingMomentsNoEvolve(tCurr, rkIdx, species)
      end
   end
end

function FluidSpecies:alloc(nRkDup)
  
   -- Allocate fields needed in RK update.
   self.moments = {}
   for i = 1, nRkDup do
      self.moments[i] = self:allocVectorMoment(self.nMoments)
      self.moments[i]:clear(0.0)
   end
   self:setActiveRKidx(1)

   -- Create Adios object for field I/O.
   self.momIo = AdiosCartFieldIo {
      elemType = self.moments[1]:elemType(),
      metaData = {polyOrder = self.basis:polyOrder(),
                  basisType = self.basis:id(),
                  charge    = self.charge,
                  mass      = self.mass,
                  grid      = GKYL_OUT_PREFIX .. "_" .. self.name .. "_grid.bp",},
   }

   self.couplingMoments = self:allocVectorMoment(self.nMoments)

   self.unitField = self:allocMoment()
   local projectUnity = Updater.ProjectOnBasis {
      onGrid   = self.grid,   evaluate = function(t, zn) return 1.0 end,
      basis    = self.basis,  onGhosts = true,
   }
   projectUnity:advance(0.0, {}, {self.unitField})

   self.noJacMom = self:allocVectorMoment(self.nMoments)   -- Moments without Jacobian.

   self.flucMom = (self.fluctuationBCs or self.perturbedDiagnostics) and self:allocVectorMoment(self.nMoments) or nil   -- Fluctuation.

   if self.positivityRescale then self.momPos = self:allocVectorMoment(self.nMoments) end

   -- Array with one component per cell to store cflRate in each cell.
   self.cflRateByCell = self:allocCartField(self.grid, 1, {1,1})
   self.cflRatePtr    = self.cflRateByCell:get(1)
   self.cflRateIdxr   = self.cflRateByCell:genIndexer()
   self.dtGlobal      = ffi.new("double[2]")
   self.dtGlobal[0], self.dtGlobal[1] = 1.0, 1.0   -- Temporary value (so diagnostics at t=0 aren't inf).

   -- Create a table of flags to indicate whether primitive.
   -- At first we consider 3 flags: 
   --   1) self primitive moments (e.g. uParSelf, TparSelf, TperpSelf).
   --   2) cross primitive moments (e.g. uParCross, TparCross, TperpCross).
   --   3) spatially varying cross-species collisionality (varNuXCross).
   self.momentFlags = {}
   self.momentFlags[1] = false
   -- The 2nd and 3rd entries need a table to store
   -- a flag for each pair of species colliding.
   self.momentFlags[2], self.momentFlags[3] = {}, {}
end

-- Note: do not call applyBc here. It is called later in initialization sequence.
function FluidSpecies:initDist(extField, species)

   if self.randomseed then
      math.randomseed(self.randomseed)
   else
      math.randomseed(47*Mpi.Comm_rank(Mpi.COMM_WORLD)+os.time())
   end

   local syncPeriodicDirs = true
   if self.fluctuationBCs then syncPeriodicDirs = false end

   local initCnt, backgroundCnt = 0, 0
   local scaleInitWithSourcePower = false
   for nm, pr in pairs(self.projections) do
      pr:fullInit(self)
      pr:advance(0.0, {extField}, {self.moments[2]})
      if string.find(nm,"init") then
         self.moments[1]:accumulate(1.0, self.moments[2])
         initCnt = initCnt + 1
         if pr.scaleWithSourcePower then scaleInitWithSourcePower = true end
      end
      if string.find(nm,"background") then
         if not self.momBackground then self.momBackground = self:allocVectorMoment(self.nMoments) end
         self.momBackground:accumulate(1.0, self.moments[2])
         self.momBackground:sync(syncPeriodicDirs)
         backgroundCnt = backgroundCnt + 1
      end
   end

   if scaleInitWithSourcePower then
      -- MF 2021/05/27: This assumes there's only one source object per species in the input file.
      for _, src in lume.orderedIter(self.sources) do self.moments[1]:scale(src.powerScalingFac) end
   end

   assert(initCnt>0, string.format("FluidSpecies: Species '%s' not initialized!", self.name))
   if self.momBackground then
      if backgroundCnt == 0 then self.momBackground:copy(self.moments[1]) end
      self.momBackground:write(string.format("%s_background_%d.bp", self.name, self.diagIoFrame), 0., self.diagIoFrame, true)
   end

   if self.fluctuationBCs then
      assert(backgroundCnt > 0, "FluidSpecies: must specify an initial background distribution with 'background' in order to use fluctuation-only BCs")
   end

   self.moments[2]:clear(0.0)

   self.posRescaler:advance(0.0, {self.moments[1]}, {self.moments[1]}, false)

   -- Compute the initial coupling moments.
   self:calcCouplingMomentsEvolve(0.0, 1, species)
end

function FluidSpecies:setActiveRKidx(rkIdx) self.activeRKidx = rkIdx end

function FluidSpecies:rkStepperFields() return self.moments end

function FluidSpecies:calcCouplingMomentsNoEvolve(tCurr, rkIdx)
   self:calcCouplingMomentsEvolve(tCurr, rkIdx)
end
function FluidSpecies:calcCouplingMoments(tCurr, rkIdx, species)
  self.calcCouplingMomentsFunc(tCurr, rkIdx, species)
end

function FluidSpecies:getMoments(rkIdx)
   return rkIdx and self:rkStepperFields()[rkIdx] or self:rkStepperFields()[self.activeRKidx]
end

function FluidSpecies:getFlucMom() return self.flucMom end

function FluidSpecies:copyRk(outIdx, aIdx)
   self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx])

   for _, bc in pairs(self.nonPeriodicBCs) do bc:copyBoundaryFluxField(aIdx, outIdx) end
end

function FluidSpecies:combineRk(outIdx, a, aIdx, ...)
   local args  = {...} -- Package up rest of args as table.
   local nFlds = #args/2
   self:rkStepperFields()[outIdx]:combine(a, self:rkStepperFields()[aIdx])
   for i = 1, nFlds do -- Accumulate rest of the fields.
      self:rkStepperFields()[outIdx]:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]])
   end

   for _, bc in pairs(self.nonPeriodicBCs) do
      bc:combineBoundaryFluxField(outIdx, a, aIdx, ...)
   end
end

function FluidSpecies:suggestDtDontEvolve() return GKYL_MAX_DOUBLE end
function FluidSpecies:suggestDtEvolve()
   local dtSuggested = math.min(self.cfl/self.cflRateByCell:reduce('max')[1], GKYL_MAX_DOUBLE)

   -- If dtSuggested == GKYL_MAX_DOUBLE, it is likely because of NaNs.
   -- If so, return 0 so that no timestep is taken, and we will abort the simulation.
   if dtSuggested == GKYL_MAX_DOUBLE then dtSuggested = 0.0 end

   return dtSuggested
end
function FluidSpecies:suggestDt() return self.suggestDtFunc() end

function FluidSpecies:setDtGlobal(dtGlobal) self.dtGlobal[0] = dtGlobal end

function FluidSpecies:clearCFL()
   -- Clear cflRateByCell for next cfl calculation.
   self.cflRateByCell:clear(0.0)
end

function FluidSpecies:clearMomentFlags(species)
   -- Clear the momentFlags table to indicate that primitive moments
   -- (and other quantities) need to be computed again.
   for iF = 1,1 do self.momentFlags[iF] = false end
   for sN, _ in lume.orderedIter(species) do
      if sN ~= self.name then
         self.momentFlags[2][sN], self.momentFlags[3][sN] = false, false
      end
   end
end

function FluidSpecies:advance(tCurr, species, emIn, inIdx, outIdx)
   -- This only clears the RHS. The :advance method should otherwise be
   -- defined in the child species.
   self:setActiveRKidx(inIdx)
   local momRhsOut = self:rkStepperFields()[outIdx]

   momRhsOut:clear(0.0)
end

function FluidSpecies:applyBcIdx(tCurr, field, externalField, inIdx, outIdx, isFirstRk)
  self.posRescalerDiffAdv(tCurr, outIdx, true, isFirstRk)

  self:applyBc(tCurr, field, externalField, inIdx, outIdx)

  self.checkPositivity(tCurr, outIdx)
end

function FluidSpecies:applyBcDontEvolve(tCurr, field, externalField, inIdx, outIdx) end
function FluidSpecies:applyBcEvolve(tCurr, field, externalField, inIdx, outIdx) 
   local tmStart = Time.clock()

   local momIn = self:rkStepperFields()[outIdx]   -- momIn is the set of evolved moments.

   self.minusBackgroundMom(momIn)

   -- Apply non-periodic BCs (to only fluctuations if fluctuation BCs).
   for _, bc in lume.orderedIter(self.nonPeriodicBCs) do bc:advance(tCurr, self, field, externalField, inIdx, outIdx) end

   -- Apply periodic BCs (to only fluctuations if fluctuation BCs)
   momIn:sync()

   self.calcFullMom(momIn, false)   -- Update ghosts, w/o enforcing periodicity.

   self.bcTime = self.bcTime + Time.clock()-tmStart
end
function FluidSpecies:applyBc(tCurr, field, externalField, inIdx, outIdx)
   self.applyBcFunc(tCurr, field, externalField, inIdx, outIdx)
end

function FluidSpecies:createDiagnostics(field)  -- More sophisticated/extensive diagnostics go in Species/Diagnostics.
   if self.tbl.diagnostics then   -- Create this species' diagnostics.
      self.diagnostics[self.name] = DiagsApp{implementation = FluidDiags()}
      self.diagnostics[self.name]:fullInit(self, field, self)
   end

   for srcNm, src in lume.orderedIter(self.sources) do
      self.diagnostics[self.name..srcNm] = src:createDiagnostics(self, field)
   end

   for bcNm, bc in lume.orderedIter(self.nonPeriodicBCs) do
      self.diagnostics[self.name..bcNm] = bc:createDiagnostics(self, field)
   end
   lume.setOrder(self.diagnostics)

   -- Many diagnostics require dividing by the Jacobian (if present).
   -- Predefine the function that does that.
   self.calcNoJacMom = self.jacobInv
      and function(tm, rkIdx) self.weakMultiply:advance(tm, {self:getMoments(rkIdx), self.jacobInv}, {self.noJacMom}) end
      or function(tm, rkIdx) self.noJacMom:copy(self:getMoments(rkIdx)) end
end

function FluidSpecies:getNoJacMoments() return self.noJacMom end

function FluidSpecies:write(tm, field, force)
   if self.evolve or force then

      -- Calculate fluctuations (if perturbed diagnostics are requested) and put them in self.flucMom.
      self.calcDeltaMom(self:rkStepperFields()[1])

      for _, dOb in lume.orderedIter(self.diagnostics) do
         dOb:resetState(tm)   -- Reset booleans indicating if diagnostic has been computed.
      end

      self.calcNoJacMom(tm, 1)  -- Many diagnostics do not include Jacobian factor.

      for _, bc in pairs(self.nonPeriodicBCs) do
         bc:computeBoundaryFluxRate(self.dtGlobal[0])
      end

      local tmStart = Time.clock()
      -- Compute integrated diagnostics.
      if self.calcIntQuantTrigger(tm) then
         for _, dOb in lume.orderedIter(self.diagnostics) do
            dOb:calcIntegratedDiagnostics(tm, self)   -- Compute integrated diagnostics (this species' and other objects').
         end
      end
      self.integratedMomentsTime = self.integratedMomentsTime + Time.clock() - tmStart
      
      -- Only write stuff if triggered.
      if self.diagIoTrigger(tm) or force then
         local momIn = self:rkStepperFields()[1]

         self.momIo:write(momIn, string.format("%s_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame, self.writeGhost)
         self.writeFluctuation(tm, self.diagIoFrame, momIn)

         for _, src in lume.orderedIter(self.sources) do
            src:write(tm, self.diagIoFrame, self)  -- Allow sources to write (aside from diagnostics).
         end

         for _, dOb in lume.orderedIter(self.diagnostics) do
            dOb:calcGridDiagnostics(tm, self)   -- Compute grid diagnostics (this species' and other objects').
         end
          
         for _, dOb in lume.orderedIter(self.diagnostics) do   -- Write grid and integrated diagnostics.
            dOb:write(tm, self.diagIoFrame)
         end

         for _, c in pairs(self.collisions) do
            c:write(tm, self.diagIoFrame)  -- Allow collisions to write (aside from diagnostics).
         end

         self.posRescalerDiffWrite(tm, self.diagIoFrame)  -- Write positivity diagnostics.

         self.diagIoFrame = self.diagIoFrame+1
      end
   else
      -- If not evolving species, don't write anything except initial conditions.
      if self.diagIoFrame == 0 then
         local momIn = self:rkStepperFields()[1]
         self.momIo:write(momIn, string.format("%s_%d.bp", self.name, 0), tm, 0, self.writeGhost)
      end
      self.diagIoFrame = self.diagIoFrame+1
   end
end

function FluidSpecies:writeRestart(tm)
   local writeGhost = false
   if self.hasSheathBCs or self.fluctuationBCs then writeGhost = true end

   self.momIo:write(self.moments[1], string.format("%s_restart.bp", self.name), tm, self.diagIoFrame, writeGhost)

   for _, dOb in lume.orderedIter(self.diagnostics) do   -- Write restart diagnostics.
      dOb:writeRestart(tm, self.diagIoFrame, self.dynVecRestartFrame)
   end

   self.dynVecRestartFrame = self.dynVecRestartFrame + 1
end

function FluidSpecies:readRestart(field, externalField)
   local readGhost = false
   if self.hasSheathBCs or self.fluctuationBCs then readGhost = true end

   local tm, fr = self.momIo:read(self.moments[1], string.format("%s_restart.bp", self.name), readGhost)
   self.diagIoFrame = fr -- Reset internal frame counter.

   -- Set ghost cells.
   self.moments[1]:sync()

   if not self.hasSheathBCs and not self.fluctuationBCs then
      self:applyBc(tm, field, externalField, 1, 1)
   end

   for _, dOb in lume.orderedIter(self.diagnostics) do   -- Read grid and integrated diagnostics.
      _, _ = dOb:readRestart()
   end
   
   -- Iterate triggers.
   self.diagIoTrigger(tm)

   return tm
end

-- Timers.
function FluidSpecies:totalSolverTime()
   if self.solver then return self.solver.totalTime end
   return 0
end
function FluidSpecies:totalBcTime() return self.bcTime end
function FluidSpecies:momCalcTime() return 0 end
function FluidSpecies:intMomCalcTime() return self.integratedMomentsTime end

return FluidSpecies
