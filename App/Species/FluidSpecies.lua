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
local SpeciesBase      = require "App.Species.SpeciesBase"
local FluidDiags       = require "App.Species.Diagnostics.FluidDiagnostics"
local Time             = require "Lib.Time"
local Updater          = require "Updater"
local ffi              = require "ffi"
local xsys             = require "xsys"
local lume             = require "Lib.lume"

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
   local tbl = self.tbl -- Previously stored table.

   self.charge = tbl.charge or 1.0
   self.mass   = tbl.mass or 1.0

   self.evolve              = xsys.pickBool(tbl.evolve, true) -- Default: evolve species.
   self.evolveCollisionless = xsys.pickBool(tbl.evolveCollisionless, self.evolve)
   self.evolveCollisions    = xsys.pickBool(tbl.evolveCollisions, self.evolve)
   self.evolveSources       = xsys.pickBool(tbl.evolveSources, self.evolve)

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

   -- Get a random seed for random initial conditions.
   self.randomseed = tbl.randomseed

   self.diagnostics = {}  -- Table in which we'll place diagnostic objects.

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
   if type(tbl.source) == "function" then
      self.projections["source"] = Projection.FluidProjection.FunctionProjection {
         func = function(t, zn) return tbl.source(t, zn, self) end,
      }
   elseif type(tbl.source) == "string" then
      self.projections["source"] = Projection.FluidProjection.ReadInput {
         inputFile = tbl.source,
      }
   end

   -- Create a keys metatable in self.projections so we always loop in the same order (better for I/O).
   local projections_keys = {}
   for k in pairs(self.projections) do table.insert(projections_keys, k) end
   table.sort(projections_keys)
   setmetatable(self.projections, projections_keys)

   self.deltaF         = xsys.pickBool(appTbl.deltaF, false)
   self.fluctuationBCs = xsys.pickBool(tbl.fluctuationBCs, false)
   if self.deltaF then self.fluctuationBCs = true end

   self.hasNonPeriodicBc   = false -- To indicate if we have non-periodic BCs.
   self.boundaryConditions = {}    -- List of BCs to apply.
   self.bcx, self.bcy, self.bcz = { }, { }, { }

   -- Read in boundary conditions.
   -- Check to see if bc type is good is now done in createBc.
   if tbl.bcx then
      self.bcx[1], self.bcx[2] = tbl.bcx[1], tbl.bcx[2]
      if self.bcx[1] == nil or self.bcx[2] == nil then assert(false, "FluidSpecies: unsupported BC type") end
      self.hasNonPeriodicBc = true
   end
   if tbl.bcy then
      self.bcy[1], self.bcy[2] = tbl.bcy[1], tbl.bcy[2]
      if self.bcy[1] == nil or self.bcy[2] == nil then assert(false, "FluidSpecies: unsupported BC type") end
      self.hasNonPeriodicBc = true
   end
   if tbl.bcz then
      self.bcz[1], self.bcz[2] = tbl.bcz[1], tbl.bcz[2]
      if self.bcz[1] == nil or self.bcz[2] == nil then assert(false, "FluidSpecies: unsupported BC type") end
      self.hasNonPeriodicBc = true
   end
   
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

   self.useShared         = xsys.pickBool(appTbl.useShared, false)
   self.positivity        = xsys.pickBool(tbl.applyPositivity, false)
   self.positivityDiffuse = xsys.pickBool(tbl.positivityDiffuse, self.positivity)
   self.positivityRescale = xsys.pickBool(tbl.positivityRescale, false)

   self.ioMethod           = "MPI"
   self.diagIoFrame        = 0 -- Frame number for diagnostics.
   self.dynVecRestartFrame = 0 -- Frame number of restarts (for DynVectors only).
   self.cfl      = 0.1
   self.nMoments = 1   -- Default to a single moment.
   self.nGhost   = 1   -- Default is 1 ghost-cell in each direction.

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
end

function FluidSpecies:setIoMethod(ioMethod) self.ioMethod = ioMethod end

function FluidSpecies:setConfBasis(basis)
   self.basis = basis
   for _, c in pairs(self.collisions) do c:setConfBasis(basis) end
end
function FluidSpecies:setConfGrid(cgrid)
   self.grid = cgrid
   self.ndim = self.grid:ndim()
   for _, c in pairs(self.collisions) do c:setConfGrid(cgrid) end
end

function FluidSpecies:createGrid(cgrid)
   self.grid = cgrid
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
      onGrid        = grid,
      numComponents = nComp,
      ghost         = ghosts,
      metaData      = metaData,
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

function FluidSpecies:bcAbsorbFunc(dir, tm, idxIn, fIn, fOut)
   -- Note that for bcAbsorb there is no operation on fIn,
   -- so skinLoop (which determines indexing of fIn) does not matter.
   for i = 1, self.nMoments*self.basis:numBasis() do fOut[i] = 0.0 end
end

function FluidSpecies:bcCopyFunc(dir, tm, idxIn, fIn, fOut)
   for i = 1, self.nMoments*self.basis:numBasis() do fOut[i] = fIn[i] end
end

-- Function to construct a BC updater.
function FluidSpecies:makeBcUpdater(dir, edge, bcList, skinLoop, hasExtFld)
   return Updater.Bc {
      onGrid             = self.grid,
      cdim               = self.ndim,
      dir                = dir,
      edge               = edge,
      boundaryConditions = bcList,
      skinLoop           = skinLoop,
      hasExtFld          = hasExtFld,
   }
end

function FluidSpecies:createBCs()
   -- Create a table to store auxiliary values needed by BCs
   -- and provided by the user in the input file.
   self.auxBCvalues = {}

   -- Functions to make life easier while reading in BCs to apply.
   -- Note: appendBoundaryConditions defined in sub-classes.
   local function handleBc(dir, bc, isPeriodic)
      table.insert(self.auxBCvalues,{nil,nil})
      
      local dirNames = {"x", "y", "z"}
      if (isPeriodic) then
         assert(bc==nil or (bc[1]==nil and bc[2]==nil), "Boundary conditions supplied in periodic direction ".. dirNames[dir]..".")
      end

      if bc[1] then
         self:appendBoundaryConditions(dir, 'lower', bc[1])
         if type(bc[1]) == "table" then self.auxBCvalues[dir][1] = bc[1][2] end
      else
         assert(isPeriodic, "Invalid lower boundary condition in non-periodic direction ".. dirNames[dir]..".")
      end

      if bc[2] then
         self:appendBoundaryConditions(dir, 'upper', bc[2])
         if type(bc[2]) == "table" then self.auxBCvalues[dir][2] = bc[2][2] end
      else
         assert(isPeriodic, "Invalid upper boundary condition in non-periodic direction ".. dirNames[dir]..".")
      end
   end

   local isPeriodic = {false, false, false}
   for _, dir in ipairs(self.grid:getPeriodicDirs()) do isPeriodic[dir] = true end

   -- Add various BCs to list of BCs to apply.
   local bc = {self.bcx, self.bcy, self.bcz}
   for d = 1, self.ndim do handleBc(d, bc[d], isPeriodic[d]) end
end

function FluidSpecies:createSolver(externalField)
   if externalField then
      -- Set up Jacobian.
      self.jacobFunc = externalField.jacobGeoFunc
      self.jacob     = externalField.geo.jacobGeo
      self.jacobInv  = externalField.geo.jacobGeoInv
   end

   -- Create solvers for collisions (diffusion).
   for _, c in pairs(self.collisions) do
      c:createSolver(externalField)
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

   -- Operators needed for time-dependent calculation and diagnostics.
   self.weakMultiply = Updater.CartFieldBinOp {
      onGrid    = self.grid,
      weakBasis = self.basis,
      operation = "Multiply",
      onGhosts  = true,
   }
   self.weakDivide = Updater.CartFieldBinOp {
      onGrid    = self.grid,
      weakBasis = self.basis,
      operation = "Divide",
      onGhosts  = true,
   }
   self.volIntegral = {
      comps1 = Updater.CartFieldIntegratedQuantCalc {
         onGrid        = self.grid,
         basis         = self.basis,
         numComponents = 1,
         quantity      = "V"
      },
      compsN = Updater.CartFieldIntegratedQuantCalc {
         onGrid        = self.grid,
         basis         = self.basis,
         numComponents = self.nMoments,
         quantity      = "V"
      }
   }
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
      metaData = {polyOrder = self.basis:polyOrder(),
                  basisType = self.basis:id(),
                  charge    = self.charge,
                  mass      = self.mass,
                  grid      = GKYL_OUT_PREFIX .. "_" .. self.name .. "_grid.bp",},
   }

   self.couplingMoments = self:allocVectorMoment(self.nMoments)

   self.unitField = self:allocMoment()
   local projectUnity = Updater.ProjectOnBasis {
      onGrid   = self.grid,
      basis    = self.basis,
      evaluate = function(t, zn) return 1.0 end,
      onGhosts = true,
   }
   projectUnity:advance(0.0, {}, {self.unitField})

   self.noJacMom = self:allocVectorMoment(self.nMoments)   -- Moments without Jacobian.

   if self.positivity then self.momPos = self:allocVectorMoment(self.nMoments) end

   -- Array with one component per cell to store cflRate in each cell.
   self.cflRateByCell = self:allocCartField(self.grid, 1, {1,1})
   self.cflRatePtr    = self.cflRateByCell:get(1)
   self.cflRateIdxr   = self.cflRateByCell:genIndexer()
   self.dtGlobal      = ffi.new("double[2]")

   self:createBCs()

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
function FluidSpecies:initDist(extField)
   if self.randomseed then
      math.randomseed(self.randomseed)
   else
      math.randomseed(47*Mpi.Comm_rank(Mpi.COMM_WORLD)+os.time())
   end

   local syncPeriodicDirs = true
   if self.fluctuationBCs then syncPeriodicDirs = false end

   local initCnt, backgroundCnt = 0, 0
   for nm, pr in pairs(self.projections) do
      pr:fullInit(self)
      pr:advance(0.0, {extField}, {self.moments[2]})
      -- This barrier is needed as when using MPI-SHM some
      -- processes will get to accumulate before projection is finished.
      Mpi.Barrier(self.grid:commSet().sharedComm)
      if string.find(nm,"init") then
         self.moments[1]:accumulate(1.0, self.moments[2])
         initCnt = initCnt + 1
      end
      if string.find(nm,"background") then
         if not self.momBackground then self.momBackground = self:allocVectorMoment(self.nMoments) end
         self.momBackground:accumulate(1.0, self.moments[2])
         self.momBackground:sync(syncPeriodicDirs)
         backgroundCnt = backgroundCnt + 1
      end
      if string.find(nm,"source") then
         if not self.mSource then self.mSource = self:allocVectorMoment(self.nMoments) end
         self.mSource:accumulate(1.0, self.moments[2])
         if self.positivityRescale then
            self.posRescaler:advance(0.0, {self.mSource}, {self.mSource})
         end
      end
   end
   assert(initCnt>0, string.format("FluidSpecies: Species '%s' not initialized!", self.name))
   if self.momBackground and backgroundCnt == 0 then self.momBackground:copy(self.moments[1]) end

   if self.fluctuationBCs then
      assert(backgroundCnt > 0, "FluidSpecies: must specify an initial background distribution with 'background' in order to use fluctuation-only BCs")
   end

   self.moments[2]:clear(0.0)

   self:setActiveRKidx(1)

   if self.positivityRescale or self.positivityDiffuse then
      self.posRescaler:advance(0.0, {self.moments[1]}, {self.moments[1]}, false)
   end
end

function FluidSpecies:setActiveRKidx(rkIdx)
   self.activeRKidx = rkIdx
end

function FluidSpecies:rkStepperFields() return self.moments end

function FluidSpecies:getMoments(rkIdx)
   if rkIdx == nil then
      return self:rkStepperFields()[self.activeRKidx]
   else
      return self:rkStepperFields()[rkIdx]
   end
end

function FluidSpecies:copyRk(outIdx, aIdx)
   self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx])
end

function FluidSpecies:combineRk(outIdx, a, aIdx, ...)
   local args = {...} -- Package up rest of args as table.
   local nFlds = #args/2
   self:rkStepperFields()[outIdx]:combine(a, self:rkStepperFields()[aIdx])
   for i = 1, nFlds do -- Accumulate rest of the fields.
      self:rkStepperFields()[outIdx]:accumulate(args[2*i-1], self:rkStepperFields()[args[2*i]])
   end
end

function FluidSpecies:suggestDt()
   if not self.evolve then return GKYL_MAX_DOUBLE end

   local dtSuggested = math.min(self.cfl/self.cflRateByCell:reduce('max')[1], GKYL_MAX_DOUBLE)

   -- If dtSuggested == GKYL_MAX_DOUBLE, it is likely because of NaNs.
   -- If so, return 0 so that no timestep is taken, and we will abort the simulation.
   if dtSuggested == GKYL_MAX_DOUBLE then dtSuggested = 0.0 end

   return dtSuggested
end

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
  self:applyBc(tCurr, self:rkStepperFields()[idx])
  if self.positivity then
     self:checkPositivity(tCurr, idx)
  end
end

function FluidSpecies:applyBc(tCurr, momIn)
   -- momIn is the set of evolved moments.
   local tmStart = Time.clock()

   if self.evolve then
      if self.fluctuationBCs then
         -- If fluctuation-only BCs, subtract off background before applying BCs.
         momIn:accumulate(-1.0, self.momBackground)
      end

      -- Apply non-periodic BCs (to only fluctuations if fluctuation BCs).
      if self.hasNonPeriodicBc then
         for _, bc in ipairs(self.boundaryConditions) do
            bc:advance(tCurr, {}, {momIn})
         end
      end

      -- Apply periodic BCs (to only fluctuations if fluctuation BCs)
      momIn:sync()

      if self.fluctuationBCs then
         -- Put back together total distribution
         momIn:accumulate(1.0, self.momBackground)

         -- Update ghosts in total distribution, without enforcing periodicity.
         momIn:sync(false)
      end
   end

   self.bcTime = self.bcTime + Time.clock()-tmStart
end

function FluidSpecies:createDiagnostics()  -- More sophisticated/extensive diagnostics go in Species/Diagnostics.
   -- Create this species' diagnostics.
   self.diagnostics[self.name] = FluidDiags{}
   self.diagnostics[self.name]:fullInit(self)

   -- Many diagnostics require dividing by the Jacobian (if present).
   -- Predefine the function that does that.
   self.calcNoJacMom = self.jacobInv
      and function(tm, rkIdx) self.weakMultiply:advance(tm, {self:getMoments(rkIdx), self.jacobInv}, {self.noJacMom}) end
      or function(tm, rkIdx) self.noJacMom:copy(self:getMoments(rkIdx)) end
end

function FluidSpecies:getNoJacMoments()
   return self.noJacMom
end

function FluidSpecies:write(tm, force)
   if self.evolve or force then

      for _, dOb in pairs(self.diagnostics) do
         dOb:resetState(tm)   -- Reset booleans indicating if diagnostic has been computed.
      end
      self.calcNoJacMom(tm, 1)  -- Many diagnostics do not include Jacobian factor.

      local tmStart = Time.clock()
      -- Compute integrated diagnostics.
      if self.calcIntQuantTrigger(tm) then
         for _, dOb in pairs(self.diagnostics) do
            dOb:calcIntegratedDiagnostics(tm, self)   -- Compute this species' integrated diagnostics.
         end
      end
      self.integratedMomentsTime = self.integratedMomentsTime + Time.clock() - tmStart
      
      -- Only write stuff if triggered.
      if self.diagIoTrigger(tm) or force then
         local momIn = self:rkStepperFields()[1]
	 self.momIo:write(momIn, string.format("%s_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame)
         if self.momBackground then
            if tm == 0.0 then
               self.momBackground:write(string.format("%s_background_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame, true)
            end
            momIn:accumulate(-1, self.momBackground)
            self.momIo:write(momIn, string.format("%s_fluctuation_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame)
            momIn:accumulate(1, self.momBackground)
         end
         if tm == 0.0 and self.mSource then
            self.momIo:write(self.mSource, string.format("%s_mSource_0.bp", self.name), tm, self.diagIoFrame)
         end

         for _, dOb in pairs(self.diagnostics) do
            dOb:calcFieldDiagnostics(tm, self)   -- Compute this species' field diagnostics.
         end

         -- Write this species' field and integrated diagnostics.
         for _, dOb in pairs(self.diagnostics) do
            dOb:write(tm, self.diagIoFrame)
         end

         if self.evolveCollisions then  -- Write collision's diagnostics.
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
      if self.diagIoFrame == 0 then
         local momIn = self:rkStepperFields()[1]
         self.momIo:write(momIn, string.format("%s_%d.bp", self.name, 0), tm, 0)
      end
      self.diagIoFrame = self.diagIoFrame+1
   end
end

function FluidSpecies:writeRestart(tm)
   local writeGhost = false
   if self.hasSheathBCs or self.fluctuationBCs then writeGhost = true end

   self.momIo:write(self.moments[1], string.format("%s_restart.bp", self.name), tm, self.diagIoFrame, writeGhost)

   -- Write this species' restart diagnostics.
   for _, dOb in pairs(self.diagnostics) do
      dOb:writeRestart(tm, self.diagIoFrame, self.dynVecRestartFrame)
   end

   self.dynVecRestartFrame = self.dynVecRestartFrame + 1
end

function FluidSpecies:readRestart()
   local readGhost = false
   if self.hasSheathBCs or self.fluctuationBCs then readGhost = true end

   local tm, fr = self.momIo:read(self.moments[1], string.format("%s_restart.bp", self.name), readGhost)
   self.diagIoFrame = fr -- Reset internal frame counter.

   -- Set ghost cells.
   self.moments[1]:sync()

   if not self.hasSheathBCs and not self.fluctuationBCs then
      self:applyBc(tm, self.moments[1])
   end

   -- Read this species field and integrated diagnostics.
   for _, dOb in pairs(self.diagnostics) do
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
function FluidSpecies:totalBcTime()
   return self.bcTime
end
function FluidSpecies:momCalcTime()
   return 0
end
function FluidSpecies:intMomCalcTime()
   return self.integratedMomentsTime
end

return FluidSpecies
