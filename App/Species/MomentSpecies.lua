-- Gkyl ------------------------------------------------------------------------
--
-- Species object constructed from moment equations.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid = require "Grid"
local DataStruct = require "DataStruct"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Lin = require "Lib.Linalg"
local LinearTrigger = require "Lib.LinearTrigger"
local Mpi = require "Comm.Mpi"
local DecompRegionCalc = require "Lib.CartDecomp"
local Proto = require "Lib.Proto"
local Projection = require "App.Projection"
local ProjectionBase = require "App.Projection.ProjectionBase"
local SpeciesBase = require "App.Species.SpeciesBase"
local Time = require "Lib.Time"
local Updater = require "Updater"
local xsys = require "xsys"
local ffi = require "ffi"
local Euler = require "Eq.Euler"
local TenMoment = require "Eq.TenMoment"

-- Species object treated as moment equations.
local MomentSpecies = Proto(SpeciesBase)

-- Add constants to object indicate various supported boundary conditions.
local SP_BC_OPEN = 1
local SP_BC_COPY = SP_BC_OPEN
local SP_BC_WALL = 4
-- The following two are not yet available for the MomentSpecies.
-- local SP_BC_DIRICHLET = 5    -- Specify the value (currently only for diffusion term).
-- local SP_BC_NEUMANN = 6    -- Specify the derivative (currently only for diffusion term).
local SP_BC_AXIS = 7

MomentSpecies.bcOpen = SP_BC_OPEN   -- Open BCs.
MomentSpecies.bcCopy = SP_BC_COPY   -- Copy BCs.
MomentSpecies.bcWall = SP_BC_WALL   -- Wall BCs.
MomentSpecies.bcAxis = SP_BC_AXIS   -- Axis BCs

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function MomentSpecies:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function MomentSpecies:fullInit(appTbl)
   local tbl = self.tbl -- Previously store table.

   self.cfl =  0.1
   self.charge = tbl.charge and tbl.charge or 1.0
   self.mass = tbl.mass and tbl.mass or 1.0
   self.ioMethod = "MPI"

   self.evolve = xsys.pickBool(tbl.evolve, true) -- Default: evolve species.
   self.forceWrite = xsys.pickBool(tbl.forceWrite, false) -- Default: don't write species if not evolved.

   -- Create triggers to write diagnostics.
   if tbl.nDiagnosticFrame then
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, tbl.nDiagnosticFrame)
   else
      self.diagIoTrigger = LinearTrigger(0, appTbl.tEnd, appTbl.nFrame)
   end

   self.diagIoFrame = 0 -- Frame number for diagnostics.
   self.dynVecRestartFrame = 0 -- Frame number of restarts (for DynVectors only)

   -- For storing integrated moments.
   self.integratedMoments = nil -- Allocated in alloc() method.

   -- Store initial condition function.
   self.initFunc = tbl.init

   self.equation = tbl.equation -- Equation system to evolve.
   self.nMoments = tbl.equation:numEquations()
   self.nGhost = 2     -- We need two ghost-cells.

   self.hasNonPeriodicBc = false -- To indicate if we have non-periodic BCs.
   self.bcx, self.bcy, self.bcz = { }, { }, { }

   -- Read in boundary conditions.
   -- Check to see if bc type is good is now done in createBc.
   if tbl.bcx then
      self.bcx[1], self.bcx[2] = tbl.bcx[1], tbl.bcx[2]
      self.hasNonPeriodicBc = true
   end
   if tbl.bcy then
      self.bcy[1], self.bcy[2] = tbl.bcy[1], tbl.bcy[2]
      self.hasNonPeriodicBc = true
   end
   if tbl.bcz then
      self.bcz[1], self.bcz[2] = tbl.bcz[1], tbl.bcz[2]
      self.hasNonPeriodicBc = true
   end

   self.ssBc = {}
   if tbl.ssBc then
      self.ssBc[1] = tbl.ssBc[1]
   end

   self.boundaryConditions = { } -- List of Bcs to apply.
   self.ssBoundaryConditions = { } -- List of stair-stepped Bcs to apply.

   self.bcTime = 0.0 -- Timer for BCs.

   -- Initialization.
   self.projections = {}
   for nm, val in pairs(tbl) do
      if ProjectionBase.is(val) then
         self.projections[nm] = val
      end
   end
   -- It is possible to use the keyword 'source' to specify a
   -- function directly without using a Projection object.
   if type(tbl.init) == "function" then
      self.projections["init"] = Projection.MomentProjection.FunctionProjection {
         func = function (t, zn)
            return tbl.init(t, zn, self)
         end,
      }
   elseif type(tbl.init) == "string" then
      -- Specify the suffix of the file with the initial condition (including the extension).
      -- The prefix is assumed to be the name of the input file.
      self.projections["init"] = Projection.MomentProjection.ReadInput {
         inputFile = tbl.init,
      }
   end
   if type(tbl.source) == "function" then
      self.projections["source"] = Projection.MomentProjection.FunctionProjection {
         func = function (t, zn)
            return tbl.source(t, zn, self)
         end,
      }
   elseif type(tbl.source) == "string" then
      -- Specify the suffix of the file with the source (including the extension).
      -- The prefix is assumed to be the name of the input file.
      self.projections["source"] = Projection.MomentProjection.ReadInput {
         inputFile = tbl.source,
      }
   end

   self.useShared = xsys.pickBool(appTbl.useShared, false)

   self.tCurr = 0.0

   self.limiter = tbl.limiter and tbl.limiter or "monotonized-centered"
   self.hyperSlvr = {} -- List of solvers.

   -- Invariant (positivity-preserving) equation system to evolve.
   self.equationInv = tbl.equationInv
   self.hyperSlvrInv = {} -- List of solvers.
   self.limiterInv = tbl.limiterInv and tbl.limiterInv or "zero"
   -- Always use invariant eqn.
   self.forceInv = tbl.forceInv and (self.equationInv ~= nil)
   -- Use invariant eqn. in next step; could change during run.
   self.tryInv = false

   self._myIsInv = ffi.new("int[2]")
   self._isInv = ffi.new("int[2]")

   self._hasSsBnd = xsys.pickBool(tbl.hasSsBnd, false)
   self._inOutFunc = tbl.inOutFunc

   self.basis = nil -- Will be set later
end

function MomentSpecies:getCharge() return self.charge end
function MomentSpecies:getMass() return self.mass end
function MomentSpecies:getEvolve() return self.evolve end

function MomentSpecies:setName(nm)
   self.name = nm
end
function MomentSpecies:setCfl(cfl)
   self.cfl = cfl
end
function MomentSpecies:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end
function MomentSpecies:setConfBasis(basis)
   self.basis = basis
end
function MomentSpecies:setConfGrid(cgrid)
   self.grid = cgrid
   self.ndim = self.grid:ndim()
end

function MomentSpecies:createGrid(cgrid)
   self.cdim = cgrid:ndim()
   self.ndim = self.cdim

   -- Create decomposition.
   local decompCuts = {}
   for d = 1, self.cdim do table.insert(decompCuts, cgrid:cuts(d)) end
   self.decomp = DecompRegionCalc.CartProd {
      cuts = decompCuts,
      useShared = self.useShared,
   }

   -- Create computational domain.
   local lower, upper, cells = {}, {}, {}
   for d = 1, self.cdim do
      table.insert(lower, cgrid:lower(d))
      table.insert(upper, cgrid:upper(d))
      table.insert(cells, cgrid:numCells(d))
   end
   self.grid = Grid.RectCart {
      lower = lower,
      upper = upper,
      cells = cells,
      periodicDirs = cgrid:getPeriodicDirs(),
      decomposition = self.decomp,
   }
end

function MomentSpecies:allocMoment()
   local m = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis(),
      ghost = {self.nGhost, self.nGhost}
   }
   return m
end
function MomentSpecies:allocVectorMoment(numComp)
   local m = DataStruct.Field {
      onGrid = self.grid,
      numComponents = self.basis:numBasis()*numComp,
      ghost = {self.nGhost, self.nGhost}
   }
   return m
end

function MomentSpecies:allocMomCouplingFields()
   return {self:allocVectorMoment(self.nMoments)}
end

function MomentSpecies:bcCopyFunc(dir, tm, idxIn, fIn, fOut)
   for i = 1, self.nMoments*self.basis:numBasis() do
      fOut[i] = fIn[i]
   end
end

-- Function to construct a BC updater.
function MomentSpecies:makeBcUpdater(dir, edge, bcList, skinLoop, hasExtFld)

   return Updater.Bc {
      onGrid = self.grid,
      boundaryConditions = bcList,
      dir = dir,
      edge = edge,
      skinLoop = skinLoop,
      cdim = self.cdim,
      vdim = self.vdim,
      hasExtFld = hasExtFld,
   }
end

-- Function to construct a stair-stepped BC updater.
function MomentSpecies:makeSsBcUpdater(dir, inOut, bcList)
   return Updater.StairSteppedBc {
      onGrid = self.grid,
      inOut = inOut,
      boundaryConditions = bcList,
      dir = dir,
   }
end

function MomentSpecies:createBCs()
   -- Create a table to store auxiliary values needed by BCs
   -- and provided by the user in the input file.
   self.auxBCvalues = {}

   -- Functions to make life easier while reading in BCs to apply.
   -- Note: appendBoundaryConditions defined in sub-classes.
   local function handleBc(dir, bc, isPeriodic)
      table.insert(self.auxBCvalues,{nil,nil})

      local dirNames = {"X", "Y", "Z"}
      if (isPeriodic) then
         assert(bc==nil or (bc[1]==nil and bc[2]==nil),
                "Boundary conditions supplied in periodic direction "..
                dirNames[dir]..".")
      end

      if bc[1] then
         self:appendBoundaryConditions(dir, 'lower', bc[1])
         if type(bc[1]) == "table" then
            self.auxBCvalues[dir][1] = bc[1][2]
         end
      else
         assert(isPeriodic,
                "Invalid lower boundary condition in non-periodic direction "..
                dirNames[dir]..".")
      end

      if bc[2] then
         self:appendBoundaryConditions(dir, 'upper', bc[2])
         if type(bc[2]) == "table" then
            self.auxBCvalues[dir][2] = bc[2][2]
         end
      else
         assert(isPeriodic,
                "Invalid upper boundary condition in non-periodic direction "..
                dirNames[dir]..".")
      end
   end

   -- Add various BCs to list of BCs to apply.
   local isPeriodic = {false, false, false}
   for _,dir in ipairs(self.grid:getPeriodicDirs()) do
      isPeriodic[dir] = true
   end
   local bc = {self.bcx, self.bcy, self.bcz}
   for d = 1, self.cdim do
     handleBc(d, bc[d], isPeriodic[d])
  end
end

function MomentSpecies:createSolver(hasE, hasB)
   if self._hasSsBnd then
      self._inOut = DataStruct.Field {
         onGrid = self.grid,
         numComponents = 1,
         ghost = {2, 2}
      }
      local project = Updater.ProjectOnBasis {
         onGrid = self.grid,
         basis = self.basis,
         evaluate = self._inOutFunc,
         onGhosts = true,
      }
      project:advance(0.0, {}, {self._inOut})
      self.momIo:write(self._inOut, string.format("%s_inOut.bp", self.name), 0, 0)
   end

   local ndim = self.grid:ndim()
   for d = 1, ndim do
      self.hyperSlvr[d] = Updater.WavePropagation {
         onGrid = self.grid,
         equation = self.equation,
         limiter = self.limiter,
         cfl = self.cfl,
         updateDirections = {d},
         hasSsBnd = self._hasSsBnd,
         inOut = self._inOut,
      }
      if self.equationInv ~= nil then
         self.hyperSlvrInv[d] = Updater.WavePropagation {
            onGrid = self.grid,
            equation = self.equationInv,
            limiter = self.limiterInv,
            cfl = self.cfl,
            updateDirections = {d},
            hasSsBnd = self._hasSsBnd,
            inOut = self._inOut
         }
      end
   end

   -- This perhaps should go to createBCs but at that point _inOut is not
   -- created yet.
   if (self._hasSsBnd) then
      local function handleSsBc(dir, bcList)
         for _, bc in ipairs(bcList) do
            if bc then
               self:appendSsBoundaryConditions(dir, self._inOut, bc)
            end
         end
      end

      for d = 1, ndim do
         handleSsBc(d, self.ssBc)
      end
   end

end

function MomentSpecies:alloc(nRkDup)
   -- Allocate fields needed in RK update.
   self.moments = {}
   for i = 1, nRkDup do
      self.moments[i] = self:allocVectorMoment(self.nMoments)
      self.moments[i]:clear(0.0)
   end
   -- Create Adios object for field I/O.
   self.momIo = AdiosCartFieldIo {
      elemType = self.moments[1]:elemType(),
      method = self.ioMethod,
      metaData = {
         polyOrder = self.basis:polyOrder(),
         basisType = self.basis:id(),
         charge = self.charge,
         mass = self.mass,
      },
   }
   self.couplingMoments = self:allocVectorMoment(self.nMoments)
   self.integratedMoments = DataStruct.DynVector { numComponents = self.nMoments }

   -- Array with one component per cell to store cflRate in each cell.
   self.cflRateByCell = DataStruct.Field {
      onGrid = self.grid,
      numComponents = 1,
      ghost = {1, 1},
   }
   self.cflRateByCell:clear(0.0)
   self.cflRatePtr = self.cflRateByCell:get(1)
   self.cflRateIdxr = self.cflRateByCell:genIndexer()
   self.dtGlobal = ffi.new("double[2]")

   self:createBCs()
end

function MomentSpecies:initDist(extField)

   local initCnt = 0
   for nm, pr in pairs(self.projections) do
      pr:fullInit(self)
      pr:advance(0.0, {}, {self.moments[2]})
      -- This barrier is needed as when using MPI-SHM some
      -- processes will get to accumulate before projection is finished.
      Mpi.Barrier(self.grid:commSet().sharedComm)
      if nm == "init" then
         self.moments[1]:accumulate(1.0, self.moments[2])
         initCnt = initCnt + 1
      end
      if nm == "source" then
         if not self.mSource then
            self.mSource = self:allocVectorMoment(self.nMoments)
         end
         self.mSource:accumulate(1.0, self.moments[2])
      end
   end
   assert(initCnt > 0,
          string.format("MomentSpecies: Species '%s' not initialized!", self.name))
   self.moments[2]:clear(0.0)

end

function MomentSpecies:rkStepperFields()
   return self.moments
end

-- For timestepping.
function MomentSpecies:copyRk(outIdx, aIdx)
   self:rkStepperFields()[outIdx]:copy(self:rkStepperFields()[aIdx])
end

function MomentSpecies:setDtGlobal(dtGlobal)
   self.dtGlobal[0] = dtGlobal
end

function MomentSpecies:appendBoundaryConditions(dir, edge, bcType)
   local function bcCopyFunc(...) return self:bcCopyFunc(...) end

   if bcType == SP_BC_COPY then
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, edge, { bcCopyFunc }))
   elseif bcType == SP_BC_WALL then
      -- FIXME better to define and use self.equation.bcWall.
      local bcWall
      if self.nMoments == 5 then
         bcWall = Euler.bcWall
      elseif self.nMoments == 10 then
         bcWall = TenMoment.bcWall
      else
         assert(false, "MomentSpecies: bcWall not provided by the equation!")
      end
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, bcWall))
   elseif bcType == SP_BC_AXIS then
      -- FIXME: see comments above
      local bcAxis
      if self.nMoments == 5 then
         bcAxis = Euler.bcAxis
      else
         assert(false, "MomentSpecies: bcAxis not provided by the equation!")
      end
      table.insert(self.boundaryConditions, self:makeBcUpdater(dir, edge, bcAxis))
   elseif type(bcType) == "table" then
      -- bcType can be literally a list of functions.
      table.insert(self.boundaryConditions,
		   self:makeBcUpdater(dir, edge, bcType ))
   else
      assert(false, "MomentSpecies: Unsupported BC type!")
   end
end

-- TODO: merge into appendBoundaryConditions.
function MomentSpecies:appendSsBoundaryConditions(dir, inOut, bcType)
   local function bcCopyFunc(...) return self:bcCopyFunc(...) end

   if bcType == SP_BC_COPY then
      table.insert(self.ssBoundaryConditions,
		   self:makeSsBcUpdater(dir, inOut, { bcCopyFunc }))
   elseif bcType == SP_BC_WALL then
      -- FIXME better to define and use self.equation.bcWall.
      local bcWall
      if self.nMoments == 5 then
        bcWall = Euler.bcWall
      elseif self.nMoments == 10 then
        bcWall = TenMoment.bcWall
      else
        assert(false, "MomentSpecies: bcWall not provided by the equation!")
      end
      table.insert(self.ssBoundaryConditions,
		   self:makeSsBcUpdater(dir, inOut, bcWall))
   elseif type(bcType) == "table" then
      -- bcType can be literally a list of functions.
      table.insert(self.ssBoundaryConditions,
		   self:makeSsBcUpdater(dir, inOut, bcType ))
   else
      assert(false, "MomentSpecies: Unsupported BC type!")
   end
end

function MomentSpecies:updateInDirection(dir, tCurr, dt, fIn, fOut, tryInv)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   local tryInv_next = false
   if self.evolve then
      self:applyBc(tCurr, fIn, dir)
      assert(self:checkInv(fIn))
      if self.forceInv or tryInv then
         self.hyperSlvrInv[dir]:setDtAndCflRate(dt, nil)
         status, dtSuggested = self.hyperSlvrInv[dir]:advance(tCurr, {fIn}, {fOut})
         -- If result is OK, do not try to use invariant eqn. in next step.
         tryInv_next = not status
         if status and (not self:checkInv(fOut)) then
            assert(false, "** Invalid output using Lax flux!")
         end
      else
         self.hyperSlvr[dir]:setDtAndCflRate(dt, nil)
         status, dtSuggested = self.hyperSlvr[dir]:advance(tCurr, {fIn}, {fOut})
         tryInv_next = status and not self:checkInv(fOut)
      end
   else
      fOut:copy(fIn)
   end
   return status, dtSuggested, tryInv_next
end

function MomentSpecies:applyBcIdx(tCurr, idx, isFirstRk)
  for dir = 1, self.ndim do
     self:applyBc(tCurr, self:rkStepperFields()[idx], dir)
  end
end

function MomentSpecies:applyBc(tCurr, fIn, dir)
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

function MomentSpecies:createDiagnostics()
   -- Create updater to compute volume-integrated moments.
   self.intMom2Calc = Updater.CartFieldIntegratedQuantCalc {
      onGrid = self.grid,
      basis = self.basis,
      numComponents = self.nMoments,
      quantity = "V"
   }
end

function MomentSpecies:write(tm, force)
   if self.evolve or self.forceWrite then
      -- Compute integrated diagnostics.
      self.intMom2Calc:advance(tm, { self.moments[1] }, { self.integratedMoments })

      -- Only write stuff if triggered.
      if self.diagIoTrigger(tm) or force then
         self.momIo:write(
            self.moments[1], string.format("%s_%d.bp", self.name, self.diagIoFrame), tm, self.diagIoFrame)
         self.integratedMoments:write(
            string.format("%s_intMom.bp", self.name), tm, self.diagIoFrame)

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

function MomentSpecies:writeRestart(tm)
   self.momIo:write(
      self.moments[1], string.format("%s_restart.bp", self.name), tm, self.diagIoFrame)

   -- Write restart files for integrated moments. Note: these are only needed for the rare case that the
   -- restart write frequency is higher than the normal write frequency from nFrame.
   -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
   self.integratedMoments:write(
      string.format("%s_intMom_restart.bp", self.name), tm, self.dynVecRestartFrame, false, false)
   self.dynVecRestartFrame = self.dynVecRestartFrame + 1
end

function MomentSpecies:readRestart()
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

function MomentSpecies:totalSolverTime()
   local tm = 0.0
   for d = 1, self.grid:ndim() do
      tm = tm + self.hyperSlvr[d].totalTime
   end
   return tm
end
function MomentSpecies:totalBcTime()
   return self.bcTime
end
function MomentSpecies:momCalcTime()
   return 0
end
function MomentSpecies:intMomCalcTime()
   return self.intMom2Calc.totalTime
end

function MomentSpecies:checkInv(fIn)
   local fInIndexer = fIn:genIndexer()
   local fInPtr = fIn:get(1)

   local isInv = true
   local localRange = fIn:localRange()   
   for idx in localRange:rowMajorIter() do
      self.grid:setIndex(idx)

      fIn:fill(fInIndexer(idx), fInPtr)
      if not self.equation:isPositive(fInPtr) then
        isInv = false
        break
      end
  end

  self._myIsInv[0] = isInv and 1 or 0
  Mpi.Allreduce(self._myIsInv, self._isInv, 1, Mpi.INT, Mpi.LAND, self.grid:commSet().comm)
  return self._isInv[0] == 1 and true or false
end

return MomentSpecies
