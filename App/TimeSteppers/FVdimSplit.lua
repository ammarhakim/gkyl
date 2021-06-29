-- Gkyl ------------------------------------------------------------------------
--
-- Function to advance solution using FV dimensionally split scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local lume = require "Lib.lume"
local Logger = require "Lib.Logger"

local FVdimSplit = Proto()

-- Store table passed to it and defer construction.
function FVdimSplit:init(tbl)
   self.tbl = tbl

   self.cflFrac = 1.0
   self.numFields = 3
   self.numStates = 7
end

function FVdimSplit:createSolver(appStatus, stepperFuncs, appsIn)

   self.copy         = stepperFuncs[2]

   self.species       = appsIn[1]
   self.field         = appsIn[2]
   self.externalField = appsIn[3]
   self.sources       = appsIn[4]

   self.stepperTime = 0

   -- Stepper status.
   self.stepStatus = {success   = true,  dt_actual    = 0.,
                      nextState = 0,     dt_suggested = 0.,}

   self.appStatus = appStatus

   -- tryInv contains indicators for each species
   -- whether the domain-invariant equation should be used in the next step;
   -- they might be changed during fvDimSplit calls.
   self.tryInv = {}
   for _, s in pairs(self.species) do self.tryInv[s] = false end

   -- Save the dimensionality for later. Fetch it from the species or the field.
   local specDim
   for _, s in pairs(self.species) do specDim = s.grid:ndim(); break end
   self.cdim = specDim or self.field.grid:ndim()

   self.log = Logger { logToFile = true }
end

function FVdimSplit:updateInDirection(dir, tCurr, dt, tryInvIn)
   local success, dt_suggested = true, GKYL_MAX_DOUBLE
   local fIdx = { {1,2}, {2,1}, {1,2} } -- For indexing inp/out fields.

   local tryInv_next = {}
   -- Update species.
   for _, s in lume.orderedIter(self.species) do
      local vars = s:rkStepperFields()
      local inp, out = vars[fIdx[dir][1]], vars[fIdx[dir][2]]
      local my_success, my_dt_suggested, myTryInv = s:updateInDirection(
         dir, tCurr, dt, inp, out, tryInvIn[s])
      tryInv_next[s] = myTryInv
      success = success and my_success
      dt_suggested = math.min(dt_suggested, my_dt_suggested)
   end
   do
      -- Update field.
      local vars = self.field:rkStepperFields()
      local inp, out = vars[fIdx[dir][1]], vars[fIdx[dir][2]]
      local my_success, my_dt_suggested = self.field:updateInDirection(dir, tCurr, dt, inp, out)
      success = success and my_success
      dt_suggested = math.min(dt_suggested, my_dt_suggested)
   end

   return success, dt_suggested, tryInv_next
end

function FVdimSplit:updateFluidSource(dataIdx, tCurr, dtIn)
   -- Update fluid sources.

   -- Make list of species data to operate on.
   local speciesVar = {}
   for nm, s in lume.orderedIter(self.species) do
      speciesVar[nm] = s:rkStepperFields()[dataIdx]
   end
   -- Field data to operate on.
   local fieldVar = self.field:rkStepperFields()[dataIdx]

   -- Expose freely-available array space space. Useful for storing
   -- intermediate quantities like spatial gradient.
   local bufIdx = (dataIdx==1) and 2 or 1
   local speciesBuf = {}
   for nm, s in lume.orderedIter(self.species) do
      speciesBuf[nm] = s:rkStepperFields()[bufIdx]
   end
   -- Field data to operate on.
   local fieldBuf = self.field:rkStepperFields()[bufIdx]

   local success, dt_suggested = true, GKYL_MAX_DOUBLE
   -- Update fluid source terms.
   for _, flSrc in lume.orderedIter(self.sources) do
      local my_success, my_dt_suggested = flSrc:updateFluidSource(
         tCurr, dtIn, speciesVar, fieldVar, speciesBuf, fieldBuf, self.species,
         self.field, self.externalField.em)
      success = success and my_success
      dt_suggested = math.min(dt_suggested, my_dt_suggested)
   end

   return success, dt_suggested
end

function FVdimSplit:advance(tCurr, dtIn)

   local dt = dtIn

   local success, dt_suggested = true, GKYL_MAX_DOUBLE
   local fIdx = { {1,2}, {2,1}, {1,2} } -- For indexing inp/out fields.

   -- Copy in case we need to take this step again.
   self.copy(3, 1)

   -- Update fluid source by half time-step.
   do
      local my_success, my_dt_suggested = self:updateFluidSource(1, tCurr, dt/2)
      success = success and my_success
      dt_suggested = math.min(dt_suggested, my_dt_suggested)
   end

   -- Update solution in each direction.
   local isInv = true
   for d = 1, self.cdim do
      local my_success, my_dt_suggested, myTryInv = self:updateInDirection(d, tCurr, dt, self.tryInv)
      success = success and my_success
      dt_suggested = math.min(dt_suggested, my_dt_suggested)
      if not success then
         self.log(" ** Time step too large! Aborting this step!")
         break
      else
         -- If an updated species is invalid, plan to use lax flux for THIS
         -- species in the re-taken step.
         for nm, s in lume.orderedIter(self.species) do
            if myTryInv[s] then
               isInv = false
               self.tryInv[s] = true
               self.log(string.format(
                           "\n ** Invalid values in %s; Will re-update using Lax flux!\n", nm))
            end
         end
         -- Break the loop if any species is invalid.
         if not isInv then break end
      end
   end
   -- Is all species is valid, do not use lax in the next step.
   if isInv then
      for _, s in lume.orderedIter(self.species) do self.tryInv[s] = false end
   end

   -- Update fluid source by half time-step.
   if success and isInv then
      local my_success, my_dt_suggested
      if fIdx[self.cdim][2] == 2 then
         my_success, my_dt_suggested = self:updateFluidSource(2, tCurr, dt/2)
      else
         my_success, my_dt_suggested = self:updateFluidSource(1, tCurr, dt/2)
      end
      success = success and my_success
      dt_suggested = math.min(dt_suggested, my_dt_suggested)
   end

   if not (success and isInv) then
      self.copy(1, 3) -- Restore old solution in case of failure.
   else
      -- If solution not already in field[1], copy for use in next
      -- time-step.
      if fIdx[self.cdim][2] == 2 then self.copy(1, 2) end
   end

   self.appStatus.success = success
   if not isInv then self.appStatus.success = false end

   self.stepStatus.dt_actual    = dt
   self.stepStatus.dt_suggested = dt_suggested

   return self.stepStatus
end

return FVdimSplit
