-- Gkyl ------------------------------------------------------------------------
--
-- Diagnostics for FluidSpecies.
-- 
-- We assume that integrated diagnostics can depend on field diagnostics, but
-- not the other way around.
-- 
-- Supported diagnostic are defined as functions at the bottom of the file.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto      = require "Lib.Proto"
local FluidDiags = require "App.Species.Diagnostics.FluidDiagnostics"
local DataStruct = require "DataStruct"
local Time       = require "Lib.Time"
local Updater    = require "Updater"
local lume       = require "Lib.lume"

-- The first entry is the diagnostic name. The second is the other diagnostics it depends on.
local allowedFieldDiags = {{"MomSq",{}},{"M0",{}},{"M1",{}},{"M2",{}}, {"M2flow",{"M1","upar"}},
                           {"upar",{}}, {"Tpar",{}}, {"Tperp",{}}, {"Temp",{}}, {"ppar",{}}, {"pperp",{}}}
local allowedIntegratedDiags = {{"intMom",{}},{"intM0",{"M0"}}, {"intM1",{"M1"}}, {"intM2",{"M2"}}, {"intM2flow",{"M2flow"}}}

local function contains(tbl, elem) return lume.any(tbl, function(e) return e[1]==elem end) end
local function isFieldDiagAllowed(nm) return contains(allowedFieldDiags, nm) end
local function isIntegratedDiagAllowed(nm) return contains(allowedIntegratedDiags, nm) end

local GyrofluidDiags = Proto(FluidDiags)  -- GyrofluidDiags is a child of FluidDiagnostics.

function GyrofluidDiags:init(mySpecies)

   self.name     = mySpecies.name
   self.grid     = mySpecies.grid
   self.basis    = mySpecies.basis
   self.nMoments = mySpecies.nMoments

   self.weakMultiply = mySpecies.weakMultiply
   self.weakDivide   = mySpecies.weakDivide
   
   self.fieldDiags, self.intDiags = {}, {}  -- Tables to store field and integrated diagnostics.

   -- Sort requested diagnostics into field and integrated diagnostics.
   for _, nm in ipairs(mySpecies.diagnostics) do
      if isFieldDiagAllowed(nm) then
         self.fieldDiags[nm] = GyrofluidDiags._diagImp[nm]
      elseif isIntegratedDiagAllowed(nm) then
         self.intDiags[nm] = GyrofluidDiags._diagImp[nm]
      else
         assert(false, string.format("FluidDiagnostics: %s is not an allowed diagnostic.",nm))
      end
   end

   local function orgDiagnostics(allowedDiags, diagsTbl, otherAllowedDiags, otherDiagsTbl) 
      -- Organize diagnostics. Insert other diagnostics needed. Set metatable for orderedIter.
      assert((otherAllowedDiags and otherDiagsTbl) or (otherAllowedDiags==nil and otherDiagsTbl==nil),
             "GyrofluidDiags:orgDiagnostics: must pass both otherAllowedDiags and otherDiagsTbl, or neither.")
      local reorderOther = false

      -- Check for dependencies and insert other necessary diagnostics.
      -- Note: this can result in extra diagnostics written out.
      for nm, diag in pairs(diagsTbl) do
         local info, _ = lume.match(allowedDiags, function(e) return e[1]==nm end)
         local depends = info and info[2] or {}
         for depI, depNm in ipairs(depends) do
            if lume.find(diagsTbl, depNm)==nil then
               diagsTbl[depNm] = GyrofluidDiags._diagImp[depNm]
            end
            -- Also check the other diagsTbl if present.
            if otherDiagsTbl then
               if lume.find(otherDiagsTbl, depNm)==nil then
                  reorderOther = true
                  otherDiagsTbl[depNm] = GyrofluidDiags._diagImp[depNm]
               end
            end
         end
      end
      local function sortDiagsTbl(allowed, diags)
         -- Create a keys list so we always loop in the same order.
         local diags_keys = {}
         for nm in pairs(diags) do table.insert(diags_keys, nm) end
         local diagSort = function(d1,d2)
            -- Check if d1 is in the list of dependencies of d2. If so it must go before d2.
            local info, _  = lume.match(allowed, function(e) return e[1]==d2 end)
            local depends2 = info[2]
            if lume.find(depends2, d1) then return true else return false end 
         end
         table.sort(diags_keys, diagSort)
         setmetatable(diags, diags_keys)
      end
      sortDiagsTbl(allowedDiags, diagsTbl)
      -- Also sort other diagsTbl if something was added.
      if otherDiagsTbl and reorderOther then sortDiagsTbl(otherAllowedDiags, otherDiagsTbl) end
   end
   orgDiagnostics(allowedFieldDiags, self.fieldDiags) 
   orgDiagnostics(allowedIntegratedDiags, self.intDiags, allowedFieldDiags, self.fieldDiags) 

   -- Initialize diagnostics.
   for _, diag in lume.orderedIter(self.fieldDiags) do diag.init(self, mySpecies) end
   for _, diag in lume.orderedIter(self.intDiags) do diag.init(self, mySpecies) end

   self.spentTime = { fieldDiags=0., intDiags=0. }
end

function GyrofluidDiags:calcIntegratedDiagnostics(tm, mySpecies)
   -- Compute integrated diagnostics.
   local tmStart = Time.clock()

   local spec = mySpecies
   for _, diag in lume.orderedIter(self.intDiags) do diag.calc(tm, spec) end

   self.spentTime.intDiags = self.spentTime.intDiags + Time.clock() - tmStart
end

function GyrofluidDiags:calcFieldDiagnostics(tm, mySpecies)
   -- Compute field diagnostics.
   local tmStart = Time.clock()

   local spec = mySpecies
   for _, diag in lume.orderedIter(self.fieldDiags) do diag.calc(tm, spec) end

   self.spentTime.fieldDiags = self.spentTime.fieldDiags + Time.clock() - tmStart
end

function GyrofluidDiags:writeFieldDiagnostics(tm, fr)
   for nm, diag in lume.orderedIter(self.fieldDiags) do
      diag.field:write(string.format("%s_%s_%d.bp", self.name, nm, fr), tm, fr)
   end
end

function GyrofluidDiags:writeIntegratedDiagnostics(tm, fr)
   for nm, diag in lume.orderedIter(self.intDiags) do
      diag.field:write(string.format("%s_%s.bp", self.name, nm), tm, fr)
   end
end

function GyrofluidDiags:write(tm, fr)
   self:writeFieldDiagnostics(tm, fr)
   self:writeIntegratedDiagnostics(tm, fr)
end

function GyrofluidDiags:writeRestartFieldDiagnostics(tm, fr)
   for nm, diag in lume.orderedIter(self.fieldDiags) do
      diag.field:write(string.format("%s_%s_restart.bp", self.name, nm), tm, fr, false)
   end
end

function GyrofluidDiags:writeRestartIntegratedDiagnostics(tm, fr)
   for nm, diag in lume.orderedIter(self.intDiags) do
      diag.field:write(string.format("%s_%s_restart.bp", self.name, nm), tm, fr, false, false)
   end
end

function GyrofluidDiags:writeRestart(tm, fieldFr, intFr)
   self:writeRestartFieldDiagnostics(tm, fieldFr)

   -- Write restart files for integrated moments. Note: these are only needed for the rare case that the
   -- restart write frequency is higher than the normal write frequency from nFrame.
   -- (the first "false" prevents flushing of data after write, the second "false" prevents appending)
   self:writeRestartIntegratedDiagnostics(tm, intFr)
end

function GyrofluidDiags:readRestart()
   local tm, fr
   for nm, diag in lume.orderedIter(self.fieldDiags) do
      tm, fr = diag.field:read(string.format("%s_%s_restart.bp", self.name, nm))
   end

   for nm, diag in lume.orderedIter(self.intDiags) do
      tm, fr = diag.field:read(string.format("%s_%s_restart.bp", self.name, nm))
   end
   return tm, fr
end

-----------------------------------
-- Place diagnostics below
--
GyrofluidDiags._diagImp = {}

-- ~~~~ The moments squared ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["MomSq"] = {}
GyrofluidDiags._diagImp["MomSq"].init = function(self, specIn)
   GyrofluidDiags._diagImp["MomSq"].field = specIn:allocMoment()
end
GyrofluidDiags._diagImp["MomSq"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   specIn.weakMultiply:advance(tm, {momIn, momIn}, {GyrofluidDiags._diagImp["MomSq"].field})
end

-- ~~~~ Number density ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["M0"] = {}
GyrofluidDiags._diagImp["M0"].init = function(self, specIn)
   GyrofluidDiags._diagImp["M0"].field = specIn:allocMoment()
end
GyrofluidDiags._diagImp["M0"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   GyrofluidDiags._diagImp["M0"].field:combineOffset(1./specIn:getMass(), momIn, specIn:getMomOff(1))
end

-- ~~~~ Momentum density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["M1"] = {}
GyrofluidDiags._diagImp["M1"].init = function(self, specIn)
   GyrofluidDiags._diagImp["M1"].field = specIn:allocMoment()
end
GyrofluidDiags._diagImp["M1"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   GyrofluidDiags._diagImp["M1"].field:combineOffset(1./specIn:getMass(), momIn, specIn:getMomOff(2))
end

-- ~~~~ Energy density (divided by mass/2) ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["M2"] = {}
GyrofluidDiags._diagImp["M2"].init = function(self, specIn)
   GyrofluidDiags._diagImp["M2"].field = specIn:allocMoment()
end
GyrofluidDiags._diagImp["M2"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   GyrofluidDiags._diagImp["M2"].field:combineOffset(2./specIn:getMass(), momIn, specIn:getMomOff(3))
end

-- ~~~~ Flow energy density (divided by mass) ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["M2flow"] = {}
GyrofluidDiags._diagImp["M2flow"].init = function(self, specIn)
   GyrofluidDiags._diagImp["M2flow"].field = specIn:allocMoment()
end
GyrofluidDiags._diagImp["M2flow"].calc = function(tm, specIn)
   specIn.weakMultiply:advance(tm, {GyrofluidDiags._diagImp["upar"].field,GyrofluidDiags._diagImp["M1"].field},
                               {GyrofluidDiags._diagImp["M2flow"].field})
end

-- ~~~~ Parallel flow speed ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["upar"] = {}
GyrofluidDiags._diagImp["upar"].init = function(self, specIn)
   GyrofluidDiags._diagImp["upar"].field = specIn:allocMoment()
   GyrofluidDiags._diagImp["upar"].fieldAux = {mJacM0=specIn:allocMoment(),mJacM1=specIn:allocMoment()}
end
GyrofluidDiags._diagImp["upar"].calc = function(tm, specIn)
   local momIn = specIn:rkStepperFields()[1]
   specIn:uParCalc(tm, momIn, GyrofluidDiags._diagImp["upar"].fieldAux.mJacM0,
                              GyrofluidDiags._diagImp["upar"].fieldAux.mJacM1,
                              GyrofluidDiags._diagImp["upar"].field)
end

-- ~~~~ Moments integrated over the domain ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["intMom"] = {}
GyrofluidDiags._diagImp["intMom"].init = function(self, specIn)
   GyrofluidDiags._diagImp["intMom"].field    = DataStruct.DynVector { numComponents = self.nMoments }
   GyrofluidDiags._diagImp["intMom"].updaters = Updater.CartFieldIntegratedQuantCalc {
                                                   onGrid        = self.grid,
                                                   basis         = self.basis,
                                                   numComponents = self.nMoments,
                                                   quantity      = "V"
                                                }
end
GyrofluidDiags._diagImp["intMom"].calc = function(tm, specIn)
   GyrofluidDiags._diagImp["intMom"].updaters:advance(tm, {specIn:rkStepperFields()[1]}, {GyrofluidDiags._diagImp["intMom"].field})
end

-- ~~~~ Integrated number density ~~~~~~~~~~~~~~~~~~~~~~
GyrofluidDiags._diagImp["intM0"] = {}
GyrofluidDiags._diagImp["intM0"].init = function(self, specIn)
   GyrofluidDiags._diagImp["intM0"].field    = DataStruct.DynVector { numComponents = 1 }
   GyrofluidDiags._diagImp["intM0"].updaters = Updater.CartFieldIntegratedQuantCalc {
                                                   onGrid        = self.grid,
                                                   basis         = self.basis,
                                                   numComponents = self.nMoments,
                                                   quantity      = "V"
                                                }
end
GyrofluidDiags._diagImp["intM0"].calc = function(tm, specIn)
   GyrofluidDiags._diagImp["intM0"].updaters:advance(tm, {GyrofluidDiags._diagImp["M0"].field}, {GyrofluidDiags._diagImp["intM0"].field})
end

return GyrofluidDiags
