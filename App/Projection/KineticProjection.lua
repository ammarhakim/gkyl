-- Gkyl ------------------------------------------------------------------------
--
-- App support code: KineticProjection object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ProjectionBase = require "App.Projection.ProjectionBase"
local Proto = require "Lib.Proto"
local Updater = require "Updater"
local xsys = require "xsys"
--local Time = require "Lib.Time"

local KineticProjection = Proto(ProjectionBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below
function KineticProjection:init(tbl)
   self.tbl = tbl
end

function KineticProjection:fullInit(species)
   self.species = species

   self.phaseBasis = species.basis
   self.phaseGrid = species.grid
   self.confBasis = species.confBasis
   self.confGrid = species.confGrid

   self.numVelDims = self.phaseGrid:ndim() - self.confGrid:ndim()

   self.isInit = xsys.pickBool(self.tbl.isInit, true)
   self.isBackground = xsys.pickBool(self.tbl.isBackground, false)
   self.isSource = xsys.pickBool(self.tbl.isSource, false)
   if self.isBackground or self.isSource then self.isInit = false end

   self.exactScaleM0 = xsys.pickBool(self.tbl.exactScaleM0, true)
   self.exactLagFixM012 = xsys.pickBool(self.tbl.exactLagFixM012, false)
end

function KineticProjection:run(t, distf)
   self.project:advance(t, 0.0, {}, {distf})
end


----------------------------------------------------------------------
local FunctionProjection = Proto(KineticProjection)

function FunctionProjection:fullInit(species)
   FunctionProjection.super.fullInit(self, species)

   local func = self.tbl.func
   if not func then
      func = self.tbl[1]
   end
   assert(func, "FunctionProjection: Must specify the function")
   assert(type(func) == "function",
	  "The input must be a table containing function")
   self.project = Updater.ProjectOnBasis {
      onGrid = self.phaseGrid,
      basis = self.phaseBasis,
      evaluate = func,
      projectOnGhosts = true
   }
end


----------------------------------------------------------------------
local MaxwellianProjection = Proto(KineticProjection)

function MaxwellianProjection:fullInit(species)
   MaxwellianProjection.super.fullInit(self, species)

   local tbl = self.tbl
   self.density = assert(tbl.density, "Maxwellian: must specify 'density'")
   self.driftSpeed = tbl.driftSpeed or function (t, zn) return nil end
   self.temperature = assert(tbl.temperature,
			     "Maxwellian: must specify 'temperature'")

   -- check for constants instead of functions
   if type(self.density) ~= "function" then
      self.density = function (t, zn) return tbl.density end
   end
   if type(self.driftSpeed) ~= "function" then
      self.driftSpeed = function (t, zn) return tbl.driftSpeed end
   end
   if type(self.temperature) ~= "function" then
      self.temperature = function (t, zn) return tbl.temperature end
   end

   local func = function (t, zn)
      return species:Maxwellian(zn, self.density(t, zn, species),
				self.temperature(t, zn, species),
				self.driftSpeed(t, zn, species))
   end

   self.project = Updater.ProjectOnBasis {
      onGrid = self.phaseGrid,
      basis = self.phaseBasis,
      evaluate = func,
      projectOnGhosts = true
   }
end

function MaxwellianProjection:scaleDensity(distf)
   local M0e, M0 = self.species:allocMoment(), self.species:allocMoment()
   local M0mod = self.species:allocMoment()

   self.species.numDensityCalc:advance(0.0, 0.0, {distf}, {M0})
   local func = function (t, zn)
      return self.density(t, zn, self.species)
   end
   local project = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = func,
      projectOnGhosts = true,
   }
   project:advance(0.0, 0.0, {}, {M0e})

   local weakDivision = Updater.CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
      onGhosts = true,
   }
   local weakMultiplication = Updater.CartFieldBinOp {
      onGrid = self.phaseGrid,
      weakBasis = self.phaseBasis,
      fieldBasis = self.confBasis,
      operation = "Multiply",
      onGhosts = true,
   }

   -- calculate M0mod = M0e / M0
   weakDivision:advance(0.0, 0.0, {M0, M0e}, {M0mod})
   -- calculate distff = M0mod * distf
   weakMultiplication:advance(0.0, 0.0, {M0mod, distf}, {distf})
end

function MaxwellianProjection:run(t, distf)
   self.project:advance(t, 0.0, {}, {distf})
   if self.exactScaleM0 then
      self:scaleDensity(distf)
   end
   assert(self.exactLagFixM012 == false, "MaxwellianProjection: Specialized version of 'MaxwellianProjection' is required. Use 'VlasovProjection.MaxwellianProjection' or 'GkProjection.MaxwellianProjection'")
end


----------------------------------------------------------------------
return {
   FunctionProjection = FunctionProjection,
   MaxwellianProjection = MaxwellianProjection,
}
