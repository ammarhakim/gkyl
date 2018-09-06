-- Gkyl ------------------------------------------------------------------------
--
-- App support code: KineticProjection object
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto = require "Lib.Proto"
local ProjectionBase = require "App.Projection.ProjectionBase"
local Updater = require "Updater"
local DataStruct = require "DataStruct"
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
end

function KineticProjection:setName(nm)
   self.name = nm
end

function KineticProjection:set(t, distf)
   self.project:advance(t, 0.0, {}, {distf})
end


----------------------------------------------------------------------
local FunctionProjection = Proto(KineticProjection)

function FunctionProjection:fullInit(species)
   FunctionProjection.super.fullInit(self, species)

   local func = assert(self.tbl[1],
		      "FunctionProjection: Must specify the function")
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
   self.drift = tbl.drift or function (t, zn) return nil end
   self.temperature = assert(tbl.temperature,
			     "Maxwellian: must specify 'temperature'")
   -- check for constants instead of functions
   if type(self.density) ~= "function" then
      self.density = function (t, zn) return tbl.density end
   end
   if type(self.drift) ~= "function" then
      if type(tbl.drift) == "number" then 
	 assert(self.numVelDims == 1,
		string.format("KineticProjection: Number of 'drift' components (1) must correspond to the number of velocity dimensions (%d)", self.numVelDims))
	 self.drift = function (t, zn) return {tbl.drift} end
      else 
	 assert(#tbl.drift == self.numVelDims,
		string.format("KineticProjection: Number of 'drift' components (%d) must correspond to the number of velocity dimensions (%d)", #tbl.drift, self.numVelDims))
	 self.drift = function (t, zn) return tbl.drift end
      end
   end
   if type(self.temperature) ~= "function" then
      self.temperature = function (t, zn) return tbl.temperature end
   end

   self.scaleM0 = xsys.pickBool(tbl.scaleM0, true)
   self.lagFixM012 = xsys.pickBool(tbl.lagFixM012, false)

   local func = function (t, zn)
      return species:Maxwellian(zn, self.density(t, zn, species),
				self.temperature(t, zn, species),
				self.drift(t, zn, species))
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

function MaxwellianProjection:lagrangeFix(distf)
   local M0, dM0 = self.species:allocMoment(), self.species:allocMoment()
   local M1 = self.species:allocVectorMoment(self.numVelDims)
   local dM1 = self.species:allocVectorMoment(self.numVelDims)
   local M2, dM2 = self.species:allocMoment(), self.species:allocMoment()

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
   project:advance(0.0, 0.0, {}, {dM0})
   dM0:accumulate(-1.0, M0)

   self.species.momDensityCalc:advance(0.0, 0.0, {distf}, {M1})
   func = function (t, zn)
      return self.density(t, zn, self.species) * 
	 self.drift(t, zn, self.species)[1]
   end
   project = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = func,
      projectOnGhosts = true,
   }
   project:advance(0.0, 0.0, {}, {dM1})
   dM1:accumulate(-1.0, M1)

   self.species.ptclEnergyCalc:advance(0.0, 0.0, {distf}, {M1})
   func = function (t, zn)
      return self.density(t, zn, self.species) * 
	 ( self.drift(t, zn, self.species)[1] * self.drift(t, zn, self.species)[1] +
	 self.temperature(t, zn, self.species) / self.species.mass )
   end
   project = Updater.ProjectOnBasis {
      onGrid = self.confGrid,
      basis = self.confBasis,
      evaluate = func,
      projectOnGhosts = true,
   }
   project:advance(0.0, 0.0, {}, {dM2})
   dM2:accumulate(-1.0, M2)

   local lagFix = Updater.LagrangeFix {
      onGrid = self.phaseGrid,
      phaseGrid = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confGrid = self.confGrid,
      confBasis = self.confBasis,
   }
   lagFix:advance(0.0, 0.0, {dM0, dM1, dM2}, {distf})
end

function MaxwellianProjection:set(t, distf)
   self.project:advance(t, 0.0, {}, {distf})
   if self.scaleM0 then
      self:scaleDensity(distf)
   end
   if self.lagFixM012 then
      self:lagrangeFix(distf)
   end
end


----------------------------------------------------------------------
return {
   FunctionProjection = FunctionProjection,
   MaxwellianProjection = MaxwellianProjection,
}
