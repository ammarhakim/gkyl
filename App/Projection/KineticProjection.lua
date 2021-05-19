-- Gkyl ------------------------------------------------------------------------
--
-- App support code: KineticProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ProjectionBase   = require "App.Projection.ProjectionBase"
local Proto            = require "Lib.Proto"
local Updater          = require "Updater"
local xsys             = require "xsys"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"

-- Shell class for kinetic projections.
local KineticProjection = Proto(ProjectionBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function KineticProjection:init(tbl)
   self.tbl = tbl
end

function KineticProjection:fullInit(species)
   self.species = species

   self.phaseBasis = species.basis
   self.phaseGrid  = species.grid
   self.confBasis  = species.confBasis
   self.confGrid   = species.confGrid
   self.mass       = species:getMass()
   self.charge     = species:getCharge()

   self.cdim = self.confGrid:ndim()
   self.vdim = self.phaseGrid:ndim() - self.confGrid:ndim()

   self.vDegFreedom = species.vDegFreedom   -- Only defined in GkSpecies.

   self.fromFile = self.tbl.fromFile

   self.exactScaleM0    = xsys.pickBool(self.tbl.exactScaleM0, true)
   self.exactScaleM012  = xsys.pickBool(self.tbl.exactScaleM012, false)
   if self.exactScaleM012 then self.exactScaleM0 = false end
   self.exactLagFixM012 = xsys.pickBool(self.tbl.exactLagFixM012, false)

   self.power = self.tbl.power
   self.scaleWithSourcePower = xsys.pickBool(self.tbl.scaleWithSourcePower, false)

   self.weakMultiplyConfPhase = Updater.CartFieldBinOp {
      onGrid     = self.phaseGrid,
      weakBasis  = self.phaseBasis,
      fieldBasis = self.confBasis,
      operation  = "Multiply",
      onGhosts   = true,
   }
end

----------------------------------------------------------------------
-- Base class for projection of arbitrary function. No re-scaling.
local FunctionProjection = Proto(KineticProjection)

function FunctionProjection:fullInit(species)
   FunctionProjection.super.fullInit(self, species)

   local func = self.tbl.func
   if not func then func = self.tbl[1] end

   assert(func, "FunctionProjection: Must specify the function")
   assert(type(func) == "function", "The input must be a table containing function")

   if self.fromFile then
      self.ioMethod  = "MPI"
      self.writeGhost = true
      self.fieldIo = AdiosCartFieldIo {
         elemType   = species.distf[1]:elemType(),
         method     = self.ioMethod,
         writeGhost = self.writeGhost,
         metaData   = {polyOrder = self.phaseBasis:polyOrder(),
                       basisType = self.phaseBasis:id()}
      }
   else
      self.project = Updater.ProjectOnBasis {
         onGrid   = self.phaseGrid,
         basis    = self.phaseBasis,
         evaluate = func,
         onGhosts = true
      }
   end
   self.initFunc = func
end

function FunctionProjection:advance(t, inFlds, outFlds)
   local distf = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
      self.project:advance(t, {}, {distf})
   end
end

----------------------------------------------------------------------
-- Base class for projection of Maxwellian function. Option for only density re-scaling.
local MaxwellianProjection = Proto(KineticProjection)

function MaxwellianProjection:fullInit(species)
   MaxwellianProjection.super.fullInit(self, species)

   local tbl = self.tbl
   self.density     = assert(tbl.density, "Maxwellian: must specify 'density'")
   self.driftSpeed  = tbl.driftSpeed or function(t, zn)
      if self.vDegFreedom then return 0.   -- Gyrokinetics. Vlasov doesn't define vDegFreedom. 
      elseif self.vdim==1 then return {0.}
      elseif self.vdim==2 then return {0., 0.}
      elseif self.vdim==3 then return {0., 0., 0.} end
   end
   self.temperature = assert(tbl.temperature,
			     "Maxwellian: must specify 'temperature'")

   -- Check for constants instead of functions.
   if type(self.density) ~= "function" then
      self.density = function (t, zn) return tbl.density end
   end
   if type(self.driftSpeed) ~= "function" then
      self.driftSpeed = function (t, zn) return tbl.driftSpeed end
   end
   if type(self.temperature) ~= "function" then
      self.temperature = function (t, zn) return tbl.temperature end
   end

   self.initFunc = function (t, zn)
      return species:Maxwellian(zn, self.density(t, zn, species),
				self.temperature(t, zn, species),
				self.driftSpeed(t, zn, species))
   end

   if self.fromFile then
      self.ioMethod  = "MPI"
      self.writeGhost = true
      self.fieldIo = AdiosCartFieldIo {
         elemType  = species.distf[1]:elemType(),
         method    = self.ioMethod,
         writeGhost = self.writeGhost,
         metaData  = {polyOrder = self.phaseBasis:polyOrder(),
                      basisType = self.phaseBasis:id(),
                      charge    = self.charge,
                      mass      = self.mass,},
      }
   else
      self.project = Updater.ProjectOnBasis {
         onGrid   = self.phaseGrid,
         basis    = self.phaseBasis,
         evaluate = self.initFunc,
         onGhosts = true
      }
   end
end

function MaxwellianProjection:scaleDensity(distf)
   local M0e, M0 = self.species:allocMoment(), self.species:allocMoment()
   local M0mod   = self.species:allocMoment()

   self.species.numDensityCalc:advance(0.0, {distf}, {M0})
   local func = function (t, zn)
      return self.density(t, zn, self.species)
   end
   local project = Updater.ProjectOnBasis {
      onGrid   = self.confGrid,
      basis    = self.confBasis,
      evaluate = func,
      onGhosts = true,
   }
   project:advance(0.0, {}, {M0e})

   local weakDivision = Updater.CartFieldBinOp {
      onGrid    = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
      onGhosts  = true,
   }

   -- Calculate M0mod = M0e / M0.
   weakDivision:advance(0.0, {M0, M0e}, {M0mod})
   -- Calculate distff = M0mod * distf.
   self.weakMultiplyConfPhase:advance(0.0, {M0mod, distf}, {distf})
end


function MaxwellianProjection:advance(t, inFlds, outFlds)
   local distf = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
      self.project:advance(t, {}, {distf})
   end
   if self.exactScaleM0 then
      self:scaleDensity(distf)
   end
   assert(self.exactLagFixM012 == false, "MaxwellianProjection: Specialized version of 'MaxwellianProjection' is required. Use 'VlasovProjection.MaxwellianProjection' or 'GkProjection.MaxwellianProjection'")
end


----------------------------------------------------------------------
return {
   FunctionProjection   = FunctionProjection,
   MaxwellianProjection = MaxwellianProjection,
}
