-- Gkyl ------------------------------------------------------------------------
--
-- App support code: VlasovProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto            = require "Lib.Proto"
local ProjectionBase   = require "App.Projection.ProjectionBase"
local xsys             = require "xsys"
local Updater          = require "Updater"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"

-- Shell class for Vlasov projections.
local VlasovProjection = Proto(ProjectionBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VlasovProjection:init(tbl) self.tbl = tbl end

function VlasovProjection:fullInit(mySpecies)
   self.mass   = mySpecies:getMass()
   self.charge = mySpecies:getCharge()

   self.fromFile = self.tbl.fromFile

   self.exactScaleM0    = xsys.pickBool(self.tbl.exactScaleM0, true)
   self.exactScaleM012  = xsys.pickBool(self.tbl.exactScaleM012, false)
   if self.exactScaleM012 then self.exactScaleM0 = false end

   self.power = self.tbl.power
   self.scaleWithSourcePower = xsys.pickBool(self.tbl.scaleWithSourcePower, false)
end

function VlasovProjection:fullInit(mySpecies)
   self.phaseBasis = mySpecies.basis
   self.phaseGrid  = mySpecies.grid
   self.confBasis  = mySpecies.confBasis
   self.confGrid   = mySpecies.confGrid

   self.cdim = self.confGrid:ndim()
   self.vdim = self.phaseGrid:ndim() - self.confGrid:ndim()

   self.weakMultiplyConfPhase = Updater.CartFieldBinOp {
      weakBasis  = self.phaseBasis,  operation = "Multiply",
      fieldBasis = self.confBasis,   onGhosts  = true,
   }
end

----------------------------------------------------------------------
local FunctionProjection = Proto(VlasovProjection)

function FunctionProjection:createSolver(mySpecies)
   FunctionProjection.super.createSolver(self, mySpecies)

   self.initFunc = self.tbl.func
   if not self.initFunc then self.initFunc = self.tbl[1] end

   assert(self.initFunc, "FunctionProjection: Must specify the function")
   assert(type(self.initFunc) == "function", "The input must be a table containing function")

   if self.fromFile then
      self.ioMethod  = "MPI"
      self.writeGhost = true
      self.fieldIo = AdiosCartFieldIo {
         elemType   = mySpecies.distf[1]:elemType(),
         method     = self.ioMethod,
         writeGhost = self.writeGhost,
         metaData   = {polyOrder = self.phaseBasis:polyOrder(),
                       basisType = self.phaseBasis:id()}
      }
   else
      self.project = Updater.ProjectOnBasis {
         onGrid = self.phaseGrid,   evaluate = self.initFunc,
         basis  = self.phaseBasis,  onGhosts = true
      }
   end
end

function FunctionProjection:scaleDensity(mySpecies, distf, currentM0, targetM0)
   local M0mod = mySpecies:allocMoment()

   local weakDivision = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,         operation = "Divide",
      onRange   = M0mod:localExtRange(),  onGhosts  = true,
   }

   -- Calculate M0mod = targetM0 / currentM0.
   weakDivision:advance(0.0, {currentM0, targetM0}, {M0mod})
   -- Calculate distff = M0mod * distf.
   self.weakMultiplyConfPhase:advance(0.0, {M0mod, distf}, {distf})
end

function FunctionProjection:advance(t, inFlds, outFlds)
   local mySpecies = inFlds[1]
   local distf     = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
      self.project:advance(t, {}, {distf})
   end
end

----------------------------------------------------------------------
local MaxwellianProjection = Proto(VlasovProjection)

function MaxwellianProjection:createSolver(mySpecies)
   MaxwellianProjection.super.createSolver(self, mySpecies)

   local tbl = self.tbl
   self.density     = assert(tbl.density, "MaxwellianProjection: must specify 'density'")
   self.driftSpeed  = tbl.driftSpeed or (
      self.vdim==1 and function(t, zn) return 0. end or (
      self.vdim==2 and function(t, zn) return 0., 0. end or (
      self.vdim==3 and function(t, zn) return 0., 0., 0. end or
      assert(false, "VlasovProjection.MaxwellianProjection: wrong vdim") ) ) )
   self.temperature = assert(tbl.temperature, "MaxwellianProjection: must specify 'temperature'")

   -- Check for constants instead of functions.
   if type(self.density) == "number" then
      self.density = function (t, zn) return tbl.density end
   end
   if type(self.driftSpeed) == "number" then
      assert(self.vdim == 1, "MaxwellianProjection: driftSpeed must return vdim values.")
      self.driftSpeed = function (t, zn) return tbl.driftSpeed end
   elseif type(self.driftSpeed) == "table" then
      assert(#tbl.driftSpeed == self.vdim, "MaxwellianProjection: driftSpeed must return vdim values.")
      self.driftSpeed =
         self.vdim==1 and function (t, zn) return tbl.driftSpeed[1] end or (
         self.vdim==2 and function (t, zn) return tbl.driftSpeed[1], tbl.driftSpeed[2] end or (
         self.vdim==3 and function (t, zn) return tbl.driftSpeed[1], tbl.driftSpeed[2], tbl.driftSpeed[3] end or 
         assert(false, "VlasovProjection.MaxwellianProjection: wrong vdim") ) )
   end
   if type(self.temperature) == "number" then
      self.temperature = function (t, zn) return tbl.temperature end
   end

   self.initFunc =
      self.vdim==1 and function (t, zn)
         local ux = self.driftSpeed(t, zn)
         return mySpecies:Maxwellian(zn, self.density(t, zn), {ux}, self.temperature(t, zn))
      end or (
      self.vdim==2 and function (t, zn)
         local ux, uy = self.driftSpeed(t, zn)
         return mySpecies:Maxwellian(zn, self.density(t, zn), {ux,uy}, self.temperature(t, zn))
      end or (
      self.vdim==3 and function (t, zn)
         local ux, uy, uz = self.driftSpeed(t, zn)
         return mySpecies:Maxwellian(zn, self.density(t, zn), {ux,uy,uz}, self.temperature(t, zn))
      end or assert(false, "VlasovProjection.MaxwellianProjection: wrong vdim")
   ) )

   if self.fromFile then
      self.ioMethod  = "MPI"
      self.writeGhost = true
      self.fieldIo = AdiosCartFieldIo {
         elemType  = mySpecies.distf[1]:elemType(),
         method    = self.ioMethod,
         writeGhost = self.writeGhost,
         metaData  = {polyOrder = self.phaseBasis:polyOrder(),
                      basisType = self.phaseBasis:id(),
                      charge    = self.charge,
                      mass      = self.mass,},
      }
   else
      self.project = Updater.ProjectOnBasis {
         onGrid = self.phaseGrid,   evaluate = self.initFunc,
         basis  = self.phaseBasis,  onGhosts = true
      }
   end
end

function MaxwellianProjection:scaleDensity(mySpecies, distf)
   local M0e, M0 = mySpecies:allocMoment(), mySpecies:allocMoment()
   local M0mod   = mySpecies:allocMoment()

   mySpecies.numDensityCalc:advance(0.0, {distf}, {M0})
   local project = Updater.ProjectOnBasis {
      onGrid = self.confGrid,   evaluate = function (t, zn) return self.density(t, zn) end,
      basis  = self.confBasis,  onGhosts = true,
   }
   project:advance(0.0, {}, {M0e})

   local weakDivision = Updater.CartFieldBinOp {
      weakBasis = self.confBasis,       operation = "Divide",
      onRange   = M0e:localExtRange(),  onGhosts  = true,
   }

   -- Calculate M0mod = M0e / M0.
   weakDivision:advance(0.0, {M0, M0e}, {M0mod})
   -- Calculate distff = M0mod * distf.
   self.weakMultiplyConfPhase:advance(0.0, {M0mod, distf}, {distf})
end

function MaxwellianProjection:advance(t, inFlds, outFlds)
   local mySpecies = inFlds[1]
   local distf     = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(distf, self.fromFile)
   else
      self.project:advance(t, {}, {distf})
   end
   if self.exactScaleM0 then
      self:scaleDensity(mySpecies, distf)
   end
end

----------------------------------------------------------------------
return {
   FunctionProjection   = FunctionProjection,
   MaxwellianProjection = MaxwellianProjection,
}
