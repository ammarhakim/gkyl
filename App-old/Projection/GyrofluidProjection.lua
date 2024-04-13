-- Gkyl ------------------------------------------------------------------------
--
-- App support code: GyrofluidProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Proto      = require "Lib.Proto"
local Updater    = require "Updater"
local xsys       = require "xsys"
local DataStruct = require "DataStruct"
local FluidProjectionParent    = require("App.Projection.FluidProjection").FluidProjection
local FunctionProjectionParent = require("App.Projection.FluidProjection").FunctionProjection

--------------------------------------------------------------------------------
-- Gyrofluid-specific FunctionProjection includes Jacobian factors in initFunc.
local FunctionProjection = Proto(FunctionProjectionParent)

function FunctionProjection:advance(time, inFlds, outFlds)
   -- Run the FluidProjection :FunctionProjection to perform the projection.
   FunctionProjection.super.advance(self, time, inFlds, outFlds)

   local extField = inFlds[1]
   local momOut   = outFlds[1]

   local jacob = extField.geo.jacobGeo
   if jacob then self.weakMultiply:advance(0, {momOut, jacob}, {momOut}) end
end

--------------------------------------------------------------------------------
-- Gyrofluid-specific projection object that takes in the density, parallel flow speed,
-- parallel temperature and perpendicular temperature, and computes the gyrofluid moments
-- mass density, parallel momentum density, energy density and perpendicular pressure
-- (divided by bmag), all multiplied by the Jacobian.

local GyrofluidProjection = Proto(FluidProjectionParent)

function GyrofluidProjection:fullInit(mySpecies)
   GyrofluidProjection.super.fullInit(self, mySpecies)

   local tbl = self.tbl
   self.density    = assert(tbl.density, "GyrofluidProjection: must specify 'density'")
   self.driftSpeed = tbl.driftSpeed or function(t, zn) return 0. end
   self.Tpar       = assert(tbl.parallelTemperature, "GyrofluidProjection: must specify 'parallelTemperature'")
   self.Tperp      = assert(tbl.perpendicularTemperature, "GyrofluidProjection: must specify 'perpendicularTemperature'")

   -- Check for constants instead of functions.
   if type(self.density) ~= "function" then
      self.density = function (t, zn) return tbl.density end
   end
   if type(self.driftSpeed) ~= "function" then
      self.driftSpeed = function (t, zn) return tbl.driftSpeed end
   end
   if type(self.Tpar) ~= "function" then
      self.Tpar = function (t, zn) return tbl.parallelTemperature end
   end
   if type(self.Tperp) ~= "function" then
      self.Tperp = function (t, zn) return tbl.perpendicularTemperature end
   end

   self.charge, self.mass = mySpecies.charge, mySpecies.mass

   if self.fromFile then
      self.writeSkin = true
      self.fieldIo = AdiosCartFieldIo {
         elemType   = mySpecies.moments[1]:elemType(),
         writeGhost = self.writeGhost,
         metaData   = {polyOrder = self.basis:polyOrder(),
                       basisType = self.basis:id(),
                       charge    = self.charge,
                       mass      = self.mass,},
      }
   else
      self.project = Updater.ProjectOnBasis {
         onGrid = self.grid,   evaluate = function(t, xn) return 0. end,   -- Set later.
         basis  = self.basis,  onGhosts = true,
      }
      self.weakDivide = Updater.CartFieldBinOp {
         weakBasis = self.basis,  operation = "Divide",
         onRange   = mySpecies.moments[1]:localExtRange(),  onGhosts = true,
      }
      -- Will also need some temporary fields. Can use those declared in GyrofluidSpecies.
      self.jacM0      = mySpecies.jacM0
      self.mJacM0     = mySpecies.mJacM0
      self.mJacM1     = mySpecies.mJacM1
      self.mJacM2     = mySpecies.mJacM2
      self.jacM2perp  = mySpecies.jacM2perp
      self.uParSelf   = mySpecies.uParSelf
      self.TperpSelf  = mySpecies.TperpSelf
      self.TparSelf   = mySpecies.TparSelf 
      self.pPerpJac   = mySpecies.pPerpJac
      self.pParJac    = mySpecies.pParJac
      self.mJacM2flow = mySpecies.mJacM2flow
      -- Useful to have a unit field.
      self.unitFld = mySpecies:allocMoment()
      self.project:setFunc(function(t,xn) return 1. end)
      self.project:advance(t, {}, {self.unitFld})

      self.momOff = mySpecies:getMomOff()
   end
end

function GyrofluidProjection:advance(tm, inFlds, outFlds)
   local extField = inFlds[1]
   local momOut   = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(momOut, self.fromFile)
   else
      local bmag  = extField.geo.bmag
      local jacob = extField.geo.jacobGeo or self.unitFld

      -- Compute the zeroth moment (number density), times mass, times Jacobian.
      self.project:setFunc(self.density)
      self.project:advance(t, {}, {self.jacM0})
      self.weakMultiply:advance(0, {self.jacM0, jacob}, {self.jacM0})
      self.mJacM0:combine(self.mass, self.jacM0)

      -- Compute the 1st moment times mass (momentum density), times Jacobian.
      self.project:setFunc(self.driftSpeed)
      self.project:advance(t, {}, {self.uParSelf})
      self.weakMultiply:advance(0, {self.mJacM0, self.uParSelf}, {self.mJacM1})

      -- Compute the perpendicular pressure divided by B, times Jacobian.
      self.project:setFunc(self.Tperp)
      self.project:advance(t, {}, {self.TperpSelf})
      self.weakMultiply:advance(0, {self.TperpSelf, self.jacM0}, {self.pPerpJac})
      self.weakDivide:advance(0, {bmag, self.pPerpJac}, {self.jacM2perp})

      -- Compute the parallel pressure, times Jacobian.
      self.project:setFunc(self.Tpar)
      self.project:advance(t, {}, {self.TparSelf})
      self.weakMultiply:advance(0, {self.TparSelf, self.jacM0}, {self.pParJac})

      -- Calculate the parallel flow energy density, times Jacobian.
      self.weakMultiply:advance(0, {self.uParSelf, self.mJacM1}, {self.mJacM2flow})

      -- Calculate the total particle kinetic energy density (times Jacobian).
      self.mJacM2:combine(0.5, self.pParJac, 0.5, self.mJacM2flow, 1., self.pPerpJac)

      -- Combine the moments into a single field.
      momOut:combineOffset(1., self.mJacM0, self.momOff[1], 
                           1., self.mJacM1, self.momOff[2],
                           1., self.mJacM2, self.momOff[3],
                           1., self.jacM2perp, self.momOff[4])

   end
end

----------------------------------------------------------------------
return {
   FunctionProjection  = FunctionProjection,
   GyrofluidProjection = GyrofluidProjection,
}
