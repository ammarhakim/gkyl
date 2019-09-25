-- Gkyl ------------------------------------------------------------------------
--
-- App support code: FluidProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ProjectionBase = require "App.Projection.ProjectionBase"
local Proto          = require "Lib.Proto"
local Updater        = require "Updater"
local xsys           = require "xsys"

-- Shell class for fluid projections.
local FluidProjection = Proto(ProjectionBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function FluidProjection:init(tbl)
   self.tbl = tbl
end

function FluidProjection:fullInit(species)
   self.species = species

   self.confBasis  = species.confBasis
   self.confGrid   = species.grid

   self.cdim = self.confGrid:ndim()

   self.isInit       = xsys.pickBool(self.tbl.isInit, true)
   self.isSource     = xsys.pickBool(self.tbl.isSource, false)
   if self.isSource then self.isInit = false end
end

----------------------------------------------------------------------
-- Base class for projection of arbitrary function. No re-scaling.
local FunctionProjection = Proto(FluidProjection)

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
      onGrid          = self.confGrid,
      basis           = self.confBasis,
      evaluate        = func,
      projectOnGhosts = true
   }
   self.initFunc = func
end

function FunctionProjection:run(t, distf)
   self.project:advance(t, {}, {distf})
end

----------------------------------------------------------------------
return {
   FunctionProjection = FunctionProjection,
}
