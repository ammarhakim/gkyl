-- Gkyl ------------------------------------------------------------------------
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Infrastructure loads
local DecompRegionCalc = require "Lib.CartDecomp"
local LinearTrigger = require "Lib.LinearTrigger"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local date = require "xsys.date"
local lfs = require "lfs"
local lume = require "Lib.lume"
local xsys = require "xsys"
local ffi = require "ffi"

-- Euler equation initialization
local Euler = function(tbl)
   return G0.Moments.Eq.Euler.new(tbl)
end

-- IsoEuler equation initialization
local IsoEuler = function(tbl)
   return G0.Moments.Eq.IsoEuler.new(tbl)
end

-- TenMoment equation initialization
local TenMoment = function(tbl)
   return G0.Moments.Eq.TenMoment.new(tbl)
end

-- Moment species initialization
local Species = function(tbl)
   return G0.Moments.Species.new(tbl)
end

-- Moment field initialization
local Field = function(tbl)
   return G0.Moments.Field.new(tbl)
end

-- top-level moments App
local App = Proto()

function App:init(tbl)
   self.g0App = G0.Moments.App.new(tbl)
end

function App:run(numSteps)
   self.g0App:run(numSteps)
end

return {
   App = App,
   Species = Species,
   Field = Field,

   -- various boundary conditions for species
   SpeciesBc = {
      bcCopy = G0.SpeciesBc.bcCopy,
      bcWall = G0.SpeciesBc.bcWall,
      bcReflect = G0.SpeciesBc.bcReflect,
      bcAbsorb = G0.SpeciesBc.bcAbsorb,
      bcNoSlip = G0.SpeciesBc.bcNoSlip,
      bcWedge = G0.SpeciesBc.bcWedge,
      bcFunc = G0.SpeciesBc.bcFunc,
      bcFixedFunc = G0.SpeciesBc.bcFixedFunc,
      bcZeroFlux = G0.SpeciesBc.bcZeroFlux,
   },

   -- various boundary conditions for fields
   FieldBc = {
      bcCopy = G0.FieldBc.bcCopy,
      bcWall = G0.FieldBc.bcWall,
      bcPECWall = G0.FieldBc.bcPECWall,
      bcSymWall = G0.FieldBc.bcSymWall,
      bcReservoir = G0.FieldBc.bcReservoir,
      bcWedge = G0.FieldBc.bcWedge,
      bcFunc = G0.FieldBc.bcFunc,
   },

   -- supported equation systems
   Eq = {
      Euler = Euler,
      IsoEuler = IsoEuler,
      TenMoment = TenMoment
   }
}
