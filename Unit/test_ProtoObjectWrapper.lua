-- Gkyl ------------------------------------------------------------------------
--
-- Test for prototype wrapper
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Proto = require "Lib.Proto"
local ProtoObjectWrapper = require "Lib.ProtoObjectWrapper"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local Particle = Proto()
   function Particle:init(tbl)
      self.x = tbl.x
      self.y = tbl.y
   end

   -- create wrapper
   local ParticleWrapper = ProtoObjectWrapper(Particle)
   --  partial creation of object
   local PartialParticle = ParticleWrapper { x = 1.0 }

   -- finish creation of particle
   local p = PartialParticle { y = 2.0 }

   assert_equal(1.0, p.x, "Checking value of x")
   assert_equal(2.0, p.y, "Checking value of y")

   local q = PartialParticle { y = 3.0 }

   assert_equal(1.0, q.x, "Checking value of x")
   assert_equal(3.0, q.y, "Checking value of y")
end

test_1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
