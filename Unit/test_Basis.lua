-- Gkyl ------------------------------------------------------------------------
--
-- Test for basis functions
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local math = require "math"
local ffi = require "ffi"
local Unit = require "Unit"
local Basis = require "Basis"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_2d_m_p1()
   local basis = Basis.CartModalMaxOrder {
      ndim = 1, polyOrder = 2,
   }

   assert_equal(1, basis:ndim(), "Checking NDIM")
   assert_equal(2, basis:polyOrder(), "Checking polynomial order")
   assert_equal(3, basis:numBasis(), "Checking number of basis functions")
   assert_equal(1, basis:numSurfBasis(), "Checking number of surface basis functions")

   -- coordinates and values
   local z, bvalues = Lin.Vec(basis:ndim()), Lin.Vec(basis:numBasis())

   z[1] = 0.0
   basis:evalBasis(z, bvalues)
   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
   assert_equal(0.0, bvalues[2], "Checking values 2")
   assert_equal(-math.sqrt(5/2^3), bvalues[3], "Checking values 3")

   z[1] = -1.0
   basis:evalBasis(z, bvalues)
   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
   assert_equal(-math.sqrt(3/2), bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(5/2), bvalues[3], "Checking values 3")

   z[1] = 1.0
   basis:evalBasis(z, bvalues)
   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
   assert_equal(math.sqrt(3/2), bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(5/2), bvalues[3], "Checking values 3")
end

function test_2d_m_p2()

   local basis = Basis.CartModalMaxOrder {
      ndim = 2, polyOrder = 2,
   }

   assert_equal(2, basis:ndim(), "Checking NDIM")
   assert_equal(2, basis:polyOrder(), "Checking polynomial order")
   assert_equal(6, basis:numBasis(), "Checking number of basis functions")
   assert_equal(3, basis:numSurfBasis(), "Checking number of surface basis functions")

   -- coordinates and values
   local z, bvalues = Lin.Vec(basis:ndim()), Lin.Vec(basis:numBasis())

   z[1], z[2] = 0.0, 0.0
   basis:evalBasis(z, bvalues)

   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(0.0, bvalues[2], "Checking values 2")
   assert_equal(0.0, bvalues[3], "Checking values 3")
   assert_equal(0.0, bvalues[4], "Checking values 4")
   assert_equal(-math.sqrt(5)/4, bvalues[5], "Checking values 5")
   assert_equal(-math.sqrt(5)/4, bvalues[6], "Checking values 6")

   z[1], z[2] = -1.0, -1.0
   basis:evalBasis(z, bvalues)

   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(-math.sqrt(3)/2, bvalues[2], "Checking values 2")
   assert_equal(-math.sqrt(3)/2, bvalues[3], "Checking values 3")
   assert_equal(3.0/2.0, bvalues[4], "Checking values 4")
   assert_equal(math.sqrt(5)/2, bvalues[5], "Checking values 5")
   assert_equal(math.sqrt(5)/2, bvalues[6], "Checking values 6")

   z[1], z[2] = 1.0, 1.0
   basis:evalBasis(z, bvalues)

   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(math.sqrt(3)/2, bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(3)/2, bvalues[3], "Checking values 3")
   assert_equal(3.0/2.0, bvalues[4], "Checking values 4")
   assert_equal(math.sqrt(5)/2, bvalues[5], "Checking values 5")
   assert_equal(math.sqrt(5)/2, bvalues[6], "Checking values 6")   
end

function test_2d_s_p1()
   local basis = Basis.CartModalSerendipity {
      ndim = 1, polyOrder = 2,
   }

   assert_equal(1, basis:ndim(), "Checking NDIM")
   assert_equal(2, basis:polyOrder(), "Checking polynomial order")
   assert_equal(3, basis:numBasis(), "Checking number of basis functions")
   assert_equal(1, basis:numSurfBasis(), "Checking number of surface basis functions")

   -- coordinates and values
   local z, bvalues = Lin.Vec(basis:ndim()), Lin.Vec(basis:numBasis())

   z[1] = 0.0
   basis:evalBasis(z, bvalues)
   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
   assert_equal(0.0, bvalues[2], "Checking values 2")
   assert_equal(-math.sqrt(5/2^3), bvalues[3], "Checking values 3")

   z[1] = -1.0
   basis:evalBasis(z, bvalues)
   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
   assert_equal(-math.sqrt(3/2), bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(5/2), bvalues[3], "Checking values 3")

   z[1] = 1.0
   basis:evalBasis(z, bvalues)
   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
   assert_equal(math.sqrt(3/2), bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(5/2), bvalues[3], "Checking values 3")
end

function test_2d_s_p2()

   local basis = Basis.CartModalSerendipity {
      ndim = 2, polyOrder = 2,
   }

   assert_equal(2, basis:ndim(), "Checking NDIM")
   assert_equal(2, basis:polyOrder(), "Checking polynomial order")
   assert_equal(8, basis:numBasis(), "Checking number of basis functions")
   assert_equal(3, basis:numSurfBasis(), "Checking number of surface basis functions")

   -- coordinates and values
   local z, bvalues = Lin.Vec(basis:ndim()), Lin.Vec(basis:numBasis())

   z[1], z[2] = 0.0, 0.0
   basis:evalBasis(z, bvalues)

   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(0.0, bvalues[2], "Checking values 2")
   assert_equal(0.0, bvalues[3], "Checking values 3")
   assert_equal(0.0, bvalues[4], "Checking values 4")
   assert_equal(-math.sqrt(5)/4, bvalues[5], "Checking values 5")
   assert_equal(-math.sqrt(5)/4, bvalues[6], "Checking values 6")
   assert_equal(0.0, bvalues[7], "Checking values 7")
   assert_equal(0.0, bvalues[8], "Checking values 8")

   z[1], z[2] = -1.0, -1.0
   basis:evalBasis(z, bvalues)

   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(-math.sqrt(3)/2, bvalues[2], "Checking values 2")
   assert_equal(-math.sqrt(3)/2, bvalues[3], "Checking values 3")
   assert_equal(3.0/2.0, bvalues[4], "Checking values 4")
   assert_equal(math.sqrt(5)/2, bvalues[5], "Checking values 5")
   assert_equal(math.sqrt(5)/2, bvalues[6], "Checking values 6")
   assert_equal(-math.sqrt(15)/2, bvalues[7], "Checking values 7")
   assert_equal(-math.sqrt(15)/2, bvalues[8], "Checking values 8")

   z[1], z[2] = 1.0, 1.0
   basis:evalBasis(z, bvalues)

   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(math.sqrt(3)/2, bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(3)/2, bvalues[3], "Checking values 3")
   assert_equal(3.0/2.0, bvalues[4], "Checking values 4")
   assert_equal(math.sqrt(5)/2, bvalues[5], "Checking values 5")
   assert_equal(math.sqrt(5)/2, bvalues[6], "Checking values 6")
   assert_equal(math.sqrt(15)/2, bvalues[7], "Checking values 7")
   assert_equal(math.sqrt(15)/2, bvalues[8], "Checking values 8")
end

test_2d_m_p1()
test_2d_m_p2()

test_2d_s_p1()
test_2d_s_p2()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
