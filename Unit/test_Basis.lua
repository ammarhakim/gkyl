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

---- Disabling tests of maximal order basis (not kept in g0).
--function test_2d_m_p1()
--   local basis = Basis.CartModalMaxOrder {
--      ndim = 1, polyOrder = 2,
--   }
--
--   assert_equal(1, basis:ndim(), "Checking NDIM")
--   assert_equal(2, basis:polyOrder(), "Checking polynomial order")
--   assert_equal(3, basis:numBasis(), "Checking number of basis functions")
--
--   -- coordinates and values
--   local z, bvalues = Lin.Vec(basis:ndim()), Lin.Vec(basis:numBasis())
--
--   z[1] = 0.0
--   basis:evalBasis(z:data(), bvalues:data())
--   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
--   assert_equal(0.0, bvalues[2], "Checking values 2")
--   assert_equal(-math.sqrt(5/2^3), bvalues[3], "Checking values 3")
--
--   z[1] = -1.0
--   basis:evalBasis(z:data(), bvalues:data())
--   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
--   assert_equal(-math.sqrt(3/2), bvalues[2], "Checking values 2")
--   assert_equal(math.sqrt(5/2), bvalues[3], "Checking values 3")
--
--   z[1] = 1.0
--   basis:evalBasis(z:data(), bvalues:data())
--   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
--   assert_equal(math.sqrt(3/2), bvalues[2], "Checking values 2")
--   assert_equal(math.sqrt(5/2), bvalues[3], "Checking values 3")
--
--   -- surface expansions
--   local vol = {1.0, 0.5, 0.25}
--   local surf = {0.0} -- to store surface expansion
--end
--
--function test_2d_m_p2()
--
--   local basis = Basis.CartModalMaxOrder {
--      ndim = 2, polyOrder = 2,
--   }
--
--   assert_equal(2, basis:ndim(), "Checking NDIM")
--   assert_equal(2, basis:polyOrder(), "Checking polynomial order")
--   assert_equal(6, basis:numBasis(), "Checking number of basis functions")
--
--   -- coordinates and values
--   local z, bvalues = Lin.Vec(basis:ndim()), Lin.Vec(basis:numBasis())
--
--   z[1], z[2] = 0.0, 0.0
--   basis:evalBasis(z:data(), bvalues:data())
--
--   assert_equal(0.5, bvalues[1], "Checking values 1")
--   assert_equal(0.0, bvalues[2], "Checking values 2")
--   assert_equal(0.0, bvalues[3], "Checking values 3")
--   assert_equal(0.0, bvalues[4], "Checking values 4")
--   assert_equal(-math.sqrt(5)/4, bvalues[5], "Checking values 5")
--   assert_equal(-math.sqrt(5)/4, bvalues[6], "Checking values 6")
--
--   z[1], z[2] = -1.0, -1.0
--   basis:evalBasis(z:data(), bvalues:data())
--
--   assert_equal(0.5, bvalues[1], "Checking values 1")
--   assert_equal(-math.sqrt(3)/2, bvalues[2], "Checking values 2")
--   assert_equal(-math.sqrt(3)/2, bvalues[3], "Checking values 3")
--   assert_equal(3.0/2.0, bvalues[4], "Checking values 4")
--   assert_equal(math.sqrt(5)/2, bvalues[5], "Checking values 5")
--   assert_equal(math.sqrt(5)/2, bvalues[6], "Checking values 6")
--
--   z[1], z[2] = 1.0, 1.0
--   basis:evalBasis(z:data(), bvalues:data())
--
--   assert_equal(0.5, bvalues[1], "Checking values 1")
--   assert_equal(math.sqrt(3)/2, bvalues[2], "Checking values 2")
--   assert_equal(math.sqrt(3)/2, bvalues[3], "Checking values 3")
--   assert_equal(3.0/2.0, bvalues[4], "Checking values 4")
--   assert_equal(math.sqrt(5)/2, bvalues[5], "Checking values 5")
--   assert_equal(math.sqrt(5)/2, bvalues[6], "Checking values 6")
--
--   -- surface expansions
--   local vol = {1.0, 0.5, 0.25, 0.1, 0.1, 0.1}
--   local surf = {0.0, 0.0} -- to store surface expansion
--end
--
--function test_fs_m_2d_p1()
--   local basis = Basis.CartModalMaxOrder { ndim = 2, polyOrder = 1 }
--   local fIn, fOut = Lin.Vec(basis:numBasis()), Lin.Vec(basis:numBasis())
--
--   for i = 1, #fIn do fIn[i] = 1 end
--
--   basis:flipSign(1, fIn:data(), fOut:data())
--   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
--   assert_equal(-fIn[2], fOut[2], "Checking if flip-sign worked")
--   assert_equal(fIn[3], fOut[3], "Checking if flip-sign worked")   
--
--   basis:flipSign(2, fIn:data(), fOut:data())
--   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
--   assert_equal(fIn[2], fOut[2], "Checking if flip-sign worked")
--   assert_equal(-fIn[3], fOut[3], "Checking if flip-sign worked")
--end
--
--function test_fs_m_2d_p2()
--   local basis = Basis.CartModalMaxOrder { ndim = 2, polyOrder = 2 }
--   local fIn, fOut = Lin.Vec(basis:numBasis()), Lin.Vec(basis:numBasis())
--
--   for i = 1, #fIn do fIn[i] = 1 end
--
--   basis:flipSign(1, fIn:data(), fOut:data())
--   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
--   assert_equal(-fIn[2], fOut[2], "Checking if flip-sign worked")
--   assert_equal(fIn[3], fOut[3], "Checking if flip-sign worked")   
--   assert_equal(-fIn[4], fOut[4], "Checking if flip-sign worked")
--   assert_equal(fIn[5], fOut[5], "Checking if flip-sign worked")
--   assert_equal(fIn[6], fOut[6], "Checking if flip-sign worked")
--
--   basis:flipSign(2, fIn:data(), fOut:data())
--   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
--   assert_equal(fIn[2], fOut[2], "Checking if flip-sign worked")
--   assert_equal(-fIn[3], fOut[3], "Checking if flip-sign worked")   
--   assert_equal(-fIn[4], fOut[4], "Checking if flip-sign worked")
--   assert_equal(fIn[5], fOut[5], "Checking if flip-sign worked")
--   assert_equal(fIn[6], fOut[6], "Checking if flip-sign worked")
--end

function test_2d_s_p1()
   local basis = Basis.CartModalSerendipity {
      ndim = 1, polyOrder = 2,
   }

   assert_equal(1, basis:ndim(), "Checking NDIM")
   assert_equal(2, basis:polyOrder(), "Checking polynomial order")
   assert_equal(3, basis:numBasis(), "Checking number of basis functions")

   -- coordinates and values
   local z, bvalues = Lin.Vec(basis:ndim()), Lin.Vec(basis:numBasis())

   z[1] = 0.0
   basis:evalBasis(z:data(), bvalues:data())
   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
   assert_equal(0.0, bvalues[2], "Checking values 2")
   assert_equal(-math.sqrt(5/2^3), bvalues[3], "Checking values 3")

   z[1] = -1.0
   basis:evalBasis(z:data(), bvalues:data())
   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
   assert_equal(-math.sqrt(3/2), bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(5/2), bvalues[3], "Checking values 3")

   z[1] = 1.0
   basis:evalBasis(z:data(), bvalues:data())
   assert_equal(1/math.sqrt(2), bvalues[1], "Checking values 1")
   assert_equal(math.sqrt(3/2), bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(5/2), bvalues[3], "Checking values 3")
end

function test_1x2v_hybrid_p1()
   local basis = Basis.CartModalHybrid {
      cdim = 1,  vdim = 2,
   }

   assert_equal(3, basis:ndim(), "Checking NDIM")
   assert_equal(1, basis:polyOrder(), "Checking polynomial order")
   assert_equal(16, basis:numBasis(), "Checking number of basis functions")

   -- coordinates and values
   local z, bvalues = Lin.Vec(basis:ndim()), Lin.Vec(basis:numBasis())

   z[1], z[2] , z[3] = 0.0, 0.0, 0.0
   basis:evalBasis(z:data(), bvalues:data())
   assert_equal(1./math.sqrt(8), bvalues[1], "Checking values 1")
   assert_equal(0.0, bvalues[2], "Checking values 2")
   assert_equal(0.0, bvalues[3], "Checking values 3")
   assert_equal(0.0, bvalues[4], "Checking values 4")
   assert_equal(0.0, bvalues[5], "Checking values 5")
   assert_equal(0.0, bvalues[6], "Checking values 6")
   assert_equal(0.0, bvalues[7], "Checking values 7")
   assert_equal(0.0, bvalues[8], "Checking values 8")
   assert_equal(-math.sqrt(5)/math.sqrt(32.), bvalues[9], "Checking values 9")
   assert_equal(0.0, bvalues[10], "Checking values 10")
   assert_equal(0.0, bvalues[11], "Checking values 11")
   assert_equal(0.0, bvalues[12], "Checking values 12")
   assert_equal(-math.sqrt(5)/math.sqrt(32.), bvalues[9], "Checking values 13")
   assert_equal(0.0, bvalues[14], "Checking values 14")
   assert_equal(0.0, bvalues[15], "Checking values 15")
   assert_equal(0.0, bvalues[16], "Checking values 16")

   z[1], z[2], z[3] = -1.0, -1.0, -1.0
   basis:evalBasis(z:data(), bvalues:data())
   assert_equal(1./math.sqrt(8), bvalues[1], "Checking values 1")
   assert_equal(-math.sqrt(3)/math.sqrt(8), bvalues[2], "Checking values 2")
   assert_equal(-math.sqrt(3)/math.sqrt(8), bvalues[3], "Checking values 3")
   assert_equal(-math.sqrt(3)/math.sqrt(8), bvalues[4], "Checking values 4")
   assert_equal(3/math.sqrt(8), bvalues[5], "Checking values 5")
   assert_equal(3/math.sqrt(8), bvalues[6], "Checking values 6")
   assert_equal(3/math.sqrt(8), bvalues[7], "Checking values 7")
   assert_equal(-math.sqrt(27)/math.sqrt(8), bvalues[8], "Checking values 8")
   assert_equal(math.sqrt(5)/math.sqrt(8), bvalues[9], "Checking values 9")
   assert_equal(-math.sqrt(15)/math.sqrt(8), bvalues[10], "Checking values 10")
   assert_equal(-math.sqrt(15)/math.sqrt(8), bvalues[11], "Checking values 11")
   assert_equal(3.*math.sqrt(5)/math.sqrt(8), bvalues[12], "Checking values 12")
   assert_equal(math.sqrt(5)/math.sqrt(8), bvalues[13], "Checking values 13")
   assert_equal(-math.sqrt(15)/math.sqrt(8), bvalues[14], "Checking values 14")
   assert_equal(-math.sqrt(15)/math.sqrt(8), bvalues[15], "Checking values 15")
   assert_equal(3.*math.sqrt(5)/math.sqrt(8), bvalues[16], "Checking values 16")

   z[1], z[2], z[3] = 1.0, 1.0, 1.0
   basis:evalBasis(z:data(), bvalues:data())
   assert_equal(1./math.sqrt(8), bvalues[1], "Checking values 1")
   assert_equal(math.sqrt(3)/math.sqrt(8), bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(3)/math.sqrt(8), bvalues[3], "Checking values 3")
   assert_equal(math.sqrt(3)/math.sqrt(8), bvalues[4], "Checking values 4")
   assert_equal(3/math.sqrt(8), bvalues[5], "Checking values 5")
   assert_equal(3/math.sqrt(8), bvalues[6], "Checking values 6")
   assert_equal(3/math.sqrt(8), bvalues[7], "Checking values 7")
   assert_equal(math.sqrt(27)/math.sqrt(8), bvalues[8], "Checking values 8")
   assert_equal(math.sqrt(5)/math.sqrt(8), bvalues[9], "Checking values 9")
   assert_equal(math.sqrt(15)/math.sqrt(8), bvalues[10], "Checking values 10")
   assert_equal(math.sqrt(15)/math.sqrt(8), bvalues[11], "Checking values 11")
   assert_equal(3.*math.sqrt(5)/math.sqrt(8), bvalues[12], "Checking values 12")
   assert_equal(math.sqrt(5)/math.sqrt(8), bvalues[13], "Checking values 13")
   assert_equal(math.sqrt(15)/math.sqrt(8), bvalues[14], "Checking values 14")
   assert_equal(math.sqrt(15)/math.sqrt(8), bvalues[15], "Checking values 15")
   assert_equal(3.*math.sqrt(5)/math.sqrt(8), bvalues[16], "Checking values 16")
end

function test_1x1v_gkhybrid_p1()
   local basis = Basis.CartModalGkHybrid {
      cdim = 1,  vdim = 1,
   }

   assert_equal(2, basis:ndim(), "Checking NDIM")
   assert_equal(1, basis:polyOrder(), "Checking polynomial order")
   assert_equal(6, basis:numBasis(), "Checking number of basis functions")

   -- coordinates and values
   local z, bvalues = Lin.Vec(basis:ndim()), Lin.Vec(basis:numBasis())

   z[1], z[2] = 0.0, 0.0
   basis:evalBasis(z:data(), bvalues:data())
   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(0.0, bvalues[2], "Checking values 2")
   assert_equal(0.0, bvalues[3], "Checking values 3")
   assert_equal(0.0, bvalues[4], "Checking values 4")
   assert_equal(-math.sqrt(5)/4, bvalues[5], "Checking values 5")
   assert_equal(0.0, bvalues[6], "Checking values 6")

   z[1], z[2] = -1.0, -1.0
   basis:evalBasis(z:data(), bvalues:data())
   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(-math.sqrt(3)/2, bvalues[2], "Checking values 2")
   assert_equal(-math.sqrt(3)/2, bvalues[3], "Checking values 3")
   assert_equal(3.0/2.0, bvalues[4], "Checking values 4")
   assert_equal(math.sqrt(5)/2, bvalues[5], "Checking values 5")
   assert_equal(-math.sqrt(15)/2, bvalues[6], "Checking values 7")

   z[1], z[2] = 1.0, 1.0
   basis:evalBasis(z:data(), bvalues:data())
   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(math.sqrt(3)/2, bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(3)/2, bvalues[3], "Checking values 3")
   assert_equal(3.0/2.0, bvalues[4], "Checking values 4")
   assert_equal(math.sqrt(5)/2, bvalues[5], "Checking values 5")
   assert_equal(math.sqrt(15)/2, bvalues[6], "Checking values 7")
end

function test_2d_s_p2()

   local basis = Basis.CartModalSerendipity {
      ndim = 2, polyOrder = 2,
   }

   assert_equal(2, basis:ndim(), "Checking NDIM")
   assert_equal(2, basis:polyOrder(), "Checking polynomial order")
   assert_equal(8, basis:numBasis(), "Checking number of basis functions")

   -- coordinates and values
   local z, bvalues = Lin.Vec(basis:ndim()), Lin.Vec(basis:numBasis())

   z[1], z[2] = 0.0, 0.0
   basis:evalBasis(z:data(), bvalues:data())

   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(0.0, bvalues[2], "Checking values 2")
   assert_equal(0.0, bvalues[3], "Checking values 3")
   assert_equal(0.0, bvalues[4], "Checking values 4")
   assert_equal(-math.sqrt(5)/4, bvalues[5], "Checking values 5")
   assert_equal(-math.sqrt(5)/4, bvalues[6], "Checking values 6")
   assert_equal(0.0, bvalues[7], "Checking values 7")
   assert_equal(0.0, bvalues[8], "Checking values 8")

   z[1], z[2] = -1.0, -1.0
   basis:evalBasis(z:data(), bvalues:data())

   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(-math.sqrt(3)/2, bvalues[2], "Checking values 2")
   assert_equal(-math.sqrt(3)/2, bvalues[3], "Checking values 3")
   assert_equal(3.0/2.0, bvalues[4], "Checking values 4")
   assert_equal(math.sqrt(5)/2, bvalues[5], "Checking values 5")
   assert_equal(math.sqrt(5)/2, bvalues[6], "Checking values 6")
   assert_equal(-math.sqrt(15)/2, bvalues[7], "Checking values 7")
   assert_equal(-math.sqrt(15)/2, bvalues[8], "Checking values 8")

   z[1], z[2] = 1.0, 1.0
   basis:evalBasis(z:data(), bvalues:data())

   assert_equal(0.5, bvalues[1], "Checking values 1")
   assert_equal(math.sqrt(3)/2, bvalues[2], "Checking values 2")
   assert_equal(math.sqrt(3)/2, bvalues[3], "Checking values 3")
   assert_equal(3.0/2.0, bvalues[4], "Checking values 4")
   assert_equal(math.sqrt(5)/2, bvalues[5], "Checking values 5")
   assert_equal(math.sqrt(5)/2, bvalues[6], "Checking values 6")
   assert_equal(math.sqrt(15)/2, bvalues[7], "Checking values 7")
   assert_equal(math.sqrt(15)/2, bvalues[8], "Checking values 8")
end

function test_fs_1d_p1()
   local basis = Basis.CartModalSerendipity { ndim = 1, polyOrder = 1 }
   local fIn, fOut = Lin.Vec(basis:numBasis()), Lin.Vec(basis:numBasis())

   for i = 1, #fIn do fIn[i] = 1 end

   basis:flipSign(1, fIn:data(), fOut:data())
   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
   assert_equal(-fIn[2], fOut[2], "Checking if flip-sign worked")
  
end

function test_fs_1d_p2()
   local basis = Basis.CartModalSerendipity { ndim = 1, polyOrder = 2 }
   local fIn, fOut = Lin.Vec(basis:numBasis()), Lin.Vec(basis:numBasis())

   for i = 1, #fIn do fIn[i] = 1 end

   basis:flipSign(1, fIn:data(), fOut:data())
   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
   assert_equal(-fIn[2], fOut[2], "Checking if flip-sign worked")
   assert_equal(fIn[3], fOut[3], "Checking if flip-sign worked")
  
end

function test_fs_s_2d_p1()
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 1 }
   local fIn, fOut = Lin.Vec(basis:numBasis()), Lin.Vec(basis:numBasis())

   for i = 1, #fIn do fIn[i] = 1 end

   basis:flipSign(1, fIn:data(), fOut:data())
   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
   assert_equal(-fIn[2], fOut[2], "Checking if flip-sign worked")
   assert_equal(fIn[3], fOut[3], "Checking if flip-sign worked")   
   assert_equal(-fIn[4], fOut[4], "Checking if flip-sign worked")

   basis:flipSign(2, fIn:data(), fOut:data())
   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
   assert_equal(fIn[2], fOut[2], "Checking if flip-sign worked")
   assert_equal(-fIn[3], fOut[3], "Checking if flip-sign worked")
   assert_equal(-fIn[4], fOut[4], "Checking if flip-sign worked")
end

function test_fs_s_2d_p2()
   local basis = Basis.CartModalSerendipity { ndim = 2, polyOrder = 2 }
   local fIn, fOut = Lin.Vec(basis:numBasis()), Lin.Vec(basis:numBasis())

   for i = 1, #fIn do fIn[i] = 1 end

   basis:flipSign(1, fIn:data(), fOut:data())
   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
   assert_equal(-fIn[2], fOut[2], "Checking if flip-sign worked")
   assert_equal(fIn[3], fOut[3], "Checking if flip-sign worked")   
   assert_equal(-fIn[4], fOut[4], "Checking if flip-sign worked")
   assert_equal(fIn[5], fOut[5], "Checking if flip-sign worked")
   assert_equal(fIn[6], fOut[6], "Checking if flip-sign worked")
   assert_equal(fIn[7], fOut[7], "Checking if flip-sign worked")
   assert_equal(-fIn[8], fOut[8], "Checking if flip-sign worked")

   basis:flipSign(2, fIn:data(), fOut:data())
   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
   assert_equal(fIn[2], fOut[2], "Checking if flip-sign worked")
   assert_equal(-fIn[3], fOut[3], "Checking if flip-sign worked")   
   assert_equal(-fIn[4], fOut[4], "Checking if flip-sign worked")
   assert_equal(fIn[5], fOut[5], "Checking if flip-sign worked")
   assert_equal(fIn[6], fOut[6], "Checking if flip-sign worked")
   assert_equal(-fIn[7], fOut[7], "Checking if flip-sign worked")
   assert_equal(fIn[8], fOut[8], "Checking if flip-sign worked")
end

function test_fs_hybrid_1x2v_p1()
   local basis = Basis.CartModalHybrid { cdim = 1, vdim = 2, }
   local fIn, fOut = Lin.Vec(basis:numBasis()), Lin.Vec(basis:numBasis())

   for i = 1, #fIn do fIn[i] = 1 end

   basis:flipSign(1, fIn:data(), fOut:data())
   assert_equal( fIn[1], fOut[1], "Checking if flip-sign worked")
   assert_equal(-fIn[2], fOut[2], "Checking if flip-sign worked")
   assert_equal( fIn[3], fOut[3], "Checking if flip-sign worked")   
   assert_equal( fIn[4], fOut[4], "Checking if flip-sign worked")
   assert_equal(-fIn[5], fOut[5], "Checking if flip-sign worked")   
   assert_equal(-fIn[6], fOut[6], "Checking if flip-sign worked")
   assert_equal( fIn[7], fOut[7], "Checking if flip-sign worked")
   assert_equal(-fIn[8], fOut[8], "Checking if flip-sign worked")
   assert_equal( fIn[9], fOut[9], "Checking if flip-sign worked")   
   assert_equal(-fIn[10], fOut[10], "Checking if flip-sign worked")
   assert_equal( fIn[11], fOut[11], "Checking if flip-sign worked")   
   assert_equal(-fIn[12], fOut[12], "Checking if flip-sign worked")
   assert_equal( fIn[13], fOut[13], "Checking if flip-sign worked")
   assert_equal(-fIn[14], fOut[14], "Checking if flip-sign worked")
   assert_equal( fIn[15], fOut[15], "Checking if flip-sign worked")
   assert_equal(-fIn[16], fOut[16], "Checking if flip-sign worked")

   basis:flipSign(2, fIn:data(), fOut:data())
   assert_equal( fIn[1], fOut[1], "Checking if flip-sign worked")
   assert_equal( fIn[2], fOut[2], "Checking if flip-sign worked")
   assert_equal(-fIn[3], fOut[3], "Checking if flip-sign worked")
   assert_equal( fIn[4], fOut[4], "Checking if flip-sign worked")
   assert_equal(-fIn[5], fOut[5], "Checking if flip-sign worked")
   assert_equal( fIn[6], fOut[6], "Checking if flip-sign worked")
   assert_equal(-fIn[7], fOut[7], "Checking if flip-sign worked")
   assert_equal(-fIn[8], fOut[8], "Checking if flip-sign worked")
   assert_equal( fIn[9], fOut[9], "Checking if flip-sign worked")
   assert_equal( fIn[10], fOut[10], "Checking if flip-sign worked")
   assert_equal( fIn[11], fOut[11], "Checking if flip-sign worked")
   assert_equal( fIn[12], fOut[12], "Checking if flip-sign worked")
   assert_equal( fIn[13], fOut[13], "Checking if flip-sign worked")
   assert_equal( fIn[14], fOut[14], "Checking if flip-sign worked")
   assert_equal(-fIn[15], fOut[15], "Checking if flip-sign worked")
   assert_equal(-fIn[16], fOut[16], "Checking if flip-sign worked")
end

function test_fs_gkhybrid_1x1v_p1()
   local basis = Basis.CartModalGkHybrid {
      cdim = 1,  vdim = 1,
   }
   local fIn, fOut = Lin.Vec(basis:numBasis()), Lin.Vec(basis:numBasis())

   for i = 1, #fIn do fIn[i] = 1 end

   basis:flipSign(1, fIn:data(), fOut:data())
   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
   assert_equal(-fIn[2], fOut[2], "Checking if flip-sign worked")
   assert_equal(fIn[3], fOut[3], "Checking if flip-sign worked")   
   assert_equal(-fIn[4], fOut[4], "Checking if flip-sign worked")
   assert_equal(fIn[5], fOut[5], "Checking if flip-sign worked")   
   assert_equal(-fIn[6], fOut[6], "Checking if flip-sign worked")

   basis:flipSign(2, fIn:data(), fOut:data())
   assert_equal(fIn[1], fOut[1], "Checking if flip-sign worked")
   assert_equal(fIn[2], fOut[2], "Checking if flip-sign worked")
   assert_equal(-fIn[3], fOut[3], "Checking if flip-sign worked")
   assert_equal(-fIn[4], fOut[4], "Checking if flip-sign worked")
   assert_equal(fIn[5], fOut[5], "Checking if flip-sign worked")
   assert_equal(fIn[6], fOut[6], "Checking if flip-sign worked")
end

---- Disabling tests of maximal order basis (not kept in g0).
--test_2d_m_p1()
--test_2d_m_p2()
--
--test_fs_m_2d_p1()
--test_fs_m_2d_p2()

test_2d_s_p1()
test_2d_s_p2()

test_fs_1d_p1()
test_fs_1d_p2()

test_fs_s_2d_p1()
test_fs_s_2d_p2()

test_1x2v_hybrid_p1()
test_fs_hybrid_1x2v_p1()

test_1x1v_gkhybrid_p1()
test_fs_gkhybrid_1x1v_p1()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
