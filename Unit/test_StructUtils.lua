-- Gkyl ------------------------------------------------------------------------
--
-- Test for range objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local StructUtils = require "Lib.StructUtils"
local Range = require "Lib.Range"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

ffi.cdef [[
  typedef struct {
    double x, y, z;
  } Particle_t;
 
  /* To illustrate behaviour of pointers */
  typedef struct {
    Particle_t *p1, *p2;
  } TwoParticles_t;

  void unit_showParticle(Particle_t *ptcl);
]]

-- constructors for structs
local Particle = StructUtils.Struct("Particle_t")
local TwoParticles = StructUtils.Struct("TwoParticles_t")

function test_1()
   local p1 = Particle {1.5, 2.5, 3.5 }
   assert_equal(1.5, p1.x, "Checking struct x")
   assert_equal(2.5, p1.y, "Checking struct y")
   assert_equal(3.5, p1.z, "Checking struct z")

   local p2 = p1:clone()
   assert_equal(1.5, p2.x, "Checking struct clone x")
   assert_equal(2.5, p2.y, "Checking struct clone y")
   assert_equal(3.5, p2.z, "Checking struct clone z")

   p2.x = 10.5
   assert_equal(10.5, p2.x, "Checking struct clone x")

   assert_equal(1.5, p1.x, "Checking struct x")

   p1:copy(p2)
   assert_equal(10.5, p1.x, "Checking struct x")
   assert_equal(2.5, p1.y, "Checking struct y")
   assert_equal(3.5, p1.z, "Checking struct z")
end

function test_2()
   local pa = Particle {1.5, 2.5, 3.5 }
   local pb = Particle {10.5, 20.5, 30.5 }

   local twoP = TwoParticles { pa, pb }

   assert_equal(1.5, twoP.p1.x, "Checking pointers")
   assert_equal(2.5, twoP.p1.y, "Checking pointers")
   assert_equal(3.5, twoP.p1.z, "Checking pointers")

   assert_equal(10.5, twoP.p2.x, "Checking pointers")
   assert_equal(20.5, twoP.p2.y, "Checking pointers")
   assert_equal(30.5, twoP.p2.z, "Checking pointers")

   local cloneTwoP = twoP:clone()
   assert_equal(1.5, cloneTwoP.p1.x, "Checking pointers")
   assert_equal(2.5, cloneTwoP.p1.y, "Checking pointers")
   assert_equal(3.5, cloneTwoP.p1.z, "Checking pointers")

   assert_equal(10.5, cloneTwoP.p2.x, "Checking pointers")
   assert_equal(20.5, cloneTwoP.p2.y, "Checking pointers")
   assert_equal(30.5, cloneTwoP.p2.z, "Checking pointers")

   -- now change pa and pb
   pa.x, pa.y, pa.z = 100.5, 200.5, 300.5
   pb.x, pb.y, pb.z = 1000.5, 2000.5, 3000.5

   -- change is reflected in both twoP as well as cloneTwoP
   assert_equal(100.5, twoP.p1.x, "Checking pointers")
   assert_equal(200.5, twoP.p1.y, "Checking pointers")
   assert_equal(300.5, twoP.p1.z, "Checking pointers")

   assert_equal(1000.5, twoP.p2.x, "Checking pointers")
   assert_equal(2000.5, twoP.p2.y, "Checking pointers")
   assert_equal(3000.5, twoP.p2.z, "Checking pointers")

   assert_equal(100.5, cloneTwoP.p1.x, "Checking pointers")
   assert_equal(200.5, cloneTwoP.p1.y, "Checking pointers")
   assert_equal(300.5, cloneTwoP.p1.z, "Checking pointers")

   assert_equal(1000.5, cloneTwoP.p2.x, "Checking pointers")
   assert_equal(2000.5, cloneTwoP.p2.y, "Checking pointers")
   assert_equal(3000.5, cloneTwoP.p2.z, "Checking pointers")
end

function test_3()
   if not GKYL_HAVE_CUDA then return end
   local pa = Particle {1.5, 2.5, 3.5 }
   local paDev = pa:cloneOnDevice()
   ffi.C.unit_showParticle(paDev)
end

test_1()
test_2()
if GKYL_HAVE_CUDA then
   test_3()
end

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
