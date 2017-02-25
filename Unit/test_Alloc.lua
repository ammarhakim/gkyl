-- Gkyl ------------------------------------------------------------------------
--
-- Test for memory allocators
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local Unit = require "Unit"
local Alloc = require "Lib.Alloc"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local len = 100
   local da = Alloc.Double(len)
   
   assert_equal(len, da:size(), "Checking size")
   da:fill(1.0)
   for i = 1, da:size() do
      assert_equal(1.0, da[i], "Checking fill")
   end
   assert_equal(1.0, da:last(), "Checking last added element")

   for i = 1, da:size() do
      da[i] = (i+0.5)*0.1
   end
   for i = 1, da:size() do
      assert_equal((i+0.5)*0.1, da[i], "Checking [] operator")
   end
   assert_equal((da:size()+0.5)*0.1, da:last(), "Checking last added element")   
end

function test_2()
   local len = 100
   local eulerAlloc = Alloc.createAllocator("struct {double rho, rhou, E;}")
   local eulerFld = eulerAlloc(len)

   assert_equal(len, eulerFld:size(), "Checking size")

   for i = 1, eulerFld:size() do
      eulerFld[i].rho = (i+0.5)*0.1
      eulerFld[i].rhou = (i+0.5)*0.2
      eulerFld[i].E = (i+0.5)*0.3
   end
   do
      local v = eulerFld:last()
      local i = eulerFld:size()
      assert_equal((i+0.5)*0.1, v.rho, "Checking [] operator")
      assert_equal((i+0.5)*0.2, v.rhou, "Checking [] operator")
      assert_equal((i+0.5)*0.3, v.E, "Checking [] operator")
   end

   for i = 1, eulerFld:size() do
      assert_equal((i+0.5)*0.1, eulerFld[i].rho, "Checking [] operator")
      assert_equal((i+0.5)*0.2, eulerFld[i].rhou, "Checking [] operator")
      assert_equal((i+0.5)*0.3, eulerFld[i].E, "Checking [] operator")
   end   
end

function test_2b()
   local len = 100
   local eulerAlloc = Alloc.createAllocator("struct {double rho, rhou, E;}")
   local eulerFld = eulerAlloc(len)

   assert_equal(len, eulerFld:size(), "Checking size")

   eFld = ffi.new(eulerFld:elemType())
   for i = 1, eulerFld:size() do
      eFld.rho = (i+0.5)*0.1
      eFld.rhou = (i+0.5)*0.2
      eFld.E = (i+0.5)*0.3
      eulerFld[i] = eFld
   end
   do
      local v = eulerFld:last()
      local i = eulerFld:size()
      assert_equal((i+0.5)*0.1, v.rho, "Checking [] operator")
      assert_equal((i+0.5)*0.2, v.rhou, "Checking [] operator")
      assert_equal((i+0.5)*0.3, v.E, "Checking [] operator")
   end

   for i = 1, eulerFld:size() do
      assert_equal((i+0.5)*0.1, eulerFld[i].rho, "Checking [] operator")
      assert_equal((i+0.5)*0.2, eulerFld[i].rhou, "Checking [] operator")
      assert_equal((i+0.5)*0.3, eulerFld[i].E, "Checking [] operator")
   end   
end

function test_3()
   local da = Alloc.Double()
   assert_equal(0, da:size(), "Testing size")
   da:expand(10)
   assert_equal(10, da:size(), "Testing size")
   local oldCap = da:capacity()
   da:expand(da:capacity()+1)
   assert_equal(oldCap+1, da:size(), "Testing size")
end

function test_4()
   local da = Alloc.Double()
   assert_equal(0, da:size(), "Testing size")

   da:push(1)
   assert_equal(1, da:size(), "Testing size")
   assert_equal(1, da[1], "Testing push")
   assert_equal(1, da:last(), "Testing last")

   for i = 2, 100 do da:push(i) end
   assert_equal(100, da:size(), "Testing size")
   
   for i = 1, da:size() do
      assert_equal(i, da[i], "Testing post-push")
   end

   da:clear()
   assert_equal(0, da:size(), "Testing size")
end

function test_5()
   local eulerAlloc = Alloc.createAllocator("struct {double rho, rhou, E;}")
   local eulerFld = eulerAlloc()

   assert_equal(0, eulerFld:size(), "Checking size")
   local eFld = ffi.new(eulerFld:elemType())
   for i = 1, 100 do
      eFld.rho = (i+0.5)*1
      eFld.rhou = (i+0.5)*2
      eFld.E = (i+0.5)*3
      eulerFld:push(eFld)
   end

   assert_equal(100, eulerFld:size(), "Checking size")
   for i = 1, eulerFld:size() do
      assert_equal((i+0.5)*1, eulerFld[i].rho, "Checking rho")
      assert_equal((i+0.5)*2, eulerFld[i].rhou, "Checking rhou")
      assert_equal((i+0.5)*3, eulerFld[i].E, "Checking E")
   end   
   
end

-- run tests
test_1()
test_2()
test_2b()
test_3()
test_4()
test_5()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
