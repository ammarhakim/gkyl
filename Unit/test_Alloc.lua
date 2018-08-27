-- Gkyl ------------------------------------------------------------------------
--
-- Test for memory allocators
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local complex = require "sci.complex"
local ffi = require "ffi"
local Unit = require "Unit"
local Alloc = require "Lib.Alloc"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

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

   eulerFld:expand(0) -- zap it, without deallocating
   eulerFld:push(ffi.new(eulerFld:elemType(), {rho = 10.5, rhou = 11.5, E = 12.5}))
   assert_equal(10.5, eulerFld[1].rho, "Checking rho")
   assert_equal(11.5, eulerFld[1].rhou, "Checking rhou")
   assert_equal(12.5, eulerFld[1].E, "Checking E")
   
end

function test_6()
   local dv = Alloc.Double()

   for i = 1, 10 do
      dv:push(i+0.5)
   end
   local l = dv:popLast()
   assert_equal(10.5, l, "Checking pop-ed value")
   local nl = dv:last()
   assert_equal(9.5, nl, "Checking new last value")
   assert_equal(9, dv:size(), "Checking size")
end

function test_7()
   local allocator = Alloc.createAllocator("struct {int i ,j; double val; }")
   local tmpData = ffi.new("struct {int i ,j; double val; }")
   tmpData.i = 1
   tmpData.j = 2
   tmpData.val = 100.5

   local data = allocator(10)
   data:fill(tmpData)

   for i = 1, data:size() do
      local v = data[i]
      assert_equal(tmpData.i, v.i, "i")
      assert_equal(tmpData.j, v.j, "j")
      assert_equal(tmpData.val, v.val, "val") 
   end
end

function test_8()
   local carr = Alloc.Complex(10)

   for i = 0, carr:size() do
      carr[i] = complex.new(i, 0.1)
   end

   for i = 0, carr:size() do
      assert_equal(math.sqrt(i^2+0.1^2), complex.abs(carr[i]), "Complex amplitude")
   end
end

-- run tests
test_1()
test_2()
test_2b()
test_3()
test_4()
test_5()
test_6()
test_7()
test_8()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
