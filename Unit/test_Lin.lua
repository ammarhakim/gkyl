-- Gkyl
-- ------------------------------------------------------------------------
--
-- Test for linear algebra objects
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Lin = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   local v = Lin.Vec(3)
   assert_equal(3, #v, "Checking length of vector")
   -- set values
   for i= 1, #v do
      v[i] = (i+0.5)*0.1
   end
   -- test them
   for i= 1, #v do
      assert_equal((i+0.5)*0.1, v[i], "Checking vector content")
   end

   local vcopy = v:copy()
   -- test copy
   for i = 1, #vcopy do
      assert_equal((i+0.5)*0.1, vcopy[i], "Checking vector copy")
   end
end

function test_2()
   local EulerVec = Lin.new_vec_ct(ffi.typeof("struct {double rho, rhou, E;}"))
   local v = EulerVec(4)
   assert_equal(4, #v, "Checking length of vector")
   -- set values
   for i = 1, #v do
      v[i].rho = i+0.0
      v[i].rhou = i+1.0
      v[i].E = i+2.0 
   end
   -- test them
   for i = 1, #v do
      assert_equal(i+0.0, v[i].rho, "Checking vector of struct contents")
      assert_equal(i+1.0, v[i].rhou, "Checking vector of struct contents")
      assert_equal(i+2.0, v[i].E, "Checking vector of struct contents")
   end

   local vcopy = v:copy()
   -- test copy
   for i = 1, #vcopy do
      assert_equal(i+0.0, vcopy[i].rho, "Checking vector of struct contents copy")
      assert_equal(i+1.0, vcopy[i].rhou, "Checking vector of struct contents copy")
      assert_equal(i+2.0, vcopy[i].E, "Checking vector of struct contents copy")
   end
end

function test_3()
   local m = Lin.Mat(3, 4)
   assert_equal(3, m:numRows(), "Num rows")
   assert_equal(4, m:numCols(), "Num cols")

   -- set values
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do   
	 m[i][j] = j+(i-1)*m:numCols()
      end
   end

   -- check them
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do
	 assert_equal(j+(i-1)*m:numCols(), m[i][j], "Checking matrix values")
      end
   end

   -- check raw data
   local count, data = 1, m:data()
   for i = 1, m:numCols()*m:numRows() do
      assert_equal(count, data[i-1], "Checking linear data")
      count = count+1
   end
end

function test_4()
   local m = Lin.Mat(3, 4)

   local count, w = 1, m[1]
   for i = 1, m:numCols() do
      w[i] = count
      count = count+1
   end

   w = m[2]
   for i = 1, m:numCols() do
      w[i] = count
      count = count+1
   end

   w = m[3]
   for i = 1, m:numCols() do
      w[i] = count
      count = count+1
   end

   -- check them
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do
	 assert_equal(j+(i-1)*m:numCols(), m[i][j], "Checking matrix values")
      end
   end
end

function test_5()
   local _data = ffi.new("double[?]", 3*4)
   local m = Lin.Mat(3, 4, _data)

   -- set values
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do   
	 m[i][j] = j+(i-1)*m:numCols()
      end
   end

   -- check them
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do
	 assert_equal(j+(i-1)*m:numCols(), m[i][j], "Checking matrix values")
      end
   end

   -- check raw data
   local count, data = 1, m:data()
   for i = 1, m:numCols()*m:numRows() do
      assert_equal(count, data[i-1], "Checking linear data")
      count = count+1
   end   
end

function test_6()
   local nelem, stride = 10, 4
   local data = ffi.new("double[?]", nelem*stride)
   local chunck = ffi.new("double[4]", {1, 2, 3, 4})

   local count = 0
   for i = 1, nelem do
      ffi.copy(data+count*stride, chunck, ffi.sizeof("double")*stride)
      count = count+1
   end

   -- check this
   count = 0
   for i = 1, nelem do
      assert_equal(1, data[count*stride+0], "Checking strided copy")
      assert_equal(2, data[count*stride+1], "Checking strided copy")
      assert_equal(3, data[count*stride+2], "Checking strided copy")
      assert_equal(4, data[count*stride+3], "Checking strided copy")
      count = count + 1
   end
end

function test_7()
   local m = Lin.Mat(3, 4)

   -- set values
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do   
	 m[i][j] = j+(i-1)*m:numCols()
      end
   end

   -- check them
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do
	 assert_equal(j+(i-1)*m:numCols(), m:g(i,j), "Checking matrix values")
      end
   end

   -- set values
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do   
	 m:s(i, j, j+(i-1)*m:numCols())
      end
   end

   -- check them
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do
	 assert_equal(j+(i-1)*m:numCols(), m:g(i,j), "Checking matrix values")
      end
   end   
end

-- Run tests
test_1()
test_2()
test_3()
test_4()
test_5()
test_6()
test_7()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
