-- Gkyl ------------------------------------------------------------------------
--
-- Test for prototype based object system
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local Unit = require "Unit"
local Proto = require "Lib.Proto"

local assert_equal = Unit.assert_equal
local stats = Unit.stats

function test_1()
   -- define a object and its methods
   local Particle = Proto()
   Particle.w = 1.5
   
   function Particle:init(x,y)
      self.x, self.y = x, y
   end
   function Particle:setw(w)
      self.w = w
   end
   function Particle:getx() return self.x end
   function Particle:gety() return self.y end
   function Particle:getw() return self.w end

   local p = Particle(2.5, 3.5)
   assert_equal(2.5, p:getx(), "Checking object")
   assert_equal(3.5, p:gety(), "Checking object")
   assert_equal(1.5, p:getw(), "Checking object")

   p:setw(10.5)
   assert_equal(10.5, p:getw(), "Checking object")

   -- derived object
   local SpecialParticle = Proto(Particle)

   function SpecialParticle:init(x, y, m)
      SpecialParticle.super.init(self, x, y)
      self.m = m
   end
   function SpecialParticle:distance()
      return math.sqrt(self.x^2 + self.y^2)
   end

   local sp = SpecialParticle(5, 10, 101.5)
   assert_equal(5.0, sp:getx(), "Checking object")
   assert_equal(10.0, sp:gety(), "Checking object")
   assert_equal(1.5, sp:getw(), "Checking object")
   assert_equal(101.5, sp.m, "Checking object")
   assert_equal(math.sqrt(125), sp:distance(), "Checking object")
end

function test_2()
   local Vars = Proto()

   function Vars:init(a,b,c,d)
      self.vals = {a,b,c,d}
   end
   function Vars:get(k)
      return self.vals[k]
   end

   Vars.__index = function(self, key)
      -- check base-class dispatch first
      local metavalue = getmetatable(self)[key]
      if metavalue then return metavalue end
      if type(key) == 'number' then
	 return self.vals[key]
      end
      return rawget(self, key)      
   end
   Vars.__newindex = function(self, k, v)
      if type(k) == "number" then
	 self.vals[k] = v
      end
      rawset(self, k, v)
   end

   local v = Vars(1,2,3,4)
   for i = 1, 4 do
      assert_equal(i, v[i], "Checking index access")
   end
   for i = 1, 4 do
      assert_equal(i, v:get(i), "Checking index access")
   end
   
   v[1] = 10.5
   assert_equal(10.5, v[1], "Checking reset index access")

end

function test_3()
   local Rational = Proto()
   function Rational:init(num, den)
      self.num, self.den = num, den
   end
   Rational.__add = function (a,b)
      local n1, d1 = a.num, a.den
      local n2, d2 = b.num, b.den
      return Rational(n1*d2+n2*d1, d1*d2)
   end

   local r1, r2 = Rational(2,1), Rational(3,1)
   local r = r1+r2

   assert_equal(r.num, 5, "Checking rational add method")
   assert_equal(r.den, 1, "Checking rational add method")
end

function test_4()
   local Base = Proto()
   function Base:init(nm)
      self.name = nm
   end
   function Base:subAdvance(t)
   end
   
   function Base:advance(t)
      print(string.format("Calling derived advance from %s", self.name))
      self:subAdvance(t)
      print("Done!")
   end

   local Wave = Proto(Base)
   function Wave:subAdvance(t)
      print("Wave advance called with ", t)
   end

   local w = Wave("Wave")
   w:advance(1.5)
end

test_1()
test_2()
test_3()
test_4()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
