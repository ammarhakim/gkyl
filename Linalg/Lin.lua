-- Gkyl ------------------------------------------------------------------------
--
-- Vectors, matrices and linear algebra
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Vector: this is a simple slice of memory with length stored. Can
-- create vector of structures (structs can not be VLS or VLA)

-- Allocate vector with type "ct" and size "n"
local function vec_alloc(ct, n)
  local v = new(ct, n)
  v._n, v._p = n, v._a
  return v
end

-- copy vector "x" into vector "y"
local function vec_memcpy(y, x)
  copy(y._p, x._p, sizeof(y:elemType())*y._n)
end

-- Vector meta-object
local function new_vec_ct(elct)
   local vec_mf = {
      elemType = function(self)	 
	 return elct
      end,
      data = function(self)
	 return self._p
      end,
      copy = function(self)
	 local v = vec_alloc(self, self._n)
	 vec_memcpy(v, self)
	 return v
      end,
   }
  local vec_mt = {
    __new = function(ct, n)
      return vec_alloc(ct, n)
    end,
    __len = function(self)
      return self._n
    end,
    __index = function(self, k)
       if type(k) == "number" then
	  return self._p[k]
       else
	  return vec_mf[k]
       end
    end,
    __newindex = function(self, i, v)
      self._p[i] = v 
    end,
  }
  local ct = typeof("struct { int32_t _n; $* _p; $ _a[?]; }", elct, elct)
  return metatype(ct, vec_mt)
end

return {
   new_vec_ct = new_vec_ct,   
   vec = new_vec_ct(typeof("double")),
}
