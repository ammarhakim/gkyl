-- Gkyl ------------------------------------------------------------------------
--
-- Decompose N elements equitably amongst a given number of pieces
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local _M = {}

-- LinearDecomp ---------------------------------------------------------------
--
--  Decompose N elements equitably amongst a given number of pieces
--------------------------------------------------------------------------------

local LinearDecomp = {}
-- constructor to make object that stores decomposed ranges
function LinearDecomp:new(tbl)
   local self = setmetatable({}, LinearDecomp)

   local domSize, numSplit = tbl.domSize, tbl.numSplit
   self._domSize, self._numSplit  = tbl.domSize, tbl.numSplit

   -- calculate local starting index and shape
   local baseShape = math.floor(domSize/numSplit)
   local remCells = domSize % numSplit -- extra cells left over   
   self._shape, self._start = {}, {}
   for c = 1, numSplit do
      self._shape[c] = baseShape + (remCells>0 and 1 or 0) -- add extra cell if any remain
      remCells = remCells-1
   end
   self._start[1] = 1
   for c = 2, numSplit do
      self._start[c] = self._start[c-1]+self._shape[c-1]
   end
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(LinearDecomp, { __call = function (self, o) return self.new(self, o) end })

LinearDecomp.__index = {
   domSize = function (self)
      return self._domSize
   end,
   numSplit = function (self)
      return self._numSplit
   end,
   lower = function (self, n)
      return self._start[n]
   end,
   shape = function (self, n)
      return self._shape[n]
   end,
   upper = function (self, n)
      return self._start[n]+self._shape[n]-1
   end,
}

return {
   LinearDecomp = LinearDecomp
}
