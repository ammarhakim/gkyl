-- Gkyl ------------------------------------------------------------------------
--
-- Decompose N elements equitably amongst a given number of "threads"
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- gkyl libraries
local Range = require "Lib.Range"
local Lin = require "Lib.Linalg"

-- LinearDecomp ----------------------------------------------------------------
--
--  Decompose N elements equitably amongst a given number of "threads"
--------------------------------------------------------------------------------

local LinearDecomp = {}
-- constructor to make object that stores decomposed ranges
function LinearDecomp:new(tbl)
   local self = setmetatable({}, LinearDecomp)

   local domSize, numSplit = tbl.domSize, tbl.numSplit
   self._domSize, self._numSplit  = tbl.domSize, tbl.numSplit

   -- calculate starting index and shape for each thread
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

-- LinearDecompRange -----------------------------------------------------------
--
--  Decompose a range object using a linear decomposition
--------------------------------------------------------------------------------

local LinearDecompRange = {}
-- constructor to make object that stores decomposed ranges
function LinearDecompRange:new(tbl)
   local self = setmetatable({}, LinearDecompRange)

   local r = tbl.range -- range to split   
   -- create linear decomp object
   self._linearDecomp = LinearDecomp {
      domSize = r:volume(), numSplit = tbl.numSplit
   }

   -- inverse indexers to map linear locations into ndim() indices
   local colInvIndexer, rowInvIndexer = Range.makeColMajorInvIndexer(r), Range.makeRowMajorInvIndexer(r)

   -- construct start indices for both column and row major layouts
   self._colStartIdx = {}
   for d = 1, tbl.numSplit do
      local idx = Lin.IntVec(r:ndim())
      colInvIndexer(self._linearDecomp:lower(d), idx)
      self._colStartIdx[d] = idx
   end

   self._rowStartIdx = {}
   for d = 1, tbl.numSplit do
      local idx = Lin.IntVec(r:ndim())
      rowInvIndexer(self._linearDecomp:lower(d), idx)
      self._rowStartIdx[d] = idx
   end   

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(LinearDecompRange, { __call = function (self, o) return self.new(self, o) end })

LinearDecompRange.__index = {
   domSize = function (self)
      return self._linearDecomp:domSize()
   end,
   numSplit = function (self)
      return self._linearDecomp:numSplit()
   end,
   shape = function (self, n)
      return self._linearDecomp:shape(n)
   end,
   rowStartIndex = function (self, n)
      return self._rowStartIdx[n]
   end,
   colStartIndex = function (self, n)
      return self._colStartIdx[n]
   end,   
}

return {
   LinearDecomp = LinearDecomp,
   LinearDecompRange = LinearDecompRange
}
