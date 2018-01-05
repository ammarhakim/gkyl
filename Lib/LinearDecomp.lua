-- Gkyl ------------------------------------------------------------------------
--
-- Decompose N elements equitably amongst a given number of "threads"
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- gkyl libraries
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"

-- LinearDecomp ----------------------------------------------------------------
--
--  Decompose N elements equitably amongst a given number of "threads"
--------------------------------------------------------------------------------

local LinearDecomp = Proto()

-- constructor to make object that stores decomposed ranges
function LinearDecomp:init(tbl)
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
end

function LinearDecomp:domSize() return self._domSize end
function LinearDecomp:numSplit() return self._numSplit end
function LinearDecomp:lower(n) return self._start[n] end
function LinearDecomp:upper(n) return self._start[n]+self._shape[n]-1 end
function LinearDecomp:shape(n) return self._shape[n] end

-- LinearDecompRange -----------------------------------------------------------
--
--  Decompose a range object using a linear decomposition
--------------------------------------------------------------------------------

local LinearDecompRange = Proto()

function LinearDecompRange:init(tbl)
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
end

function LinearDecompRange:domSize() return self._linearDecomp:domSize() end
function LinearDecompRange:numSplit() return self._linearDecomp:numSplit() end
function LinearDecompRange:lower(n) return self._linearDecomp:lower(n) end
function LinearDecompRange:upper(n) return self._linearDecomp:upper(n) end
function LinearDecompRange:shape(n) return self._linearDecomp:shape(n) end
function LinearDecompRange:rowStartIndex(n) return self._rowStartIdx[n] end
function LinearDecompRange:colStartIndex(n) return self._colStartIdx[n] end

return {
   LinearDecomp = LinearDecomp,
   LinearDecompRange = LinearDecompRange
}
