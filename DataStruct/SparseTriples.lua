-- Gkyl ------------------------------------------------------------------------
--
-- Triples of the form (i,j,val). Mainly for use in filling up sparse
-- matrices for linear system solvers.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")


-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Range = require "Lib.Range"

-- SparseTriples ---------------------------------------------------------------
--
-- Triples of the form (i,j,val). Mainly for use in filling up sparse
-- matrices for linear system solvers.
--------------------------------------------------------------------------------

local SparseTriples = {}
function SparseTriples:new(tbl)
   local self = setmetatable({}, SparseTriples)

   self._nRows = tbl[1]
   self._nCols = tbl[2]

   local allocator = Alloc.createAllocator("struct { int i ,j; double val; }")
   self._triples = allocator() -- triples
   self._tmpTriple = new("struct { int i ,j; double val; }") -- temp storage for single entry

   return self
end
setmetatable(SparseTriples, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
SparseTriples.__index = {
   numRows = function(self)
      return self._nRows
   end,
   numColumns = function (self)
      return self._nCols
   end,
   size = function (self)
      return self._triples:size()
   end,
   append = function(self, i, j, val)
      self._tmpTriple.i, self._tmpTriple.j, self._tmpTriple.val = i, j, val
      self._triples:push(self._tmpTriple)
   end,
   g = function (self, n)
      return self._triples[n]
   end,
}

return {SparseTriples = SparseTriples}


