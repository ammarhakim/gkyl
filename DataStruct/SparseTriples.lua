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
local Proto = require "Lib.Proto"

-- SparseTriples ---------------------------------------------------------------
--
-- Triples of the form (i,j,val). Mainly for use in filling up sparse
-- matrices for linear system solvers.
--------------------------------------------------------------------------------

local SparseTriples = Proto()

function SparseTriples:init(tbl)
   self._nRows = tbl[1]
   self._nCols = tbl[2]

   local allocator = Alloc.createAllocator("struct { int i ,j; double val; }")
   self._triples = allocator() -- triples
   self._tmpTriple = new("struct { int i ,j; double val; }") -- temp storage for single entry
end

function SparseTriples:numRows() return self._nRows end
function SparseTriples:numColumns() return self._nCols end
function SparseTriples:size() return self._triples:size() end
function SparseTriples:get(k) return self._triples[k] end

function SparseTriples:append(i, j, val)
   self._tmpTriple.i, self._tmpTriple.j, self._tmpTriple.val = i, j, val
   self._triples:push(self._tmpTriple)
end

return {SparseTriples = SparseTriples}


