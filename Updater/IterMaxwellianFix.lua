-- Gkyl ------------------------------------------------------------------------
--
-- Iterative fix for Maxwell projection
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local DataStruct = require "DataStruct"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Logger = require "Lib.Logger"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"
local xsys = require "xsys"

local ffi  = require "ffi"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Iterative Poisson solver updater object
local IterMaxwellianFix = Proto(UpdaterBase)

-- constructor
function IterMaxwellianFix:init(tbl)
   IterMaxwellianFix.super.init(self, tbl) -- setup base object
   
   -- read data from input file
   self.onGrid = assert(
      tbl.onGrid, "Updater.IterMaxwellianFix: Must provide grid object using 'onGrid'")
   self.basis = assert(
      tbl.basis, "Updater.IterMaxwellianFix: Must specify basis functions to use using 'basis'")

   local polyOrder = self.basis:polyOrder()

end

-- advance method
function IterMaxwellianFix:_advance(tCurr, inFld, outFld)

   local grid, basis = self.onGrid, self.basis
   local polyOrder = basis:numBasis()

   local fIn = assert(inFld[1], "IterMaxwellianFix.advance: Must-specify an input distribution function")
   local fOut = assert(outFld[1], "IterMaxwellianFix.advance: Must-specify an output distribution function")

end

return IterMaxwellianFix
