-- Gkyl ------------------------------------------------------------------------
--
-- Updater to project function on basis functions. Uses Gaussian
-- quadrature. The projection is not exact unless the function being
-- projected is a polynomial of order less than the basis function
-- polyOrder.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Base = require "Updater.Base"
local GaussQuadRules = require "Lib.GaussQuadRules"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Projection  updater object
local ProjectOnBasis = {}

function ProjectOnBasis:new(tbl)
   local self = setmetatable({}, ProjectOnBasis)
   Base.setup(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid, "Updater.ProjectOnBasis: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.onGrid, "Updater.Provide: Must specify basis functions to use using 'basis'")

   local N = self._basis:polyOrder()+1 -- number of quadrature points in each direction
   -- 1D weights and ordinates
   local ordinates, weights = GaussQuadRules.ordinates[N], GaussQuadRules.weights[N]

   local l, u = {}, {}
   for d = 0, self._basis:ndim() do l[i], u[i] = 1, N end
   -- construct weights and ordinates for integration in multiple dimensions
   local nr = Range.Range(l, u)
   for idx in nr:colMajorIter() do
   end

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(ProjectOnBasis, { __call = function (self, o) return self.new(self, o) end })

-- advance method
local function advance(self, tCurr, dt, inFld, outFld)
   return true, GKYL_MAX_DOUBLE
end

-- Methods in updater
ProjectOnBasis.__index = { advance = Base.advanceFuncWrap(advance) }

return {
   ProjectOnBasis = ProjectOnBasis
}
