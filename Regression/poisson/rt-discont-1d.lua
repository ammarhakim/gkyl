-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"

local polyOrder = 1
local lower = {0}
local upper = {1}
local cells = {24}
local periodicDirs = {}

local grid = Grid.RectCart {
   lower = {lower[1]},
   upper = {upper[1]},
   cells = {cells[1]},
   periodicDirs = periodicDirs,
}
local basis = Basis.CartModalSerendipity {
   ndim = grid:ndim(),
   polyOrder = polyOrder
}

local function getField()
   return DataStruct.Field {
      onGrid = grid,
      numComponents = basis:numBasis(),
      ghost = {1, 1},
      metaData = {
	 polyOrder = basis:polyOrder(),
	 basisType = basis:id(),
      },
   }
end
local fIn = getField()
local fOut = getField()
local fExact = getField()

-- Initial conditions from:
-- http://ammar-hakim.org/sj/je/je11/je11-fem-poisson.html
local aD = 2
local initDistD = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x = z[1]
      return 1-aD*x^2
   end,
}
local cD = aD/12.0 - 0.5
local initExactD = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x = z[1]
      return x^2/2-aD*x^4/12+cD*x + 1
   end,
}
local discontPoissonD = Updater.DiscontPoisson {
   onGrid = grid,
   basis = basis,
   bcLower = { {T="D", V=1.0} },
   bcUpper = { {T="D", V=1.0} },
}

local aN = 5
local initDistN = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x = z[1]
      return 1-aN*x^2
   end,
}
local cN = aN/12.0 - 0.5
local initExactN = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x = z[1]
      return x^2/2-aN*x^4/12+cN
   end,
}
local discontPoissonN = Updater.DiscontPoisson {
   onGrid = grid,
   basis = basis,
   bcLower = { {T="N", V=0.0} },
   bcUpper = { {T="D", V=0.0} },
}

initDistD:advance(0.0, {}, {fIn})
--initExactD:advance(0.0, {}, {fExact})
discontPoissonD:advance(0.0, {fIn}, {fOut})
fIn:write('fDIn.bp', 0.0, 0)
--fExact:write('fDExact.bp', 0.0, 0)
fOut:write('fDOut.bp', 0.0, 0)

initDistN:advance(0.0, {}, {fIn})
--initExactN:advance(0.0, {}, {fExact})
discontPoissonN:advance(0.0, {fIn}, {fOut})
fIn:write('fNIn.bp', 0.0, 0)
--fExact:write('fNExact.bp', 0.0, 0)
fOut:write('fNOut.bp', 0.0, 0)
