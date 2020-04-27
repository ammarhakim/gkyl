-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"

local polyOrder = 1
local lower = {0, 0}
local upper = {1, 1}
local cells = {16, 16}
local periodicDirs = {}

local grid = Grid.RectCart {
   lower = {lower[1], lower[2]},
   upper = {upper[1], upper[2]},
   cells = {cells[1], cells[2]},
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
local initDist = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x, y = z[1], z[2]
      local a, b = 2, 5
      local c1, d0 = 0, 0
      local c0 = a/12 - 1/2
      local d1 = b/12 - 1/2
      local t1 = (1-a*x^2)*(-b*y^4/12 + y^2/2 + d0*y + d1)
      local t2 = (1-b*y^2)*(-a*x^4/12 + x^2/2 + c0*x + c1)
      return -t1-t2
   end,
}
local initExact = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x, y = z[1], z[2]
      local a, b = 2, 5
      local c1, d0 = 0, 0
      local c0 = a/12 - 1/2
      local d1 = b/12 - 1/2
      local t1 = x^2/2 - a*x^4/12 + c0*x + c1
      local t2 = y^2/2 - b*y^4/12 + d0*y + d1
      return t1*t2
   end,
}

local discontPoisson = Updater.DiscontPoisson {
   onGrid = grid,
   basis = basis,
   bcLower = { { T = "D", V = 0.0 }, { T = "N", V = 0.0 } },
   bcUpper = { { T = "D", V = 0.0 }, { T = "D", V = 0.0 } },
}

initDist:advance(0.0, {}, {fIn})
--initExact:advance(0.0, {}, {fExact})
discontPoisson:advance(0.0, {fIn}, {fOut})

fIn:write('fIn.bp', 0.0, 0)
--fExact:write('fExact.bp', 0.0, 0)
fOut:write('fOut.bp', 0.0, 0)
