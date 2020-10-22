-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"

local polyOrder = 1
local lower = {-1, -1}
local upper = {1, 1}
local cells = {32, 32}
local periodicDirs = {1,2}

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
-- http://ammar-hakim.org/sj/je/je1/je1-periodic-poisson.html
local initSource = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x, y = z[1], z[2]
      return -math.exp(-10*(2*x*x+4*x*y+5*y*y))
   end,
}
-- local initExact = Updater.ProjectOnBasis {
--    onGrid = grid,
--    basis = basis,
--    evaluate = function(t, z)
--       local x, y = z[1], z[2]
--       local a, b = 2, 5
--       local c1, d0 = 0, 0
--       local c0 = a/12 - 1/2
--       local d1 = b/12 - 1/2
--       local t1 = x^2/2 - a*x^4/12 + c0*x + c1
--       local t2 = y^2/2 - b*y^4/12 + d0*y + d1
--       return t1*t2
--    end,
-- }

local discontPoisson = Updater.DiscontPoisson {
   onGrid = grid,
   basis = basis,
   bcLower = { { }, { } },
   bcUpper = { { }, { } },
}

initSource:advance(0.0, {}, {fIn})
--initExact:advance(0.0, {}, {fExact})
discontPoisson:advance(0.0, {fIn}, {fOut})

fIn:write('fIn.bp', 0.0, 0)
--fExact:write('fExact.bp', 0.0, 0)
fOut:write('fOut.bp', 0.0, 0)
