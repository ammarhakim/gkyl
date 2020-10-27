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
local upper = {2*math.pi, 2*math.pi}
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
   evaluate = function(t, xn)
      local x, y = xn[1], xn[2]
      local amn = {{0,10,0}, {10,0,0}, {10,0,0}}
      local bmn = {{0,10,0}, {10,0,0}, {10,0,0}}
      local t1, t2 = 0.0, 0.0
      local f = 0.0
      for m = 0,2 do
	 for n = 0,2 do
	    t1 = amn[m+1][n+1]*math.cos(m*x)*math.cos(n*y)
	    t2 = bmn[m+1][n+1]*math.sin(m*x)*math.sin(n*y)
	    f = f + -(m*m+n*n)*(t1+t2)
	 end
      end
      return -f/50.0
   end,
}

local iterPoisson = Updater.IterPoisson {
   onGrid = grid,
   basis = basis,
}

initSource:advance(0.0, {}, {fIn})
iterPoisson:advance(0.0, {fIn}, {fOut})

fIn:write('fIn.bp', 0.0, 0)
fOut:write('fOut.bp', 0.0, 0)
