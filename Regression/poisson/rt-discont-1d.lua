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

local function maxwellian (vx,vy,ux,uy,vtx,vty)
   return  1/(2*math.pi)*math.exp(-(vx-ux)^2/(2*vtx^2))*math.exp(-(vy-uy)^2/(2*vty^2))
end

-- Initial conditions from:
-- http://ammar-hakim.org/sj/je/je11/je11-fem-poisson.html
local a = 2
local initDist = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x = z[1]
      return 1-a*x^2
   end,
}
local c = a/12.0 - 0.5
local initExact = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   evaluate = function(t, z)
      local x = z[1]
      return x^2/2-a*x^4/12+c*x
   end,
}

local discontPoisson = Updater.DiscontPoisson {
   onGrid = grid,
   basis = basis,
   bcLeft = { T = "D", V = 0.0 },
   bcRight = { T = "D", V = 0.0 },
}

initDist:advance(0.0, {}, {fIn})
initExact:advance(0.0, {}, {fExact})
discontPoisson:advance(0.0, {fIn}, {fOut})

fIn:write('fIn.bp', 0.0, 0)
fExact:write('fExact.bp', 0.0, 0)
fOut:write('fOut.bp', 0.0, 0)
