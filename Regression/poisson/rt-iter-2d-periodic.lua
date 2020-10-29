-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Time = require "Lib.Time"

local polyOrder = 1
local lower = {0, 0}
local upper = {2*math.pi, 2*math.pi}
local cells = {32, 32}
local periodicDirs = {1, 2}

local grid = Grid.RectCart {
   lower = lower,
   upper = upper,
   cells = cells,
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

-- Initial conditions from
-- http://ammar-hakim.org/sj/je/je11/je11-fem-poisson.html#convergence-of-2d-solver
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

   -- there parameters will eventually be replaced by internal
   -- heuristics
   
   errEps = 1e-8, -- maximum residual error
   factor = 120, -- factor over explicit scheme
   extraStages = 2, -- extra stages
   cflFrac = 1.0, -- CFL frac for internal iterations
   stepper = 'RKL1', -- stepper to use 'RKL1' or 'RKL2'
   extrapolateInterval = 2, -- extrapolate every these many steps
}

initSource:advance(0.0, {}, {fIn})

local tmStart = Time.clock()
iterPoisson:advance(0.0, {fIn}, {fOut})
print(string.format("Simulation took %g", Time.clock()-tmStart))

fIn:write('fIn.bp', 0.0, 0)
fOut:write('fOut.bp', 0.0, 0)

iterPoisson:writeDiagnostics()
