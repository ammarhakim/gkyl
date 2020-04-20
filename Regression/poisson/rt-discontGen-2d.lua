-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"

local x0, y0 = 0, 0
local bx = function(x, y) return -(y-y0)/math.sqrt((x-x0)^2+(y-y0)^2) end
local by = function(x, y) return (x-x0)/math.sqrt((x-x0)^2+(y-y0)^2) end

local zet = 1e9 -- Dpar/Dperp
local Dpar = 1.0 -- parallel heat conduction
local Dperp = Dpar/zet -- perpendicular heat conduction
-- diffusion coefficients
local DxxFn = function(t, z)
   local x, y = z[1], z[2]
   return Dpar*bx(x,y)^2 + Dperp*by(x,y)^2
end
local DyyFn = function(t, z)
   local x, y = z[1], z[2]
   return Dperp*bx(x,y)^2 + Dpar*by(x,y)^2
end
local DxyFn = function(t, z)
   local x, y = z[1], z[2]
   return (Dpar-Dperp)*bx(x,y)*by(x,y)
end

local numQuad = 7
local solFn = function(t, z)
   local x, y = z[1], z[2]
   return x*(x+0.5)*(x-0.5)*y*(y+0.5)*(y-0.5)
end
local srcFn = function(t, z)
   local x, y = z[1], z[2]
   return (-1.0*((x-0.5)*x*(x+0.5)*y*(y+0.5)+(x-0.5)*x*(x+0.5)*(y-0.5)*(y+0.5)+(x-0.5)*x*(x+0.5)*(y-0.5)*y)*((-(2.0*(y-1.0*y0)*(y0-1.0*y)^2)/(((y-1.0*y0)^2+(x-1.0*x0)^2)^2*zet))-(2.0*(y0-1.0*y))/(((y-1.0*y0)^2+(x-1.0*x0)^2)*zet)-(2.0*(x-1.0*x0)^2*(y-1.0*y0))/((y-1.0*y0)^2+(x-1.0*x0)^2)^2))-1.0*(2.0*(x-0.5)*x*(x+0.5)*(y+0.5)+2.0*(x-0.5)*x*(x+0.5)*y+2.0*(x-0.5)*x*(x+0.5)*(y-0.5))*((y0-1.0*y)^2/(((y-1.0*y0)^2+(x-1.0*x0)^2)*zet)+(x-1.0*x0)^2/((y-1.0*y0)^2+(x-1.0*x0)^2))-1.0*(x*(x+0.5)*(y-0.5)*y*(y+0.5)+(x-0.5)*(x+0.5)*(y-0.5)*y*(y+0.5)+(x-0.5)*x*(y-0.5)*y*(y+0.5))*((2.0*(x-1.0*x0))/(((y-1.0*y0)^2+(x-1.0*x0)^2)*zet)-(2.0*(x-1.0*x0)^3)/(((y-1.0*y0)^2+(x-1.0*x0)^2)^2*zet)-(2.0*(x-1.0*x0)*(y0-1.0*y)^2)/((y-1.0*y0)^2+(x-1.0*x0)^2)^2)-1.0*(2.0*(x+0.5)*(y-0.5)*y*(y+0.5)+2.0*x*(y-0.5)*y*(y+0.5)+2.0*(x-0.5)*(y-0.5)*y*(y+0.5))*((x-1.0*x0)^2/(((y-1.0*y0)^2+(x-1.0*x0)^2)*zet)+(y0-1.0*y)^2/((y-1.0*y0)^2+(x-1.0*x0)^2))+(2.0*(x-1.0*x0)*(x*(x+0.5)*(y-0.5)*y*(y+0.5)+(x-0.5)*(x+0.5)*(y-0.5)*y*(y+0.5)+(x-0.5)*x*(y-0.5)*y*(y+0.5))*(y-1.0*y0)*(y0-1.0*y)*(1.0-1.0/zet))/((y-1.0*y0)^2+(x-1.0*x0)^2)^2-(1.0*((x-0.5)*x*(x+0.5)*y*(y+0.5)+(x-0.5)*x*(x+0.5)*(y-0.5)*(y+0.5)+(x-0.5)*x*(x+0.5)*(y-0.5)*y)*(y0-1.0*y)*(1.0-1.0/zet))/((y-1.0*y0)^2+(x-1.0*x0)^2)-(2.0*(x-1.0*x0)*(x*(x+0.5)*y*(y+0.5)+(x-0.5)*(x+0.5)*y*(y+0.5)+(x-0.5)*x*y*(y+0.5)+x*(x+0.5)*(y-0.5)*(y+0.5)+(x-0.5)*(x+0.5)*(y-0.5)*(y+0.5)+(x-0.5)*x*(y-0.5)*(y+0.5)+x*(x+0.5)*(y-0.5)*y+(x-0.5)*(x+0.5)*(y-0.5)*y+(x-0.5)*x*(y-0.5)*y)*(y0-1.0*y)*(1.0-1.0/zet))/((y-1.0*y0)^2+(x-1.0*x0)^2)+(2.0*(x-1.0*x0)^2*((x-0.5)*x*(x+0.5)*y*(y+0.5)+(x-0.5)*x*(x+0.5)*(y-0.5)*(y+0.5)+(x-0.5)*x*(x+0.5)*(y-0.5)*y)*(y0-1.0*y)*(1.0-1.0/zet))/((y-1.0*y0)^2+(x-1.0*x0)^2)^2+((x-1.0*x0)*(x*(x+0.5)*(y-0.5)*y*(y+0.5)+(x-0.5)*(x+0.5)*(y-0.5)*y*(y+0.5)+(x-0.5)*x*(y-0.5)*y*(y+0.5))*(1.0-1.0/zet))/((y-1.0*y0)^2+(x-1.0*x0)^2)
end

local grid = Grid.RectCart {
   lower = {-0.5, -0.5},
   upper = {0.5, 0.5},
   cells = {16, 16},
   periodicDirs = {}
}
local basis = Basis.CartModalSerendipity {
   ndim = grid:ndim(),
   polyOrder = 1,
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

local src = getField()
local initSource = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = numQuad,
   evaluate = srcFn,
}
initSource:advance(0.0, {}, {src})

local solExact = getField()
local initSol = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = numQuad,
   evaluate = solFn,
}
initSol:advance(0.0, {}, {solExact})

local Dxx = getField()
local initDxx = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = numQuad,
   evaluate = DxxFn,
   projectOnGhosts = true,
}
initDxx:advance(0.0, {}, {Dxx})

local Dyy = getField()
local initDyy = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = numQuad,
   evaluate = DyyFn,
   projectOnGhosts = true,
}
initDyy:advance(0.0, {}, {Dyy})

local Dxy = getField()
local initDxy = Updater.ProjectOnBasis {
   onGrid = grid,
   basis = basis,
   numQuad = numQuad,
   evaluate = DxyFn,
   projectOnGhosts = true,
}
initDxy:advance(0.0, {}, {Dxy})

local solSim = getField()
local discontPoisson = Updater.DiscontGenPoisson {
   onGrid = grid,
   basis = basis,
   Dxx = Dxx,
   Dyy = Dyy,
   Dxy = Dxy,
   bcLower = { {D=1, N=0, val=0.0}, {D=1, N=0, val=0.0} },
   bcUpper = { {D=1, N=0, val=0.0}, {D=1, N=0, val=0.0} },
   writeMatrix = false,
}
discontPoisson:advance(0.0, {src}, {solSim})

src:write(string.format('src.bp'), 0.0, 0)
solExact:write(string.format('solExact.bp'), 0.0, 0)
solSim:write(string.format('solSim.bp'), 0.0, 0)
Dxx:write(string.format('Dxx.bp'), 0.0, 0)
Dyy:write(string.format('Dyy.bp'), 0.0, 0)
Dxy:write(string.format('Dxy.bp'), 0.0, 0)
