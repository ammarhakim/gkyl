-- Gkyl ------------------------------------------------------------------------
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Time = require "Lib.Time"
local Basis = require "Basis"
local Updater = require "Updater"


polyOrder = 1 -- polynomial order
VDIM = 1 -- velocity dimensions
nMom = VDIM -- number of momentum component
nPrs = VDIM*(VDIM+1)/2 -- number of pressure tensor component

local phaseGrid = Grid.RectCart {
   lower = {-1.0, -6.0},
   upper = {1.0, 6.0},
   cells = {32, 8},
}
local confGrid = Grid.RectCart {
   lower = { phaseGrid:lower(1) },
   upper = { phaseGrid:upper(1) },
   cells = { phaseGrid:numCells(1) },
}

-- basis functions
local phaseBasis = Basis.CartModalSerendipity { ndim = phaseGrid:ndim(), polyOrder = polyOrder }
local confBasis = Basis.CartModalSerendipity { ndim = confGrid:ndim(), polyOrder = polyOrder }

-- fields
local distf = DataStruct.Field {
   onGrid = phaseGrid,
   numComponents = phaseBasis:numBasis(),
   ghost = {0, 0},
}

local numDensity = DataStruct.Field {
   onGrid = confGrid,
   numComponents = confBasis:numBasis(),
   ghost = {0, 0},
}
local momentum = DataStruct.Field {
   onGrid = confGrid,
   numComponents = nMom*confBasis:numBasis(),
   ghost = {0, 0},
}
local pressureTensor = DataStruct.Field {
   onGrid = confGrid,
   numComponents = nPrs*confBasis:numBasis(),
   ghost = {0, 0},
}
local ptclEnergy = DataStruct.Field {
   onGrid = confGrid,
   numComponents = confBasis:numBasis(),
   ghost = {0, 0},
}

--------------
-- Updaters --
--------------

numDensityCalc = Updater.DistFuncMomentCalc {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   confBasis = confBasis,
   moment = "M0",
}
momentumCalc = Updater.DistFuncMomentCalc {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   confBasis = confBasis,
   moment = "M1i",
}		
pressureTensorCalc = Updater.DistFuncMomentCalc {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   confBasis = confBasis,
   moment = "M2ij",
}
ptclEnergyCalc = Updater.DistFuncMomentCalc {
   onGrid = phaseGrid,
   phaseBasis = phaseBasis,
   confBasis = confBasis,
   moment = "M2",
}

-- initial condition to apply
function maxwellian(x, vx)
   local Pi = math.pi
   local n = 1.0*math.sin(2*Pi*x)
   local ux = 0.1*math.cos(2*Pi*x)
   local Txx = 0.75 + 0.25*math.cos(2*Pi*x)
   
   local u2 = (vx-ux)^2/(2*Txx)
   return n/math.sqrt(2*Pi*Txx)*math.exp(-u2)
end

-- updater to initialize distribution function
local project = Updater.ProjectOnBasis {
   onGrid = phaseGrid,
   basis = phaseBasis,
   evaluate = function (t, xn)
      return maxwellian(xn[1], xn[2])
   end
}
project:advance(0.0, {}, {distf})
distf:write("distf.bp", 0.0)

local tStart = Time.clock()
-- compute moments
numDensityCalc:advance(0.0, {distf}, {numDensity})
momentumCalc:advance(0.0, {distf}, {momentum})
pressureTensorCalc:advance(0.0, {distf}, {pressureTensor})
ptclEnergyCalc:advance(0.0, {distf}, {ptclEnergy})
local tEnd = Time.clock()
print("Moment calculations took", tEnd-tStart)

-- write data
numDensity:write("numDensity.bp", 0.0)
momentum:write("momentum.bp", 0.0)
pressureTensor:write("pressureTensor.bp", 0.0)
ptclEnergy:write("ptclEnergy.bp", 0.0)
