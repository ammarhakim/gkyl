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


polyOrder = 2 -- polynomial order
VDIM = 3 -- velocity dimensions
nMom = VDIM -- number of momentum component
nPrs = VDIM*(VDIM+1)/2 -- number of pressure tensor component

local phaseGrid = Grid.RectCart {
   lower = {-1.0, -1.0, -1.0, -6.0, -6.0, -6.0},
   upper = {1.0, 1.0, 1.0, 6.0, 6.0, 6.0},
   cells = {4, 4, 4, 8, 8, 8},
}
local confGrid = Grid.RectCart {
   lower = { phaseGrid:lower(1), phaseGrid:lower(2) , phaseGrid:lower(3) },
   upper = { phaseGrid:upper(1), phaseGrid:upper(2) , phaseGrid:upper(3) },
   cells = { phaseGrid:numCells(1), phaseGrid:numCells(2) , phaseGrid:numCells(3) },
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
function maxwellian(x,y,z,vx,vy,vz)
   local Pi = math.pi   
   local n = 1.0*math.sin(2*Pi*x)*math.sin(2*Pi*y)*math.sin(2*Pi*z)
   local ux = 0.1*math.cos(2*Pi*x)*math.cos(2*Pi*y)*math.cos(2*Pi*z)
   local uy = 0.2*math.sin(2*Pi*x)*math.sin(2*Pi*y)*math.sin(2*Pi*z)
   local uz = 0.1*math.cos(2*Pi*x)*math.cos(2*Pi*y)*math.cos(2*Pi*z)

   local Txx = 0.75 + 0.25*math.cos(2*Pi*x)*math.cos(2*Pi*y)*math.cos(2*Pi*z)
   local Tyy = 0.75 + 0.25*math.sin(2*Pi*x)*math.sin(2*Pi*y)*math.sin(2*Pi*z)
   local Tzz = 0.75 + 0.1*math.sin(2*Pi*x)*math.sin(2*Pi*y)*math.sin(2*Pi*z)
   local Txy = 0.5 + 0.1*math.sin(2*Pi*x)*math.sin(2*Pi*y)*math.sin(2*Pi*z)
   local Txz = 0.25 + 0.1*math.sin(2*Pi*x)*math.sin(2*Pi*y)*math.sin(2*Pi*z)
   local Tyz = 0.125 + 0.1*math.sin(2*Pi*x)*math.sin(2*Pi*y)*math.sin(2*Pi*z)

   local cx = vx-ux
   local cy = vy-uy
   local cz = vz-uz

   local detT = Txx*(Tyy*Tzz-Tyz^2)-Txy*(Txy*Tzz-Txz*Tyz)+Txz*(Txy*Tyz-Txz*Tyy)
   local u2 = cx*(cx*(Tyy*Tzz-Tyz^2)+cy*(Txz*Tyz-Txy*Tzz)+cz*(Txy*Tyz-Txz*Tyy))+cy*(cx*(Txz*Tyz-Txy*Tzz)+cy*(Txx*Tzz-Txz^2)+cz*(Txy*Txz-Txx*Tyz))+cz*(cx*(Txy*Tyz-Txz*Tyy)+cy*(Txy*Txz-Txx*Tyz)+cz*(Txx*Tyy-Txy^2))
   u2 = u2/(2*detT)

   return n/(math.pow(2*Pi, VDIM/2)*math.sqrt(detT))*math.exp(-u2)
end

-- updater to initialize distribution function
local project = Updater.ProjectOnBasis {
   onGrid = phaseGrid,
   basis = phaseBasis,
   evaluate = function (t, xn)
      return maxwellian(xn[1], xn[2], xn[3], xn[4], xn[5], xn[6])
   end
}

local tStart = Time.clock()
project:advance(0.0, {}, {distf})
distf:write("distf.bp", 0.0)
local tEnd = Time.clock()

local tStart = Time.clock()
-- compute moments
numDensityCalc:advance(0.0, {distf}, {numDensity})
momentumCalc:advance(0.0, {distf}, {momentum})
pressureTensorCalc:advance(0.0, {distf}, {pressureTensor})
ptclEnergyCalc:advance(0.0, {distf}, {ptclEnergy})
local tEnd = Time.clock()

-- write data
numDensity:write("numDensity.bp", 0.0)
momentum:write("momentum.bp", 0.0)
pressureTensor:write("pressureTensor.bp", 0.0)
ptclEnergy:write("ptclEnergy.bp", 0.0)
