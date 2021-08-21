-- Gkyl ------------------------------------------------------------------------
--
-- Test the twist-shift updater for twist-shift BCs in gyrokinetics.
--
-- This tests the 3x2v twist shift.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"

local vt   = 1.0   -- Thermal speed.
local mass = 1.0
local B0   = 1.0   -- Magnetic field magnitude. 

local polyOrder       = 1
local lower           = {-2.0, -1.50, -3.0, -5.0*vt, 0.}
local upper           = { 2.0,  1.50,  3.0,  5.0*vt, mass*((5.0*vt)^2)/(2.0*B0)}
local numCells        = {10, 10, 4, 4, 4}
local periodicDirs    = {2}
local yShiftPolyOrder = 1

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
      metaData      = {polyOrder = basis:polyOrder(), basisType = basis:id()},
   }
   fld:clear(0.)
   return fld
end

local wrapNum = function (val, lims, pickUpper)
   -- Wrap a number to range [lower,upper]. If pickUpper=true, output upper when
   -- val is a multiple of upper. Otherwise multiples of upper wrap to lower.
   local lower, upper = lims.lo, lims.up
   local L        = upper - lower
   local disp     = (val - lower) % L
   local newCoord = lower + (L + disp) % L
   local eps      = 1.e-12
   if ( (lower-eps < newCoord and newCoord < lower + eps) or
        (upper-eps < newCoord and newCoord < upper + eps) ) then
      if pickUpper then 
         return upper
      else
         return lower
      end
   else
      return newCoord
   end
end

local grid = Grid.RectCart {
   lower = lower,  cells        = numCells,
   upper = upper,  periodicDirs = periodicDirs,
}
local confGrid = Grid.RectCart {
   lower = {lower[1],lower[2],lower[3]},  cells        = {numCells[1],numCells[2],numCells[3]},
   upper = {upper[1],upper[2],upper[3]},  periodicDirs = periodicDirs,
}
local basis         = Basis.CartModalSerendipity { ndim = grid:ndim(), polyOrder = polyOrder }
local confBasis     = Basis.CartModalSerendipity { ndim = confGrid:ndim(), polyOrder = polyOrder }
local fldDo         = createField(grid, basis)
local fldDoShifted  = createField(grid, basis)
local fldDoNoGhosts = createField(grid, basis)
local fldTar        = createField(grid, basis)

-- Create a 1D grid and project the function that determines the shift.
-- In flux-tube gyrokinetics this shift is a function of the magnetic
-- safety profile, something like yShift = L_z*C_y(x)*q(x).
local polyOrder1D = polyOrder
local grid1D      = Grid.RectCart {
   lower = {lower[1]},  cells        = {numCells[1]},
   upper = {upper[1]},  periodicDirs = {},
}
local basis1D = Basis.CartModalSerendipity { ndim = grid1D:ndim(), polyOrder = polyOrder1D }

-- Function describing the shift in y of the BC.
-- It has to be everywhere >0 or everywhere <0 (it cannot be zero, or too close to it).
local yShiftFunc = function(t, xn)
                      local x = xn[1]
--                      return 1./(1.+0.25*x)
                      return -0.*x+0.97
                   end

-- Donor field function.
local fldDoFunc = function(t, xn)
   local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
   local vExp  = -(vpar^2 + 2.*math.abs(mu)*B0/mass)/(2.0*(vt^2))
   local vFunc = (1./math.sqrt((2.*math.pi*(vt^2))^3)) * math.exp(vExp)

   local muX, muY   = 0., 0.
   local sigX, sigY = 0.3, 0.3
   return math.exp(-((x-muX)^2)/(2.*(sigX^2))-((y-muY)^2)/(2.*(sigY^2)))*vFunc
end
local fldDoFuncUpperOnly = function(t, xn)
   local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
   if ((z > grid:upper(3)) and 
      (vpar > grid:lower(4)) and (vpar < grid:upper(4)) and 
      (mu > grid:lower(5)) and (mu < grid:upper(5))) then
      return fldDoFunc(t, xn)
   else
      return 0.
   end
end

-- Shifted donor field function.
local fldDoShiftedFunc = function(t, xn)
   local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
   local vExp  = -(vpar^2 + 2.*math.abs(mu)*B0/mass)/(2.0*(vt^2))
   local vFunc = (1./math.sqrt((2.*math.pi*(vt^2))^3)) * math.exp(vExp)

   local muX, muY   = 0., 0.
   local sigX, sigY = 0.3, 0.3
   return math.exp(-((x-muX)^2)/(2.*(sigX^2))
                   -((wrapNum(y-yShiftFunc(0,xn),{lo=grid:lower(2),up=grid:upper(2)},true)-muY)^2)/(2.*(sigY^2)))*vFunc
end

-- Projection updaters.
local project = Updater.ProjectOnBasis {
   onGrid = grid,   evaluate = function(t, xn) return 1. end,
   basis  = basis,
}
local projectOnGhosts = Updater.ProjectOnBasis {
   onGrid = grid,   evaluate = function(t, xn) return 1. end,
   basis  = basis,  onGhosts = true,
}
local projectBmag = Updater.EvalOnNodes {
   onGrid = confGrid,   evaluate = function(t, xn) return B0 end,
   basis  = confBasis,  onGhosts = true,
}

-- Project the magnetic field amplitude.
local bmag = createField(confGrid, confBasis)
projectBmag:advance(0., {}, {bmag})

-- Project donor field function onto basis.
projectOnGhosts:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
projectOnGhosts:advance(0., {}, {fldDo})
fldDo:write("fldDo.bp")
-- Project shifted donor field function onto basis.
projectOnGhosts:setFunc(function(t,xn) return fldDoShiftedFunc(t,xn) end)
projectOnGhosts:advance(0., {}, {fldDoShifted})
fldDoShifted:write("fldDoShifted.bp")

-- Project the donor field function into the target field but not
-- in the ghost cells. Twist-shift will fill the ghost cells.
project:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
project:advance(0., {}, {fldDoNoGhosts})
fldTar:copy(fldDoNoGhosts)

-- Compute the velocity moments of the donor field in the upper ghost cells only,
-- to be compared with the moments of the target field in the lower ghost cells.
projectOnGhosts:setFunc(function(t,xn) return fldDoFuncUpperOnly(t,xn) end)
projectOnGhosts:advance(0., {}, {fldDo})
fldDo:write("fldDo_upperGhost.bp", 0, 0, true)
local m0Do = createField(confGrid, confBasis)
local m1Do = createField(confGrid, confBasis)
local m2Do = createField(confGrid, confBasis)
local numDensityCalcOnGhosts = Updater.DistFuncMomentCalc {
   onGrid     = grid,       moment   = "GkM0", -- GkM0 = < f >
   phaseBasis = basis,      gkfacs   = {mass, bmag},
   confBasis  = confBasis,  onGhosts = true,
}
local momDensityCalcOnGhosts = Updater.DistFuncMomentCalc {
   onGrid     = grid,       moment   = "GkM1", -- GkM1 = < v_parallel f >
   phaseBasis = basis,      gkfacs   = {mass, bmag},
   confBasis  = confBasis,  onGhosts = true,
}
local ptclEnergyCalcOnGhosts = Updater.DistFuncMomentCalc {
   onGrid     = grid,       moment   = "GkM2", -- GkM2 = < (v_parallel^2 + 2*mu*B/m) f >
   phaseBasis = basis,      gkfacs   = {mass, bmag},
   confBasis  = confBasis,  onGhosts = true,
}
numDensityCalcOnGhosts:advance(0., {fldDo}, {m0Do})
momDensityCalcOnGhosts:advance(0., {fldDo}, {m1Do})
ptclEnergyCalcOnGhosts:advance(0., {fldDo}, {m2Do})
m0Do:write("fldDo_M0.bp",0., 0, true)
m1Do:write("fldDo_M1.bp",0., 0, true)
m2Do:write("fldDo_M2.bp",0., 0, true)

-- Compute integrated moments.
local intQuantOnGhosts = Updater.CartFieldIntegratedQuantCalc {
   onGrid   = confGrid,   numComponents = 1,
   basis    = confBasis,  quantity      = "V",
   onGhosts = true,
}
local intM0Do = DataStruct.DynVector { numComponents = 1, }
local intM1Do = DataStruct.DynVector { numComponents = 1, }
local intM2Do = DataStruct.DynVector { numComponents = 1, }
intQuantOnGhosts:advance(0., {m0Do}, {intM0Do})
intQuantOnGhosts:advance(0., {m1Do}, {intM1Do})
intQuantOnGhosts:advance(0., {m2Do}, {intM2Do})
intM0Do:write("fldDo_intM0.bp",0., 0)
intM1Do:write("fldDo_intM1.bp",0., 0)
intM2Do:write("fldDo_intM2.bp",0., 0)

-- Compute the integral of the donor field.
local intQuantPhase = Updater.CartFieldIntegratedQuantCalc {
   onGrid   = grid,   numComponents = 1,
   basis    = basis,  quantity      = "V",
   onGhosts = true,
}
local intFldDo = DataStruct.DynVector { numComponents = 1, }
intQuantPhase:advance(0., {fldDo}, {intFldDo})
intFldDo:write("intFldDo.bp", 0., 0)

-- Apply the shift to the target field.
local twistShiftUpd = Updater.TwistShift {
   onGrid    = grid,       yShiftFunc      = yShiftFunc, 
   basis     = basis,      yShiftPolyOrder = yShiftPolyOrder, 
   confBasis = confBasis,  edge            = "lower",
}

local t1 = os.clock()
twistShiftUpd:_advance(0., {}, {fldTar})
local t2 = os.clock()
io.write("Total test time: ", t2-t1, " s\n")

fldTar:write("fldTar.bp", 0., 0, true)

-- Compute the velocity moments of the target field.
-- Subtract the donor field in the inner cells so we only compute
-- ghost cell moments.
fldTar:accumulate(-1., fldDoNoGhosts)
fldTar:write("fldTar_lowerGhost.bp", 0, 0, true)
local m0Tar = createField(confGrid, confBasis)
local m1Tar = createField(confGrid, confBasis)
local m2Tar = createField(confGrid, confBasis)
numDensityCalcOnGhosts:advance(0., {fldTar}, {m0Tar})
momDensityCalcOnGhosts:advance(0., {fldTar}, {m1Tar})
ptclEnergyCalcOnGhosts:advance(0., {fldTar}, {m2Tar})
m0Tar:write("fldTar_M0.bp",0., 0, true)
m1Tar:write("fldTar_M1.bp",0., 0, true)
m2Tar:write("fldTar_M2.bp",0., 0, true)

-- Compute integrated moments.
local intM0Tar = DataStruct.DynVector { numComponents = 1, }
local intM1Tar = DataStruct.DynVector { numComponents = 1, }
local intM2Tar = DataStruct.DynVector { numComponents = 1, }
intQuantOnGhosts:advance(0., {m0Tar}, {intM0Tar})
intQuantOnGhosts:advance(0., {m1Tar}, {intM1Tar})
intQuantOnGhosts:advance(0., {m2Tar}, {intM2Tar})
intM0Tar:write("fldTar_intM0.bp",0., 0, true)
intM1Tar:write("fldTar_intM1.bp",0., 0, true)
intM2Tar:write("fldTar_intM2.bp",0., 0, true)

-- Compute the integral of the target field.
local intFldTar = DataStruct.DynVector { numComponents = 1, }
intQuantPhase:advance(0., {fldTar}, {intFldTar})
intFldTar:write("intFldTar.bp", 0., 0)
