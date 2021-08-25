-- Gkyl ------------------------------------------------------------------------
--
-- Test the twist-shift updater for twist-shift BCs in gyrokinetics.
--
-- This tests the 3x2v twist shift by taking a field that is only non-zero in
-- the inner cells, and using the twistShift update to populate its z ghost
-- cells.
--
-- In order to check that the moments of that ghost cell are the same as those
-- of the skin cell it originated from, we need to create a field that projects
-- the donor function only in the ghost cell or only in the skin cell.
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
local dy              = (upper[2]-lower[2])/numCells[2]
local dz              = (upper[3]-lower[3])/numCells[3]

-- Function describing the shift in y of the BC.
-- It has to be everywhere >0 or everywhere <0 (it cannot be zero, or too close to it).
local yShiftFunc = function(t, xn)
                      local x = xn[1]
--                      return 1./(1.+0.25*x)
                      return -0.3*x+0.97
--                      return dy/2.
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
-- ....................... END OF USER INPUTS (maybe) ........................... --


-- Donor field function that only allows projection in the z skin cell.
local fldDoFuncZskinOnly = {
   lower = function(t, xn)
      local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
      if ((   x > lower[1]) and (   x < upper[1])) and 
         ((   y > lower[2]) and (   y < upper[2])) and 
         ((   z > lower[3]) and (   z < lower[3]+dz)) and 
         ((vpar > lower[4]) and (vpar < upper[4])) and 
         (  (mu > lower[5]) and (  mu < upper[5])) then
         return fldDoFunc(t, xn)
      else
         return 0.
      end
   end,
   upper = function(t, xn)
      local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
      if ((   x > lower[1])    and (   x < upper[1])) and 
         ((   y > lower[2])    and (   y < upper[2])) and 
         ((   z > upper[3]-dz) and (   z < upper[3])) and 
         ((vpar > lower[4])    and (vpar < upper[4])) and 
         (  (mu > lower[5])    and (  mu < upper[5])) then
         return fldDoFunc(t, xn)
      else
         return 0.
      end
   end
}

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

-- Shifted donor field function.
local fldDoShiftedFunc = function(t, xn)
   local x, y, z, vpar, mu = xn[1], xn[2], xn[3], xn[4], xn[5]
   local yS = wrapNum(y-yShiftFunc(0,xn),{lo=lower[2],up=upper[2]},true)
   return fldDoFunc(t, {x, yS, z, vpar, mu})
end

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

-- Project donor field function onto basis (including ghost cells).
projectOnGhosts:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
projectOnGhosts:advance(0., {}, {fldDo})
fldDo:write("fldDo.bp")
-- Project shifted donor field function onto basis.
projectOnGhosts:setFunc(function(t,xn) return fldDoShiftedFunc(t,xn) end)
projectOnGhosts:advance(0., {}, {fldDoShifted})
fldDoShifted:write("fldDoShifted.bp")

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
local intQuantOnGhosts = Updater.CartFieldIntegratedQuantCalc {
   onGrid   = confGrid,   numComponents = 1,
   basis    = confBasis,  quantity      = "V",
   onGhosts = true,
}

local m0Do  = createField(confGrid, confBasis)
local m1Do  = createField(confGrid, confBasis)
local m2Do  = createField(confGrid, confBasis)
local m0Tar = createField(confGrid, confBasis)
local m1Tar = createField(confGrid, confBasis)
local m2Tar = createField(confGrid, confBasis)

local intM0Do, intM1Do, intM2Do    = {}, {}, {}
local intM0Tar, intM1Tar, intM2Tar = {}, {}, {}
local intFldDo  = {}
local intFldTar = {}

local intQuantPhase = Updater.CartFieldIntegratedQuantCalc {
   onGrid   = grid,   numComponents = 1,
   basis    = basis,  quantity      = "V",
   onGhosts = true,
}

local twistShiftUpd = {}

local edges     = {"lower","upper"}
local dualEdges = {"upper","lower"}

local shiftFuncs = {lower=yShiftFunc, upper=function(t,xn) return -yShiftFunc(t,xn) end}

for i, edge in ipairs(edges) do

   -- Project the donor field function into the target field but not
   -- in the ghost cells. Twist-shift will fill the ghost cells.
   project:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
   project:advance(0., {}, {fldDoNoGhosts})
   fldTar:copy(fldDoNoGhosts)

   local dualEdge = dualEdges[i] 
   -- Compute the velocity moments of the donor field in the upper ghost cells only,
   -- to be compared with the moments of the target field in the lower ghost cells.
   projectOnGhosts:setFunc(function(t,xn) return fldDoFuncZskinOnly[dualEdge](t,xn) end)
   projectOnGhosts:advance(0., {}, {fldDo})
   fldDo:write("fldDo_"..dualEdge..".bp", 0, 0)
   numDensityCalcOnGhosts:advance(0., {fldDo}, {m0Do})
   momDensityCalcOnGhosts:advance(0., {fldDo}, {m1Do})
   ptclEnergyCalcOnGhosts:advance(0., {fldDo}, {m2Do})
   m0Do:write("fldDo_M0_"..dualEdge..".bp",0., 0)
   m1Do:write("fldDo_M1_"..dualEdge..".bp",0., 0)
   m2Do:write("fldDo_M2_"..dualEdge..".bp",0., 0)
   
   -- Compute integrated moments.
   intM0Do[edge] = DataStruct.DynVector { numComponents = 1, }
   intM1Do[edge] = DataStruct.DynVector { numComponents = 1, }
   intM2Do[edge] = DataStruct.DynVector { numComponents = 1, }
   intQuantOnGhosts:advance(0., {m0Do}, {intM0Do[edge]})
   intQuantOnGhosts:advance(0., {m1Do}, {intM1Do[edge]})
   intQuantOnGhosts:advance(0., {m2Do}, {intM2Do[edge]})
   intM0Do[edge]:write("fldDo_intM0_"..dualEdge..".bp",0., 0)
   intM1Do[edge]:write("fldDo_intM1_"..dualEdge..".bp",0., 0)
   intM2Do[edge]:write("fldDo_intM2_"..dualEdge..".bp",0., 0)
   
   -- Compute the integral of the donor field.
   intFldDo[edge] = DataStruct.DynVector { numComponents = 1, }
   intQuantPhase:advance(0., {fldDo}, {intFldDo[edge]})
   intFldDo[edge]:write("fldDo_intV_"..dualEdge..".bp", 0., 0)
   
   -- Apply the shift to the target field.
   twistShiftUpd[edge] = Updater.TwistShiftBC {
      onGrid    = grid,       yShiftFunc      = shiftFuncs[edge], 
      basis     = basis,      yShiftPolyOrder = yShiftPolyOrder, 
      confBasis = confBasis,  edge            = edge,
   }
   
   local t1 = os.clock()
   twistShiftUpd[edge]:advance(0., {}, {fldTar})
   local t2 = os.clock()
   io.write(edge.." shift time: ", t2-t1, " s\n")
   
   -- Compute the velocity moments of the target field.
   -- Subtract the donor field in the inner cells so we only compute
   -- ghost cell moments.
   fldTar:accumulate(-1., fldDoNoGhosts)
   fldTar:write("fldTar_"..edge..".bp", 0., 0, true)
   numDensityCalcOnGhosts:advance(0., {fldTar}, {m0Tar})
   momDensityCalcOnGhosts:advance(0., {fldTar}, {m1Tar})
   ptclEnergyCalcOnGhosts:advance(0., {fldTar}, {m2Tar})
   m0Tar:write("fldTar_M0_"..edge..".bp",0., 0, true)
   m1Tar:write("fldTar_M1_"..edge..".bp",0., 0, true)
   m2Tar:write("fldTar_M2_"..edge..".bp",0., 0, true)
   
   -- Compute integrated moments.
   intM0Tar[edge] = DataStruct.DynVector { numComponents = 1, }
   intM1Tar[edge] = DataStruct.DynVector { numComponents = 1, }
   intM2Tar[edge] = DataStruct.DynVector { numComponents = 1, }
   intQuantOnGhosts:advance(0., {m0Tar}, {intM0Tar[edge]})
   intQuantOnGhosts:advance(0., {m1Tar}, {intM1Tar[edge]})
   intQuantOnGhosts:advance(0., {m2Tar}, {intM2Tar[edge]})
   intM0Tar[edge]:write("fldTar_intM0_"..edge..".bp",0., 0, true)
   intM1Tar[edge]:write("fldTar_intM1_"..edge..".bp",0., 0, true)
   intM2Tar[edge]:write("fldTar_intM2_"..edge..".bp",0., 0, true)
   
   -- Compute the integral of the target field.
   intFldTar[edge] = DataStruct.DynVector { numComponents = 1, }
   intQuantPhase:advance(0., {fldTar}, {intFldTar[edge]})
   intFldTar[edge]:write("fldTar_intV_"..edge..".bp", 0., 0)

end
