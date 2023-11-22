-- Gkyl ------------------------------------------------------------------------
--
-- Test the interpolation needed for twist-shift BCs in gyrokinetics.
--
-- Create two fields on a 2D grid, the donor field and the target field. Then
-- the donor field with a shift in y (that may be a function of x) assuming
-- periodicity in y to obtain the target field.
-- 
-- This test uses a rectangular donor field and a quadratic shift.
--
-- Some pgkyl commands to check the results:
--   pgkyl rt-twistShift-2x_fldDo.bp -l '$f(x,y)$' rt-twistShift-2x_fldTar.bp -l '$g(x,y)=f(x,y-\mathcal{S}(x))$' rt-twistShift-2x_fldDoBack.bp -l '$g(x,y+\mathcal{S}(x))$' pl
--   pgkyl rt-twistShift-2x_fldDo.bp -l '$f(x,y)$' rt-twistShift-2x_fldTar.bp -l '$g(x,y)=f(x,y-\mathcal{S}(x))$' rt-twistShift-2x_fldDoBack.bp -l '$g(x,y+\mathcal{S}(x))$' interp pl -b
--   pgkyl rt-twistShift-2x_fldDo_intV.bp rt-twistShift-2x_fldTar_intV.bp ev 'f[0] f[1] -' pr
--   pgkyl rt-twistShift-2x_fldDo_intV.bp rt-twistShift-2x_fldDoBack_intV.bp ev 'f[0] f[1] -' pr
-- 
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"
local Mpi        = require "Comm.Mpi"
local Adios      = require "Io.Adios"

-- The following domain and shift are similar to what is used in a CBC ITG benchmark.
local L_ref  = 1.67
local a      = 0.36*L_ref
local x0, y0 = 0.5*a, 0.

local Lx = 0.026648936170213
local Ly = 0.13395216335732
local Lz = 2*math.pi

function qprofile(r)
   return 2.52*(r/a)^2 - 0.16*(r/a) + 0.86
end

local q0 = qprofile(x0)
--print("x0/q0*Lz = ",x0/q0*Lz)

local polyOrder       = 1
local lower           = {x0-Lx/2, y0-Ly/2}
local upper           = {x0+Lx/2, y0+Ly/2}
local cells           = {{10,10}, {20,20}, {40, 40}, {16,64}, {32, 128}, {64, 256}}  -- MF 2021/09/15: 10x10, 20x20 and 40x40 are failing.
local periodicDirs    = {1,2}
local yShiftPolyOrder = polyOrder  -- Poly order for discretizing the shift function, = polyOrder, 1 if polyOrder=0, or 1 or 2 if polyOrder=1.

-- Function describing the shift in y of the BC.
-- It has to be everywhere >0 or everywhere <0 (it cannot be zero, or too close to it).
local yShiftFunc = function(t, xn)
                      local x = xn[1]
                      return (x0/q0)*qprofile(xn[1])*Lz
                   end
local function box(xn)
   local rx2, ry2 = (xn[1]-x0)^2, (xn[2]-y0)^2
   if rx2 < (Lx/4)^2 and ry2 < (Ly/4)^2 then
      return 1.
   end
   return 1.0e-10
end
-- Donor field function.
local fldDoFunc = function(t, xn)
   return box(xn)
end

-- ....................... END OF USER INPUTS (maybe) ........................... --

GKYL_ADIOS2_MPI = GKYL_ADIOS2_MPI or Adios.init_mpi(Mpi.COMM_WORLD)
local yShiftBackFunc = function(t, xn) return -yShiftFunc(t, xn) end

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

-- Shifted donor field function (should be the same as fldDoFunc but shifted in y).
local fldDoShiftedFunc = function(t, xn)
   local x, y = xn[1], xn[2]
   local yS   = wrapNum(y-yShiftFunc(0,xn),{lo=lower[2],up=upper[2]},true)
   return box({x,yS})
end

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
      metaData      = {polyOrder = basis:polyOrder(), basisType = basis:id()},
   }
   return fld
end

local basis = Basis.CartModalSerendipity { ndim = #lower, polyOrder = polyOrder }

local grids = {}
local fldDos, fldDoShifteds, fldTars = {}, {}, {}
local projects, intQuants = {}, {}
local intFldDos, intFldTars, intFldDoBacks = {}, {}, {}
local twistShiftUpds, twistShiftBackUpds = {}, {}

for gI, numCells in ipairs(cells) do 

   grids[gI] = Grid.RectCart {
      lower = lower,  cells        = numCells,
      upper = upper,  periodicDirs = periodicDirs,
   }
   local grid = grids[gI]

   local function fileName(varNm)
      return string.format(varNm .. "_Nx%dNy%dP%d_yShP%d.bp",grid:numCells(1),grid:numCells(2),polyOrder,yShiftPolyOrder)
   end
   
   fldDos[gI]        = createField(grid, basis)
   fldDoShifteds[gI] = createField(grid, basis)
   fldTars[gI]       = createField(grid, basis)
   local fldDo, fldDoShifted, fldTar = fldDos[gI], fldDoShifteds[gI], fldTars[gI]
   
   projects[gI] = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = function(t, xn) return 1. end,  -- Set later.
   }
   local project = projects[gI]
   -- Project donor field function onto basis.
   project:setFunc(function(t,xn) return fldDoFunc(t,xn) end)
   project:advance(0., {}, {fldDo})
   fldDo:write(fileName("fldDo"))
   -- Project shifted donor field function onto basis.
   project:setFunc(function(t,xn) return fldDoShiftedFunc(t,xn) end)
   project:advance(0., {}, {fldDoShifted})
   fldDoShifted:write(fileName("fldDoShifted"))
   
   -- Compute the integral of the donor field.
   intQuants[gI] = Updater.CartFieldIntegrate {
      onGrid = grid,  basis = basis,
   }
   local intQuant = intQuants[gI]
   intFldDos[gI] = DataStruct.DynVector { numComponents = 1, }
   local intFldDo = intFldDos[gI]
   intQuant:advance(0., {fldDo}, {intFldDo})
   intFldDo:write(fileName("fldDo_intV"), 0., 0)
   
   twistShiftUpds[gI] = Updater.TwistShiftBC {
      onGrid = grid,   yShiftFunc      = yShiftFunc, 
      basis  = basis,  yShiftPolyOrder = yShiftPolyOrder, 
   }
   local twistShiftUpd = twistShiftUpds[gI]
   
   -- Apply the forward shift.
   local t1 = os.clock()
   twistShiftUpd:_advance2x(0., {fldDo}, {fldTar})
   local t2 = os.clock()
   --io.write("Forward shift time: ", t2-t1, " s\n")
   
   fldTar:write(fileName("fldTar"))
   
   -- Compute the integral of the target field.
   intFldTars[gI] = DataStruct.DynVector { numComponents = 1, }
   local intFldTar = intFldTars[gI] 
   intQuant:advance(0., {fldTar}, {intFldTar})
   intFldTar:write(fileName("fldTar_intV"), 0., 0)
   
   
   -- ............... SHIFT BACK .................. --
   twistShiftBackUpds[gI] = Updater.TwistShiftBC {
      onGrid = grid,   yShiftFunc      = yShiftBackFunc, 
      basis  = basis,  yShiftPolyOrder = yShiftPolyOrder, 
   }
   local twistShiftBackUpd = twistShiftBackUpds[gI]
   
   fldDo:clear(0.)
   local t1 = os.clock()
   -- Apply the backward shift.
   twistShiftBackUpd:_advance2x(0., {fldTar}, {fldDo})
   local t2 = os.clock()
   --io.write("Backward shift time: ", t2-t1, " s\n")
   
   fldDo:write(fileName("fldDoBack"))
   
   -- Compute the integral of the new target field.
   intFldDoBacks[gI] = DataStruct.DynVector { numComponents = 1, }
   local intFldDoBack = intFldDoBacks[gI]
   intQuant:advance(0., {fldDo}, {intFldDoBack})
   intFldDoBack:write(fileName("fldDoBack_intV"), 0., 0)

end

if GKYL_ADIOS2_MPI then Adios.finalize(GKYL_ADIOS2_MPI);  GKYL_ADIOS2_MPI = nil end
