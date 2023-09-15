-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to fix GK maxwellian with iteration method, 1x1v.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit = require "Unit"
local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"
local Time = require "Lib.Time"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats = Unit.stats

local function createGrid(lo, up, nCells, pDirs)
   pDirs = pDirs or {}
   local gridOut = Grid.RectCart {
      lower = lo,
      cells = nCells,
      upper = up,
      periodicDirs = pDirs,
   }
   return gridOut
end

local function createBasis(dim, pOrder, bKind)
   bKind = bKind or "Ser"
   local basis
   if (bKind=="Ser") then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Max") then
      basis = Basis.CartModalOrder { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Tensor") then
      basis = Basis.CartModalTensor { ndim = dim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

local function createField(grid, basis, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid = grid,
      numComponents = basis:numBasis()*vComp,
      ghost = {1, 1},
      metaData = {
         polyOrder  = basis:polyOrder(),
         basisType  = basis:id(),
      },
   }
   return fld
end

function testGK_1x1v()
   local n0 = 1.0
   local u0 = 0.50
   local vt0 = 1.0
   local B0 = 10.0
   local mass = 1.0

   local lower = {-0.50, -6.0*vt0}
   local upper = {0.50,   6.0*vt0}
   local numCells = {16, 16}

   for polyOrder = 1, 2 do
      local phaseGrid = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid = createGrid({lower[1]}, {upper[1]}, {numCells[1]})
      local confBasis = createBasis(confGrid:ndim(), polyOrder)

      local fM          = createField(phaseGrid, phaseBasis)
      local fOut        = createField(phaseGrid, phaseBasis)
      local moms_in     = createField(confGrid, confBasis, 3)
      local moms_out    = createField(confGrid, confBasis, 3)
      local jacobTot    = createField(confGrid, confBasis)
      local bmag        = createField(confGrid, confBasis)
   
      local eps_n = 0.0
      local eps_u = 0.0
      local eps_vt = 0.0
      local m0Func = function (t, xn) return n0*(1.0+eps_n*math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return u0*(1.0+eps_u*math.cos(2.*math.pi*xn[1])) end
      local vtSqFunc = function (t, xn) return vt0^2 end
      local momsFunc = function (t, xn)
         local m0, udrift, vtsq = m0Func(t, xn), uDriftFunc(t,xn), vtSqFunc(t,xn)
         local m1 = m0*udrift
         local m2 = udrift*m1 + m0*vtsq
         return m0, m1, m2
      end
      local bmagFunc = function (t, xn) return B0 end
      local jacobTotFunc = function (t, xn) return 1.0 end
   
      local projMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }  
      local projPhaseScalar = Updater.ProjectOnBasis {
         onGrid = phaseGrid,
         basis = phaseBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projConfScalar = Updater.ProjectOnBasis {
         onGrid = confGrid,
         basis = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      projMoms:advance(0.0, {}, {moms_in}) 
      projConfScalar:setFunc(function(t,xn) return bmagFunc(t,xn) end)
      projConfScalar:advance(0.0, {}, {bmag}) 
      projConfScalar:setFunc(function(t,xn) return jacobTotFunc(t,xn) end)
      projConfScalar:advance(0.0, {}, {jacobTot}) 
   
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid = phaseGrid,
         confGrid = confGrid,
         phaseBasis = phaseBasis,
         confBasis = confBasis,
         mass = mass,
      }
      maxwellian:advance(0.0, {moms_in, bmag, jacobTot}, {fM})
   
      local iterFix = Updater.CorrectMaxwellian {
         onGrid = phaseGrid,
         phaseBasis = phaseBasis,
         confGrid = confGrid,
         confBasis = confBasis,
         mass = mass,
         bmag = bmag,
         jacobTot = jacobTot,
         iter_max = 100,
         err_max = 1e-14,
         useDevice = GKYL_USE_GPU,
      }
      iterFix:advance(0.0, {fM, moms_in}, {fOut})
 
      local calcMoms = Updater.DistFuncMomentCalc {
         onGrid = phaseGrid,
         confBasis = confBasis,
         phaseBasis = phaseBasis,
         moment = "GkThreeMoments",
         gkfacs = {mass, bmag},
      }  
      calcMoms:advance(0.0, {fOut}, {moms_out})

      if GKYL_USE_GPU then
         moms_out:copyDeviceToHost()
      end

      local indexer    = moms_in:genIndexer()
      local localRange = moms_in:localRange()
      local momsInPtr  = moms_in:get(1)
      local momsOutPtr    = moms_out:get(1)
      for idx in localRange:rowMajorIter() do
         momsInPtr = moms_in:get(indexer(idx))
         momsOutPtr = moms_out:get(indexer(idx))
         for k = 1, confBasis:numBasis()*3 do
            assert_close(momsInPtr[1], momsOutPtr[1], 1.e-14, "Checking mom CorrectMaxwellian 1x1v.")
         end
      end
   end
end

-- Run gyrokinetic tests.
testGK_1x1v()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
