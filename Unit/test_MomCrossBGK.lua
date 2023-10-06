-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to compute the cross moments for the BGK collision operator.
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
   local ms = 1.0
   local mr = 2.0
   local beta = 0.0
   local nu_sr = 10.0
   local nu_rs = 5.0

   local ns = 10.0
   local nr = 10.0
   local us = 0.0
   local ur = 1.0
   local vts = 10.0
   local vtr = 6.0

   local lower = {-0.50, -5.0*vts}
   local upper = {0.50,   5.0*vts}
   local numCells = {8, 8}

   local polyOrder = 1
   local phaseGrid = createGrid(lower, upper, numCells)
   local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
   local confGrid = createGrid({lower[1]}, {upper[1]}, {numCells[1]})
   local confBasis = createBasis(confGrid:ndim(), polyOrder)

   local moms_s = createField(confGrid, confBasis, 3)
   local moms_r = createField(confGrid, confBasis, 3)
   local moms_tar = createField(confGrid, confBasis, 6)
   local moms_cross = createField(confGrid, confBasis, 6)
   
   local eps_n = 0.2
   local eps_u = 0.2
   local eps_vt = 0.2
   local m0sFunc = function (t, xn) return ns*(1.0+eps_n*math.cos(2.*math.pi*xn[1])) end
   local m0rFunc = function (t, xn) return nr*(1.0+eps_n*math.cos(2.*math.pi*xn[1])) end
   local usFunc = function (t, xn) return us*(1.0+eps_u*math.cos(2.*math.pi*xn[1])) end
   local urFunc = function (t, xn) return ur*(1.0+eps_u*math.cos(2.*math.pi*xn[1])) end
   local vtSqsFunc = function (t, xn) return (vts*(1.0+eps_vt*math.cos(2.*math.pi*xn[1])))^2 end
   local vtSqrFunc = function (t, xn) return (vtr*(1.0+eps_vt*math.cos(2.*math.pi*xn[1])))^2 end
   local moms_sFunc = function (t, xn)
      local m0, udrift, vtsq = m0sFunc(t,xn), usFunc(t,xn), vtSqsFunc(t,xn)
      local m1 = m0*udrift
      local m2 = udrift*m1 + m0*vtsq
      return m0, m1, m2
   end
   local moms_rFunc = function (t, xn)
      local m0, udrift, vtsq = m0rFunc(t,xn), urFunc(t,xn), vtSqrFunc(t,xn)
      local m1 = m0*udrift
      local m2 = udrift*m1 + m0*vtsq
      return m0, m1, m2
   end

   local alphaFunc = function (t, xn)
      local m0s, m0r = m0sFunc(t,xn), m0rFunc(t,xn)
      return (1+beta)/(ms+mr)*2*(ms*nu_sr*m0s*mr*nu_rs*m0r)/(ms*nu_sr*m0s+mr*nu_rs*m0r)
   end
   local moms_tarFunc = function (t, xn)
      local alpha = alphaFunc(t,xn)
      local m0s, udrifts, vtsqs = m0sFunc(t,xn), usFunc(t,xn), vtSqsFunc(t,xn)
      local m0r, udriftr, vtsqr = m0rFunc(t,xn), urFunc(t,xn), vtSqrFunc(t,xn)
      local m1s = m0s*udrifts
      local m2s = udrifts*m1s + m0s*vtsqs  
      local m1r = m0r*udriftr
      local m2r = udriftr*m1r + m0r*vtsqr

      local m0sr, m0rs = m0s, m0r
      local m1sr = m1s + alpha*(ms+mr)/(2*ms*nu_sr*m0s*m0r)*(m0s*m1r-m0r*m1s)
      local m1rs = m1r - alpha*(ms+mr)/(2*mr*nu_rs*m0s*m0r)*(m0s*m1r-m0r*m1s)
      local m2sr = m2s + alpha/(ms*nu_sr*m0s*m0r)*(mr*m0s*m2r-ms*m0r*m2s+(ms-mr)*m1s*m1r)
      local m2rs = m2r - alpha/(mr*nu_rs*m0s*m0r)*(mr*m0s*m2r-ms*m0r*m2s+(ms-mr)*m1s*m1r)
      return m0sr, m0rs, m1sr, m1rs, m2sr, m2rs
   end
   
   local projMomsSingle = Updater.ProjectOnBasis {
      onGrid = confGrid,
      basis = confBasis,
      evaluate = function (t, xn) return 1.0 end   -- Set later.
   }
   local projMomsDouble = Updater.ProjectOnBasis {
      onGrid = confGrid,
      basis = confBasis,
      evaluate = function (t, xn) return 1.0 end   -- Set later.
   }
   projMomsSingle:setFunc(function(t,xn) return moms_sFunc(t,xn) end)
   projMomsSingle:advance(0.0, {}, {moms_s})  
   projMomsSingle:setFunc(function(t,xn) return moms_rFunc(t,xn) end)
   projMomsSingle:advance(0.0, {}, {moms_r})  
   projMomsDouble:setFunc(function(t,xn) return moms_tarFunc(t,xn) end)
   projMomsDouble:advance(0.0, {}, {moms_tar})  
 
   local MomCrossBGK = Updater.MomCrossBGK {
      phaseBasis = phaseBasis,
      confBasis = confBasis,
      beta = beta,
      mass_s = ms,
      mass_r = mr,
      nu_sr = nu_sr,
      nu_rs = nu_rs,
      --useDevice = GKYL_USE_GPU,
   }
   MomCrossBGK:advance(0.0, {moms_s, moms_r}, {moms_cross})
 
   if GKYL_USE_GPU then
      moms_cross:copyDeviceToHost()
   end

   local indexer    = moms_cross:genIndexer()
   local localRange = moms_cross:localRange()
   local momsTarPtr  = moms_tar:get(1)
   local momsCrossPtr  = moms_cross:get(1)
   for idx in localRange:rowMajorIter() do
      momsTarPtr = moms_tar:get(indexer(idx))
      momsCrossPtr = moms_cross:get(indexer(idx))
      for k = 1, confBasis:numBasis()*6 do
         assert_close(momsTarPtr[1], momsCrossPtr[1], 1.e-14, "Checking cross moments 1x1v.")
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
