-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to compute moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"
local Lin        = require "Lib.Linalg"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats        = Unit.stats

local function createGrid(lo, up, nCells, pDirs)
   pDirs = pDirs or {}
   local gridOut = Grid.RectCart {
      lower        = lo,
      upper        = up,
      cells        = nCells,
      periodicDirs = pDirs,
   }
   return gridOut
end

local function createBasis(dim, pOrder, bKind)
   local basis
   if (bKind=="Ser") then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Max") then
      basis = Basis.CartModalMaxOrder { ndim = dim, polyOrder = pOrder }
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
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
      metaData = {
         polyOrder = basis:polyOrder(),
         basisType = basis:id()
      }
   }
   fld:clear(0.0)
   return fld
end

function test_vm_ser_1x1v(polyOrder)
   -- Note that for partial moments the v=0 point has to lie on a cell boundary.
   local lower     = {0.0, -4.0}
   local upper     = {1.0, 12.0}
   local numCells  = {4, 16}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(lower, upper, numCells)
   local confGrid  = createGrid({phaseGrid:lower(1)}, {phaseGrid:upper(1)}, {phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder, "Ser")
   local confBasis  = createBasis(confGrid:ndim(), polyOrder, "Ser")
   -- Fields.
   local distf      = createField(phaseGrid, phaseBasis)
   local numDensity = createField(confGrid, confBasis)
   local momDensity = createField(confGrid, confBasis, phaseGrid:ndim()-confGrid:ndim())

   -- Updater to initialize distribution function.
   local project = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,
      basis    = phaseBasis,
      evaluate = function (t, xn) return 1 end   -- Set below.
   }
   project:setFunc(function (t, xn) return 1/(phaseGrid:upper(2)-phaseGrid:lower(2)) end)
   project:advance(0.0, {}, {distf})

   -- Moment updaters.
   local calcNumDensity = Updater.DistFuncMomentCalc {
      advanceArgs = {{distf}, {numDensity}},
      onGrid      = phaseGrid,
      phaseBasis  = phaseBasis,
      confBasis   = confBasis,
      moment      = "M0",
   }
   local calcMomDensity = Updater.DistFuncMomentCalc {
      advanceArgs = {{distf}, {momDensity}},
      onGrid      = phaseGrid,
      phaseBasis  = phaseBasis,
      confBasis   = confBasis,
      moment      = "M1i",
   }
   local calcNumDensityPvx = Updater.DistFuncMomentCalc {
      advanceArgs = {{distf}, {numDensity}},
      onGrid      = phaseGrid,
      phaseBasis  = phaseBasis,
      confBasis   = confBasis,
      moment      = "M0Pvx",
   }
   local calcMomDensityPvx = Updater.DistFuncMomentCalc {
      advanceArgs = {{distf}, {numDensity}},
      onGrid      = phaseGrid,
      phaseBasis  = phaseBasis,
      confBasis   = confBasis,
      moment      = "M1iPvx",
   }

   -- Check M0, number density.
   calcNumDensity:advance(0.0, {distf}, {numDensity})
   local momIdxr = numDensity:genIndexer()
   local nItr    = numDensity:get(momIdxr( {1} ))
   assert_equal(1, nItr[1]/math.sqrt(2), "Checking M0")

   -- Check M0Pvx, number density of positive particles only.
   project:setFunc(function (t, xn) return 1/(phaseGrid:upper(2)-0.0) end)
   project:advance(0.0, {}, {distf})
   calcNumDensityPvx:advance(0.0, {distf}, {numDensity})
   nItr    = numDensity:get(momIdxr( {1} ))
   assert_equal(1, nItr[1]/math.sqrt(2), "Checking M0Pvx")
   
   -- Check M1i, momentum density.
   project:setFunc(function (t, xn) 
         return 2./(phaseGrid:upper(2)^2-phaseGrid:lower(2)^2)
      end)
   project:advance(0.0, {}, {distf})
   calcMomDensity:advance(0.0, {distf}, {momDensity})
--   momDensity:write("momDensity.bp",0.0)
   local momIdxr = momDensity:genIndexer()
   local nItr    = momDensity:get(momIdxr( {1} ))
   assert_equal(1, nItr[1]/math.sqrt(2), "Checking M1i")

   -- Check M1iPvx, momentum density of positive particles only.
   project:setFunc(function (t, xn)
         return 2./(phaseGrid:upper(2)^2-0.0^2)
      end)
   project:advance(0.0, {}, {distf})
   calcMomDensityPvx:advance(0.0, {distf}, {momDensity})
--   momDensity:write("momDensityPvx.bp",0.0)
   local momItr = momDensity:get(momIdxr( {1} ))
   assert_equal(1, momItr[1]/math.sqrt(2), "Checking M1iPvx")
end

function bmag_1x(t, xn)
   local x = xn[1]
   return math.cos((2.*math.pi/(2.*2.*math.pi))*x)
end
function distf_1x1v(t, xn)
   local x, vpar = xn[1], xn[2]
   local bmag = bmag_1x(t, xn)
   return bmag*(x^2)*(vpar-0.5)^2
end
function distf_1x2v(t, xn)
   local x, vpar, mu = xn[1], xn[2], xn[3]
   local bmag = bmag_1x(t, xn)
   return bmag*(x^2)*(vpar-0.5)^2
end

function test_gk_ser(lower, upper, numCells, polyOrder, cdim)
   local mass = 1.0
   local pdim = #lower
   local vdim = pdim - cdim
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(lower, upper, numCells)
   local confLower, confUpper, confNumCells = {}, {}, {}
   for d = 1, cdim do
      confLower[d], confUpper[d], confNumCells[d] = lower[d], upper[d], numCells[d]
   end
   local confGrid = createGrid(confLower, confUpper, confNumCells)
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder, "Ser")
   local confBasis  = createBasis(confGrid:ndim(), polyOrder, "Ser")
   -- Fields.
   local distf      = createField(phaseGrid, phaseBasis)
   local numDensity = createField(confGrid, confBasis)
   local momDensity = createField(confGrid, confBasis)
   local keDensity  = createField(confGrid, confBasis)
   local bmag       = createField(confGrid, confBasis)

   -- Updater to initialize distribution function.
   local phaseProj = Updater.ProjectOnBasis {
      onGrid   = phaseGrid,  basis = phaseBasis,
      evaluate = function (t, xn) return 1 end   -- Set below.
   }
   local distfFunc
   if pdim == 2     then distfFunc = distf_1x1v
   elseif pdim == 3 then distfFunc = distf_1x2v
   elseif pdim == 4 then distfFunc = distf_2x2v
   elseif pdim == 5 then distfFunc = distf_3x2v end
   phaseProj:setFunc(distfFunc)
   phaseProj:advance(0.0, {}, {distf})

   -- Project the magnetic field magnitude onto basis.
   local confProj = Updater.ProjectOnBasis {
      onGrid   = confGrid,  basis = confBasis,
      evaluate = function (t, xn) return 1 end   -- Set below.
   }
   local distfFunc
   if cdim == 1     then distfFunc = bmag_1x
   elseif cdim == 2 then distfFunc = bmag_2x
   elseif cdim == 3 then distfFunc = bmag_3x end
   confProj:setFunc(bmag_1x)
   confProj:advance(0.0, {}, {bmag})

   -- Moment updaters.
   local calcNumDensity = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,   confBasis = confBasis,
      phaseBasis = phaseBasis,  moment    = "GkM0",
      gkfacs     = {mass, bmag},
   }
   local calcMomDensity = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,   confBasis = confBasis,
      phaseBasis = phaseBasis,  moment    = "GkM1",
      gkfacs     = {mass, bmag},
   }
   local calcKeDensity = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,   confBasis = confBasis,
      phaseBasis = phaseBasis,  moment    = "GkM2",
      gkfacs     = {mass, bmag},
   }

   -- Correct answers.
   local m0Correct, m1Correct, m2Correct
   if pdim == 2 then
      -- 1x1v
      if polyOrder == 1 then
         m0Correct = {
            { 1.52537436025689e+01,  3.57234787161818e+00},
            { 6.08286277145113e+00, -5.10949091314880e+00},
            { 6.08286277145113e+00,  5.10949091314880e+00},
            { 1.52537436025689e+01, -3.57234787161818e+00},
         }
         m1Correct = {
            {-1.28452577705844e+01, -3.00829294452057e+00},
            {-5.12241075490621e+00,  4.30272919002004e+00},
            {-5.12241075490621e+00, -4.30272919002004e+00},
            {-1.28452577705844e+01,  3.00829294452057e+00},
         }
         m2Correct = {
            {3.31835825740096e+01,  7.77142344001148e+00},
            {1.32328944501744e+01, -1.11153837408851e+01},
            {1.32328944501744e+01,  1.11153837408851e+01},
            {3.31835825740096e+01, -7.77142344001148e+00},
         }
      elseif polyOrder == 2 then
         m0Correct = {
            {1.526837339934706e+01,   3.951518219554417e+00,  -3.363344534446567e+00},
            {6.052510010088350e+00,  -4.868229034295940e+00,   8.480389975048731e-01},
            {6.052510010088350e+00,   4.868229034295939e+00,   8.480389975048728e-01},
            {1.526837339934706e+01,  -3.951518219554418e+00,  -3.363344534446568e+00},
         }
         m1Correct = {
            {-1.285757759945016e+01,  -3.327594290151089e+00,   2.832290134270792e+00},
            {-5.096850534811242e+00,   4.099561292038686e+00,  -7.141381031619991e-01},
            {-5.096850534811242e+00,  -4.099561292038685e+00,  -7.141381031619988e-01},
            {-1.285757759945016e+01,   3.327594290151089e+00,   2.832290134270792e+00},
         }
         m2Correct = {
            {3.407258063854292e+01,   8.818124868900387e+00,  -7.316749513532877e+00 },
            {1.350665391724979e+01,  -1.086383742390252e+01,   1.844856766501832e+00 },
            {1.350665391724979e+01,   1.086383742390252e+01,   1.844856766501832e+00 },
            {3.407258063854292e+01,  -8.818124868900387e+00,  -7.316749513532877e+00 },
         }
      end
   elseif pdim == 3 then
      -- 1x2v
      if polyOrder == 1 then
         m0Correct = {
            { 191.6841953662915,   44.89144731817122},
            {  76.4395079823428,   -64.2077564653283},
            { 76.43950798234282,   64.20775646532829},
            {191.68419536629153,  -44.89144731817124},
         }
         m1Correct = {
            {-161.41826978214021,   -37.80332405740736},
            { -64.37011198513079,   54.069689655013306},
            { -64.37011198513079,   -54.06968965501329},
            {-161.41826978214021,    37.80332405740738},
         }
         m2Correct = {
            {578.5971330556217,   210.75432886828457},
            {292.8699479183419,  -242.1332009483718 },
            {292.8699479183419,   242.13320094837178},
            {578.5971330556217,  -210.75432886828463},
         }
      elseif polyOrder == 2 then
         m0Correct = {
            {1.918680388146181e+02,  4.965624243631351e+01, -4.226503392363477e+01},
            {7.605808393388897e+01, -6.117597028054660e+01,  1.065677233807588e+01},
            {7.605808393388898e+01,  6.117597028054661e+01,  1.065677233807589e+01},
            {1.918680388146182e+02, -4.965624243631353e+01, -4.226503392363478e+01},
         }
         m1Correct = {
            {-1.615730853175732e+02, -4.181578310426401e+01,  3.559160751463980e+01},
            {-6.404891278643279e+01,  5.151660655203925e+01, -8.974124074169159e+00},
            {-6.404891278643282e+01, -5.151660655203925e+01, -8.974124074169165e+00},
            {-1.615730853175732e+02,  4.181578310426401e+01,  3.559160751463981e+01},
         }
         m2Correct = {
            {5.924940346248225e+02,  2.106245140015062e+02, -1.080258558862571e+02},
            {2.957803339447422e+02, -2.297438543767285e+02,  2.952994007812927e+01},
            {2.957803339447422e+02,  2.297438543767285e+02,  2.952994007812929e+01},
            {5.924940346248227e+02, -2.106245140015063e+02, -1.080258558862570e+02},
         }
      end
   end

   -- Check M0, number density.
   calcNumDensity:advance(0.0, {distf}, {numDensity})
--   numDensity:write("numDensity.bp",0.0)
   local momIdxr = numDensity:genIndexer()
   local m0Ptr = numDensity:get(1)
   for i = 1, confNumCells[1] do
      numDensity:fill(momIdxr({i}), m0Ptr)
      for k = 1, confBasis:numBasis() do
         assert_close(m0Correct[i][k], m0Ptr:data()[k-1], 1.e-13, "Checking M0")
      end
   end

   -- Check M1, momentum density.
   calcMomDensity:advance(0.0, {distf}, {momDensity})
--   momDensity:write("momDensity.bp",0.0)
   local momIdxr = momDensity:genIndexer()
   local m1Ptr = momDensity:get(1)
   for i = 1, confNumCells[1] do
      momDensity:fill(momIdxr({i}), m1Ptr)
      for k = 1, confBasis:numBasis() do
         assert_close(m1Correct[i][k], m1Ptr:data()[k-1], 1.e-13, "Checking M1")
      end
   end

   -- Check M2, momentum density.
   calcKeDensity:advance(0.0, {distf}, {keDensity})
--   keDensity:write("keDensity.bp",0.0)
   local momIdxr = keDensity:genIndexer()
   local m2Ptr = keDensity:get(1)
   for i = 1, confNumCells[1] do
      keDensity:fill(momIdxr({i}), m2Ptr)
      for k = 1, confBasis:numBasis() do
         assert_close(m2Correct[i][k], m2Ptr:data()[k-1], 1.e-13, "Checking M1")
      end
   end

end

test_vm_ser_1x1v(2)

local lower, upper, numcells, polyOrder, cdim
-- Test 1x1v p=1
lower    = {-math.pi, -2.0}
upper    = { math.pi,  2.0}
numCells = {4, 2}
polyOrder, cdim = 1, 1 
test_gk_ser(lower, upper, numCells, polyOrder, cdim)
-- Test 1x1v p=2
polyOrder = 2
test_gk_ser(lower, upper, numCells, polyOrder, cdim)
-- Test 1x2v p=1
lower    = {-math.pi, -2.0, 0.0}
upper    = { math.pi,  2.0, 2.0}
numCells = {4, 2, 2}
polyOrder, cdim = 1, 1 
test_gk_ser(lower, upper, numCells, polyOrder, cdim)
-- Test 1x2v p=2
polyOrder = 2
test_gk_ser(lower, upper, numCells, polyOrder, cdim)

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
