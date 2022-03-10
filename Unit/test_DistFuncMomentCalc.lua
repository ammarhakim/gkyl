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
function bmag_2x(t, xn)
   local x, y = xn[1], xn[2]
   return math.cos((2.*math.pi/(2.*2.*math.pi))*x)*math.exp(-(y^2)/(2.*(math.pi/3.)^2))
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
function distf_2x2v(t, xn)
   local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
   local bmag = bmag_2x(t, xn)
   return bmag*(x^2+y^2)*(vpar-0.5)^2
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
   local bmagFunc
   if cdim == 1     then bmagFunc = bmag_1x
   elseif cdim == 2 then bmagFunc = bmag_2x
   elseif cdim == 3 then bmagFunc = bmag_3x end
   confProj:setFunc(bmagFunc)
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
   elseif pdim == 4 then
      -- 2x2v
      if polyOrder == 1 then
         m0Correct = {
            5.674373976691270e+01,  2.201144875962325e+01,  3.651966323006221e+01,  1.314016650990054e+01,
            2.219569636062632e+02,  6.028640616624954e+01,  4.206235100555175e+01,  3.552942178215428e+00,
            2.219569636062633e+02,  6.028640616624952e+01, -4.206235100555176e+01, -3.552942178215411e+00,
            5.674373976691270e+01,  2.201144875962326e+01, -3.651966323006221e+01, -1.314016650990054e+01,
            7.709667278852976e+01, -3.719961743450785e+00,  4.321049456829391e+01, -4.192642488321794e+00,
            1.403753088009537e+02, -5.979190987388544e+01, -2.255457479022396e+01, -2.512741619288884e+01,
            1.403753088009537e+02, -5.979190987388544e+01,  2.255457479022397e+01,  2.512741619288884e+01,
            7.709667278852976e+01, -3.719961743450784e+00, -4.321049456829390e+01,  4.192642488321796e+00,
            7.709667278852976e+01,  3.719961743450789e+00,  4.321049456829391e+01,  4.192642488321798e+00,
            1.403753088009537e+02,  5.979190987388544e+01, -2.255457479022396e+01,  2.512741619288883e+01,
            1.403753088009538e+02,  5.979190987388544e+01,  2.255457479022397e+01, -2.512741619288886e+01,
            7.709667278852974e+01,  3.719961743450781e+00, -4.321049456829391e+01, -4.192642488321795e+00,
            5.674373976691270e+01, -2.201144875962325e+01,  3.651966323006221e+01, -1.314016650990054e+01,
            2.219569636062632e+02, -6.028640616624953e+01,  4.206235100555175e+01, -3.552942178215412e+00,
            2.219569636062632e+02, -6.028640616624953e+01, -4.206235100555175e+01,  3.552942178215411e+00,
            5.674373976691268e+01, -2.201144875962326e+01, -3.651966323006221e+01,  1.314016650990054e+01
         };
         m1Correct = {
           -4.778420190897913e+01, -1.853595685020905e+01, -3.075340061478923e+01, -1.106540337675835e+01,
           -1.869111272473795e+02, -5.076749992947330e+01, -3.542092716256990e+01, -2.991951307970891e+00,
           -1.869111272473796e+02, -5.076749992947327e+01,  3.542092716256991e+01,  2.991951307970870e+00,
           -4.778420190897913e+01, -1.853595685020905e+01,  3.075340061478923e+01,  1.106540337675835e+01,
           -6.492351392718295e+01,  3.132599362905920e+00, -3.638778489961594e+01,  3.530646305955197e+00,
           -1.182107863586978e+02,  5.035108199906144e+01,  1.899332613913597e+01,  2.115992942559060e+01,
           -1.182107863586978e+02,  5.035108199906144e+01, -1.899332613913598e+01, -2.115992942559060e+01,
           -6.492351392718297e+01,  3.132599362905925e+00,  3.638778489961592e+01, -3.530646305955196e+00,
           -6.492351392718297e+01, -3.132599362905928e+00, -3.638778489961594e+01, -3.530646305955196e+00,
           -1.182107863586979e+02, -5.035108199906144e+01,  1.899332613913597e+01, -2.115992942559060e+01,
           -1.182107863586979e+02, -5.035108199906144e+01, -1.899332613913597e+01,  2.115992942559062e+01,
           -6.492351392718295e+01, -3.132599362905919e+00,  3.638778489961592e+01,  3.530646305955196e+00,
           -4.778420190897911e+01,  1.853595685020907e+01, -3.075340061478922e+01,  1.106540337675835e+01,
           -1.869111272473795e+02,  5.076749992947325e+01, -3.542092716256990e+01,  2.991951307970878e+00,
           -1.869111272473796e+02,  5.076749992947332e+01,  3.542092716256988e+01, -2.991951307970872e+00,
           -4.778420190897911e+01,  1.853595685020906e+01,  3.075340061478923e+01, -1.106540337675835e+01
         };
         m2Correct = {
            1.317741859760154e+02,  5.432242387904527e+01,  8.726478357864879e+01,  3.461290078036008e+01,
            6.282598132762612e+02,  2.349967248260282e+02,  1.585668164431439e+02,  5.344698656702761e+01,
            6.282598132762612e+02,  2.349967248260282e+02, -1.585668164431439e+02, -5.344698656702760e+01,
            1.317741859760154e+02,  5.432242387904527e+01, -8.726478357864879e+01, -3.461290078036008e+01,
            1.892086113026537e+02, -7.382118864615809e+00,  1.138208348715406e+02, -8.593088099826790e+00,
            4.706732142611326e+02, -2.016751811606539e+02, -2.651705819234549e+01, -1.091874194410921e+02,
            4.706732142611326e+02, -2.016751811606539e+02,  2.651705819234548e+01,  1.091874194410921e+02,
            1.892086113026537e+02, -7.382118864615817e+00, -1.138208348715406e+02,  8.593088099826790e+00,
            1.892086113026537e+02,  7.382118864615816e+00,  1.138208348715406e+02,  8.593088099826790e+00,
            4.706732142611326e+02,  2.016751811606539e+02, -2.651705819234548e+01,  1.091874194410921e+02,
            4.706732142611327e+02,  2.016751811606539e+02,  2.651705819234547e+01, -1.091874194410922e+02,
            1.892086113026537e+02,  7.382118864615806e+00, -1.138208348715406e+02, -8.593088099826787e+00,
            1.317741859760154e+02, -5.432242387904528e+01,  8.726478357864879e+01, -3.461290078036008e+01,
            6.282598132762611e+02, -2.349967248260282e+02,  1.585668164431439e+02, -5.344698656702759e+01,
            6.282598132762612e+02, -2.349967248260283e+02, -1.585668164431438e+02,  5.344698656702758e+01,
            1.317741859760154e+02, -5.432242387904529e+01, -8.726478357864879e+01,  3.461290078036008e+01
         };
      elseif polyOrder == 2 then
         m0Correct = {
           5.648341168299724e+01, 2.250387513900885e+01, 3.652502205171847e+01, 1.335144639473130e+01, -7.166399139585262e+00, 7.761081311319963e+00,
           -5.437466693389656e+00, 2.161243931947423e+00, 2.224236975980227e+02, 6.521607861574371e+01, 4.094337722353598e+01, 5.582134934124797e+00,
           -4.387646236656694e+01, -8.976042695293010e+00, -1.237376484382488e+01, -2.753498001080819e+00, 2.224236975980228e+02, 6.521607861574371e+01,
           -4.094337722353598e+01, -5.582134934124794e+00, -4.387646236656691e+01, -8.976042695293000e+00, 1.237376484382486e+01, -2.753498001080819e+00,
           5.648341168299725e+01, 2.250387513900884e+01, -3.652502205171847e+01, -1.335144639473129e+01, -7.166399139585268e+00, 7.761081311319959e+00,
           5.437466693389650e+00, 2.161243931947424e+00, 7.646435675172015e+01, -3.325734694746438e+00, 4.121233955674857e+01, -4.386344478763552e+00,
           1.342837549799004e-01, 4.123263549958735e+00, 5.440652066518922e-01, -2.190342820971765e+00, 1.406410562141511e+02, -5.667019945920631e+01,
           -1.815281995189654e+01, -2.239122325134061e+01, 9.440015133641253e+00, -6.509924435714448e+00, 4.183511326759515e+00, 2.060415698733552e+00,
           1.406410562141512e+02, -5.667019945920630e+01, 1.815281995189653e+01, 2.239122325134061e+01, 9.440015133641275e+00, -6.509924435714441e+00,
           -4.183511326759525e+00, 2.060415698733553e+00, 7.646435675172015e+01, -3.325734694746449e+00, -4.121233955674858e+01, 4.386344478763548e+00,
           1.342837549798988e-01, 4.123263549958738e+00, -5.440652066518882e-01, -2.190342820971765e+00, 7.646435675172016e+01, 3.325734694746450e+00,
           4.121233955674858e+01, 4.386344478763550e+00, 1.342837549798992e-01, 4.123263549958738e+00, 5.440652066518876e-01, 2.190342820971762e+00,
           1.406410562141513e+02, 5.667019945920629e+01, -1.815281995189653e+01, 2.239122325134060e+01, 9.440015133641266e+00, -6.509924435714446e+00,
           4.183511326759531e+00, -2.060415698733550e+00, 1.406410562141513e+02, 5.667019945920630e+01, 1.815281995189654e+01, -2.239122325134062e+01,
           9.440015133641259e+00, -6.509924435714444e+00, -4.183511326759524e+00, -2.060415698733558e+00, 7.646435675172017e+01, 3.325734694746441e+00,
           -4.121233955674857e+01, -4.386344478763546e+00, 1.342837549798972e-01, 4.123263549958737e+00, -5.440652066518941e-01, 2.190342820971768e+00,
           5.648341168299724e+01, -2.250387513900884e+01, 3.652502205171847e+01, -1.335144639473129e+01, -7.166399139585268e+00, 7.761081311319965e+00,
           -5.437466693389656e+00, -2.161243931947421e+00, 2.224236975980227e+02, -6.521607861574371e+01, 4.094337722353598e+01, -5.582134934124787e+00,
           -4.387646236656693e+01, -8.976042695293009e+00, -1.237376484382487e+01, 2.753498001080815e+00, 2.224236975980228e+02, -6.521607861574373e+01,
           -4.094337722353598e+01, 5.582134934124794e+00, -4.387646236656688e+01, -8.976042695292994e+00, 1.237376484382487e+01, 2.753498001080822e+00,
           5.648341168299726e+01, -2.250387513900885e+01, -3.652502205171848e+01, 1.335144639473129e+01, -7.166399139585263e+00, 7.761081311319963e+00,
           5.437466693389656e+00, -2.161243931947425e+00,
         }
         m1Correct = {
           -4.756497825936609e+01, -1.895063169600745e+01, -3.075791330671029e+01, -1.124332327977373e+01, 6.034862433334959e+00, -6.535647420058911e+00,
           4.578919320749178e+00, -1.819994890060987e+00, -1.873041663983350e+02, -5.491880304483679e+01, -3.447863345139874e+01, -4.700745207684037e+00,
           3.694859988763532e+01, 7.558772796036234e+00, 1.042001250006305e+01, 2.318735158804893e+00, -1.873041663983350e+02, -5.491880304483681e+01,
           3.447863345139872e+01, 4.700745207684039e+00, 3.694859988763526e+01, 7.558772796036222e+00, -1.042001250006304e+01, 2.318735158804894e+00,
           -4.756497825936611e+01, -1.895063169600745e+01, 3.075791330671029e+01, 1.124332327977372e+01, 6.034862433334967e+00, -6.535647420058911e+00,
           -4.578919320749177e+00, -1.819994890060993e+00, -6.439103726460642e+01, 2.800618690312792e+00, -3.470512804778827e+01, 3.693763771590361e+00,
           -1.130810568251719e-01, -3.472221936807352e+00, -4.581601740226504e-01, 1.844499217660437e+00, -1.184345736540220e+02, 4.772227322880535e+01,
           1.528658522264971e+01, 1.885576694849735e+01, -7.949486428329483e+00, 5.482041630075337e+00, -3.522956906744856e+00, -1.735086904196671e+00,
           -1.184345736540220e+02, 4.772227322880529e+01, -1.528658522264972e+01, -1.885576694849735e+01, -7.949486428329486e+00, 5.482041630075328e+00,
           3.522956906744861e+00, -1.735086904196669e+00, -6.439103726460642e+01, 2.800618690312799e+00, 3.470512804778828e+01, -3.693763771590359e+00,
           -1.130810568251768e-01, -3.472221936807356e+00, 4.581601740226479e-01, 1.844499217660435e+00, -6.439103726460642e+01, -2.800618690312798e+00,
           -3.470512804778826e+01, -3.693763771590358e+00, -1.130810568251744e-01, -3.472221936807353e+00, -4.581601740226480e-01, -1.844499217660433e+00,
           -1.184345736540221e+02, -4.772227322880531e+01, 1.528658522264971e+01, -1.885576694849735e+01, -7.949486428329481e+00, 5.482041630075334e+00,
           -3.522956906744867e+00, 1.735086904196668e+00, -1.184345736540222e+02, -4.772227322880531e+01, -1.528658522264971e+01, 1.885576694849737e+01,
           -7.949486428329483e+00, 5.482041630075329e+00, 3.522956906744854e+00, 1.735086904196673e+00, -6.439103726460645e+01, -2.800618690312791e+00,
           3.470512804778826e+01, 3.693763771590357e+00, -1.130810568251665e-01, -3.472221936807356e+00, 4.581601740226530e-01, -1.844499217660438e+00,
           -4.756497825936611e+01, 1.895063169600745e+01, -3.075791330671029e+01, 1.124332327977372e+01, 6.034862433334961e+00, -6.535647420058913e+00,
           4.578919320749178e+00, 1.819994890060988e+00, -1.873041663983349e+02, 5.491880304483681e+01, -3.447863345139872e+01, 4.700745207684028e+00,
           3.694859988763531e+01, 7.558772796036230e+00, 1.042001250006304e+01, -2.318735158804888e+00, -1.873041663983350e+02, 5.491880304483683e+01,
           3.447863345139870e+01, -4.700745207684037e+00, 3.694859988763528e+01, 7.558772796036221e+00, -1.042001250006305e+01, -2.318735158804896e+00,
           -4.756497825936611e+01, 1.895063169600746e+01, 3.075791330671029e+01, -1.124332327977372e+01, 6.034862433334961e+00, -6.535647420058913e+00,
           -4.578919320749179e+00, 1.819994890060991e+00,
         }
         m2Correct = {
           1.346810209205576e+02, 5.636184644070686e+01, 9.081689031496188e+01, 3.638800942322894e+01, -1.531915442293607e+01, 2.162068038060000e+01,
           -1.154920978322388e+01, 8.056050757541808e+00, 6.437526165222573e+02, 2.383699952856996e+02, 1.543494567891024e+02, 4.915181136302768e+01,
           -1.058058544975506e+02, -2.407700215608105e+01, -3.508892496659295e+01, -1.032850505731126e+01, 6.437526165222575e+02, 2.383699952856996e+02,
           -1.543494567891024e+02, -4.915181136302770e+01, -1.058058544975505e+02, -2.407700215608102e+01, 3.508892496659290e+01, -1.032850505731125e+01,
           1.346810209205576e+02, 5.636184644070686e+01, -9.081689031496188e+01, -3.638800942322894e+01, -1.531915442293607e+01, 2.162068038059999e+01,
           1.154920978322387e+01, 8.056050757541811e+00, 1.921238511758901e+02, -6.770111068912961e+00, 1.140587007866537e+02, -9.496106164226617e+00,
           -2.226967471797892e-01, 1.915556929571275e+01, 6.944928133740057e-01, -5.022783835435162e+00, 4.823471699659661e+02, -1.893777713359310e+02,
           -1.637529707819997e+01, -9.468924055391847e+01, 2.310532509390621e+01, -3.594756731257781e+01, 1.302051527839764e+01, 3.653556204232846e-01,
           4.823471699659663e+02, -1.893777713359309e+02, 1.637529707819996e+01, 9.468924055391845e+01, 2.310532509390625e+01, -3.594756731257780e+01,
           -1.302051527839767e+01, 3.653556204232817e-01, 1.921238511758900e+02, -6.770111068912972e+00, -1.140587007866538e+02, 9.496106164226617e+00,
           -2.226967471797914e-01, 1.915556929571276e+01, -6.944928133739985e-01, -5.022783835435162e+00, 1.921238511758901e+02, 6.770111068912968e+00,
           1.140587007866538e+02, 9.496106164226621e+00, -2.226967471797906e-01, 1.915556929571276e+01, 6.944928133739963e-01, 5.022783835435160e+00,
           4.823471699659665e+02, 1.893777713359309e+02, -1.637529707819991e+01, 9.468924055391845e+01, 2.310532509390623e+01, -3.594756731257780e+01,
           1.302051527839768e+01, -3.653556204232722e-01, 4.823471699659665e+02, 1.893777713359309e+02, 1.637529707819994e+01, -9.468924055391849e+01,
           2.310532509390620e+01, -3.594756731257780e+01, -1.302051527839765e+01, -3.653556204232893e-01, 1.921238511758901e+02, 6.770111068912962e+00,
           -1.140587007866538e+02, -9.496106164226616e+00, -2.226967471797981e-01, 1.915556929571276e+01, -6.944928133740064e-01, 5.022783835435166e+00,
           1.346810209205576e+02, -5.636184644070686e+01, 9.081689031496187e+01, -3.638800942322894e+01, -1.531915442293607e+01, 2.162068038060000e+01,
           -1.154920978322388e+01, -8.056050757541808e+00, 6.437526165222574e+02, -2.383699952856996e+02, 1.543494567891024e+02, -4.915181136302766e+01,
           -1.058058544975506e+02, -2.407700215608102e+01, -3.508892496659293e+01, 1.032850505731125e+01, 6.437526165222575e+02, -2.383699952856996e+02,
           -1.543494567891024e+02, 4.915181136302769e+01, -1.058058544975505e+02, -2.407700215608101e+01, 3.508892496659291e+01, 1.032850505731126e+01,
           1.346810209205576e+02, -5.636184644070686e+01, -9.081689031496188e+01, 3.638800942322894e+01, -1.531915442293607e+01, 2.162068038060000e+01,
           1.154920978322388e+01, -8.056050757541813e+00,
         }
      end
   end

   if pdim < 4 then
      -- Check M0, number density.
      calcNumDensity:advance(0.0, {distf}, {numDensity})
--      numDensity:write("numDensity.bp",0.0)
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
--      momDensity:write("momDensity.bp",0.0)
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
--      keDensity:write("keDensity.bp",0.0)
      local momIdxr = keDensity:genIndexer()
      local m2Ptr = keDensity:get(1)
      for i = 1, confNumCells[1] do
        keDensity:fill(momIdxr({i}), m2Ptr)
        for k = 1, confBasis:numBasis() do
           assert_close(m2Correct[i][k], m2Ptr:data()[k-1], 1.e-13, "Checking M1")
         end
      end

   elseif pdim==4 then

      -- Check M0, number density.
      calcNumDensity:advance(0.0, {distf}, {numDensity})
--      numDensity:write("numDensity.bp",0.0)
      local momIdxr = numDensity:genIndexer()
      local m0Ptr = numDensity:get(1)
      for i = 0, numCells[1]-1 do
         for j = 0, numCells[2]-1 do
            local idx = {i+1, j+1} -- Account for ghost cells.
            local linidx = momIdxr(idx)
            numDensity:fill(linidx, m0Ptr)
            for k = 1, confBasis:numBasis() do
               assert_close(m0Correct[(i*numCells[2]+j)*confBasis:numBasis()+k], m0Ptr:data()[k-1], 1.e-13, "Checking M0")
            end
         end
      end

      -- Check M1, momentum density.
      calcMomDensity:advance(0.0, {distf}, {momDensity})
--      momDensity:write("momDensity.bp",0.0)
      local momIdxr = momDensity:genIndexer()
      local m1Ptr = momDensity:get(1)
      for i = 0, numCells[1]-1 do
         for j = 0, numCells[2]-1 do
            local idx = {i+1, j+1} -- Account for ghost cells.
            local linidx = momIdxr(idx)
            momDensity:fill(linidx, m1Ptr)
            for k = 1, confBasis:numBasis() do
               assert_close(m1Correct[(i*numCells[2]+j)*confBasis:numBasis()+k], m1Ptr:data()[k-1], 1.e-13, "Checking M1")
            end
         end
      end

      -- Check M2, momentum density.
      calcKeDensity:advance(0.0, {distf}, {keDensity})
--      keDensity:write("keDensity.bp",0.0)
      local momIdxr = keDensity:genIndexer()
      local m2Ptr = keDensity:get(1)
      for i = 0, numCells[1]-1 do
         for j = 0, numCells[2]-1 do
            local idx = {i+1, j+1} -- Account for ghost cells.
            local linidx = momIdxr(idx)
            keDensity:fill(linidx, m2Ptr)
            for k = 1, confBasis:numBasis() do
               assert_close(m2Correct[(i*numCells[2]+j)*confBasis:numBasis()+k], m2Ptr:data()[k-1], 1.e-12, "Checking M2")
            end
         end
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
-- Test 2x2v p=1
lower    = {-math.pi, -math.pi, -2.0, 0.0}
upper    = { math.pi,  math.pi,  2.0, 2.0}
numCells = {4, 4, 2, 2}
polyOrder, cdim = 1, 2 
test_gk_ser(lower, upper, numCells, polyOrder, cdim)
-- Test 2x2v p=2
polyOrder = 2
test_gk_ser(lower, upper, numCells, polyOrder, cdim)

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
