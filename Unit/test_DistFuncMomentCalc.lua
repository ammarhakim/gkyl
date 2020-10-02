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

function test_ser_1x1v()
   local polyOrder = 2 
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
      onGrid     = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "M0",
   }
   local calcNumDensityPvx = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "M0Pvx",
   }
   local calcMomDensity = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "M1i",
   }
   local calcMomDensityPvx = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "M1iPvx",
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
   momDensity:write("momDensity.bp",0.0)
   local momIdxr = momDensity:genIndexer()
   local nItr    = momDensity:get(momIdxr( {1} ))
   assert_equal(1, nItr[1]/math.sqrt(2), "Checking M1i")

   -- Check M1iPvx, momentum density of positive particles only.
   project:setFunc(function (t, xn)
         return 2./(phaseGrid:upper(2)^2-0.0^2)
      end)
   project:advance(0.0, {}, {distf})
   calcMomDensityPvx:advance(0.0, {distf}, {momDensity})
   momDensity:write("momDensityPvx.bp",0.0)
   local momItr = momDensity:get(momIdxr( {1} ))
   assert_equal(1, momItr[1]/math.sqrt(2), "Checking M1iPvx")
end

function test_max_1x1v()
   -- Phase-space and config-space grids.
   local phaseGrid = Grid.RectCart {
      lower = {0.0, -6.0},
      upper = {1.0, 6.0},
      cells = {1, 16},
   }
   local confGrid = Grid.RectCart {
      lower = { phaseGrid:lower(1) },
      upper = { phaseGrid:upper(1) },
      cells = { phaseGrid:numCells(1) },
   }
   -- Basis functions
   local phaseBasis = Basis.CartModalMaxOrder { ndim = 2, polyOrder = 2 }
   local confBasis  = Basis.CartModalMaxOrder { ndim = 1, polyOrder = 2 }
   -- Fields.
   local distf = DataStruct.Field {
      onGrid        = phaseGrid,
      numComponents = phaseBasis:numBasis(),
      ghost         = {0, 0},
   }
   local numDensity = DataStruct.Field {
      onGrid        = confGrid,
      numComponents = confBasis:numBasis(),
      ghost         = {0, 0},
   }

   -- Moment updaters.
   local calcNumDensity = Updater.DistFuncMomentCalc {
      onGrid     = phaseGrid,
      phaseBasis = phaseBasis,
      confBasis  = confBasis,
      moment     = "M0",
   }
   
end

test_ser_1x1v()
--test_max_1x1v()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
