-- Gkyl ------------------------------------------------------------------------
--
-- Test updater that computes sqrt(f)^q where f is a CartField using quadrature.
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
local Time       = require "Lib.Time"

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
   bKind = bKind or "Ser"
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
         polyOrder  = basis:polyOrder(),
         basisType  = basis:id(),
      }
   }
   return fld
end

function test_1x()
   local lower    = {-0.50}
   local upper    = { 0.50}
   local numCells = {8}

   -- Tolerances for polyOrder=1-3 and for sqrt(f) and sqrt(f)^3.
   local tols = {{1.e-14,1.e-14},{1.e-14,1.e-14},{1.e-10,1.e-9}}

   for polyOrder = 1, 3 do

      local grid  = createGrid(lower, upper, numCells)
      local basis = createBasis(grid:ndim(), polyOrder)

      local fld = createField(grid, basis)
      local sqrtFld   = createField(grid, basis)
      local sqrtFldR3 = createField(grid, basis)
      local sqrtFldProj   = createField(grid, basis)
      local sqrtFldR3Proj = createField(grid, basis)

      local projScalarFunc = Updater.ProjectOnBasis {
         onGrid   = grid,
         basis    = basis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }

      local fFunc = function(t, xn)  -- Function to project onto basis.
         local x = xn[1]
         return 1.+0.25*math.cos((2.*math.pi/(upper[1]-lower[1]))*x)
      end
      projScalarFunc:setFunc(fFunc)
      projScalarFunc:advance(0., {}, {fld})
      projScalarFunc:setFunc(function(t,xn) return math.sqrt(fFunc(t,xn)) end)
      projScalarFunc:advance(0., {}, {sqrtFldProj})
      projScalarFunc:setFunc(function(t,xn) return (math.sqrt(fFunc(t,xn)))^3 end)
      projScalarFunc:advance(0., {}, {sqrtFldR3Proj})

      local sqrtUpd = Updater.SqrtOnBasis {onGrid = grid,  basis = basis,}

      sqrtUpd:advance(0., {fld}, {sqrtFld})  -- Compute sqrt(f).
      sqrtUpd:advance(0., {fld, 3.0}, {sqrtFldR3})  -- Compute sqrt(f)^3.

      local indexer = fld:genIndexer()
      local sqrtFldPtr, sqrtFldR3Ptr = sqrtFld:get(1), sqrtFldR3:get(1)
      local sqrtFldProjPtr, sqrtFldR3ProjPtr = sqrtFldProj:get(1), sqrtFldR3Proj:get(1)
      local localRange = fld:localRange()
      for idx in localRange:rowMajorIter() do
         sqrtFld:fill(indexer(idx), sqrtFldPtr)
         sqrtFldR3:fill(indexer(idx), sqrtFldR3Ptr)
         sqrtFldProj:fill(indexer(idx), sqrtFldProjPtr)
         sqrtFldR3Proj:fill(indexer(idx), sqrtFldR3ProjPtr)
         for k = 1, basis:numBasis() do
            assert_close(sqrtFldPtr[1], sqrtFldProjPtr[1], tols[polyOrder][1], string.format("Checking 1x p%d sqrt(f) vs ProjectOnBasis.",polyOrder))
            assert_close(sqrtFldR3Ptr[1], sqrtFldR3ProjPtr[1], tols[polyOrder][2], string.format("Checking 1x p%d sqrt(f)^3 vs ProjectOnBasis.",polyOrder))
         end
      end
   end
end

function test_2x()
   local lower    = {-0.50, -0.50}
   local upper    = { 0.50,  0.50}
   local numCells = {8, 8}

   -- Tolerances for polyOrder=1-3 and for sqrt(f) and sqrt(f)^3.
   local tols = {{1.e-14,1.e-14},{1.e-7,1.e-7},{1.e-7,1.e-6}}

   for polyOrder = 1, 3 do

      local grid  = createGrid(lower, upper, numCells)
      local basis = createBasis(grid:ndim(), polyOrder)

      local fld = createField(grid, basis)
      local sqrtFld   = createField(grid, basis)
      local sqrtFldR3 = createField(grid, basis)
      local sqrtFldProj   = createField(grid, basis)
      local sqrtFldR3Proj = createField(grid, basis)

      local projScalarFunc = Updater.ProjectOnBasis {
         onGrid   = grid,
         basis    = basis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }

      local fFunc = function(t, xn)  -- Function to project onto basis.
         local x, y = xn[1], xn[2]
         return 1.+0.25*math.cos((2.*math.pi/(upper[1]-lower[1]))*x)
                       *math.sin(2.*(2.*math.pi/(upper[2]-lower[2]))*y)
      end
      projScalarFunc:setFunc(fFunc)
      projScalarFunc:advance(0., {}, {fld})
      projScalarFunc:setFunc(function(t,xn) return math.sqrt(fFunc(t,xn)) end)
      projScalarFunc:advance(0., {}, {sqrtFldProj})
      projScalarFunc:setFunc(function(t,xn) return (math.sqrt(fFunc(t,xn)))^3 end)
      projScalarFunc:advance(0., {}, {sqrtFldR3Proj})

      local sqrtUpd = Updater.SqrtOnBasis {onGrid = grid,  basis = basis,}

      sqrtUpd:advance(0., {fld}, {sqrtFld})  -- Compute sqrt(f).
      sqrtUpd:advance(0., {fld, 3.0}, {sqrtFldR3})  -- Compute sqrt(f)^3.

      local indexer = fld:genIndexer()
      local sqrtFldPtr, sqrtFldR3Ptr = sqrtFld:get(1), sqrtFldR3:get(1)
      local sqrtFldProjPtr, sqrtFldR3ProjPtr = sqrtFldProj:get(1), sqrtFldR3Proj:get(1)
      local localRange = fld:localRange()
      for idx in localRange:rowMajorIter() do
         sqrtFld:fill(indexer(idx), sqrtFldPtr)
         sqrtFldR3:fill(indexer(idx), sqrtFldR3Ptr)
         sqrtFldProj:fill(indexer(idx), sqrtFldProjPtr)
         sqrtFldR3Proj:fill(indexer(idx), sqrtFldR3ProjPtr)
         for k = 1, basis:numBasis() do
            assert_close(sqrtFldPtr[1], sqrtFldProjPtr[1], tols[polyOrder][1], string.format("Checking 2x p%d sqrt(f) vs ProjectOnBasis.",polyOrder))
            assert_close(sqrtFldR3Ptr[1], sqrtFldR3ProjPtr[1], tols[polyOrder][2], string.format("Checking 2x p%d sqrt(f)^3 vs ProjectOnBasis.",polyOrder))
         end
      end
   end
end

function test_3x()
   local lower    = {-0.50, -0.50, -1.0}
   local upper    = { 0.50,  0.50,  1.0}
   local numCells = {8, 8, 12}

   -- Tolerances for polyOrder=1-3 and for sqrt(f) and sqrt(f)^3.
   local tols = {{1.e-14,1.e-14},{1.e-7,1.e-6},{1.e-7,1.e-6}}

   for polyOrder = 3, 3 do

      local grid  = createGrid(lower, upper, numCells)
      local basis = createBasis(grid:ndim(), polyOrder)

      local fld = createField(grid, basis)
      local sqrtFld   = createField(grid, basis)
      local sqrtFldR3 = createField(grid, basis)
      local sqrtFldProj   = createField(grid, basis)
      local sqrtFldR3Proj = createField(grid, basis)

      local projScalarFunc = Updater.ProjectOnBasis {
         onGrid   = grid,
         basis    = basis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }

      local fFunc = function(t, xn)  -- Function to project onto basis.
         local x, y, z = xn[1], xn[2], xn[3]
         return (1.+0.25*math.cos((2.*math.pi/(upper[1]-lower[1]))*x)
                       *math.sin(2.*(2.*math.pi/(upper[2]-lower[2]))*y))*(1.+0.25*z)
      end
      projScalarFunc:setFunc(fFunc)
      projScalarFunc:advance(0., {}, {fld})
      projScalarFunc:setFunc(function(t,xn) return math.sqrt(fFunc(t,xn)) end)
      projScalarFunc:advance(0., {}, {sqrtFldProj})
      projScalarFunc:setFunc(function(t,xn) return (math.sqrt(fFunc(t,xn)))^3 end)
      projScalarFunc:advance(0., {}, {sqrtFldR3Proj})

      local sqrtUpd = Updater.SqrtOnBasis {onGrid = grid,  basis = basis,}

      sqrtUpd:advance(0., {fld}, {sqrtFld})  -- Compute sqrt(f).
      sqrtUpd:advance(0., {fld, 3.0}, {sqrtFldR3})  -- Compute sqrt(f)^3.

      local indexer = fld:genIndexer()
      local sqrtFldPtr, sqrtFldR3Ptr = sqrtFld:get(1), sqrtFldR3:get(1)
      local sqrtFldProjPtr, sqrtFldR3ProjPtr = sqrtFldProj:get(1), sqrtFldR3Proj:get(1)
      local localRange = fld:localRange()
      for idx in localRange:rowMajorIter() do
         sqrtFld:fill(indexer(idx), sqrtFldPtr)
         sqrtFldR3:fill(indexer(idx), sqrtFldR3Ptr)
         sqrtFldProj:fill(indexer(idx), sqrtFldProjPtr)
         sqrtFldR3Proj:fill(indexer(idx), sqrtFldR3ProjPtr)
         for k = 1, basis:numBasis() do
            assert_close(sqrtFldPtr[1], sqrtFldProjPtr[1], tols[polyOrder][1], string.format("Checking 3x p%d sqrt(f) vs ProjectOnBasis.",polyOrder))
            assert_close(sqrtFldR3Ptr[1], sqrtFldR3ProjPtr[1], tols[polyOrder][2], string.format("Checking 3x p%d sqrt(f)^3 vs ProjectOnBasis.",polyOrder))
         end
      end
   end
end

-- Run tests.
test_1x()
test_2x()
test_3x()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
