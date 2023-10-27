-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to compute integrated moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit       = require "Unit"
local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

local function createGrid(lo,up,nCells)
   local gridOut = Grid.RectCart {
      lower = lo,
      upper = up,
      cells = nCells,
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
      ghost         = {0, 0},
   }
   return fld
end

local function createProject(grid,basis)
   local projUp = Updater.ProjectOnBasis {
      onGrid   = grid,
      basis    = basis,
      evaluate = function(t, xn) return 1 end
   }
   return projUp
end

local function distFmoment(fIn,momOut,grid,pBasis,cBasis,mom)
   local momUp = Updater.DistFuncMomentCalc {
      advanceArgs = {{fIn}, {momOut}},
      onGrid      = grid,
      phaseBasis  = pBasis,
      confBasis   = cBasis,
      moment      = mom,
   }
   return momUp
end

local function intQuantCalc(grid,basis,op,vComp) 
   vComp = vComp or 1
   local intQuant = Updater.CartFieldIntegratedQuantCalc {
      onGrid        = grid,
      basis         = basis,
      numComponents = vComp,
      operator      = op,
   }
   return intQuant
end

local function intDGmom(gridIn,basisIn,momIn,addBasisIn)
   addBasisIn = false or addBasisIn
   local intMom
   if addBasisIn then
      intMom = Updater.IntegratedDGMoment {
         onGrid    = gridIn,
         basis     = basisIn,
         confBasis = addBasisIn,
         moment    = momIn,
      }
   else
      intMom = Updater.IntegratedDGMoment {
         onGrid = gridIn,
         basis  = basisIn,
         moment = momIn,
      }
   end
   return intMom
end

local function test_one(pOrder, basis)
   local pLower = {0.0, -6.0}
   local pUpper = {1.0,  6.0}
   local pN     = {1, 16}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1)},{phaseGrid:upper(1)},{phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local distf      = createField(phaseGrid,phaseBasis)
   local numDensity = createField(confGrid,confBasis)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/((phaseGrid:upper(1)-phaseGrid:lower(1))*(phaseGrid:upper(2)-phaseGrid:lower(2)))
      end)
   project:advance(0.0, {}, {distf})

   -- Moment updater.
   local calcNumDensity = distFmoment(distf,numDensity,phaseGrid,phaseBasis,confBasis,"M0")
   calcNumDensity:advance(0.0, {distf}, {numDensity})

   -- Compute integrated f and integrated moment with CartFieldIntegratedQuantCalc.
   local intF  = DataStruct.DynVector {numComponents = 1,}
   local intM0 = DataStruct.DynVector {numComponents = 1,}
   local intQuantP = intQuantCalc(phaseGrid,phaseBasis,"none")
   local intQuantC = intQuantCalc(confGrid,confBasis,"none")
   intQuantP:advance(0.0, {distf}, {intF})
   intQuantC:advance(0.0, {numDensity}, {intM0})
   local _, NP = intF:lastData()
   local _, NC = intM0:lastData()

   -- Compute integral of f and of zeroth moment with IntegratedDGMoment.
   local intFn   = DataStruct.DynVector {numComponents = 1,}
   local intM0n  = DataStruct.DynVector {numComponents = 1,}
   local intMomP = intDGmom(phaseGrid,phaseBasis,"one")
   local intMomC = intDGmom(confGrid,confBasis,"one")
   intMomP:advance(0.0, {distf}, {intFn})
   intMomC:advance(0.0, {numDensity}, {intM0n})
   local _, NPn = intFn:lastData()
   local _, NCn = intM0n:lastData()

   assert_equal(NP[1], NPn[1], "Checking 'one' moment")
   assert_equal(NP[1], NCn[1], "Checking 'one' moment")
end

local function test_x1(pOrder, basis)
   local cLower = {1.0}
   local cUpper = {3.0}
   local cN     = {16}
   -- Config-space grid.
   local confGrid = createGrid(cLower, cUpper, cN)
   -- Basis functions.
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local confFld = createField(confGrid,confBasis)

   -- Updater to initialize distribution function
   local project = createProject(confGrid,confBasis)
   project:setFunc(function (t, xn)
         return 1/( (0.5*(confGrid:upper(1)^2-confGrid:lower(1)^2)) )
      end
   )
   project:advance(0.0, {}, {confFld})

   -- Compute integral of x_1*f with IntegratedDGMoment.
   local intF       = DataStruct.DynVector {numComponents = 1,}
   local intMom     = intDGmom(confGrid,confBasis,"x1")
   intMom:advance(0.0, {confFld}, {intF})
   local _, intM = intF:lastData()

   assert_equal(1, intM[1], "Checking x1 moment")
end

local function test_x2v1(pOrder, basis)
   local pLower = {0.0, -2.0}
   local pUpper = {1.0,  6.0}
   local pN     = {1, 16}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1)},{phaseGrid:upper(1)},{phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local vDim       = phaseBasis:ndim()-confBasis:ndim()
   local distf      = createField(phaseGrid,phaseBasis)
   local momDensity = createField(confGrid,confBasis,vDim)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *(0.5*(phaseGrid:upper(2)^2-phaseGrid:lower(2)^2)) )
      end
   )
   project:advance(0.0, {}, {distf})

   -- Moment updater.
   local calcMomDensity = distFmoment(distf,momDensity,phaseGrid,phaseBasis,confBasis,"M1i")
   calcMomDensity:advance(0.0, {distf}, {momDensity})

   -- Compute integrated M1_1 moment with CartFieldIntegratedQuantCalc.
   local intM1_1   = DataStruct.DynVector {numComponents = vDim,}
   local intQuantC = intQuantCalc(confGrid,confBasis,"none",vDim)
   intQuantC:advance(0.0, {momDensity}, {intM1_1})
   local _, intM1  = intM1_1:lastData()

   -- Compute integral of v_1*f with IntegratedDGMoment.
   local intFn     = DataStruct.DynVector {numComponents = 1,}
   local intMomP   = intDGmom(phaseGrid,phaseBasis,"v1",confBasis)
   intMomP:advance(0.0, {distf}, {intFn})
   local _, intM1n = intFn:lastData()

   local intFnA     = DataStruct.DynVector {numComponents = 1,}
   local intMomPA   = intDGmom(phaseGrid,phaseBasis,"x2",confBasis)
   intMomPA:advance(0.0, {distf}, {intFnA})
   local _, intM1nA = intFnA:lastData()

   assert_equal(intM1[1], intM1n[1], "Checking v1 moment")
   assert_equal(intM1n[1], intM1nA[1], "Checking x2 moment")
end

local function test_x3v2(pOrder, basis)
   local pLower = {0.0, -2.0, -2.0}
   local pUpper = {1.0,  6.0, 6.0}
   local pN     = {1, 16, 16}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower,pUpper,pN)
   local confGrid  = createGrid({phaseGrid:lower(1)},{phaseGrid:upper(1)},{phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local vDim       = phaseBasis:ndim()-confBasis:ndim()
   local distf      = createField(phaseGrid,phaseBasis)
   local momDensity = createField(confGrid,confBasis,vDim)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *(phaseGrid:upper(2)-phaseGrid:lower(2)) 
                   *(0.5*(phaseGrid:upper(3)^2-phaseGrid:lower(3)^2)) )
      end
   )
   project:advance(0.0, {}, {distf})

   -- Moment updater.
   local calcMomDensity = distFmoment(distf,momDensity,phaseGrid,phaseBasis,confBasis,"M1i")
   calcMomDensity:advance(0.0, {distf}, {momDensity})

   -- Compute integrated M1_1 moment with CartFieldIntegratedQuantCalc.
   local intM1_1   = DataStruct.DynVector {numComponents = vDim,}
   local intQuantC = intQuantCalc(confGrid,confBasis,"none",vDim)
   intQuantC:advance(0.0, {momDensity}, {intM1_1})
   local _, intM1  = intM1_1:lastData()

   -- Compute integral of v_2*f with IntegratedDGMoment.
   local intFn     = DataStruct.DynVector {numComponents = 1,}
   local intMomP   = intDGmom(phaseGrid,phaseBasis,"v2",confBasis)
   intMomP:advance(0.0, {distf}, {intFn})
   local _, intM1n = intFn:lastData()

   local intFnA     = DataStruct.DynVector {numComponents = 1,}
   local intMomPA   = intDGmom(phaseGrid,phaseBasis,"x3",confBasis)
   intMomPA:advance(0.0, {distf}, {intFnA})
   local _, intM1nA = intFnA:lastData()

   assert_equal(intM1[2], intM1n[1], "Checking v2 moment")
   assert_equal(intM1n[1], intM1nA[1], "Checking x3 moment")
end

local function test_x4v3(pOrder, basis)
   local pLower = {0.0, -2.0, -2.0, -2.0}
   local pUpper = {1.0,  6.0, 6.0, 6.0}
   local pN     = {1, 8, 8, 8}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1)},{phaseGrid:upper(1)},{phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local vDim       = phaseBasis:ndim()-confBasis:ndim()
   local distf      = createField(phaseGrid,phaseBasis)
   local momDensity = createField(confGrid,confBasis,vDim)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *(phaseGrid:upper(2)-phaseGrid:lower(2))
                   *(phaseGrid:upper(3)-phaseGrid:lower(3))
                   *(0.5*(phaseGrid:upper(4)^2-phaseGrid:lower(4)^2)) )
      end
   )
   project:advance(0.0, {}, {distf})

   -- Moment updater.
   local calcMomDensity = distFmoment(distf,momDensity,phaseGrid,phaseBasis,confBasis,"M1i")
   calcMomDensity:advance(0.0, {distf}, {momDensity})

   -- Compute integrated M1_1 moment with CartFieldIntegratedQuantCalc.
   local intM1_1   = DataStruct.DynVector {numComponents = vDim,}
   local intQuantC = intQuantCalc(confGrid,confBasis,"none",vDim)
   intQuantC:advance(0.0, {momDensity}, {intM1_1})
   local _, intM1  = intM1_1:lastData()

   -- Compute integral of v_3*f with IntegratedDGMoment.
   local intFn     = DataStruct.DynVector {numComponents = 1,}
   local intMomP   = intDGmom(phaseGrid,phaseBasis,"v3",confBasis)
   intMomP:advance(0.0, {distf}, {intFn})
   local _, intM1n = intFn:lastData()

   local intFnA     = DataStruct.DynVector {numComponents = 1,}
   local intMomPA   = intDGmom(phaseGrid,phaseBasis,"x4",confBasis)
   intMomPA:advance(0.0, {distf}, {intFnA})
   local _, intM1nA = intFnA:lastData()

   assert_equal(intM1[3], intM1n[1], "Checking v3 moment")
   assert_equal(intM1n[1], intM1nA[1], "Checking x4 moment")
end

local function test_x5v3(pOrder, basis)
   local pLower = {0.0, 0.0, -2.0, -2.0, -2.0}
   local pUpper = {1.0, 1.0, 6.0, 6.0, 6.0}
   local pN     = {2, 2, 8, 8, 8}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1), phaseGrid:lower(2)},
                                {phaseGrid:upper(1), phaseGrid:upper(2)},
                                {phaseGrid:numCells(1), phaseGrid:numCells(2)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local vDim       = phaseBasis:ndim()-confBasis:ndim()
   local distf      = createField(phaseGrid,phaseBasis)
   local momDensity = createField(confGrid,confBasis,vDim)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *(phaseGrid:upper(2)-phaseGrid:lower(2))
                   *(phaseGrid:upper(3)-phaseGrid:lower(3))
                   *(phaseGrid:upper(4)-phaseGrid:lower(4))
                   *(0.5*(phaseGrid:upper(5)^2-phaseGrid:lower(5)^2)) )
      end
   )
   project:advance(0.0, {}, {distf})

   -- Moment updater.
   local calcMomDensity = distFmoment(distf,momDensity,phaseGrid,phaseBasis,confBasis,"M1i")
   calcMomDensity:advance(0.0, {distf}, {momDensity})

   -- Compute integrated M1_1 moment with CartFieldIntegratedQuantCalc.
   local intM1_1   = DataStruct.DynVector {numComponents = vDim,}
   local intQuantC = intQuantCalc(confGrid,confBasis,"none",vDim)
   intQuantC:advance(0.0, {momDensity}, {intM1_1})
   local _, intM1  = intM1_1:lastData()

   -- Compute integral of v_3*f with IntegratedDGMoment.
   local intFn     = DataStruct.DynVector {numComponents = 1,}
   local intMomP   = intDGmom(phaseGrid,phaseBasis,"v3",confBasis)
   intMomP:advance(0.0, {distf}, {intFn})
   local _, intM1n = intFn:lastData()

   local intFnA     = DataStruct.DynVector {numComponents = 1,}
   local intMomPA   = intDGmom(phaseGrid,phaseBasis,"x5",confBasis)
   intMomPA:advance(0.0, {distf}, {intFnA})
   local _, intM1nA = intFnA:lastData()

   assert_equal(intM1[3], intM1n[1], "Checking v3 moment")
   assert_equal(intM1n[1], intM1nA[1], "Checking x5 moment")
end

local function test_x6v3(pOrder, basis)
   local pLower = {0.0, 0.0, 0.0, -2.0, -2.0, -2.0}
   local pUpper = {1.0, 1.0, 1.0, 6.0, 6.0, 6.0}
   local pN     = {1, 1, 1, 8, 8, 8}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1),   phaseGrid:lower(2),   phaseGrid:lower(3)},   
                                {phaseGrid:upper(1),   phaseGrid:upper(2),   phaseGrid:upper(3)},   
                                {phaseGrid:numCells(1),phaseGrid:numCells(2),phaseGrid:numCells(3)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local vDim       = phaseBasis:ndim()-confBasis:ndim()
   local distf      = createField(phaseGrid,phaseBasis)
   local momDensity = createField(confGrid,confBasis,vDim)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *(phaseGrid:upper(2)-phaseGrid:lower(2))
                   *(phaseGrid:upper(3)-phaseGrid:lower(3))
                   *(phaseGrid:upper(4)-phaseGrid:lower(4))
                   *(phaseGrid:upper(5)-phaseGrid:lower(5))
                   *(0.5*(phaseGrid:upper(6)^2-phaseGrid:lower(6)^2)) )
      end
   )
   project:advance(0.0, {}, {distf})

   -- Compute integral of v_3*f with IntegratedDGMoment.
   local intFn     = DataStruct.DynVector {numComponents = 1,}
   local intMomP   = intDGmom(phaseGrid,phaseBasis,"v3",confBasis)
   intMomP:advance(0.0, {distf}, {intFn})
   local _, intM1n = intFn:lastData()

   -- Compute integral of x_6*f with IntegratedDGMoment.
   local intFnA     = DataStruct.DynVector {numComponents = 1,}
   local intMomPA   = intDGmom(phaseGrid,phaseBasis,"x6",confBasis)
   intMomPA:advance(0.0, {distf}, {intFnA})
   local _, intM1nA = intFnA:lastData()

   assert_equal(1, intM1n[1], "Checking v3 moment")
   assert_equal(intM1n[1], intM1nA[1], "Checking x6 moment")
end

local function test_x1SqxSq(pOrder, basis)
   local cLower = {1.0}
   local cUpper = {3.0}
   local cN     = {16}
   -- Config-space grid.
   local confGrid = createGrid(cLower, cUpper, cN)
   -- Basis functions.
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local confFld = createField(confGrid,confBasis)

   -- Updater to initialize distribution function
   local project = createProject(confGrid,confBasis)
   project:setFunc(function (t, xn)
         return 1/( ((1/3)*(confGrid:upper(1)^3-confGrid:lower(1)^3)) )
      end
   )
   project:advance(0.0, {}, {confFld})

   -- Compute integral of xSq*f with IntegratedDGMoment.
   local intxSqF = DataStruct.DynVector {numComponents = 1,}
   local intMom  = intDGmom(confGrid,confBasis,"xSq")
   intMom:advance(0.0, {confFld}, {intxSqF})
   local _, intM = intxSqF:lastData()

   -- Compute integral of x_1^2*f with IntegratedDGMoment.
   local intx1SqF   = DataStruct.DynVector {numComponents = 1,}
   local intMom1    = intDGmom(confGrid,confBasis,"x1Sq")
   intMom1:advance(0.0, {confFld}, {intx1SqF})
   local _, intM1 = intx1SqF:lastData()

   assert_equal(1, intM1[1], "Checking x1Sq moment")
   assert_equal(intM1[1], intM[1], "Checking xSq moment")
end

local function test_x2SqxSq(pOrder, basis)
   local cLower = {1.0, 1.0}
   local cUpper = {3.0, 3.0}
   local cN     = {16, 16}
   -- Config-space grid.
   local confGrid  = createGrid(cLower, cUpper, cN)
   -- Basis functions.
   local confBasis = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local confFld = createField(confGrid,confBasis)

   -- Updater to initialize distribution function.
   local project = createProject(confGrid,confBasis)
   project:setFunc(function (t, xn)
         return 1/( (confGrid:upper(1)-confGrid:lower(1))
                   *((1/3)*(confGrid:upper(2)^3-confGrid:lower(2)^3)) )
      end
   )
   project:advance(0.0, {}, {confFld})

   -- Compute integral of xSq*f with IntegratedDGMoment.
   local intxSqF = DataStruct.DynVector {numComponents = 1,}
   local intMom  = intDGmom(confGrid,confBasis,"xSq")
   intMom:advance(0.0, {confFld}, {intxSqF})
   local _, intM = intxSqF:lastData()

   -- Compute integral of x_1^2*f with IntegratedDGMoment.
   local intx1SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom1  = intDGmom(confGrid,confBasis,"x1Sq")
   intMom1:advance(0.0, {confFld}, {intx1SqF})
   local _, intM1 = intx1SqF:lastData()

   -- Compute integral of x_2^2*f with IntegratedDGMoment.
   local intx2SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom2  = intDGmom(confGrid,confBasis,"x2Sq")
   intMom2:advance(0.0, {confFld}, {intx2SqF})
   local _, intM2 = intx2SqF:lastData()

   assert_equal(1, intM2[1], "Checking x2Sq moment")
   assert_equal(intM1[1]+intM2[1], intM[1], "Checking xSq moment")
end

local function test_x3SqxSq(pOrder, basis)
   local cLower = {1.0, 1.0, 1.0}
   local cUpper = {3.0, 3.0, 3.0}
   local cN     = {8, 8, 8}
   -- Config-space grid.
   local confGrid  = createGrid(cLower, cUpper, cN)
   -- Basis functions.
   local confBasis = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local confFld = createField(confGrid,confBasis)

   -- Updater to initialize distribution function.
   local project = createProject(confGrid,confBasis)
   project:setFunc(function (t, xn)
         return 1/( (confGrid:upper(1)-confGrid:lower(1))
                   *(confGrid:upper(2)-confGrid:lower(2))
                   *((1/3)*(confGrid:upper(3)^3-confGrid:lower(3)^3)) )
      end
   )
   project:advance(0.0, {}, {confFld})

   -- Compute integral of xSq*f with IntegratedDGMoment.
   local intxSqF = DataStruct.DynVector {numComponents = 1,}
   local intMom  = intDGmom(confGrid,confBasis,"xSq")
   intMom:advance(0.0, {confFld}, {intxSqF})
   local _, intM = intxSqF:lastData()

   -- Compute integral of x_1^2*f with IntegratedDGMoment.
   local intx1SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom1  = intDGmom(confGrid,confBasis,"x1Sq")
   intMom1:advance(0.0, {confFld}, {intx1SqF})
   local _, intM1 = intx1SqF:lastData()

   -- Compute integral of x_2^2*f with IntegratedDGMoment.
   local intx2SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom2  = intDGmom(confGrid,confBasis,"x2Sq")
   intMom2:advance(0.0, {confFld}, {intx2SqF})
   local _, intM2 = intx2SqF:lastData()

   -- Compute integral of x_3^2*f with IntegratedDGMoment.
   local intx3SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom3  = intDGmom(confGrid,confBasis,"x3Sq")
   intMom3:advance(0.0, {confFld}, {intx3SqF})
   local _, intM3 = intx3SqF:lastData()

   assert_equal(1, intM3[1], "Checking x3Sq moment")
   assert_equal(intM1[1]+intM2[1]+intM3[1], intM[1], "Checking xSq moment")
end

local function test_x4SqxSq(pOrder, basis)
   local cLower = {1.0, 1.0, 1.0, 1.0}
   local cUpper = {3.0, 3.0, 3.0, 3.0}
   local cN     = {1, 8, 8, 8}
   -- Config-space grid.
   local confGrid  = createGrid(cLower, cUpper, cN)
   -- Basis functions.
   local confBasis = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local confFld = createField(confGrid,confBasis)

   -- Updater to initialize distribution function.
   local project = createProject(confGrid,confBasis)
   project:setFunc(function (t, xn)
         return 1/( (confGrid:upper(1)-confGrid:lower(1))
                   *(confGrid:upper(2)-confGrid:lower(2))
                   *(confGrid:upper(3)-confGrid:lower(3))
                   *((1/3)*(confGrid:upper(4)^3-confGrid:lower(4)^3)) )
      end
   )
   project:advance(0.0, {}, {confFld})

   -- Compute integral of xSq*f with IntegratedDGMoment.
   local intxSqF = DataStruct.DynVector {numComponents = 1,}
   local intMom  = intDGmom(confGrid,confBasis,"xSq")
   intMom:advance(0.0, {confFld}, {intxSqF})
   local _, intM = intxSqF:lastData()

   -- Compute integral of x_1^2*f with IntegratedDGMoment.
   local intx1SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom1  = intDGmom(confGrid,confBasis,"x1Sq")
   intMom1:advance(0.0, {confFld}, {intx1SqF})
   local _, intM1 = intx1SqF:lastData()

   -- Compute integral of x_2^2*f with IntegratedDGMoment.
   local intx2SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom2  = intDGmom(confGrid,confBasis,"x2Sq")
   intMom2:advance(0.0, {confFld}, {intx2SqF})
   local _, intM2 = intx2SqF:lastData()

   -- Compute integral of x_3^2*f with IntegratedDGMoment.
   local intx3SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom3  = intDGmom(confGrid,confBasis,"x3Sq")
   intMom3:advance(0.0, {confFld}, {intx3SqF})
   local _, intM3 = intx3SqF:lastData()

   -- Compute integral of x_4^2*f with IntegratedDGMoment.
   local intx4SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom4  = intDGmom(confGrid,confBasis,"x4Sq")
   intMom4:advance(0.0, {confFld}, {intx4SqF})
   local _, intM4 = intx4SqF:lastData()

   assert_equal(1, intM4[1], "Checking x3Sq moment")
   assert_equal(intM1[1]+intM2[1]+intM3[1]+intM4[1], intM[1], "Checking xSq moment")
end

local function test_x5SqxSq(pOrder, basis)
   local cLower = {1.0, 1.0, 1.0, 1.0, 1.0}
   local cUpper = {3.0, 3.0, 3.0, 3.0, 3.0}
   local cN     = {1, 1, 8, 8, 8}
   -- Config-space grid.
   local confGrid  = createGrid(cLower, cUpper, cN)
   -- Basis functions.
   local confBasis = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local confFld = createField(confGrid,confBasis)

   -- Updater to initialize distribution function.
   local project = createProject(confGrid,confBasis)
   project:setFunc(function (t, xn)
         return 1/( (confGrid:upper(1)-confGrid:lower(1))
                   *(confGrid:upper(2)-confGrid:lower(2))
                   *(confGrid:upper(3)-confGrid:lower(3))
                   *(confGrid:upper(4)-confGrid:lower(4))
                   *((1/3)*(confGrid:upper(5)^3-confGrid:lower(5)^3)) )
      end
   )
   project:advance(0.0, {}, {confFld})

   -- Compute integral of xSq*f with IntegratedDGMoment.
   local intxSqF = DataStruct.DynVector {numComponents = 1,}
   local intMom  = intDGmom(confGrid,confBasis,"xSq")
   intMom:advance(0.0, {confFld}, {intxSqF})
   local _, intM = intxSqF:lastData()

   -- Compute integral of x_1^2*f with IntegratedDGMoment.
   local intx1SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom1  = intDGmom(confGrid,confBasis,"x1Sq")
   intMom1:advance(0.0, {confFld}, {intx1SqF})
   local _, intM1 = intx1SqF:lastData()

   -- Compute integral of x_2^2*f with IntegratedDGMoment.
   local intx2SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom2  = intDGmom(confGrid,confBasis,"x2Sq")
   intMom2:advance(0.0, {confFld}, {intx2SqF})
   local _, intM2 = intx2SqF:lastData()

   -- Compute integral of x_3^2*f with IntegratedDGMoment.
   local intx3SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom3  = intDGmom(confGrid,confBasis,"x3Sq")
   intMom3:advance(0.0, {confFld}, {intx3SqF})
   local _, intM3 = intx3SqF:lastData()

   -- Compute integral of x_4^2*f with IntegratedDGMoment.
   local intx4SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom4  = intDGmom(confGrid,confBasis,"x4Sq")
   intMom4:advance(0.0, {confFld}, {intx4SqF})
   local _, intM4 = intx4SqF:lastData()

   -- Compute integral of x_5^2*f with IntegratedDGMoment.
   local intx5SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom5  = intDGmom(confGrid,confBasis,"x5Sq")
   intMom5:advance(0.0, {confFld}, {intx5SqF})
   local _, intM5 = intx5SqF:lastData()

   assert_equal(1, intM5[1], "Checking x5Sq moment")
   assert_equal(intM1[1]+intM2[1]+intM3[1]+intM4[1]+intM5[1], intM[1], "Checking xSq moment")
end

local function test_x6SqxSq(pOrder, basis)
   local cLower = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
   local cUpper = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0}
   local cN     = {1, 1, 1, 8, 8, 8}
   -- Config-space grid.
   local confGrid  = createGrid(cLower, cUpper, cN)
   -- Basis functions.
   local confBasis = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local confFld = createField(confGrid,confBasis)

   -- Updater to initialize distribution function.
   local project = createProject(confGrid,confBasis)
   project:setFunc(function (t, xn)
         return 1/( (confGrid:upper(1)-confGrid:lower(1))
                   *(confGrid:upper(2)-confGrid:lower(2))
                   *(confGrid:upper(3)-confGrid:lower(3))
                   *(confGrid:upper(4)-confGrid:lower(4))
                   *(confGrid:upper(5)-confGrid:lower(5))
                   *((1/3)*(confGrid:upper(6)^3-confGrid:lower(6)^3)) )
      end
   )
   project:advance(0.0, {}, {confFld})

   -- Compute integral of xSq*f with IntegratedDGMoment.
   local intxSqF = DataStruct.DynVector {numComponents = 1,}
   local intMom  = intDGmom(confGrid,confBasis,"xSq")
   intMom:advance(0.0, {confFld}, {intxSqF})
   local _, intM = intxSqF:lastData()

   -- Compute integral of x_1^2*f with IntegratedDGMoment.
   local intx1SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom1  = intDGmom(confGrid,confBasis,"x1Sq")
   intMom1:advance(0.0, {confFld}, {intx1SqF})
   local _, intM1 = intx1SqF:lastData()

   -- Compute integral of x_2^2*f with IntegratedDGMoment.
   local intx2SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom2  = intDGmom(confGrid,confBasis,"x2Sq")
   intMom2:advance(0.0, {confFld}, {intx2SqF})
   local _, intM2 = intx2SqF:lastData()

   -- Compute integral of x_3^2*f with IntegratedDGMoment.
   local intx3SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom3  = intDGmom(confGrid,confBasis,"x3Sq")
   intMom3:advance(0.0, {confFld}, {intx3SqF})
   local _, intM3 = intx3SqF:lastData()

   -- Compute integral of x_4^2*f with IntegratedDGMoment.
   local intx4SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom4  = intDGmom(confGrid,confBasis,"x4Sq")
   intMom4:advance(0.0, {confFld}, {intx4SqF})
   local _, intM4 = intx4SqF:lastData()

   -- Compute integral of x_5^2*f with IntegratedDGMoment.
   local intx5SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom5  = intDGmom(confGrid,confBasis,"x5Sq")
   intMom5:advance(0.0, {confFld}, {intx5SqF})
   local _, intM5 = intx5SqF:lastData()

   -- Compute integral of x_6^2*f with IntegratedDGMoment.
   local intx6SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom6  = intDGmom(confGrid,confBasis,"x6Sq")
   intMom6:advance(0.0, {confFld}, {intx6SqF})
   local _, intM6 = intx6SqF:lastData()

   assert_equal(1, intM6[1], "Checking x6Sq moment")
   assert_equal(intM1[1]+intM2[1]+intM3[1]+intM4[1]+intM5[1]+intM6[1], intM[1], "Checking xSq moment")
end

local function test_v1SqvSq(pOrder, basis)
   local pLower = {0.0, -6.0}
   local pUpper = {1.0,  6.0}
   local pN     = {1, 16}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1)},{phaseGrid:upper(1)},{phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local distf         = createField(phaseGrid,phaseBasis)
   local energyDensity = createField(confGrid,confBasis)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *((1/3)*(phaseGrid:upper(2)^3-phaseGrid:lower(2)^3)) )
      end
   )
   project:advance(0.0, {}, {distf})

   -- Moment updater.
   local calcEnergyDensity = distFmoment(distf,energyDensity,phaseGrid,phaseBasis,confBasis,"M2")
   calcEnergyDensity:advance(0.0, {distf}, {energyDensity})

   -- Compute integrated f and integrated moment with CartFieldIntegratedQuantCalc.
   local intM2     = DataStruct.DynVector {numComponents = 1,}
   local intQuantC = intQuantCalc(confGrid,confBasis,"none")
   intQuantC:advance(0.0, {energyDensity}, {intM2})
   local _, NC     = intM2:lastData()

   -- Compute integral of vSq*f with IntegratedDGMoment.
   local intFn   = DataStruct.DynVector {numComponents = 1,}
   local intMomP = intDGmom(phaseGrid,phaseBasis,"vSq",confBasis)
   intMomP:advance(0.0, {distf}, {intFn})
   local _, NPn  = intFn:lastData()

   -- Compute integral of v1Sq*f with IntegratedDGMoment.
   local intFnA   = DataStruct.DynVector {numComponents = 1,}
   local intMomPA = intDGmom(phaseGrid,phaseBasis,"v1Sq",confBasis)
   intMomPA:advance(0.0, {distf}, {intFnA})
   local _, NPnA  = intFnA:lastData()

   assert_equal(NC[1], NPn[1], "Checking vSq moment")
   assert_equal(NPn[1], NPnA[1], "Checking v1Sq moment")
end

local function test_v2SqvSq(pOrder, basis)
   local pLower = {0.0, -6.0, -6.0}
   local pUpper = {1.0,  6.0,  6.0}
   local pN     = {1, 16, 16}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1)},{phaseGrid:upper(1)},{phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local distf         = createField(phaseGrid,phaseBasis)
   local energyDensity = createField(confGrid,confBasis)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *(phaseGrid:upper(2)-phaseGrid:lower(2))
                   *((1/3)*(phaseGrid:upper(3)^3-phaseGrid:lower(3)^3)) )
      end
   )
   project:advance(0.0, {}, {distf})

   -- Compute integral of vSq*f with IntegratedDGMoment.
   local intFn   = DataStruct.DynVector {numComponents = 1,}
   local intMomP = intDGmom(phaseGrid,phaseBasis,"vSq",confBasis)
   intMomP:advance(0.0, {distf}, {intFn})
   local _, NPn  = intFn:lastData()

   -- Compute integral of v1Sq*f with IntegratedDGMoment.
   local intv1SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom1  = intDGmom(phaseGrid,phaseBasis,"v1Sq",confBasis)
   intMom1:advance(0.0, {distf}, {intv1SqF})
   local _, NP1   = intv1SqF:lastData()

   -- Compute integral of v2Sq*f with IntegratedDGMoment.
   local intv2SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom2  = intDGmom(phaseGrid,phaseBasis,"v2Sq",confBasis)
   intMom2:advance(0.0, {distf}, {intv2SqF})
   local _, NP2   = intv2SqF:lastData()

   assert_equal(1, NP2[1], "Checking v2Sq moment")
   assert_equal(NPn[1], NP1[1]+NP2[1], "Checking v1Sq moment")
end

local function test_v3SqvSq(pOrder, basis)
   local pLower = {0.0, -6.0, -6.0, -6.0}
   local pUpper = {1.0,  6.0,  6.0,  6.0}
   local pN     = {1, 8, 8, 8}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1)},{phaseGrid:upper(1)},{phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local distf         = createField(phaseGrid,phaseBasis)
   local energyDensity = createField(confGrid,confBasis)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *(phaseGrid:upper(2)-phaseGrid:lower(2))
                   *(phaseGrid:upper(3)-phaseGrid:lower(3))
                   *((1/3)*(phaseGrid:upper(4)^3-phaseGrid:lower(4)^3)) )
      end
   )
   project:advance(0.0, {}, {distf})

   -- Compute integral of vSq*f with IntegratedDGMoment.
   local intFn   = DataStruct.DynVector {numComponents = 1,}
   local intMomP = intDGmom(phaseGrid,phaseBasis,"vSq",confBasis)
   intMomP:advance(0.0, {distf}, {intFn})
   local _, NPn  = intFn:lastData()

   -- Compute integral of v1Sq*f with IntegratedDGMoment.
   local intv1SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom1  = intDGmom(phaseGrid,phaseBasis,"v1Sq",confBasis)
   intMom1:advance(0.0, {distf}, {intv1SqF})
   local _, NP1   = intv1SqF:lastData()

   -- Compute integral of v2Sq*f with IntegratedDGMoment.
   local intv2SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom2  = intDGmom(phaseGrid,phaseBasis,"v2Sq",confBasis)
   intMom2:advance(0.0, {distf}, {intv2SqF})
   local _, NP2   = intv2SqF:lastData()

   -- Compute integral of v3Sq*f with IntegratedDGMoment.
   local intv3SqF = DataStruct.DynVector {numComponents = 1,}
   local intMom3  = intDGmom(phaseGrid,phaseBasis,"v3Sq",confBasis)
   intMom3:advance(0.0, {distf}, {intv3SqF})
   local _, NP3   = intv3SqF:lastData()

   assert_equal(1, NP3[1], "Checking v3Sq moment")
   assert_equal(NPn[1], NP1[1]+NP2[1]+NP3[1], "Checking v1Sq+v2Sq moment")
end

local function test_xi(pOrder, basis)
   local pLower = {0.0, -2.0, -2.0}
   local pUpper = {1.0,  6.0,  6.0}
   local pN     = {4, 4, 4}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1)},{phaseGrid:upper(1)},{phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local vDim       = phaseBasis:ndim()-confBasis:ndim()
   local distf      = createField(phaseGrid,phaseBasis)
   local momDensity = createField(confGrid,confBasis,vDim)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *(phaseGrid:upper(2)-phaseGrid:lower(2))
                   *(phaseGrid:upper(3)-phaseGrid:lower(3))
                   *(0.5*(phaseGrid:upper(4)^2-phaseGrid:lower(4)^2)) )
      end
   )
   project:advance(0.0, {}, {distf})
   -- Compute integral of x_1*f, x_2*f and x_3*f at the same time with IntegratedDGMoment.
   local intxiF  = DataStruct.DynVector {numComponents = phaseGrid:ndim(),}
   local intMom  = intDGmom(phaseGrid,phaseBasis,"xi")
   intMom:advance(0.0, {distf}, {intxiF})
   local _, intM = intxiF:lastData()

   -- Compute integral of x_1*f with IntegratedDGMoment.
   local intx1F   = DataStruct.DynVector {numComponents = 1,}
   local intMom1  = intDGmom(phaseGrid,phaseBasis,"x1")
   intMom1:advance(0.0, {distf}, {intx1F})
   local _, intM1 = intx1F:lastData()
   -- Compute integral of x_2*f with IntegratedDGMoment.
   local intx2F   = DataStruct.DynVector {numComponents = 1,}
   local intMom2  = intDGmom(phaseGrid,phaseBasis,"x2")
   intMom2:advance(0.0, {distf}, {intx2F})
   local _, intM2 = intx2F:lastData()
   -- Compute integral of x_3*f with IntegratedDGMoment.
   local intx3F   = DataStruct.DynVector {numComponents = 1,}
   local intMom3  = intDGmom(phaseGrid,phaseBasis,"x3")
   intMom3:advance(0.0, {distf}, {intx3F})
   local _, intM3 = intx3F:lastData()

   assert_equal(intM[1], intM1[1], "Checking x1 moment")
   assert_equal(intM[2], intM2[1], "Checking x2 moment")
   assert_equal(intM[3], intM3[1], "Checking x2 moment")
end

local function test_vi(pOrder, basis)
   local pLower = {0.0, -2.0, -2.0, -2.0}
   local pUpper = {1.0,  6.0, 6.0, 6.0}
   local pN     = {1, 8, 8, 8}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1)},{phaseGrid:upper(1)},{phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local vDim       = phaseBasis:ndim()-confBasis:ndim()
   local distf      = createField(phaseGrid,phaseBasis)
   local momDensity = createField(confGrid,confBasis,vDim)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *(phaseGrid:upper(2)-phaseGrid:lower(2))
                   *(phaseGrid:upper(3)-phaseGrid:lower(3))
                   *(0.5*(phaseGrid:upper(4)^2-phaseGrid:lower(4)^2)) )
      end
   )
   project:advance(0.0, {}, {distf})
   -- Compute integral of v_1*f, v_2*f and v_3*f at the same time with IntegratedDGMoment.
   local intviF  = DataStruct.DynVector {numComponents = vDim,}
   local intMom  = intDGmom(phaseGrid,phaseBasis,"vi",confBasis)
   intMom:advance(0.0, {distf}, {intviF})
   local _, intM = intviF:lastData()

   -- Compute integral of v_1*f with IntegratedDGMoment.
   local intv1F   = DataStruct.DynVector {numComponents = 1,}
   local intMom1  = intDGmom(phaseGrid,phaseBasis,"v1",confBasis)
   intMom1:advance(0.0, {distf}, {intv1F})
   local _, intM1 = intv1F:lastData()
   -- Compute integral of v_2*f with IntegratedDGMoment.
   local intv2F   = DataStruct.DynVector {numComponents = 1,}
   local intMom2  = intDGmom(phaseGrid,phaseBasis,"v2",confBasis)
   intMom2:advance(0.0, {distf}, {intv2F})
   local _, intM2 = intv2F:lastData()
   -- Compute integral of v_3*f with IntegratedDGMoment.
   local intv3F   = DataStruct.DynVector {numComponents = 1,}
   local intMom3  = intDGmom(phaseGrid,phaseBasis,"v3",confBasis)
   intMom3:advance(0.0, {distf}, {intv3F})
   local _, intM3 = intv3F:lastData()

   assert_equal(intM[1], intM1[1], "Checking v1 moment")
   assert_equal(intM[2], intM2[1], "Checking v2 moment")
   assert_equal(intM[3], intM3[1], "Checking v2 moment")
end

local function test_intM(pOrder, basis)
   local pLower = {0.0, -2.0, -2.0, -2.0}
   local pUpper = {1.0,  6.0, 6.0, 6.0}
   local pN     = {1, 8, 8, 8}
   -- Phase-space and config-space grids.
   local phaseGrid = createGrid(pLower, pUpper, pN)
   local confGrid  = createGrid({phaseGrid:lower(1)},{phaseGrid:upper(1)},{phaseGrid:numCells(1)})
   -- Basis functions.
   local phaseBasis = createBasis(phaseGrid:ndim(), pOrder, basis)
   local confBasis  = createBasis(confGrid:ndim(), pOrder, basis)
   -- Fields.
   local vDim       = phaseBasis:ndim()-confBasis:ndim()
   local distf      = createField(phaseGrid,phaseBasis)
   local momDensity = createField(confGrid,confBasis,vDim)

   -- Updater to initialize distribution function
   local project = createProject(phaseGrid,phaseBasis)
   project:setFunc(function (t, xn)
         return 1/( (phaseGrid:upper(1)-phaseGrid:lower(1))
                   *(phaseGrid:upper(2)-phaseGrid:lower(2))
                   *(phaseGrid:upper(3)-phaseGrid:lower(3))
                   *(0.5*(phaseGrid:upper(4)^2-phaseGrid:lower(4)^2)) )
      end
   )
   project:advance(0.0, {}, {distf})
   -- Compute integral of v_1*f, v_2*f and v_3*f at the same time with IntegratedDGMoment.
   local intMF   = DataStruct.DynVector {numComponents = 2+vDim,}
   local intMom  = intDGmom(phaseGrid,phaseBasis,"intM",confBasis)
   intMom:advance(0.0, {distf}, {intMF})
   local _, intM = intMF:lastData()

   -- Compute integral of f with IntegratedDGMoment.
   local intF       = DataStruct.DynVector {numComponents = 1,}
   local intMom0    = intDGmom(phaseGrid,phaseBasis,"one",confBasis)
   intMom0:advance(0.0, {distf}, {intF})
   local _, intOneF = intF:lastData()
   -- Compute integral of v_1*f with IntegratedDGMoment.
   local intv1F   = DataStruct.DynVector {numComponents = 1,}
   local intMom1  = intDGmom(phaseGrid,phaseBasis,"v1",confBasis)
   intMom1:advance(0.0, {distf}, {intv1F})
   local _, intM1 = intv1F:lastData()
   -- Compute integral of v_2*f with IntegratedDGMoment.
   local intv2F   = DataStruct.DynVector {numComponents = 1,}
   local intMom2  = intDGmom(phaseGrid,phaseBasis,"v2",confBasis)
   intMom2:advance(0.0, {distf}, {intv2F})
   local _, intM2 = intv2F:lastData()
   -- Compute integral of v_3*f with IntegratedDGMoment.
   local intv3F   = DataStruct.DynVector {numComponents = 1,}
   local intMom3  = intDGmom(phaseGrid,phaseBasis,"v3",confBasis)
   intMom3:advance(0.0, {distf}, {intv3F})
   local _, intM3 = intv3F:lastData()
   -- Compute integral of v^2*f with IntegratedDGMoment.
   local intvSqF  = DataStruct.DynVector {numComponents = 1,}
   local intMomE  = intDGmom(phaseGrid,phaseBasis,"vSq",confBasis)
   intMomE:advance(0.0, {distf}, {intvSqF})
   local _, intME = intvSqF:lastData()

   assert_equal(intM[1], intOneF[1], "Checking 'one' moment")
   assert_equal(intM[2], intM1[1], "Checking v1 moment")
   assert_equal(intM[3], intM2[1], "Checking v2 moment")
   assert_equal(intM[4], intM3[1], "Checking v2 moment")
   assert_equal(intM[5], intME[1], "Checking vSq moment")
end

local polyOrderMax = 2
local basisType    = "Ser"

for polyOrder = 1, polyOrderMax do
  test_one( polyOrder, basisType)
  test_x1( polyOrder, basisType)
  test_x2v1( polyOrder, basisType)
  test_x3v2( polyOrder, basisType)
  test_x4v3( polyOrder, basisType)
  test_x5v3( polyOrder, basisType)
  if polyOrder == 1 then test_x6v3( polyOrder, basisType) end
  test_x1SqxSq( polyOrder, basisType)
  test_x2SqxSq( polyOrder, basisType)
  test_x3SqxSq( polyOrder, basisType)
  test_x4SqxSq( polyOrder, basisType)
  test_x5SqxSq( polyOrder, basisType)
  if polyOrder == 1 then test_x6SqxSq( polyOrder, basisType) end
  test_v1SqvSq( polyOrder, basisType)
  test_v2SqvSq( polyOrder, basisType)
  test_v3SqvSq( polyOrder, basisType)
  test_xi( polyOrder, basisType)
  test_vi( polyOrder, basisType)
  test_intM( polyOrder, basisType)
end

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
