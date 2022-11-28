-- Gkyl ------------------------------------------------------------------------
--
-- Test for updater to project on a Maxwellian onto DG basis.
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
      lower = lo,  cells        = nCells,
      upper = up,  periodicDirs = pDirs,
   }
   return gridOut
end

local function createConfBasis(dim, pOrder, bKind)
   bKind = bKind or "Ser"
   local basis
   if (bKind=="Ser") then
      basis = Basis.CartModalSerendipity { ndim = dim, polyOrder = pOrder }
   elseif (bKind=="Tensor") then
      basis = Basis.CartModalTensor { ndim = dim, polyOrder = pOrder }
   else
      assert(false,"Invalid basis")
   end
   return basis
end

local function createPhaseBasis(cdim, vdim, pOrder, bKind)
   bKind = bKind or "Ser"
   local pdim = cdim+vdim
   local basis
   if (bKind=="Ser") then
      if pOrder == 1 then
         basis = Basis.CartModalHybrid { cdim = cdim, vdim = vdim }
      else
         basis = Basis.CartModalSerendipity { ndim = pdim, polyOrder = pOrder }
      end
   elseif (bKind=="Tensor") then
      basis = Basis.CartModalTensor { ndim = pdim, polyOrder = pOrder }
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
      metaData = {polyOrder  = basis:polyOrder(),
                  basisType  = basis:id(),},
   }
   return fld
end

function test_1x1v()
   local m0     = 1.0
   local uDrift = 0.75
   local vt     = 1.0
   local vdim   = 1

   local lower     = {-0.50, -6.0*vt}
   local upper     = { 0.50,  6.0*vt}
   local numCells  = {8, 8}

   for polyOrder = 1, 3 do
      local lowerC, upperC, numCellsC = {}, {}, {}
      for d=1,1 do
         lowerC[d] = lower[d]
         upperC[d] = upper[d]
         numCellsC[d] = numCells[d]
      end
      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createPhaseBasis(#lowerC, #lower-#lowerC, polyOrder)
      local confGrid   = createGrid(lowerC, upperC, numCellsC)
      local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)
   
      local momsFld     = createField(confGrid, confBasis, vdim+2)
      local primMomsFld = createField(confGrid, confBasis, vdim+1)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)
   
      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local momsFunc   = function (t, xn)
         local m0, udrift, vtsq = m0Func(t, xn), uDriftFunc(t,xn), vtSqFunc(t,xn)
         local m1 = m0*udrift 
         local m2 = udrift*m1 + m0*vtsq
         return m0, m1, m2
      end
      local primMomsFunc = function (t, xn) return uDriftFunc(t,xn), vtSqFunc(t, xn) end
   
      local projectMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }
      local projectPrimMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return primMomsFunc(t,xn) end
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis = confBasis,
         phaseBasis = phaseBasis,
      }
      local maxwellianPrim = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis   = confBasis,
         phaseBasis = phaseBasis,  usePrimMoms = true,
      }
   
      -- Project moments and primitive moments onto configuration space basis.
      projectMoms:advance(0.0, {}, {momsFld})
      projectPrimMoms:advance(0.0, {}, {primMomsFld})

      -- Do projection.
      maxwellian:advance(0.0, {momsFld}, {fM})
      maxwellianPrim:advance(0.0, {momsFld,primMomsFld}, {distf})
   
      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking mom and prim_mom MaxwellianOnBasis 1x1v.")
         end
      end

   end
end

function testGK_1x1v()
   local m0     = 1.0
   local uDrift = 0.50
   local vt     = 1.0
   local B0     = 0.5
   local mass   = 1.0

   local lower     = {-0.50, -6.0*vt}
   local upper     = { 0.50,  6.0*vt}
   local numCells  = {16, 16}

   for polyOrder = 1, 2 do
      local lowerC, upperC, numCellsC = {}, {}, {}
      for d=1,1 do
         lowerC[d] = lower[d]
         upperC[d] = upper[d]
         numCellsC[d] = numCells[d]
      end
      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createPhaseBasis(#lowerC, #lower-#lowerC, polyOrder)
      local confGrid   = createGrid(lowerC, upperC, numCellsC)
      local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)

      local momsFld     = createField(confGrid, confBasis, 3)
      local primMomsFld = createField(confGrid, confBasis, 2)
      local jacobTotFld = createField(confGrid, confBasis)
      local bmagFld     = createField(confGrid, confBasis)
      local distf       = createField(phaseGrid, phaseBasis)
      local fM          = createField(phaseGrid, phaseBasis)

      local m0Func       = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc   = function (t, xn) return uDrift end
      local vtSqFunc     = function (t, xn) return vt^2 end
      local bmagFunc     = function (t, xn) return B0 end
      local jacobTotFunc = function (t, xn) return 1.0 end
      local momsFunc     = function (t, xn)
         local m0, udrift, vtsq = m0Func(t, xn), uDriftFunc(t,xn), vtSqFunc(t,xn)
         local m1 = m0*udrift
         local m2 = udrift*m1 + m0*vtsq
         return m0, m1, m2
      end
      local primMomsFunc = function (t, xn) return uDriftFunc(t,xn), vtSqFunc(t, xn) end

      local projectMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }
      local projectPrimMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return primMomsFunc(t,xn) end
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis = confBasis,
         phaseBasis = phaseBasis,
         mass       = mass,
      }
      local maxwellianPrim = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis   = confBasis,
         phaseBasis = phaseBasis,  usePrimMoms = true,
         mass       = mass,
      }

      -- Project moments and primitive moments onto configuration space basis.
      projectMoms:advance(0.0, {}, {momsFld})
      projectPrimMoms:advance(0.0, {}, {primMomsFld})

      -- Project magnetic field and jacobian onto configuration space basis.
      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      projectScalar:setFunc(function(t,xn) return bmagFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {bmagFld})
      projectScalar:setFunc(function(t,xn) return jacobTotFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {jacobTotFld})

      -- Do projection.
      maxwellian:advance(0.0, {momsFld,bmagFld,jacobTotFld}, {fM})
      maxwellianPrim:advance(0.0, {momsFld,primMomsFld,bmagFld,jacobTotFld}, {distf})

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking mom and prim_mom GkMaxwellianOnBasis 1x1v.")
         end
      end

--      -- Now create a gyrokinetic Maxwellian taking a 3v flow velocity.
--      local uDriftFld  = createField(confGrid, confBasis, 3)
--      local uDriftFunc = function (t, xn) return 1.5, 2.5, uDrift end
--      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
--      projectVec:advance(0.0, {}, {uDriftFld})
--      local uDriftMaxwellian = Updater.MaxwellianOnBasis {
--         onGrid         = phaseGrid,
--         phaseBasis     = phaseBasis,
--         confGrid       = confGrid,
--         confBasis      = confBasis,
--         implementation = "C",
--         mass           = mass,
--      }
--      -- Do projection.
--      local tmStart = Time.clock()
--      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
--      local tmMid = Time.clock()
--      uDriftMaxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})
--      local tmEnd = Time.clock()
--
--      -- Check projection.
--      local indexer    = distf:genIndexer()
--      local distFPtr   = distf:get(1)
--      local fMPtr      = fM:get(1)
--      local localRange = distf:localRange()
--      for idx in localRange:rowMajorIter() do
--         distfPtr = distf:get(indexer(idx))
--         fMPtr    = fM:get(indexer(idx))
--         for k = 1, phaseBasis:numBasis() do
--            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C GkMaxwellianOnBasis (with 3v uFlow) 1x1v.")
--         end
--      end

   end

end

function test_1x2v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75}
   local vt     = 1.0
   local vdim   = 2

   local lower     = {-0.50, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  6.0*vt,  6.0*vt}
   local numCells  = {8, 8, 8}

   for polyOrder = 1, 3 do
      local lowerC, upperC, numCellsC = {}, {}, {}
      for d=1,1 do
         lowerC[d] = lower[d]
         upperC[d] = upper[d]
         numCellsC[d] = numCells[d]
      end
      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createPhaseBasis(#lowerC, #lower-#lowerC, polyOrder)
      local confGrid   = createGrid(lowerC, upperC, numCellsC)
      local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)

      local momsFld     = createField(confGrid, confBasis, vdim+2)
      local primMomsFld = createField(confGrid, confBasis, vdim+1)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift[1], uDrift[2] end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local momsFunc   = function (t, xn)
         local m1, udrift = {}, {}
         for d = 1, vdim do m1[d], udrift[d] = 0., 0. end
         local m0 = m0Func(t, xn)
         udrift[1], udrift[2] = uDriftFunc(t,xn)
         local vtsq = vtSqFunc(t,xn)
         local m2 = 0.
         for d = 1, vdim do
            m1[d] = m0*udrift[d]
            m2 = m2 + udrift[d]*m1[d]
         end 
         local m2 = m2 + vdim*m0*vtsq
         return m0, m1[1], m1[2], m2
      end
      local primMomsFunc = function (t, xn)
         local udrift = {}
         for d = 1, vdim do udrift[d] = 0. end
         udrift[1], udrift[2] = uDriftFunc(t,xn)
         local vtsq = vtSqFunc(t, xn)
         return udrift[1], udrift[2], vtsq
      end

      local projectMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }
      local projectPrimMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return primMomsFunc(t,xn) end
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis = confBasis,
         phaseBasis = phaseBasis,
      }
      local maxwellianPrim = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis   = confBasis,
         phaseBasis = phaseBasis,  usePrimMoms = true,
      }

      -- Project moments and primitive moments onto configuration space basis.
      projectMoms:advance(0.0, {}, {momsFld})
      projectPrimMoms:advance(0.0, {}, {primMomsFld})

      -- Do projection.
      maxwellian:advance(0.0, {momsFld}, {fM})
      maxwellianPrim:advance(0.0, {momsFld,primMomsFld}, {distf})

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking mom and prim_mom MaxwellianOnBasis 1x2v.")
         end
      end

   end

end

function testGK_1x2v()
   local m0     = 1.0
   local uDrift = 0.50
   local vt     = 1.0
   local B0     = 0.5
   local mass   = 1.0

   local lower     = {-0.50, -6.0*vt, 0.0}
   local upper     = { 0.50,  6.0*vt, mass*(6.0*(vt)^2)/(2*B0)}
   local numCells  = {16, 16, 16}

   for polyOrder = 1, 2 do
      local lowerC, upperC, numCellsC = {}, {}, {}
      for d=1,1 do
         lowerC[d] = lower[d]
         upperC[d] = upper[d]
         numCellsC[d] = numCells[d]
      end
      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createPhaseBasis(#lowerC, #lower-#lowerC, polyOrder)
      local confGrid   = createGrid(lowerC, upperC, numCellsC)
      local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)

      local momsFld     = createField(confGrid, confBasis, 3)
      local primMomsFld = createField(confGrid, confBasis, 2)
      local jacobTotFld = createField(confGrid, confBasis)
      local bmagFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local bmagFunc   = function (t, xn) return B0 end
      local jacobTotFunc = function (t, xn) return 1.0 end
      local momsFunc     = function (t, xn)
         local m0, udrift, vtsq = m0Func(t, xn), uDriftFunc(t,xn), vtSqFunc(t,xn)
         local m1 = m0*udrift
         local m2 = udrift*m1 + 3.*m0*vtsq
         return m0, m1, m2
      end
      local primMomsFunc = function (t, xn) return uDriftFunc(t,xn), vtSqFunc(t, xn) end

      local projectMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }
      local projectPrimMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return primMomsFunc(t,xn) end
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis = confBasis,
         phaseBasis = phaseBasis,
         mass       = mass,
      }
      local maxwellianPrim = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis   = confBasis,
         phaseBasis = phaseBasis,  usePrimMoms = true,
         mass       = mass,
      }

      -- Project moments and primitive moments onto configuration space basis.
      projectMoms:advance(0.0, {}, {momsFld})
      projectPrimMoms:advance(0.0, {}, {primMomsFld})

      -- Project magnetic field and jacobian onto configuration space basis.
      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      projectScalar:setFunc(function(t,xn) return bmagFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {bmagFld})
      projectScalar:setFunc(function(t,xn) return jacobTotFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {jacobTotFld})

      -- Do projection.
      maxwellian:advance(0.0, {momsFld,bmagFld,jacobTotFld}, {fM})
      maxwellianPrim:advance(0.0, {momsFld,primMomsFld,bmagFld,jacobTotFld}, {distf})

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking mom and prim_mom GkMaxwellianOnBasis 1x2v.")
         end
      end

--      -- Now create a gyrokinetic Maxwellian taking a 3v flow velocity.
--      local uDriftFld  = createField(confGrid, confBasis, 3)
--      local uDriftFunc = function (t, xn) return 1.5, 2.5, uDrift end
--      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
--      projectVec:advance(0.0, {}, {uDriftFld})
--      local uDriftMaxwellian = Updater.MaxwellianOnBasis {
--         onGrid         = phaseGrid,
--         phaseBasis     = phaseBasis,
--         confGrid       = confGrid,
--         confBasis      = confBasis,
--         implementation = "C",
--         mass           = mass,
--      }
--      -- Do projection.
--      local tmStart = Time.clock()
--      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
--      local tmMid = Time.clock()
--      uDriftMaxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})
--      local tmEnd = Time.clock()
--
--      -- Check projection.
--      local indexer    = distf:genIndexer()
--      local distFPtr   = distf:get(1)
--      local fMPtr      = fM:get(1)
--      local localRange = distf:localRange()
--      for idx in localRange:rowMajorIter() do
--         distfPtr = distf:get(indexer(idx))
--         fMPtr    = fM:get(indexer(idx))
--         for k = 1, phaseBasis:numBasis() do
--            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C GkMaxwellianOnBasis (with 3v uFlow) 1x2v.")
--         end
--      end

   end

end

function test_1x3v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75, 0.1}
   local vt     = 1.0
   local vdim   = 3

   local lower     = {-0.50, -6.0*vt, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  6.0*vt,  6.0*vt,  6.0*vt}
   local numCells  = {6, 8, 8, 4}

   for polyOrder = 1, 3 do
      local lowerC, upperC, numCellsC = {}, {}, {}
      for d=1,1 do
         lowerC[d] = lower[d]
         upperC[d] = upper[d]
         numCellsC[d] = numCells[d]
      end
      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createPhaseBasis(#lowerC, #lower-#lowerC, polyOrder)
      local confGrid   = createGrid(lowerC, upperC, numCellsC)
      local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)

      local momsFld     = createField(confGrid, confBasis, vdim+2)
      local primMomsFld = createField(confGrid, confBasis, vdim+1)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift[1], uDrift[2], uDrift[3] end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local momsFunc   = function (t, xn)
         local m1, udrift = {}, {}
         for d = 1, vdim do m1[d], udrift[d] = 0., 0. end
         local m0 = m0Func(t, xn)
         udrift[1], udrift[2], udrift[3] = uDriftFunc(t,xn)
         local vtsq = vtSqFunc(t,xn)
         local m2 = 0.
         for d = 1, vdim do
            m1[d] = m0*udrift[d]
            m2 = m2 + udrift[d]*m1[d]
         end
         local m2 = m2 + vdim*m0*vtsq
         return m0, m1[1], m1[2], m1[3], m2
      end
      local primMomsFunc = function (t, xn)
         local udrift = {}
         for d = 1, vdim do udrift[d] = 0. end
         udrift[1], udrift[2], udrift[3] = uDriftFunc(t,xn)
         local vtsq = vtSqFunc(t, xn)
         return udrift[1], udrift[2], udrift[3], vtsq
      end

      local projectMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }
      local projectPrimMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return primMomsFunc(t,xn) end
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis = confBasis,
         phaseBasis = phaseBasis,
      }
      local maxwellianPrim = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis   = confBasis,
         phaseBasis = phaseBasis,  usePrimMoms = true,
      }

      -- Project moments and primitive moments onto configuration space basis.
      projectMoms:advance(0.0, {}, {momsFld})
      projectPrimMoms:advance(0.0, {}, {primMomsFld})

      -- Do projection.
      local tmStart = Time.clock()
      maxwellian:advance(0.0, {momsFld}, {fM})
      local tmMid = Time.clock()
      maxwellianPrim:advance(0.0, {momsFld,primMomsFld}, {distf})
      local tmEnd = Time.clock()

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking prim and prim_mom MaxwellianOnBasis 1x3v.")
         end
      end

--      -- Now create a Maxwellian with uPar and assuming zero flow velocity in vx and vy.
--      local uParFld  = createField(confGrid, confBasis, 1)
--      local uParFunc = function (t, xn) return uDrift[3] end
--      projectScalar:setFunc(function(t,xn) return uParFunc(t,xn) end)
--      projectScalar:advance(0.0, {}, {uParFld})
--      local uParMaxwellian = Updater.MaxwellianOnBasis {
--         onGrid         = phaseGrid,
--         phaseBasis     = phaseBasis,
--         confGrid       = confGrid,
--         confBasis      = confBasis,
--         implementation = "C",
--      }
--      -- Do projection.
--      local tmStart = Time.clock()
--      maxwellianLua:advance(0.0, {m0Fld,uParFld,vtSqFld}, {distf})
--      local tmMid = Time.clock()
--      uParMaxwellian:advance(0.0, {m0Fld,uParFld,vtSqFld}, {fM})
--      local tmEnd = Time.clock()
--
--      -- Check projection.
--      local indexer    = distf:genIndexer()
--      local distFPtr   = distf:get(1)
--      local fMPtr      = fM:get(1)
--      local localRange = distf:localRange()
--      for idx in localRange:rowMajorIter() do
--         distfPtr = distf:get(indexer(idx))
--         fMPtr    = fM:get(indexer(idx))
--         for k = 1, phaseBasis:numBasis() do
--            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C MaxwellianOnBasis (with uPar) 1x3v.")
--         end
--      end

   end
end

function test_2x2v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75}
   local vt     = 1.0
   local vdim   = 2

   local lower     = {-0.50, -0.50, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  0.50,  6.0*vt,  6.0*vt}
   local numCells  = {6, 12, 8, 8}

   for polyOrder = 1, 3 do
      local lowerC, upperC, numCellsC = {}, {}, {}
      for d=1,2 do
         lowerC[d] = lower[d]
         upperC[d] = upper[d]
         numCellsC[d] = numCells[d]
      end
      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createPhaseBasis(#lowerC, #lower-#lowerC, polyOrder)
      local confGrid   = createGrid(lowerC, upperC, numCellsC)
      local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)

      local momsFld     = createField(confGrid, confBasis, vdim+2)
      local primMomsFld = createField(confGrid, confBasis, vdim+1)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift[1], uDrift[2] end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local momsFunc   = function (t, xn)
         local m1, udrift = {}, {}
         for d = 1, vdim do m1[d], udrift[d] = 0., 0. end
         local m0 = m0Func(t, xn)
         udrift[1], udrift[2] = uDriftFunc(t,xn)
         local vtsq = vtSqFunc(t,xn)
         local m2 = 0.
         for d = 1, vdim do
            m1[d] = m0*udrift[d]
            m2 = m2 + udrift[d]*m1[d]
         end
         local m2 = m2 + vdim*m0*vtsq
         return m0, m1[1], m1[2], m2
      end
      local primMomsFunc = function (t, xn)
         local udrift = {}
         for d = 1, vdim do udrift[d] = 0. end
         udrift[1], udrift[2] = uDriftFunc(t,xn)
         local vtsq = vtSqFunc(t, xn)
         return udrift[1], udrift[2], vtsq
      end

      local projectMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }
      local projectPrimMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return primMomsFunc(t,xn) end
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis = confBasis,
         phaseBasis = phaseBasis,
      }
      local maxwellianPrim = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis   = confBasis,
         phaseBasis = phaseBasis,  usePrimMoms = true,
      }

      -- Project moments and primitive moments onto configuration space basis.
      projectMoms:advance(0.0, {}, {momsFld})
      projectPrimMoms:advance(0.0, {}, {primMomsFld})

      -- Do projection.
      local tmStart = Time.clock()
      maxwellian:advance(0.0, {momsFld}, {fM})
      local tmMid = Time.clock()
      maxwellianPrim:advance(0.0, {momsFld,primMomsFld}, {distf})
      local tmEnd = Time.clock()

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C MaxwellianOnBasis 2x2v.")
         end
      end

   end

end

function testGK_2x2v()
   local m0     = 1.0
   local uDrift = 0.50
   local vt     = 1.0
   local B0     = 0.5
   local mass   = 1.0

   local lower     = {-0.50, -0.50, -6.0*vt, 0.0}
   local upper     = { 0.50,  0.50,  6.0*vt, mass*(6.0*(vt)^2)/(2*B0)}
   local numCells  = {8, 8, 8, 8}

   for polyOrder = 1, 2 do
      local lowerC, upperC, numCellsC = {}, {}, {}
      for d=1,2 do
         lowerC[d] = lower[d]
         upperC[d] = upper[d]
         numCellsC[d] = numCells[d]
      end
      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createPhaseBasis(#lowerC, #lower-#lowerC, polyOrder)
      local confGrid   = createGrid(lowerC, upperC, numCellsC)
      local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)

      local momsFld     = createField(confGrid, confBasis, 3)
      local primMomsFld = createField(confGrid, confBasis, 2)
      local jacobTotFld = createField(confGrid, confBasis)
      local bmagFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local bmagFunc   = function (t, xn) return B0 end
      local jacobTotFunc = function (t, xn) return 1.0 end
      local momsFunc     = function (t, xn)
         local m0, udrift, vtsq = m0Func(t, xn), uDriftFunc(t,xn), vtSqFunc(t,xn)
         local m1 = m0*udrift
         local m2 = udrift*m1 + 3.*m0*vtsq
         return m0, m1, m2
      end
      local primMomsFunc = function (t, xn) return uDriftFunc(t,xn), vtSqFunc(t, xn) end

      local projectMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }
      local projectPrimMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return primMomsFunc(t,xn) end
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis = confBasis,
         phaseBasis = phaseBasis,
         mass       = mass,
      }
      local maxwellianPrim = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis   = confBasis,
         phaseBasis = phaseBasis,  usePrimMoms = true,
         mass       = mass,
      }

      -- Project moments and primitive moments onto configuration space basis.
      projectMoms:advance(0.0, {}, {momsFld})
      projectPrimMoms:advance(0.0, {}, {primMomsFld})

      -- Project magnetic field and jacobian onto configuration space basis.
      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      projectScalar:setFunc(function(t,xn) return bmagFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {bmagFld})
      projectScalar:setFunc(function(t,xn) return jacobTotFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {jacobTotFld})

      -- Do projection.
      maxwellian:advance(0.0, {momsFld,bmagFld,jacobTotFld}, {fM})
      maxwellianPrim:advance(0.0, {momsFld,primMomsFld,bmagFld,jacobTotFld}, {distf})


      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking mom and prim_mom GkMaxwellianOnBasis 2x2v.")
         end
      end

--      -- Now create a gyrokinetic Maxwellian taking a 3v flow velocity.
--      local uDriftFld  = createField(confGrid, confBasis, 3)
--      local uDriftFunc = function (t, xn) return 1.5, 2.5, uDrift end
--      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
--      projectVec:advance(0.0, {}, {uDriftFld})
--      local uDriftMaxwellian = Updater.MaxwellianOnBasis {
--         onGrid         = phaseGrid,
--         phaseBasis     = phaseBasis,
--         confGrid       = confGrid,
--         confBasis      = confBasis,
--         implementation = "C",
--         mass           = mass,
--      }
--      -- Do projection.
--      local tmStart = Time.clock()
--      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
--      local tmMid = Time.clock()
--      uDriftMaxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})
--      local tmEnd = Time.clock()
--
--      -- Check projection.
--      local indexer    = distf:genIndexer()
--      local distFPtr   = distf:get(1)
--      local fMPtr      = fM:get(1)
--      local localRange = distf:localRange()
--      for idx in localRange:rowMajorIter() do
--         distfPtr = distf:get(indexer(idx))
--         fMPtr    = fM:get(indexer(idx))
--         for k = 1, phaseBasis:numBasis() do
--            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C GkMaxwellianOnBasis (with 3v uFlow) 2x2v.")
--         end
--      end

   end

end

function test_2x3v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75, 0.1}
   local vt     = 1.0
   local vdim   = 3

   local lower     = {-0.50, -0.50, -6.0*vt, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  0.50,  6.0*vt,  6.0*vt,  6.0*vt}
   local numCells  = {6, 12, 8, 8, 4}

   for polyOrder = 1, 2 do
      local lowerC, upperC, numCellsC = {}, {}, {}
      for d=1,2 do
         lowerC[d] = lower[d]
         upperC[d] = upper[d]
         numCellsC[d] = numCells[d]
      end
      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createPhaseBasis(#lowerC, #lower-#lowerC, polyOrder)
      local confGrid   = createGrid(lowerC, upperC, numCellsC)
      local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)

      local momsFld     = createField(confGrid, confBasis, vdim+2)
      local primMomsFld = createField(confGrid, confBasis, vdim+1)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift[1], uDrift[2], uDrift[3] end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local momsFunc   = function (t, xn)
         local m1, udrift = {}, {}
         for d = 1, vdim do m1[d], udrift[d] = 0., 0. end
         local m0 = m0Func(t, xn)
         udrift[1], udrift[2], udrift[3] = uDriftFunc(t,xn)
         local vtsq = vtSqFunc(t,xn)
         local m2 = 0.
         for d = 1, vdim do
            m1[d] = m0*udrift[d]
            m2 = m2 + udrift[d]*m1[d]
         end
         local m2 = m2 + vdim*m0*vtsq
         return m0, m1[1], m1[2], m1[3], m2
      end
      local primMomsFunc = function (t, xn)
         local udrift = {}
         for d = 1, vdim do udrift[d] = 0. end
         udrift[1], udrift[2], udrift[3] = uDriftFunc(t,xn)
         local vtsq = vtSqFunc(t, xn)
         return udrift[1], udrift[2], udrift[3], vtsq
      end

      local projectMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }
      local projectPrimMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return primMomsFunc(t,xn) end
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis = confBasis,
         phaseBasis = phaseBasis,
      }
      local maxwellianPrim = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis   = confBasis,
         phaseBasis = phaseBasis,  usePrimMoms = true,
      }

      -- Project moments and primitive moments onto configuration space basis.
      projectMoms:advance(0.0, {}, {momsFld})
      projectPrimMoms:advance(0.0, {}, {primMomsFld})

      -- Do projection.
      local tmStart = Time.clock()
      maxwellian:advance(0.0, {momsFld}, {fM})
      local tmMid = Time.clock()
      maxwellianPrim:advance(0.0, {momsFld,primMomsFld}, {distf})
      local tmEnd = Time.clock()

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking prim and prim_mom MaxwellianOnBasis 2x3v.")
         end
      end

--      -- Now create a Maxwellian with uPar and assuming zero flow velocity in vx and vy.
--      local uParFld  = createField(confGrid, confBasis, 1)
--      local uParFunc = function (t, xn) return uDrift[3] end
--      projectScalar:setFunc(function(t,xn) return uParFunc(t,xn) end)
--      projectScalar:advance(0.0, {}, {uParFld})
--      local uParMaxwellian = Updater.MaxwellianOnBasis {
--         onGrid         = phaseGrid,
--         phaseBasis     = phaseBasis,
--         confGrid       = confGrid,
--         confBasis      = confBasis,
--         implementation = "C",
--      }
--      -- Do projection.
--      local tmStart = Time.clock()
--      maxwellianLua:advance(0.0, {m0Fld,uParFld,vtSqFld}, {distf})
--      local tmMid = Time.clock()
--      uParMaxwellian:advance(0.0, {m0Fld,uParFld,vtSqFld}, {fM})
--      local tmEnd = Time.clock()
--
--      -- Check projection.
--      local indexer    = distf:genIndexer()
--      local distFPtr   = distf:get(1)
--      local fMPtr      = fM:get(1)
--      local localRange = distf:localRange()
--      for idx in localRange:rowMajorIter() do
--         distfPtr = distf:get(indexer(idx))
--         fMPtr    = fM:get(indexer(idx))
--         for k = 1, phaseBasis:numBasis() do
--            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C MaxwellianOnBasis (with uPar) 2x3v.")
--         end
--      end

   end

end

function testGK_3x2v()
   local m0     = 1.0
   local uDrift = 2.50
   local vt     = 1.0
   local B0     = 0.5
   local mass   = 1.0

   local lower     = {-0.50, -0.50, -0.50, -6.0*vt, 0.0}
   local upper     = { 0.50,  0.50,  0.50,  6.0*vt, mass*(6.0*(vt)^2)/(2*B0)}
   local numCells  = {6, 8, 8, 8, 4}

   for polyOrder = 1, 2 do
      local lowerC, upperC, numCellsC = {}, {}, {}
      for d=1,3 do
         lowerC[d] = lower[d]
         upperC[d] = upper[d]
         numCellsC[d] = numCells[d]
      end
      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createPhaseBasis(#lowerC, #lower-#lowerC, polyOrder)
      local confGrid   = createGrid(lowerC, upperC, numCellsC)
      local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)

      local momsFld     = createField(confGrid, confBasis, 3)
      local primMomsFld = createField(confGrid, confBasis, 2)
      local jacobTotFld = createField(confGrid, confBasis)
      local bmagFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local bmagFunc   = function (t, xn) return B0 end
      local jacobTotFunc = function (t, xn) return 1.0 end
      local momsFunc     = function (t, xn)
         local m0, udrift, vtsq = m0Func(t, xn), uDriftFunc(t,xn), vtSqFunc(t,xn)
         local m1 = m0*udrift
         local m2 = udrift*m1 + 3.*m0*vtsq
         return m0, m1, m2
      end
      local primMomsFunc = function (t, xn) return uDriftFunc(t,xn), vtSqFunc(t, xn) end

      local projectMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }
      local projectPrimMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return primMomsFunc(t,xn) end
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis = confBasis,
         phaseBasis = phaseBasis,
         mass       = mass,
      }
      local maxwellianPrim = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis   = confBasis,
         phaseBasis = phaseBasis,  usePrimMoms = true,
         mass       = mass,
      }

      -- Project moments and primitive moments onto configuration space basis.
      projectMoms:advance(0.0, {}, {momsFld})
      projectPrimMoms:advance(0.0, {}, {primMomsFld})

      -- Project magnetic field and jacobian onto configuration space basis.
      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      projectScalar:setFunc(function(t,xn) return bmagFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {bmagFld})
      projectScalar:setFunc(function(t,xn) return jacobTotFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {jacobTotFld})

      -- Do projection.
      local tmStart = Time.clock()
      maxwellian:advance(0.0, {momsFld,bmagFld,jacobTotFld}, {fM})
      local tmMid = Time.clock()
      maxwellianPrim:advance(0.0, {momsFld,primMomsFld,bmagFld,jacobTotFld}, {distf})
      local tmEnd = Time.clock()

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking mom and prim_mom GkMaxwellianOnBasis 3x2v.")
         end
      end

--      -- Now create a gyrokinetic Maxwellian taking a 3v flow velocity.
--      local uDriftFld  = createField(confGrid, confBasis, 3)
--      local uDriftFunc = function (t, xn) return 0.13, 0.23, uDrift end
--      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
--      projectVec:advance(0.0, {}, {uDriftFld})
--      local uDriftMaxwellian = Updater.MaxwellianOnBasis {
--         onGrid         = phaseGrid,
--         phaseBasis     = phaseBasis,
--         confGrid       = confGrid,
--         confBasis      = confBasis,
--         implementation = "C",
--         mass           = mass,
--      }
--      -- Do projection.
--      local tmStart = Time.clock()
--      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
--      local tmMid = Time.clock()
--      uDriftMaxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})
--      local tmEnd = Time.clock()
--
--      -- Check projection.
--      local indexer    = distf:genIndexer()
--      local distFPtr   = distf:get(1)
--      local fMPtr      = fM:get(1)
--      local localRange = distf:localRange()
--      for idx in localRange:rowMajorIter() do
--         distfPtr = distf:get(indexer(idx))
--         fMPtr    = fM:get(indexer(idx))
--         for k = 1, phaseBasis:numBasis() do
--            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C GkMaxwellianOnBasis (with 3v uFlow) 3x2v.")
--         end
--      end


   end

end

function test_3x3v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75, 0.1}
   local vt     = 1.0
   local vdim   = 3

   local lower     = {-0.50, -0.50, -0.50, -6.0*vt, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  0.50,  0.50,  6.0*vt,  6.0*vt,  6.0*vt}
   local numCells  = {4, 4, 4, 4, 4, 4}

   for polyOrder = 1, 1 do
      local lowerC, upperC, numCellsC = {}, {}, {}
      for d=1,3 do
         lowerC[d] = lower[d]
         upperC[d] = upper[d]
         numCellsC[d] = numCells[d]
      end
      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createPhaseBasis(#lowerC, #lower-#lowerC, polyOrder)
      local confGrid   = createGrid(lowerC, upperC, numCellsC)
      local confBasis  = createConfBasis(confGrid:ndim(), polyOrder)

      local momsFld     = createField(confGrid, confBasis, vdim+2)
      local primMomsFld = createField(confGrid, confBasis, vdim+1)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift[1], uDrift[2], uDrift[3] end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local momsFunc   = function (t, xn)
         local m1, udrift = {}, {}
         for d = 1, vdim do m1[d], udrift[d] = 0., 0. end
         local m0 = m0Func(t, xn)
         udrift[1], udrift[2], udrift[3] = uDriftFunc(t,xn)
         local vtsq = vtSqFunc(t,xn)
         local m2 = 0.
         for d = 1, vdim do
            m1[d] = m0*udrift[d]
            m2 = m2 + udrift[d]*m1[d]
         end
         local m2 = m2 + vdim*m0*vtsq
         return m0, m1[1], m1[2], m1[3], m2
      end
      local primMomsFunc = function (t, xn)
         local udrift = {}
         for d = 1, vdim do udrift[d] = 0. end
         udrift[1], udrift[2], udrift[3] = uDriftFunc(t,xn)
         local vtsq = vtSqFunc(t, xn)
         return udrift[1], udrift[2], udrift[3], vtsq
      end

      local projectMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return momsFunc(t,xn) end
      }
      local projectPrimMoms = Updater.ProjectOnBasis {
         onGrid   = confGrid,  basis = confBasis,
         evaluate = function (t, xn) return primMomsFunc(t,xn) end
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis = confBasis,
         phaseBasis = phaseBasis,
      }
      local maxwellianPrim = Updater.MaxwellianOnBasis {
         onGrid     = phaseGrid,   confBasis   = confBasis,
         phaseBasis = phaseBasis,  usePrimMoms = true,
      }

      -- Project moments and primitive moments onto configuration space basis.
      projectMoms:advance(0.0, {}, {momsFld})
      projectPrimMoms:advance(0.0, {}, {primMomsFld})

      -- Do projection.
      local tmStart = Time.clock()
      maxwellian:advance(0.0, {momsFld}, {fM})
      local tmMid = Time.clock()
      maxwellianPrim:advance(0.0, {momsFld,primMomsFld}, {distf})
      local tmEnd = Time.clock()

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking mom and prim_mom MaxwellianOnBasis 3x3v.")
         end
      end

--      -- Now create a Maxwellian with uPar and assuming zero flow velocity in vx and vy.
--      local uParFld  = createField(confGrid, confBasis, 1)
--      local uParFunc = function (t, xn) return uDrift[3] end
--      projectScalar:setFunc(function(t,xn) return uParFunc(t,xn) end)
--      projectScalar:advance(0.0, {}, {uParFld})
--      local uParMaxwellian = Updater.MaxwellianOnBasis {
--         onGrid         = phaseGrid,
--         phaseBasis     = phaseBasis,
--         confGrid       = confGrid,
--         confBasis      = confBasis,
--         implementation = "C",
--      }
--      -- Do projection.
--      local tmStart = Time.clock()
--      maxwellianLua:advance(0.0, {m0Fld,uParFld,vtSqFld}, {distf})
--      local tmMid = Time.clock()
--      uParMaxwellian:advance(0.0, {m0Fld,uParFld,vtSqFld}, {fM})
--      local tmEnd = Time.clock()
--
--      -- Check projection.
--      local indexer    = distf:genIndexer()
--      local distFPtr   = distf:get(1)
--      local fMPtr      = fM:get(1)
--      local localRange = distf:localRange()
--      for idx in localRange:rowMajorIter() do
--         distfPtr = distf:get(indexer(idx))
--         fMPtr    = fM:get(indexer(idx))
--         for k = 1, phaseBasis:numBasis() do
--            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C MaxwellianOnBasis (with uPar) 3x3v.")
--         end
--      end
   end

end

-- Run tests.
test_1x1v()
test_1x2v()
test_1x3v()
test_2x2v()
test_2x3v()
test_3x3v()
-- Gyrokinetic tests.
testGK_1x1v()
testGK_1x2v()
testGK_2x2v()
testGK_3x2v()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
