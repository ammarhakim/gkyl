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

function test_1x1v()
   local m0     = 1.0
   local uDrift = 0.75
   local vt     = 1.0

   local lower     = {-0.50, -6.0*vt}
   local upper     = { 0.50,  6.0*vt}
   local numCells  = {8, 8}

   for polyOrder = 1, 3 do

      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid   = createGrid({lower[1]}, {upper[1]}, {numCells[1]})
      local confBasis  = createBasis(confGrid:ndim(), polyOrder)
   
      local m0Fld     = createField(confGrid, confBasis)
      local uDriftFld = createField(confGrid, confBasis)
      local vtSqFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)
   
      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local maxwellianFunc = function (t, xn)
         local x, v = xn[1], xn[2] 
         local fOut = (m0/math.sqrt(2.*math.pi*vt^2))*math.exp(-((v-uDrift)^2)/(2*(vt^2))) 
         return fOut
      end
   
      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projectVec = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local maxwellianLua = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "Lua"
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C"
      }
   
      -- Project the primitive moments onto configuration space basis.
      projectScalar:setFunc(function(t,xn) return m0Func(t,xn) end)
      projectScalar:advance(0.0, {}, {m0Fld})
      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectVec:advance(0.0, {}, {uDriftFld})
      projectScalar:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {vtSqFld})
   
      -- Do projection.
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {distf})
      maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {fM})
   
      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
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

      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid   = createGrid({lower[1]}, {upper[1]}, {numCells[1]})
      local confBasis  = createBasis(confGrid:ndim(), polyOrder)

      local m0Fld     = createField(confGrid, confBasis)
      local uDriftFld = createField(confGrid, confBasis)
      local vtSqFld   = createField(confGrid, confBasis)
      local bmagFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local bmagFunc   = function (t, xn) return B0 end
      local maxwellianFunc = function (t, xn)
         local x, vpar = xn[1], xn[2]
         local fOut = (m0Func(t,xn)/(math.sqrt(2.*math.pi*vtSqFunc(t,xn))))
                     *math.exp(-((vpar-uDriftFunc(t,xn))^2)/(2.*(vtSqFunc(t,xn)))) 
         return fOut
      end

      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projectVec = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0, 1.0, 1.0 end   -- Set later.
      }
      local maxwellianLua = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "Lua",
         mass           = mass,
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
         mass           = mass,
      }

      -- Project the primitive moments onto configuration space basis.
      projectScalar:setFunc(function(t,xn) return m0Func(t,xn) end)
      projectScalar:advance(0.0, {}, {m0Fld})
      projectScalar:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {uDriftFld})
      projectScalar:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {vtSqFld})
      projectScalar:setFunc(function(t,xn) return bmagFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {bmagFld})

      -- Do projection.
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
      maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

      -- Now create a gyrokinetic Maxwellian taking a 3v flow velocity.
      local uDriftFld  = createField(confGrid, confBasis, 3)
      local uDriftFunc = function (t, xn) return 1.5, 2.5, uDrift end
      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectVec:advance(0.0, {}, {uDriftFld})
      local uDriftMaxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
         mass           = mass,
      }
      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
      local tmMid = Time.clock()
      uDriftMaxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

   end

end

function test_1x2v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75}
   local vt     = 1.0

   local lower     = {-0.50, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  6.0*vt,  6.0*vt}
   local numCells  = {8, 8, 8}

   for polyOrder = 1, 3 do

      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid   = createGrid({lower[1]}, {upper[1]}, {numCells[1]})
      local confBasis  = createBasis(confGrid:ndim(), polyOrder)

      local m0Fld     = createField(confGrid, confBasis)
      local uDriftFld = createField(confGrid, confBasis, #uDrift)
      local vtSqFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift[1], uDrift[2] end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local maxwellianFunc = function (t, xn)
         local x, vx, vy = xn[1], xn[2], xn[3]
         local uDriftX, uDriftY = uDriftFunc(t,xn)
         local fOut = (m0Func(t,xn)/(2.*math.pi*vtSqFunc(t,xn)))*math.exp(-((vx-uDriftX)^2+(vy-uDriftY)^2)/(2.*(vtSqFunc(t,xn)))) 
         return fOut
      end

      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projectVec = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0, 1.0 end   -- Set later.
      }
      local maxwellianLua = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "Lua"
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C"
      }

      -- Project the primitive moments onto configuration space basis.
      projectScalar:setFunc(function(t,xn) return m0Func(t,xn) end)
      projectScalar:advance(0.0, {}, {m0Fld})
      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectVec:advance(0.0, {}, {uDriftFld})
      projectScalar:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {vtSqFld})

      -- Do projection.
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {distf})
      maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {fM})

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
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

      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid   = createGrid({lower[1]}, {upper[1]}, {numCells[1]})
      local confBasis  = createBasis(confGrid:ndim(), polyOrder)

      local m0Fld     = createField(confGrid, confBasis)
      local uDriftFld = createField(confGrid, confBasis)
      local vtSqFld   = createField(confGrid, confBasis)
      local bmagFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local bmagFunc   = function (t, xn) return B0 end
      local maxwellianFunc = function (t, xn)
         local x, vpar, mu = xn[1], xn[2], xn[3]
         local fOut = (m0Func(t,xn)/(math.sqrt(2.*math.pi*vtSqFunc(t,xn))^3))
                     *math.exp(-((vpar-uDriftFunc(t,xn))^2+2*mu*bmagFunc(t,xn)/mass)/(2.*(vtSqFunc(t,xn)))) 
         return fOut
      end

      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projectVec = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0, 1.0, 1.0 end   -- Set later.
      }
      local maxwellianLua = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "Lua",
         mass           = mass,
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
         mass           = mass,
      }

      -- Project the primitive moments onto configuration space basis.
      projectScalar:setFunc(function(t,xn) return m0Func(t,xn) end)
      projectScalar:advance(0.0, {}, {m0Fld})
      projectScalar:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {uDriftFld})
      projectScalar:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {vtSqFld})
      projectScalar:setFunc(function(t,xn) return bmagFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {bmagFld})

      -- Do projection.
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
      maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

      -- Now create a gyrokinetic Maxwellian taking a 3v flow velocity.
      local uDriftFld  = createField(confGrid, confBasis, 3)
      local uDriftFunc = function (t, xn) return 1.5, 2.5, uDrift end
      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectVec:advance(0.0, {}, {uDriftFld})
      local uDriftMaxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
         mass           = mass,
      }
      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
      local tmMid = Time.clock()
      uDriftMaxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

   end

end

function test_1x3v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75, 0.1}
   local vt     = 1.0

   local lower     = {-0.50, -6.0*vt, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  6.0*vt,  6.0*vt,  6.0*vt}
   local numCells  = {6, 8, 8, 4}

   for polyOrder = 1, 3 do

      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid   = createGrid({lower[1]}, {upper[1]}, {numCells[1]})
      local confBasis  = createBasis(confGrid:ndim(), polyOrder)

      local m0Fld     = createField(confGrid, confBasis)
      local uDriftFld = createField(confGrid, confBasis, #uDrift)
      local vtSqFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift[1], uDrift[2], uDrift[3] end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local maxwellianFunc = function (t, xn)
         local x, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]
         local fOut = (m0/(math.sqrt(2.*math.pi*vt^2)^3))
                     *math.exp(-((vx-uDrift[1])^2+(vy-uDrift[2])^2+(vz-uDrift[3])^2)/(2.*(vt^2))) 
         return fOut
      end

      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projectVec = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0, 1.0, 1.0 end   -- Set later.
      }
      local maxwellianLua = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "Lua"
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C"
      }

      -- Project the primitive moments onto configuration space basis.
      projectScalar:setFunc(function(t,xn) return m0Func(t,xn) end)
      projectScalar:advance(0.0, {}, {m0Fld})
      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectVec:advance(0.0, {}, {uDriftFld})
      projectScalar:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {vtSqFld})

      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {distf})
      local tmMid = Time.clock()
      maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

      -- Now create a Maxwellian with uPar and assuming zero flow velocity in vx and vy.
      local uParFld  = createField(confGrid, confBasis, 1)
      local uParFunc = function (t, xn) return uDrift[3] end
      projectScalar:setFunc(function(t,xn) return uParFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {uParFld})
      local uParMaxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
      }
      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uParFld,vtSqFld}, {distf})
      local tmMid = Time.clock()
      uParMaxwellian:advance(0.0, {m0Fld,uParFld,vtSqFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end
   end
end

function test_2x2v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75}
   local vt     = 1.0

   local lower     = {-0.50, -0.50, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  0.50,  6.0*vt,  6.0*vt}
   local numCells  = {6, 12, 8, 8}

   for polyOrder = 1, 3 do

      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid   = createGrid({lower[1],lower[2]}, {upper[1],upper[2]}, {numCells[1],numCells[2]})
      local confBasis  = createBasis(confGrid:ndim(), polyOrder)

      local m0Fld     = createField(confGrid, confBasis)
      local uDriftFld = createField(confGrid, confBasis, #uDrift)
      local vtSqFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift[1], uDrift[2] end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local maxwellianFunc = function (t, xn)
         local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
         local uDriftX, uDriftY = uDriftFunc(t,xn)
         local fOut = (m0Func(t,xn)/(2.*math.pi*vtSqFunc(t,xn)))*math.exp(-((vx-uDriftX)^2+(vy-uDriftY)^2)/(2.*(vtFunc(t,xn)))) 
         return fOut
      end

      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projectVec = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0, 1.0 end   -- Set later.
      }
      local maxwellianLua = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "Lua"
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C"
      }

      -- Project the primitive moments onto configuration space basis.
      projectScalar:setFunc(function(t,xn) return m0Func(t,xn) end)
      projectScalar:advance(0.0, {}, {m0Fld})
      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectVec:advance(0.0, {}, {uDriftFld})
      projectScalar:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {vtSqFld})

      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {distf})
      local tmMid = Time.clock()
      maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
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

      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid   = createGrid({lower[1],lower[2]}, {upper[1],upper[2]}, {numCells[1],numCells[2]})
      local confBasis  = createBasis(confGrid:ndim(), polyOrder)

      local m0Fld     = createField(confGrid, confBasis)
      local uDriftFld = createField(confGrid, confBasis)
      local vtSqFld   = createField(confGrid, confBasis)
      local bmagFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local bmagFunc   = function (t, xn) return B0 end
      local maxwellianFunc = function (t, xn)
         local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
         local fOut = (m0Func(t,xn)/(math.sqrt(2.*math.pi*vtSqFunc(t,xn))^3))
                     *math.exp(-((vpar-uDriftFunc(t,xn))^2+2*mu*bmagFunc(t,xn)/mass)/(2.*(vtSqFunc(t,xn)))) 
         return fOut
      end

      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projectVec = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0, 1.0, 1.0 end   -- Set later.
      }
      local maxwellianLua = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "Lua",
         mass           = mass,
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
         mass           = mass,
      }

      -- Project the primitive moments onto configuration space basis.
      projectScalar:setFunc(function(t,xn) return m0Func(t,xn) end)
      projectScalar:advance(0.0, {}, {m0Fld})
      projectScalar:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {uDriftFld})
      projectScalar:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {vtSqFld})
      projectScalar:setFunc(function(t,xn) return bmagFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {bmagFld})

      -- Do projection.
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
      maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

      -- Now create a gyrokinetic Maxwellian taking a 3v flow velocity.
      local uDriftFld  = createField(confGrid, confBasis, 3)
      local uDriftFunc = function (t, xn) return 1.5, 2.5, uDrift end
      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectVec:advance(0.0, {}, {uDriftFld})
      local uDriftMaxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
         mass           = mass,
      }
      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
      local tmMid = Time.clock()
      uDriftMaxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

   end

end

function test_2x3v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75, 0.1}
   local vt     = 1.0

   local lower     = {-0.50, -0.50, -6.0*vt, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  0.50,  6.0*vt,  6.0*vt,  6.0*vt}
   local numCells  = {6, 12, 8, 8, 4}

   for polyOrder = 1, 2 do

      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid   = createGrid({lower[1],lower[2]}, {upper[1],upper[2]}, {numCells[1],numCells[2]})
      local confBasis  = createBasis(confGrid:ndim(), polyOrder)

      local m0Fld     = createField(confGrid, confBasis)
      local uDriftFld = createField(confGrid, confBasis, #uDrift)
      local vtSqFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift[1], uDrift[2], uDrift[3] end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local maxwellianFunc = function (t, xn)
         local x, y, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]
         local fOut = (m0/(math.sqrt(2.*math.pi*vt^2)^3))
                     *math.exp(-((vx-uDrift[1])^2+(vy-uDrift[2])^2+(vz-uDrift[3])^2)/(2.*(vt^2))) 
         return fOut
      end

      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projectVec = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0, 1.0, 1.0 end   -- Set later.
      }
      local maxwellianLua = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "Lua"
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C"
      }

      -- Project the primitive moments onto configuration space basis.
      projectScalar:setFunc(function(t,xn) return m0Func(t,xn) end)
      projectScalar:advance(0.0, {}, {m0Fld})
      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectVec:advance(0.0, {}, {uDriftFld})
      projectScalar:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {vtSqFld})

      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {distf})
      local tmMid = Time.clock()
      maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

      -- Now create a Maxwellian with uPar and assuming zero flow velocity in vx and vy.
      local uParFld  = createField(confGrid, confBasis, 1)
      local uParFunc = function (t, xn) return uDrift[3] end
      projectScalar:setFunc(function(t,xn) return uParFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {uParFld})
      local uParMaxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
      }
      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uParFld,vtSqFld}, {distf})
      local tmMid = Time.clock()
      uParMaxwellian:advance(0.0, {m0Fld,uParFld,vtSqFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

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

   for polyOrder = 1, 1 do

      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid   = createGrid({lower[1],lower[2],lower[3]}, {upper[1],upper[2],upper[3]}, {numCells[1],numCells[2],numCells[3]})
      local confBasis  = createBasis(confGrid:ndim(), polyOrder)

      local m0Fld     = createField(confGrid, confBasis)
      local uDriftFld = createField(confGrid, confBasis)
      local vtSqFld   = createField(confGrid, confBasis)
      local bmagFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local bmagFunc   = function (t, xn) return B0 end
      local maxwellianFunc = function (t, xn)
         local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
         local fOut = (m0Func(t,xn)/(math.sqrt(2.*math.pi*vtSqFunc(t,xn))^3))
                     *math.exp(-((vpar-uDriftFunc(t,xn))^2+2*mu*bmagFunc(t,xn)/mass)/(2.*(vtSqFunc(t,xn)))) 
         return fOut
      end

      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projectVec = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0, 1.0, 1.0 end   -- Set later.
      }
      local maxwellianLua = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "Lua",
         mass           = mass,
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
         mass           = mass,
      }

      -- Project the primitive moments onto configuration space basis.
      projectScalar:setFunc(function(t,xn) return m0Func(t,xn) end)
      projectScalar:advance(0.0, {}, {m0Fld})
      projectScalar:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {uDriftFld})
      projectScalar:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {vtSqFld})
      projectScalar:setFunc(function(t,xn) return bmagFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {bmagFld})

      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
      local tmMid = Time.clock()
      maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})
      local tmEnd = Time.clock()

      -- Check projection.
      local indexer    = distf:genIndexer()
      local distFPtr   = distf:get(1)
      local fMPtr      = fM:get(1)
      local localRange = distf:localRange()
      for idx in localRange:rowMajorIter() do
         distfPtr = distf:get(indexer(idx))
         fMPtr    = fM:get(indexer(idx))
         for k = 1, 1 do --phaseBasis:numBasis() do
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

      -- Now create a gyrokinetic Maxwellian taking a 3v flow velocity.
      local uDriftFld  = createField(confGrid, confBasis, 3)
      local uDriftFunc = function (t, xn) return 0.13, 0.23, uDrift end
      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectVec:advance(0.0, {}, {uDriftFld})
      local uDriftMaxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
         mass           = mass,
      }
      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {distf})
      local tmMid = Time.clock()
      uDriftMaxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld,bmagFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end


   end

end

function test_3x3v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75, 0.1}
   local vt     = 1.0

   local lower     = {-0.50, -0.50, -0.50, -6.0*vt, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  0.50,  0.50,  6.0*vt,  6.0*vt,  6.0*vt}
   local numCells  = {4, 4, 4, 4, 4, 4}

   for polyOrder = 1, 1 do

      local phaseGrid  = createGrid(lower, upper, numCells)
      local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
      local confGrid   = createGrid({lower[1],lower[2],lower[3]}, {upper[1],upper[2],upper[3]}, {numCells[1],numCells[2],numCells[3]})
      local confBasis  = createBasis(confGrid:ndim(), polyOrder)

      local m0Fld     = createField(confGrid, confBasis)
      local uDriftFld = createField(confGrid, confBasis, #uDrift)
      local vtSqFld   = createField(confGrid, confBasis)
      local distf     = createField(phaseGrid, phaseBasis)
      local fM        = createField(phaseGrid, phaseBasis)

      local m0Func     = function (t, xn) return m0*(2.0+math.cos(2.*math.pi*xn[1])) end
      local uDriftFunc = function (t, xn) return uDrift[1], uDrift[2], uDrift[3] end
      local vtSqFunc   = function (t, xn) return vt^2 end
      local maxwellianFunc = function (t, xn)
         local x, y, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]
         local uDriftX, uDriftY, uDriftZ = uDriftFunc(t,xn)
         local fOut = (m0Func(t,xn)/(math.sqrt(2.*math.pi*vtFunc(t,xn))^3))
                     *math.exp(-((vx-uDriftX)^2+(vy-uDriftY)^2+(vz-uDriftZ)^2)/(2.*(vtFunc(t,xn)))) 
         return fOut
      end

      local projectScalar = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0 end   -- Set later.
      }
      local projectVec = Updater.ProjectOnBasis {
         onGrid   = confGrid,
         basis    = confBasis,
         evaluate = function (t, xn) return 1.0, 1.0, 1.0 end   -- Set later.
      }
      local maxwellianLua = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "Lua"
      }
      local maxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C"
      }

      -- Project the primitive moments onto configuration space basis.
      projectScalar:setFunc(function(t,xn) return m0Func(t,xn) end)
      projectScalar:advance(0.0, {}, {m0Fld})
      projectVec:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
      projectVec:advance(0.0, {}, {uDriftFld})
      projectScalar:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {vtSqFld})

      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {distf})
      local tmMid = Time.clock()
      maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end

      -- Now create a Maxwellian with uPar and assuming zero flow velocity in vx and vy.
      local uParFld  = createField(confGrid, confBasis, 1)
      local uParFunc = function (t, xn) return uDrift[3] end
      projectScalar:setFunc(function(t,xn) return uParFunc(t,xn) end)
      projectScalar:advance(0.0, {}, {uParFld})
      local uParMaxwellian = Updater.MaxwellianOnBasis {
         onGrid         = phaseGrid,
         phaseBasis     = phaseBasis,
         confGrid       = confGrid,
         confBasis      = confBasis,
         implementation = "C",
      }
      -- Do projection.
      local tmStart = Time.clock()
      maxwellianLua:advance(0.0, {m0Fld,uParFld,vtSqFld}, {distf})
      local tmMid = Time.clock()
      uParMaxwellian:advance(0.0, {m0Fld,uParFld,vtSqFld}, {fM})
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
            assert_close(distfPtr[1], fMPtr[1], 1.e-14, "Checking Lua and C implementations of MaxwellianOnBasis.")
         end
      end
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
