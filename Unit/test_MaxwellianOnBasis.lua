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
local Time       = require "Lib.Time"

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
   local polyOrder = 1

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

   local project = Updater.ProjectOnBasis {
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
   project:setFunc(function(t,xn) return m0Func(t,xn) end)
   project:advance(0.0, {}, {m0Fld})
   project:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
   project:advance(0.0, {}, {uDriftFld})
   project:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
   project:advance(0.0, {}, {vtSqFld})

   -- Do projection.
   maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {distf})
   maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {fM})

end

function test_1x2v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75}
   local vt     = 1.0

   local lower     = {-0.50, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  6.0*vt,  6.0*vt}
   local numCells  = {8, 8, 8}
   local polyOrder = 1

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
      local fOut = (m0/(2.*math.pi*vt^2))*math.exp(-((vx-uDrift[1])^2+(vy-uDrift[2])^2)/(2.*(vt^2))) 
      return fOut
   end

   local project = Updater.ProjectOnBasis {
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
   project:setFunc(function(t,xn) return m0Func(t,xn) end)
   project:advance(0.0, {}, {m0Fld})
   project:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
   project:advance(0.0, {}, {uDriftFld})
   project:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
   project:advance(0.0, {}, {vtSqFld})

   -- Do projection.
   maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {distf})
   maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {fM})

end

function test_2x2v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75}
   local vt     = 1.0

   local lower     = {-0.50, -0.50, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  0.50,  6.0*vt,  6.0*vt}
   local numCells  = {6, 12, 8, 8}
   local polyOrder = 1

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
      local fOut = (m0/(2.*math.pi*vt^2))*math.exp(-((vx-uDrift[1])^2+(vy-uDrift[2])^2)/(2.*(vt^2))) 
      return fOut
   end

   local project = Updater.ProjectOnBasis {
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
   project:setFunc(function(t,xn) return m0Func(t,xn) end)
   project:advance(0.0, {}, {m0Fld})
   project:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
   project:advance(0.0, {}, {uDriftFld})
   project:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
   project:advance(0.0, {}, {vtSqFld})

   -- Do projection.
   local tmStart = Time.clock()
   maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {distf})
   local tmMid = Time.clock()
   print(" Lua took = ", tmMid - tmStart, " s")
   maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {fM})
   print(" C took   = ", Time.clock() - tmMid, " s")

end

function test_2x3v()
   local m0     = 1.0
   local uDrift = {0.0, 0.75, 0.1}
   local vt     = 1.0

   local lower     = {-0.50, -0.50, -6.0*vt, -6.0*vt, -6.0*vt}
   local upper     = { 0.50,  0.50,  6.0*vt,  6.0*vt,  6.0*vt}
   local numCells  = {6, 12, 8, 8, 8}
   local polyOrder = 1

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

   local project = Updater.ProjectOnBasis {
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
   project:setFunc(function(t,xn) return m0Func(t,xn) end)
   project:advance(0.0, {}, {m0Fld})
   project:setFunc(function(t,xn) return uDriftFunc(t,xn) end)
   project:advance(0.0, {}, {uDriftFld})
   project:setFunc(function(t,xn) return vtSqFunc(t,xn) end)
   project:advance(0.0, {}, {vtSqFld})

   -- Do projection.
   local tmStart = Time.clock()
   maxwellianLua:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {distf})
   local tmMid = Time.clock()
   print(" Lua took = ", tmMid - tmStart, " s")
   maxwellian:advance(0.0, {m0Fld,uDriftFld,vtSqFld}, {fM})
   print(" C took   = ", Time.clock() - tmMid, " s")

end

-- Run tests.
test_1x1v()
test_1x2v()
test_2x2v()
test_2x3v()

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
