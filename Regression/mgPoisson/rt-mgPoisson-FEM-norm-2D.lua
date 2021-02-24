-- Gkyl ------------------------------------------------------------------------
--
-- Test 1D FEM norms used in the MGpoisson updater.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"
local Unit       = require "Unit"

local assert_equal = Unit.assert_equal
local stats        = Unit.stats

-- .......................... User inputs ............................ --

local pOrder    = 1                -- Polynomial order.
local basisName = "Ser"            -- Type of polynomial basis.

local numCells  = {8, 8}           -- Number of cells.
local lower     = {0.0, 0.0}       -- Lower boundary of the domain.
local upper     = {2.0, 1.0}       -- Upper boundary of the domain.

-- .................. end of user inputs (MAYBE) .................... --

local tests = {}
-- Lower and Upper x and y homogeneous Dirichlet BCs.
tests[1] = {}
tests[1]["bcStr"]        = "LxHomoDUxHomoDLyHomoDUyHomoD"
tests[1]["periodicDirs"] = {}
tests[1]["bcLower"]      = { {T="D", V=1.0}, {T="D", V=1.0} }
tests[1]["bcUpper"]      = { {T="D", V=1.0}, {T="D", V=1.0} }
-- Lower x homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[2] = {}
tests[2]["bcStr"]        = "LxHomoNUxHomoDLyHomoDUyHomoD"
tests[2]["periodicDirs"] = {}
tests[2]["bcLower"]      = { {T="N", V=0.0}, {T="D", V=1.0} }
tests[2]["bcUpper"]      = { {T="D", V=1.0}, {T="D", V=1.0} }
-- Upper x homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[3] = {}
tests[3]["bcStr"]        = "LxHomoDUxHomoNLyHomoDUyHomoD"
tests[3]["periodicDirs"] = {}
tests[3]["bcLower"]      = { {T="D", V=1.0}, {T="D", V=1.0} }
tests[3]["bcUpper"]      = { {T="N", V=0.0}, {T="D", V=1.0} }
-- Lower y homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[4] = {}
tests[4]["bcStr"]        = "LxHomoDUxHomoDLyHomoNUyHomoD"
tests[4]["periodicDirs"] = {}
tests[4]["bcLower"]      = { {T="D", V=1.0}, {T="N", V=0.0} }
tests[4]["bcUpper"]      = { {T="D", V=1.0}, {T="D", V=1.0} }
-- Upper y homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[5] = {}
tests[5]["bcStr"]        = "LxHomoDUxHomoDLyHomoDUyHomoN"
tests[5]["periodicDirs"] = {}
tests[5]["bcLower"]      = { {T="D", V=1.0}, {T="D", V=1.0} }
tests[5]["bcUpper"]      = { {T="D", V=1.0}, {T="N", V=0.0} }
-- Periodic.
tests[6] = {}
tests[6]["bcStr"]        = "xPeriodicyPeriodic"
tests[6]["periodicDirs"] = {1,2}
tests[6]["bcLower"]      = { {T="P", V=0.0}, {T="P", V=0.0} }
tests[6]["bcUpper"]      = { {T="P", V=0.0}, {T="P", V=0.0} }
-- x homogeneous Dirichlet BCs, y periodic.
tests[7] = {}
tests[7]["bcStr"]        = "LxHomoDUxHomoDyPeriodic"
tests[7]["periodicDirs"] = {2}
tests[7]["bcLower"]      = { {T="D", V=1.0}, {T="P", V=0.0} }
tests[7]["bcUpper"]      = { {T="D", V=1.0}, {T="P", V=0.0} }
-- x periodic, y homogeneous Dirichlet BCs.
tests[8] = {}
tests[8]["bcStr"]        = "xPeriodicLyHomoDUyHomoD"
tests[8]["periodicDirs"] = {1}
tests[8]["bcLower"]      = { {T="P", V=0.0}, {T="D", V=1.0} }
tests[8]["bcUpper"]      = { {T="P", V=0.0}, {T="D", V=1.0} }
-- Lower x homogeneous Dirichlet BCs, upper x homogeneous Neumann, y periodic.
tests[9] = {}
tests[9]["bcStr"]        = "LxHomoDUxHomoNyPeriodic"
tests[9]["periodicDirs"] = {2}
tests[9]["bcLower"]      = { {T="D", V=1.0}, {T="P", V=0.0} }
tests[9]["bcUpper"]      = { {T="N", V=0.0}, {T="P", V=0.0} }

-- Function to project onto the DG basis.
local function testFunction(xn, pDirs)
   local x = xn[1]
   return 1.0
end

local function createGrid(lo, up, nCells, pDirs)
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
   return fld
end

local function createProject(grid,basis)
   local projUp = Updater.EvalOnNodes {
      onGrid   = grid,
      basis    = basis,
      evaluate = function(t, xn) return 1 end,
      onGhosts = true,
   }
   return projUp
end

for iT = 1, #tests do

   -- Grids and basis.
   local grid  = createGrid(lower, upper, numCells, tests[iT]["periodicDirs"])
   local basis = createBasis(grid:ndim(), pOrder, basisName)
   
   -- Fields.
   local phiDG  = createField(grid, basis)
   local phiFEM = createField(grid, basis)
   
   -- Project the test function onto the basis.
   local project = createProject(grid, basis)
   project:setFunc(
      function (t, xn)
         local periodicDirs = tests[iT]["periodicDirs"]
         return testFunction(xn, periodicDirs)
      end) 
   project:advance(0.0, {}, {phiDG})
   phiDG:sync()
   
   -- Create the MG Poisson solver.
   local poissonSlv = Updater.MGpoisson {
      solverType = 'FEM',
      onGrid     = grid,
      basis      = basis,
      bcLower    = tests[iT]["bcLower"],
      bcUpper    = tests[iT]["bcUpper"],
   }
   
   -- Translate the DG field to an FEM field, and output the FEM field.
   poissonSlv:translateDGtoFEM(phiDG, phiFEM)

   local phiNorm = DataStruct.DynVector { numComponents = 1 }
   
   -- Compute the norm of the phi field.
   poissonSlv.l2normCalcAdv(0.0,{phiFEM}, {phiNorm})
   
   -- Check that the value matches expectation.
   local _, fldNorm = phiNorm:lastData()
   assert_equal(1.0*math.sqrt((upper[1]-lower[1])*(upper[2]-lower[2])), fldNorm[1], "Checking L2 norm")

   -- For periodic domains, also check the integral over the whole domain:
   if poissonSlv.isPeriodicDomain then
      poissonSlv.intCalcAdv(0.0,{phiFEM},{phiNorm})
      local _, fldNorm = phiNorm:lastData()
      assert_equal(1.0*(upper[1]-lower[1])*(upper[2]-lower[2]), fldNorm[1], "Checking M0 norm")
   end

   -- Write norm(s) out.
   phiNorm:write(string.format("phiL2Norm_%sP%i_Nx%iNy%i_%s.bp",basisName,pOrder,numCells[1],numCells[2],tests[iT]["bcStr"]), 0.0)

end

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end



