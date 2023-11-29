-- Gkyl ------------------------------------------------------------------------
--
-- Test whole FEM MGpoisson solver in 1D.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"

-- .......................... User inputs ............................ --

local testLaplace = false  -- Choose whether to test a Laplace solver, else solve Poisson.

local pOrder    = 1       -- Polynomial order.
local basisName = "Ser"   -- Type of polynomial basis.

local numCells  = {16}    -- Number of cells.
local lower     = {0.0}   -- Lower boundary of the domain.
local upper     = {1.0}   -- Upper boundary of the domain.

local gamma      = 1                 -- V-cycles=1, W-cycles=2.
local relaxType  = {'DampedJacobi','DampedGaussSeidel'}   -- Types of relaxations to test. 
local numRelax   = {{1,2,1000},{1,1,1000}}                  -- Number of pre,post and coarsest-grid smoothings.
local tolerance  = 1.e-6                                  -- Do cycles until reaching this relative residue norm.
local relaxOmega = {2./3., 1.0}                           -- Damping parameter.
 
-- .................. end of user inputs (MAYBE) .................... --

-- Length of the simulation domain.
local Lx = {upper[1]-lower[1]}

-- Create short names for the relaxation algorithms (for file names).
local relaxShortNm = {}
relaxShortNm["DampedJacobi"]      = "wJac"
relaxShortNm["DampedGaussSeidel"] = "wGS"

local tests = {}
-- Lower and Upper homogeneous Dirichlet BCs.
tests[1] = {}
tests[1]["a"]            = 2.0
tests[1]["c"]            = {(tests[1]["a"]/12.0-0.5), 0.0}
tests[1]["mu"]           = 0.0
tests[1]["bcStr"]        = "LxHomoDUxHomoD"
tests[1]["periodicDirs"] = {}
tests[1]["bcLower"]      = { {T="D", V=0.0} }
tests[1]["bcUpper"]      = { {T="D", V=0.0} }
tests[1]["lapInitial"]   = function (xn, lx, nCells)
                              local x    = xn[1]
                              local dkx  = math.pi/lx[1]
                              local fOut = 0
                              for k = 1, nCells[1]/2+1 do
                                 local kx = (k-1)*dkx
                                 fOut = fOut + math.sin(kx*x)
                              end
                              return fOut
                           end
-- Lower homogeneous Neumann, upper homogeneous Dirichlet.
tests[2] = {}
tests[2]["a"]            = 5.0
tests[2]["c"]            = {0.0, (tests[2]["a"]/12.0-0.5)}
tests[2]["mu"]           = 0.0
tests[2]["bcStr"]        = "LxHomoNUxHomoD"
tests[2]["periodicDirs"] = {}
tests[2]["bcLower"]      = { {T="N", V=0.0} }
tests[2]["bcUpper"]      = { {T="D", V=0.0} }
tests[2]["lapInitial"]   = function (xn, lx, nCells)
                              local x    = xn[1]
                              local dkx  = 2.*math.pi/(1.5*lx[1])
                              local fOut = 0
                              for k = 1, nCells[1]/2+1 do
                                 local kx = (2*k-1)*dkx
                                 fOut = fOut + math.cos(kx*x)
                              end
                              return fOut
                           end
-- Lower homogeneous Dirichlet, upper homogeneous Neumann.
tests[3] = {}
tests[3]["a"]            = 5.0
tests[3]["c"]            = {0.0, (tests[3]["a"]/12.0-0.5)}
tests[3]["mu"]           = 1.0
tests[3]["bcStr"]        = "LxHomoDUxHomoN"
tests[3]["periodicDirs"] = {}
tests[3]["bcLower"]      = { {T="D", V=0.0} }
tests[3]["bcUpper"]      = { {T="N", V=0.0} }
tests[3]["lapInitial"]   = function (xn, lx, nCells)
                              local x    = xn[1]
                              local dkx  = 2.*math.pi/(1.5*lx[1])
                              local fOut = 0
                              for k = 1, nCells[1]/2+1 do
                                 local kx = (2*k-1)*dkx
                                 fOut = fOut + math.cos(kx*x+math.pi/2)
                              end
                              return fOut
                           end
-- Periodic.
tests[4] = {}
tests[4]["a"]            = {1, 2, 5}
tests[4]["c"]            = {1, 2, 5}
tests[4]["mu"]           = 0.0
tests[4]["bcStr"]        = "xPeriodic"
tests[4]["periodicDirs"] = {1}
tests[4]["bcLower"]      = { {T="P", V=0.0} }
tests[4]["bcUpper"]      = { {T="P", V=0.0} }
tests[4]["lapInitial"]   = function (xn, lx, nCells)
                              local x    = xn[1]
                              local dkx  = math.pi/lx[1]
                              local fOut = 0
                              for k = 1, nCells[1]/2+1 do
                                 local kx = (k-1)*dkx
                                 fOut = fOut + math.sin(kx*x)
                              end
                              return fOut
                           end


-- Solution (phi) to be projected onto the DG basis.
local function phiFunction(xn, aIn, cIn, muIn, pDirs)
   local x = xn[1]
   if #pDirs < 1 then
      return ((x-muIn)^2)/2.0-aIn*((x-muIn)^4)/12.0+cIn[1]*x+cIn[2]
   else
      local am         = aIn
      local bm         = cIn
      local cosT, sinT = 0.0, 0.0
      local f          = 0.0
      for m = 0,2 do
         cosT = am[m+1]*math.cos(2.*math.pi*m*x)
         sinT = bm[m+1]*math.sin(2.*math.pi*m*x)
         f  = f + (cosT+sinT)
      end
      return f
   end
end

-- Right-side source (rho) to be projected onto the DG basis.
local function srcFunction(xn, aIn, cIn, muIn, pDirs)
   local x = xn[1]
   if #pDirs < 1 then
      return aIn*((x-muIn)^2)-1.0
   else
      local am         = aIn
      local bm         = cIn
      local cosT, sinT = 0.0, 0.0
      local f          = 0.0
      for m = 1,2 do
         cosT = am[m+1]*((2.*math.pi*m)^2)*math.cos(2.*math.pi*m*x)
         sinT = bm[m+1]*((2.*math.pi*m)^2)*math.sin(2.*math.pi*m*x)
         f  = f + (cosT+sinT)
      end
      return f
   end
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

local function createField(grid, basis, gam, rlxKind, numRlx, rlxOmega, tol, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
      metaData = {
         polyOrder  = basis:polyOrder(),
         basisType  = basis:id(),
         -- Metadata for post-processing during development.
         gamma      = gam,        -- V-cycles=1, W-cycles=2.
         relaxType  = rlxKind,     -- Type of relaxation algorithm.
         relaxNum   = numRlx,      -- Number of pre,post and coarsest-grid smoothings.
         relaxOmega = rlxOmega,    -- Relaxation damping parameter.
         tolerance  = tol          -- Do cycles until reaching this relative residue norm.
      }
   }
   return fld
end

local function createProject(grid, basis)
   local projUp = Updater.EvalOnNodes {
      onGrid   = grid,
      basis    = basis,
      evaluate = function(t, xn) return 1 end,
      onGhosts = true,
   }
   return projUp
end


for rlxI = 1, #relaxType do   -- Loop over types of relaxation.
   local relaxKind      = relaxType[rlxI]
   local numRelaxations = numRelax[rlxI]
   local omegaRelax     = relaxOmega[rlxI]
   local relaxStr       = relaxShortNm[relaxKind]

   for iT = 1, #tests do   -- Loop over various boundary condition tests.
   
      -- Grids and basis.
      local grid  = createGrid(lower, upper, numCells, tests[iT]["periodicDirs"])
      local basis = createBasis(grid:ndim(), pOrder, basisName)
      
      -- Fields.
      local phiDG  = createField(grid, basis, relaxKind, numRelaxations, omegaRelax)
      local rhoDG  = createField(grid, basis, relaxKind, numRelaxations, omegaRelax)
      local phiFEM = createField(grid, basis, relaxKind, numRelaxations, omegaRelax)
      local rhoFEM = createField(grid, basis, relaxKind, numRelaxations, omegaRelax)
      
      local project = createProject(grid,basis)
      if testLaplace then
         -- Seed the initial condition with Fourier modes & solve the Laplace equation.
         project:setFunc(
            function (t, xn)
               return tests[iT]["lapInitial"](xn, Lx, numCells)
            end)
         project:advance(0.0, {}, {phiDG})
         rhoDG:clear(0.0)
      else
         -- Zero initial guess and some prescribed source function (rho).
         phiDG:clear(0.0)
         project:setFunc(
            function (t, xn)
               local a  = tests[iT]["a"]
               local c  = tests[iT]["c"]
               local mu = tests[iT]["mu"]
               local periodicDirs = tests[iT]["periodicDirs"]
               return srcFunction(xn, a, c, mu, periodicDirs)
            end)
         project:advance(0.0, {}, {rhoDG})
      end
      phiDG:sync()
      rhoDG:sync()

      -- Create the MG Poisson solver.
      local poissonSlv = Updater.MGpoisson {
         solverType  = 'FEM',
         onGrid      = grid,
         basis       = basis,
         bcLower     = tests[iT]["bcLower"],
         bcUpper     = tests[iT]["bcUpper"],
         relaxType   = relaxKind,        -- DampedJacobi or DampedGaussSeidel
         relaxNum    = numRelaxations,   -- Number of pre,post and coarsest-grid smoothings.
         relaxOmega  = omegaRelax,       -- Relaxation damping parameter.
         gamma       = gamma,            -- V-cycles=1, W-cycles=2.
         tolerance   = tolerance,        -- Do cycles until reaching this relative residue norm.
      }
   
      poissonSlv:advance(0.0,{rhoDG,phiDG},{phiDG})   -- Call the MG solver.
      
      -- Output numerical solution and the realtive residue norm.
      phiDG:write(string.format("phiDG_%sP%i_Nx%i_%s_%sV%i%i_final.bp",basisName,pOrder,numCells[1],tests[iT]["bcStr"],relaxStr,numRelaxations[1],numRelaxations[2]), 0.0)
      poissonSlv.relResNorm:write(string.format("relResidue_RmsV_%sP%i_Nx%i_%s_%sV%i%i_final.bp",basisName,pOrder,numCells[1],tests[iT]["bcStr"],relaxStr,numRelaxations[1],numRelaxations[2]),0.0)
   
      -- Project the solution (phi) onto the basis and write it out (may not be needed).
      if testLaplace then
         phiFEM:clear(0.0)
      else
         project:setFunc(
            function (t, xn)
               local a  = tests[iT]["a"]
               local c  = tests[iT]["c"]
               local mu = tests[iT]["mu"]
               local periodicDirs = tests[iT]["periodicDirs"]
               return phiFunction(xn, a, c, mu, periodicDirs)
            end)
         project:advance(0.0, {}, {phiDG})
      end
      phiDG:write(string.format("phiDGanalyticSol_%sP%i_Nx%i_%s.bp",basisName,pOrder,numCells[1],tests[iT]["bcStr"]), 0.0)
   
   end

end

