-- Gkyl ------------------------------------------------------------------------
--
-- Test 2D FEM relaxation in the MGpoisson updater.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"

-- .......................... User inputs ............................ --

local testLaplace = true  -- Choose whether to test a Laplace solver, else solve Poisson.

local pOrder    = 1                -- Polynomial order.
local basisName = "Ser"            -- Type of polynomial basis.

local numCells  = {16, 16}         -- Number of cells.
local lower     = {0.0, 0.0}       -- Lower boundary of the domain.
local upper     = {1.0, 1.0}       -- Upper boundary of the domain.

local relaxType  = {'DampedJacobi','DampedGaussSeidel'}   -- Types of relaxations to test.
local numRelax   = {16, 16}                               -- Number of relaxations.
local relaxOmega = {4./5., 1.0}                           -- Damping parameter.

-- .................. end of user inputs (MAYBE) .................... --

-- Length of the simulation domain.
local Lx = {upper[1]-lower[1], upper[2]-lower[2]}

-- Create short names for the relaxation algorithms (for file names).
local relaxShortNm = {}
relaxShortNm["DampedJacobi"]      = "wJac"
relaxShortNm["DampedGaussSeidel"] = "wGS"

local tests = {}
-- Lower and Upper x and y homogeneous Dirichlet BCs.
tests[1] = {}
tests[1]["a"]            = 2.0
tests[1]["b"]            = 2.0
tests[1]["c"]            = {(tests[1]["a"]/12.0-0.5), 0.0}
tests[1]["d"]            = {(tests[1]["b"]/12.0-0.5), 0.0}
tests[1]["mu"]           = {0.0, 0.0}
tests[1]["bcStr"]        = "LxHomoDUxHomoDLyHomoDUyHomoD"
tests[1]["periodicDirs"] = {}
tests[1]["bcLower"]      = { {T="D", V=0.0}, {T="D", V=0.0} }
tests[1]["bcUpper"]      = { {T="D", V=0.0}, {T="D", V=0.0} }
tests[1]["lapInitial"]   = function (xn, lx, nCells)
                              local x, y     = xn[1], xn[2]
                              local dkx, dky = math.pi/lx[1], math.pi/lx[2]
                              local fOut     = 0.0
                              for i = 1, nCells[1]/2+1 do
                                 for j = 1, nCells[2]/2+1 do
                                    local kx, ky = (i-1)*dkx, (j-1)*dky
                                    fOut = fOut + math.sin(kx*x)*math.sin(ky*y)
                                 end
                              end
                              return fOut
                           end
-- Lower x homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[2] = {}
tests[2]["a"]            = 5.0
tests[2]["b"]            = 2.0
tests[2]["c"]            = {0.0, (tests[2]["a"]/12.0-0.5)}
tests[2]["d"]            = {(tests[2]["b"]/12.0-0.5), 0.0}
tests[2]["mu"]           = {0.0, 0.0}
tests[2]["bcStr"]        = "LxHomoNUxHomoDLyHomoDUyHomoD"
tests[2]["periodicDirs"] = {}
tests[2]["bcLower"]      = { {T="N", V=0.0}, {T="D", V=0.0} }
tests[2]["bcUpper"]      = { {T="D", V=0.0}, {T="D", V=0.0} }
tests[2]["lapInitial"]   = function (xn, lx, nCells)
                              local x, y     = xn[1], xn[2]
                              local dkx, dky = 2.*math.pi/(1.5*lx[1]), math.pi/lx[2]
                              local fOut     = 0.0
                              for i = 1, nCells[1]/2+1 do
                                 for j = 1, nCells[2]/2+1 do
                                    local kx, ky = (2*i-1)*dkx, (j-1)*dky
                                    fOut = fOut + math.cos(kx*x)*math.sin(ky*y)
                                 end
                              end
                              return fOut
                           end
-- Upper x homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[3] = {}
tests[3]["a"]            = 5.0
tests[3]["b"]            = 2.0
tests[3]["c"]            = {0.0, (tests[3]["a"]/12.0-0.5)}
tests[3]["d"]            = {(tests[3]["b"]/12.0-0.5), 0.0}
tests[3]["mu"]           = {1.0, 0.0}
tests[3]["bcStr"]        = "LxHomoDUxHomoNLyHomoDUyHomoD"
tests[3]["periodicDirs"] = {}
tests[3]["bcLower"]      = { {T="D", V=0.0}, {T="D", V=0.0} }
tests[3]["bcUpper"]      = { {T="N", V=0.0}, {T="D", V=0.0} }
tests[3]["lapInitial"]   = function (xn, lx, nCells)
                              local x, y     = xn[1], xn[2]
                              local dkx, dky = 2.*math.pi/(1.5*lx[1]), math.pi/lx[2]
                              local fOut     = 0.0
                              for i = 1, nCells[1]/2+1 do
                                 for j = 1, nCells[2]/2+1 do
                                    local kx, ky = (2*i-1)*dkx, (j-1)*dky
                                    fOut = fOut + math.cos(kx*x+math.pi/2)*math.sin(ky*y)
                                 end
                              end
                              return fOut
                           end
-- Lower y homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[4] = {}
tests[4]["a"]            = 2.0
tests[4]["b"]            = 5.0
tests[4]["c"]            = {(tests[4]["a"]/12.0-0.5), 0.0}
tests[4]["d"]            = {0.0, (tests[4]["b"]/12.0-0.5)}
tests[4]["mu"]           = {0.0, 0.0}
tests[4]["bcStr"]        = "LxHomoDUxHomoDLyHomoNUyHomoD"
tests[4]["periodicDirs"] = {}
tests[4]["bcLower"]      = { {T="D", V=0.0}, {T="N", V=0.0} }
tests[4]["bcUpper"]      = { {T="D", V=0.0}, {T="D", V=0.0} }
tests[4]["lapInitial"]   = function (xn, lx, nCells)
                              local x, y     = xn[1], xn[2]
                              local dkx, dky = math.pi/lx[1], 2.*math.pi/(1.5*lx[2])
                              local fOut     = 0.0
                              for i = 1, nCells[1]/2+1 do
                                 for j = 1, nCells[2]/2+1 do
                                    local kx, ky = (i-1)*dkx, (2*j-1)*dky
                                    fOut = fOut + math.sin(kx*x)*math.cos(ky*y)
                                 end
                              end
                              return fOut
                           end
-- Upper y homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[5] = {}
tests[5]["a"]            = 2.0
tests[5]["b"]            = 5.0
tests[5]["c"]            = {(tests[5]["a"]/12.0-0.5), 0.0}
tests[5]["d"]            = {0.0, (tests[5]["b"]/12.0-0.5)}
tests[5]["mu"]           = {0.0, 1.0}
tests[5]["bcStr"]        = "LxHomoDUxHomoDLyHomoDUyHomoN"
tests[5]["periodicDirs"] = {}
tests[5]["bcLower"]      = { {T="D", V=0.0}, {T="D", V=0.0} }
tests[5]["bcUpper"]      = { {T="D", V=0.0}, {T="N", V=0.0} }
tests[5]["lapInitial"]   = function (xn, lx, nCells)
                              local x, y     = xn[1], xn[2]
                              local dkx, dky = math.pi/lx[1], 2.*math.pi/(1.5*lx[2])
                              local fOut     = 0.0
                              for i = 1, nCells[1]/2+1 do
                                 for j = 1, nCells[2]/2+1 do
                                    local kx, ky = (i-1)*dkx, (2*j-1)*dky
                                    fOut = fOut + math.sin(kx*x)*math.cos(ky*y+math.pi/2)
                                 end
                              end
                              return fOut
                           end
-- Periodic.
tests[6] = {}
tests[6]["a"]            = {{1, 2, 5}, {1, 2, 5}, {1, 2, 5}}
tests[6]["b"]            = {{1, 2, 5}, {1, 2, 5}, {1, 2, 5}}
tests[6]["c"]            = 0.0
tests[6]["d"]            = 0.0
tests[6]["mu"]           = 0.0
tests[6]["bcStr"]        = "xPeriodicyPeriodic"
tests[6]["periodicDirs"] = {1,2}
tests[6]["bcLower"]      = { {T="P", V=0.0}, {T="P", V=0.0} }
tests[6]["bcUpper"]      = { {T="P", V=0.0}, {T="P", V=0.0} }
tests[6]["lapInitial"]   = function (xn, lx, nCells)
                              local x, y     = xn[1], xn[2]
                              local dkx, dky = math.pi/lx[1], math.pi/lx[2]
                              local fOut     = 0.0
                              for i = 1, nCells[1]/2+1 do
                                 for j = 1, nCells[2]/2+1 do
                                    local kx, ky = (i-1)*dkx, (j-1)*dky
                                    fOut = fOut + math.sin(kx*x)*math.sin(ky*y)
                                 end
                              end
                              return fOut
                           end
-- x homogeneous Dirichlet BCs, y periodic.
tests[7] = {}
tests[7]["a"]            = 2.0
tests[7]["b"]            = {1, 2, 5}
tests[7]["c"]            = {(tests[7]["a"]/12.0-0.5), 0.0}
tests[7]["d"]            = 0.0
tests[7]["mu"]           = {0.0, 0.0}
tests[7]["bcStr"]        = "LxHomoDUxHomoDyPeriodic"
tests[7]["periodicDirs"] = {2}
tests[7]["bcLower"]      = { {T="D", V=0.0}, {T="P", V=0.0} }
tests[7]["bcUpper"]      = { {T="D", V=0.0}, {T="P", V=0.0} }
tests[7]["lapInitial"]   = function (xn, lx, nCells)
                              local x, y     = xn[1], xn[2]
                              local dkx, dky = math.pi/lx[1], math.pi/lx[2]
                              local fOut     = 0.0
                              for i = 1, nCells[1]/2+1 do
                                 for j = 1, nCells[2]/2+1 do
                                    local kx, ky = (i-1)*dkx, (j-1)*dky
                                    fOut = fOut + math.sin(kx*x)*math.cos(ky*y)
                                 end
                              end
                              return fOut
                           end
-- x periodic, y homogeneous Dirichlet BCs.
tests[8] = {}
tests[8]["a"]            = {1, 2, 5}
tests[8]["b"]            = 2.0
tests[8]["c"]            = 0.0
tests[8]["d"]            = {(tests[8]["b"]/12.0-0.5), 0.0}
tests[8]["mu"]           = {0.0, 0.0}
tests[8]["bcStr"]        = "xPeriodicLyHomoDUyHomoD"
tests[8]["periodicDirs"] = {1}
tests[8]["bcLower"]      = { {T="P", V=0.0}, {T="D", V=0.0} }
tests[8]["bcUpper"]      = { {T="P", V=0.0}, {T="D", V=0.0} }
tests[8]["lapInitial"]   = function (xn, lx, nCells)
                              local x, y     = xn[1], xn[2]
                              local dkx, dky = math.pi/lx[1], math.pi/lx[2]
                              local fOut     = 0.0
                              for i = 1, nCells[1]/2+1 do
                                 for j = 1, nCells[2]/2+1 do
                                    local kx, ky = (i-1)*dkx, (j-1)*dky
                                    fOut = fOut + math.cos(kx*x)*math.sin(ky*y)
                                 end
                              end
                              return fOut
                           end
-- Lower x homogeneous Dirichlet BCs, upper x homogeneous Neumann, y periodic.
tests[9] = {}
tests[9]["a"]            = 5.0
tests[9]["b"]            = {1, 2, 5}
tests[9]["c"]            = {0.0, (tests[9]["a"]/12.0-0.5)}
tests[9]["d"]            = 0.0
tests[9]["mu"]           = {1.0, 0.0}
tests[9]["bcStr"]        = "LxHomoDUxHomoNyPeriodic"
tests[9]["periodicDirs"] = {2}
tests[9]["bcLower"]      = { {T="D", V=0.0}, {T="P", V=0.0} }
tests[9]["bcUpper"]      = { {T="N", V=0.0}, {T="P", V=0.0} }
tests[9]["lapInitial"]   = function (xn, lx, nCells)
                              local x, y     = xn[1], xn[2]
                              local dkx, dky = 2.*math.pi/(1.5*lx[1]), math.pi/lx[2]
                              local fOut     = 0.0
                              for i = 1, nCells[1]/2+1 do
                                 for j = 1, nCells[2]/2+1 do
                                    local kx, ky = (2*i-1)*dkx, (j-1)*dky
                                    fOut = fOut + math.cos(kx*x+math.pi/2)*math.cos(ky*y)
                                 end
                              end
                              return fOut
                           end

-- Solution (phi) to be projected onto the DG basis.
local function phiFunction(xIn, aIn, bIn, cIn, dIn, muIn, pDirs)
   local x  = xIn[1]
   local y  = xIn[2]
   local fx, fy
   if #pDirs < 1 then
      fx = ((x-muIn[1])^2)/2.0-aIn*((x-muIn[1])^4)/12.0+cIn[1]*x+cIn[2]
      fy = ((y-muIn[2])^2)/2.0-bIn*((y-muIn[2])^4)/12.0+dIn[1]*y+dIn[2]
   else
      if #pDirs == 2 then
         local amn        = aIn
         local bmn        = bIn
         local cosT, sinT = 0.0, 0.0
         fx, fy = 0.0, 1.0
         for m = 0,2 do
            for n = 0,2 do
               cosT = amn[m+1][n+1]*math.cos(2.*math.pi*m*x)*math.cos(2.*math.pi*n*y)
               sinT = bmn[m+1][n+1]*math.sin(2.*math.pi*m*x)*math.sin(2.*math.pi*n*y)
               fx   = fx + cosT + sinT
            end
         end
      elseif pDirs[1] == 1 then
         local amn        = aIn
         local bmn        = aIn
         local cosT, sinT = 0.0, 0.0
         fx, fy = 0.0, ((y-muIn[2])^2)/2.0-bIn*((y-muIn[2])^4)/12.0+dIn[1]*y+dIn[2]
         for m = 0,2 do
            cosT = amn[m+1]*math.cos(2.*math.pi*m*x)
            sinT = bmn[m+1]*math.sin(2.*math.pi*m*x)
            fx   = fx + cosT + sinT
         end
      elseif pDirs[1] == 2 then
         local amn        = bIn
         local bmn        = bIn
         local cosT, sinT = 0.0, 0.0
         fx, fy = ((x-muIn[1])^2)/2.0-aIn*((x-muIn[1])^4)/12.0+cIn[1]*x+cIn[2], 0.0
         for n = 0,2 do
            cosT = amn[n+1]*math.cos(2.*math.pi*n*y)
            sinT = bmn[n+1]*math.sin(2.*math.pi*n*y)
            fy   = fy + cosT + sinT
         end
      end
   end
   return fx*fy
end

-- Right-side source (rho) to be projected onto the DG basis.
local function srcFunction(xIn, aIn, bIn, cIn, dIn, muIn, pDirs)
   local x  = xIn[1]
   local y  = xIn[2]
   if #pDirs < 1 then
      local fx = ((x-muIn[1])^2)/2.0-aIn*((x-muIn[1])^4)/12.0+cIn[1]*x+cIn[2]
      local fy = ((y-muIn[2])^2)/2.0-bIn*((y-muIn[2])^4)/12.0+dIn[1]*y+dIn[2]
      return (aIn*((x-muIn[1])^2)-1.0)*fy+fx*(bIn*((y-muIn[2])^2)-1.0)
   else
      if #pDirs == 2 then
         local amn        = aIn
         local bmn        = bIn
         local cosT, sinT = 0.0, 0.0
         local f          = 0.0
         for m = 0,2 do
            for n = 0,2 do
               cosT = amn[m+1][n+1]*((2.*math.pi*m)^2+(2.*math.pi*n)^2)*math.cos(2.*math.pi*m*x)*math.cos(2.*math.pi*n*y)
               sinT = bmn[m+1][n+1]*((2.*math.pi*m)^2+(2.*math.pi*n)^2)*math.sin(2.*math.pi*m*x)*math.sin(2.*math.pi*n*y)
               f    = f + cosT + sinT
            end
         end
         return f
      elseif pDirs[1] == 1 then
         local amn         = aIn
         local bmn         = aIn
         local cosT, sinT  = 0.0, 0.0
         local fx, fxp, fy = 0.0, 0.0, ((y-muIn[2])^2)/2.0-bIn*((y-muIn[2])^4)/12.0+dIn[1]*y+dIn[2]
         for m = 0,2 do
            cosT  = amn[m+1]*math.cos(2.*math.pi*m*x)
            sinT  = bmn[m+1]*math.sin(2.*math.pi*m*x)
            cosTp = amn[m+1]*((2.*math.pi*m)^2)*math.cos(2.*math.pi*m*x)
            sinTp = bmn[m+1]*((2.*math.pi*m)^2)*math.sin(2.*math.pi*m*x)
            fx    = fx + cosT + sinT
            fxp   = fxp + cosTp + sinTp
         end
         return fxp*fy+fx*(bIn*((y-muIn[2])^2)-1.0)
      elseif pDirs[1] == 2 then
         local amn         = bIn
         local bmn         = bIn
         local cosT, sinT  = 0.0, 0.0
         local fx, fy, fyp = ((x-muIn[1])^2)/2.0-aIn*((x-muIn[1])^4)/12.0+cIn[1]*x+cIn[2], 0.0, 0.0
         for n = 0,2 do
            cosT  = amn[n+1]*math.cos(2.*math.pi*n*y)
            sinT  = bmn[n+1]*math.sin(2.*math.pi*n*y)
            cosTp = amn[n+1]*((2.*math.pi*n)^2)*math.cos(2.*math.pi*n*y)
            sinTp = bmn[n+1]*((2.*math.pi*n)^2)*math.sin(2.*math.pi*n*y)
            fy    = fy + cosT + sinT
            fyp   = fyp + cosTp + sinTp
         end
         return (aIn*((x-muIn[1])^2)-1.0)*fy+fx*fyp
      end
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

local function createField(grid, basis, rlxKind, numRlx, rlxOmega, vComp)
   vComp = vComp or 1
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis()*vComp,
      ghost         = {1, 1},
      metaData = {
         polyOrder  = basis:polyOrder(),
         basisType  = basis:id(),
         -- Metadata for post-processing during development.
         relaxType  = rlxKind,     -- Type of relaxation algorithm.
         relaxNum   = numRlx,      -- Number of pre,post and coarsest-grid smoothings.
         relaxOmega = rlxOmega,    -- Relaxation damping parameter.
      }
   }
   return fld
end

local function createProject(grid, basis)
   local projUp = Updater.EvalOnNodes {
      onGrid          = grid,
      basis           = basis,
      projectOnGhosts = true,
      evaluate        = function(t, xn) return 1 end
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
      local phiFEM = createField(grid, basis, relaxKind, numRelaxations, omegaRelax)
      local rhoDG  = createField(grid, basis, relaxKind, numRelaxations, omegaRelax)
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
               local b  = tests[iT]["b"]
               local c  = tests[iT]["c"]
               local d  = tests[iT]["d"]
               local mu = tests[iT]["mu"]
               local periodicDirs = tests[iT]["periodicDirs"]
               return srcFunction(xn, a, b, c, d, mu, periodicDirs)
            end)
         project:advance(0.0, {}, {rhoDG})
      end
      phiDG:sync()
      rhoDG:sync()
      
      -- Create the MG Poisson solver.
      local poissonSlv = Updater.MGpoisson {
         solverType     = 'FEM',
         onGrid         = grid,
         basis          = basis,
         bcLower        = tests[iT]["bcLower"],
         bcUpper        = tests[iT]["bcUpper"],
         relaxType      = relaxKind,
         relaxNum       = {1,1,1},   -- Call the relaxations in a separate loop below.
         relaxOmega     = omegaRelax,
      }
      -- Translate phi and right side source to FEM.
      poissonSlv:translateDGtoFEM(phiDG, phiFEM)
      poissonSlv:translateDGtoFEM(rhoDG, rhoFEM)

      -- Output initial FEM guess.
      phiFEM:write(string.format("phiFEM_%sP%i_Nx%iNy%i_%s_%s_relax%i.bp",basisName,pOrder,numCells[1],numCells[2],tests[iT]["bcStr"],relaxStr,0), 0.0)

      if (#tests[iT]["periodicDirs"] == grid:ndim()) then
         -- For periodic domains subtract the integral of right-side source from the right side.
         poissonSlv.intCalcAdv(0.0,{rhoFEM},{poissonSlv.dynVbuf})
         local  _, intSrc = poissonSlv.dynVbuf:lastData()
         local intSrcVol = intSrc[1]/rhoFEM:grid():gridVolume()
         poissonSlv:accumulateConst(-intSrcVol, rhoFEM)
      end

      -- Project right side source onto FEM basis.
      rhoDG:copy(rhoFEM)   -- Use rhoDG as a temporary buffer.
      poissonSlv:projectFEM(rhoDG,rhoFEM)
      
      -- Relaxations.
      for nR = 1, numRelaxations do
         poissonSlv:relax(1, phiFEM, rhoFEM)
      
         -- Output the latest iterate.
         phiFEM:write(string.format("phiFEM_%sP%i_Nx%iNy%i_%s_%s_relax%i.bp",basisName,pOrder,numCells[1],numCells[2],tests[iT]["bcStr"],relaxStr,nR), nR)
      end

      -- Project the solution (phi) onto the basis and write it out (may not be needed).
      if testLaplace then
         phiFEM:clear(0.0)
      else
         project:setFunc(
            function (t, xn)
               local a  = tests[iT]["a"]
               local b  = tests[iT]["b"]
               local c  = tests[iT]["c"]
               local d  = tests[iT]["d"]
               local mu = tests[iT]["mu"]
               local periodicDirs = tests[iT]["periodicDirs"]
               return phiFunction(xn, a, b, c, d, mu, periodicDirs)
            end)
         project:advance(0.0, {}, {phiDG})
         poissonSlv:translateDGtoFEM(phiDG, phiFEM)
      end
      phiFEM:write(string.format("phiFEManalyticSol_%sP%i_Nx%iNy%i_%s.bp",basisName,pOrder,numCells[1],numCells[2],tests[iT]["bcStr"]), 0.0)

   end
end



