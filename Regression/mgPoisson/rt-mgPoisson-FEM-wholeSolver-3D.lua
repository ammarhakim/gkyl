-- Gkyl ------------------------------------------------------------------------
--
-- Test whole FEM MGpoisson solver in 2D.
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

local numCells  = {16, 16, 16}     -- Number of cells.
local lower     = {0.0, 0.0, 0.0}   -- Lower boundary of the domain.
local upper     = {1.0, 1.0, 1.0}   -- Upper boundary of the domain.

local gamma      = 1                 -- V-cycles=1, W-cycles=2.
local relaxType  = {'DampedJacobi','DampedGaussSeidel'}   -- Types of relaxations to test. 
local numRelax   = {{1,2,1000},{1,1,1000}}                  -- Number of pre,post and coarsest-grid smoothings.
local tolerance  = 1.e-6                                  -- Do cycles until reaching this relative residue norm.
local relaxOmega = {0.9, 1.0}                           -- Damping parameter.
 
-- .................. end of user inputs (MAYBE) .................... --

-- Length of the simulation domain.
local Lx = {upper[1]-lower[1], upper[2]-lower[2], upper[3]-lower[3]}

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
tests[1]["e"]            = 2.0
tests[1]["f"]            = {(tests[1]["e"]/12.0-0.5), 0.0}
tests[1]["mu"]           = {0.0, 0.0, 0.0}
tests[1]["bcStr"]        = "LxHomoDUxHomoDLyHomoDUyHomoDLzHomoDUzHomoD"
tests[1]["periodicDirs"] = {}
tests[1]["bcLower"]      = { {T="D", V=0.0}, {T="D", V=0.0}, {T="D", V=0.0} }
tests[1]["bcUpper"]      = { {T="D", V=0.0}, {T="D", V=0.0}, {T="D", V=0.0} }
tests[1]["lapInitial"]   = function (xn, lx, nCells)
                              local x, y, z       = xn[1], xn[2], xn[3]
                              local dkx, dky, dkz = math.pi/lx[1], math.pi/lx[2], math.pi/lx[3]
                              local fOut          = 0.0
                              for i = 1, nCells[1]/2+1 do
                                 for j = 1, nCells[2]/2+1 do
                                    for k = 1, nCells[3]/2+1 do
                                       local kx, ky, kz = (i-1)*dkx, (j-1)*dky, (k-1)*dkz
                                       fOut = fOut + math.sin(kx*x)*math.sin(ky*y)*math.sin(kz*z)
                                    end
                                 end
                              end
                              return fOut
                           end

-- Solution (phi) to be projected onto the DG basis.
local function phiFunction(xIn, aIn, bIn, cIn, dIn, eIn, fIn, muIn, pDirs)
   local x, y, z = xIn[1], xIn[2], xIn[3]
   local fx, fy, fz
   if #pDirs < 1 then
      fx = ((x-muIn[1])^2)/2.0-aIn*((x-muIn[1])^4)/12.0+cIn[1]*x+cIn[2]
      fy = ((y-muIn[2])^2)/2.0-bIn*((y-muIn[2])^4)/12.0+dIn[1]*y+dIn[2]
      fz = ((z-muIn[3])^2)/2.0-eIn*((z-muIn[3])^4)/12.0+fIn[1]*z+fIn[2]
   else
      if #pDirs == 2 then
         local amn        = aIn
         local bmn        = bIn
         local cosT, sinT = 0.0, 0.0
         fx, fy, fz = 0.0, 1.0, 1.0
         for m = 0,2 do
            for n = 0,2 do
               for q = 0,2 do
                  cosT = amn[m+1][n+1][q+1]*math.cos(2.*math.pi*m*x)*math.cos(2.*math.pi*n*y)*math.cos(2.*math.pi*q*z)
                  sinT = bmn[m+1][n+1][q+1]*math.sin(2.*math.pi*m*x)*math.sin(2.*math.pi*n*y)*math.sin(2.*math.pi*q*z)
                  fx   = fx + cosT + sinT
               end
            end
         end
      elseif pDirs[1] == 1 then
         local amn        = aIn
         local bmn        = aIn
         local cosT, sinT = 0.0, 0.0
         fx = 0.0
         fy = ((y-muIn[2])^2)/2.0-bIn*((y-muIn[2])^4)/12.0+dIn[1]*y+dIn[2]
         fz = ((z-muIn[3])^2)/2.0-bIn*((z-muIn[3])^4)/12.0+eIn[1]*z+eIn[2]
         for m = 0,2 do
            cosT = amn[m+1]*math.cos(2.*math.pi*m*x)
            sinT = bmn[m+1]*math.sin(2.*math.pi*m*x)
            fx   = fx + cosT + sinT
         end
      elseif pDirs[1] == 2 then
         local amn        = bIn
         local bmn        = bIn
         local cosT, sinT = 0.0, 0.0
         fx = ((x-muIn[1])^2)/2.0-aIn*((x-muIn[1])^4)/12.0+cIn[1]*x+cIn[2]
         fy = 0.0
         fz = ((z-muIn[3])^2)/2.0-bIn*((z-muIn[3])^4)/12.0+eIn[1]*z+eIn[2]
         for n = 0,2 do
            cosT = amn[n+1]*math.cos(2.*math.pi*n*y)
            sinT = bmn[n+1]*math.sin(2.*math.pi*n*y)
            fy   = fy + cosT + sinT
         end
      elseif pDirs[1] == 3 then
         local amn        = eIn
         local bmn        = eIn
         local cosT, sinT = 0.0, 0.0
         fx = ((x-muIn[1])^2)/2.0-aIn*((x-muIn[1])^4)/12.0+cIn[1]*x+cIn[2]
         fy = ((y-muIn[2])^2)/2.0-bIn*((y-muIn[2])^4)/12.0+dIn[1]*y+dIn[2]
         fz = 0.0
         for q = 0,2 do
            cosT = amn[q+1]*math.cos(2.*math.pi*q*z)
            sinT = bmn[q+1]*math.sin(2.*math.pi*q*z)
            fz   = fz + cosT + sinT
         end
      end
   end
   return fx*fy*fz
end

-- Right-side source (rho) to be projected onto the DG basis.
local function srcFunction(xIn, aIn, bIn, cIn, dIn, eIn, fIn, muIn, pDirs)
   local x, y, z = xIn[1], xIn[2], xIn[3]
   if #pDirs < 1 then
      local fx = ((x-muIn[1])^2)/2.0-aIn*((x-muIn[1])^4)/12.0+cIn[1]*x+cIn[2]
      local fy = ((y-muIn[2])^2)/2.0-bIn*((y-muIn[2])^4)/12.0+dIn[1]*y+dIn[2]
      local fz = ((z-muIn[3])^2)/2.0-eIn*((z-muIn[3])^4)/12.0+fIn[1]*z+fIn[2]
      return (aIn*((x-muIn[1])^2)-1.0)*fy*fz+fx*(bIn*((y-muIn[2])^2)-1.0)*fz+fx*fy*(eIn*((z-muIn[3])^2)-1.0)
   else
      if #pDirs == 2 then
         local amn        = aIn
         local bmn        = bIn
         local cosT, sinT = 0.0, 0.0
         local f          = 0.0
         for m = 0,2 do
            for n = 0,2 do
               for n = 0,2 do
                  cosT = amn[m+1][n+1][q+1]*((2.*math.pi*m)^2+(2.*math.pi*n)^2+(2.*math.pi*q)^2)
                        *math.cos(2.*math.pi*m*x)*math.cos(2.*math.pi*n*y)*math.cos(2.*math.pi*q*z)
                  sinT = bmn[m+1][n+1][q+1]*((2.*math.pi*m)^2+(2.*math.pi*n)^2+(2.*math.pi*q)^2)
                        *math.sin(2.*math.pi*m*x)*math.sin(2.*math.pi*n*y)*math.sin(2.*math.pi*q*z)
                  f    = f + cosT + sinT
               end
            end
         end
         return f
      elseif pDirs[1] == 1 then
         local amn         = aIn
         local bmn         = aIn
         local cosT, sinT  = 0.0, 0.0
         local fx, fxp     = 0.0, 0.0
         local fy          = ((y-muIn[2])^2)/2.0-bIn*((y-muIn[2])^4)/12.0+dIn[1]*y+dIn[2]
         local fz          = ((z-muIn[3])^2)/2.0-eIn*((z-muIn[3])^4)/12.0+fIn[1]*z+fIn[2]
         for m = 0,2 do
            cosT  = amn[m+1]*math.cos(2.*math.pi*m*x)
            sinT  = bmn[m+1]*math.sin(2.*math.pi*m*x)
            cosTp = amn[m+1]*((2.*math.pi*m)^2)*math.cos(2.*math.pi*m*x)
            sinTp = bmn[m+1]*((2.*math.pi*m)^2)*math.sin(2.*math.pi*m*x)
            fx    = fx + cosT + sinT
            fxp   = fxp + cosTp + sinTp
         end
         return fxp*fy*fz+fx*(bIn*((y-muIn[2])^2)-1.0)*fz+fx*fy*(eIn*((z-muIn[3])^2)-1.0)
      elseif pDirs[1] == 2 then
         local amn         = bIn
         local bmn         = bIn
         local cosT, sinT  = 0.0, 0.0
         local fx          = ((x-muIn[1])^2)/2.0-aIn*((x-muIn[1])^4)/12.0+cIn[1]*x+cIn[2]
         local fy, fyp     = 0.0, 0.0
         local fz          = ((z-muIn[3])^2)/2.0-eIn*((z-muIn[3])^4)/12.0+fIn[1]*z+fIn[2]
         for n = 0,2 do
            cosT  = amn[n+1]*math.cos(2.*math.pi*n*y)
            sinT  = bmn[n+1]*math.sin(2.*math.pi*n*y)
            cosTp = amn[n+1]*((2.*math.pi*n)^2)*math.cos(2.*math.pi*n*y)
            sinTp = bmn[n+1]*((2.*math.pi*n)^2)*math.sin(2.*math.pi*n*y)
            fy    = fy + cosT + sinT
            fyp   = fyp + cosTp + sinTp
         end
         return (aIn*((x-muIn[1])^2)-1.0)*fy*fz+fx*fyp*fz+fx*fy*(eIn*((z-muIn[3])^2)-1.0)
      elseif pDirs[1] == 3 then
         local amn         = eIn
         local bmn         = eIn
         local cosT, sinT  = 0.0, 0.0
         local fx          = ((x-muIn[1])^2)/2.0-aIn*((x-muIn[1])^4)/12.0+cIn[1]*x+cIn[2]
         local fy          = ((y-muIn[2])^2)/2.0-bIn*((y-muIn[2])^4)/12.0+dIn[1]*y+dIn[2]
         local fz, fzp     = 0.0, 0.0
         for q = 0,2 do
            cosT  = amn[q+1]*math.cos(2.*math.pi*q*z)
            sinT  = bmn[q+1]*math.sin(2.*math.pi*q*z)
            cosTp = amn[q+1]*((2.*math.pi*q)^2)*math.cos(2.*math.pi*q*z)
            sinTp = bmn[q+1]*((2.*math.pi*q)^2)*math.sin(2.*math.pi*q*z)
            fz    = fz + cosT + sinT
            fzp   = fzp + cosTp + sinTp
         end
         return (aIn*((x-muIn[1])^2)-1.0)*fy*fz+fx*fyp*fz+fx*fy*fzp
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
               local b  = tests[iT]["b"]
               local c  = tests[iT]["c"]
               local d  = tests[iT]["d"]
               local e  = tests[iT]["e"]
               local f  = tests[iT]["f"]
               local mu = tests[iT]["mu"]
               local periodicDirs = tests[iT]["periodicDirs"]
               return srcFunction(xn, a, b, c, d, e, f, mu, periodicDirs)
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
      phiDG:write(string.format("phiDG_%sP%i_Nx%iNy%iNz%i_%s_%sV%i%i_final.bp",basisName,pOrder,numCells[1],numCells[2],numCells[3],tests[iT]["bcStr"],relaxStr,numRelaxations[1],numRelaxations[2]), 0.0)
      poissonSlv.relResNorm:write(string.format("relResidue_RmsV_%sP%i_Nx%iNy%iNz%i_%s_%sV%i%i_final.bp",basisName,pOrder,numCells[1],numCells[2],numCells[3],tests[iT]["bcStr"],relaxStr,numRelaxations[1],numRelaxations[2]),0.0)
   
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
               local e  = tests[iT]["e"]
               local f  = tests[iT]["f"]
               local mu = tests[iT]["mu"]
               local periodicDirs = tests[iT]["periodicDirs"]
               return phiFunction(xn, a, b, c, d, e, f, mu, periodicDirs)
            end)
         project:advance(0.0, {}, {phiDG})
      end
      phiDG:write(string.format("phiDGanalyticSol_%sP%i_Nx%iNy%iNz%i_%s.bp",basisName,pOrder,numCells[1],numCells[2],numCells[3],tests[iT]["bcStr"]), 0.0)
   
   end

end

