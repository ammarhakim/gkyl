-- Gkyl ------------------------------------------------------------------------
--
-- Test 2D FEM inter-grid operators in the MGpoisson updater.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Grid       = require "Grid"
local DataStruct = require "DataStruct"
local Basis      = require "Basis"
local Updater    = require "Updater"

-- .......................... User inputs ............................ --

local pOrder    = 1       -- Polynomial order.
local basisName = "Ser"   -- Type of polynomial basis.
 
local numCells  = {{8, 8}, {16, 16}}   -- Number of cells.
local lower     = {0.0, 0.0}           -- Lower boundary of the domain.
local upper     = {1.0, 1.0}           -- Upper boundary of the domain.

-- .................. end of user inputs (MAYBE) .................... --

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
-- Lower x homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[2] = {}
tests[2]["a"]            = 5.0
tests[2]["b"]            = 2.0
tests[2]["c"]            = {0.0, (tests[1]["a"]/12.0-0.5)}
tests[2]["d"]            = {(tests[1]["b"]/12.0-0.5), 0.0}
tests[2]["mu"]           = {0.0, 0.0}
tests[2]["bcStr"]        = "LxHomoNUxHomoDLyHomoDUyHomoD"
tests[2]["periodicDirs"] = {}
tests[2]["bcLower"]      = { {T="N", V=0.0}, {T="D", V=0.0} }
tests[2]["bcUpper"]      = { {T="D", V=0.0}, {T="D", V=0.0} }
-- Upper x homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[3] = {}
tests[3]["a"]            = 5.0
tests[3]["b"]            = 2.0
tests[3]["c"]            = {0.0, (tests[1]["a"]/12.0-0.5)}
tests[3]["d"]            = {(tests[1]["b"]/12.0-0.5), 0.0}
tests[3]["mu"]           = {1.0, 0.0}
tests[3]["bcStr"]        = "LxHomoDUxHomoNLyHomoDUyHomoD"
tests[3]["periodicDirs"] = {}
tests[3]["bcLower"]      = { {T="D", V=0.0}, {T="D", V=0.0} }
tests[3]["bcUpper"]      = { {T="N", V=0.0}, {T="D", V=0.0} }
-- Lower y homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[4] = {}
tests[4]["a"]            = 2.0
tests[4]["b"]            = 5.0
tests[4]["c"]            = {(tests[1]["a"]/12.0-0.5), 0.0}
tests[4]["d"]            = {0.0, (tests[1]["b"]/12.0-0.5)}
tests[4]["mu"]           = {0.0, 0.0}
tests[4]["bcStr"]        = "LxHomoDUxHomoDLyHomoNUyHomoD"
tests[4]["periodicDirs"] = {}
tests[4]["bcLower"]      = { {T="D", V=0.0}, {T="N", V=0.0} }
tests[4]["bcUpper"]      = { {T="D", V=0.0}, {T="D", V=0.0} }
-- Upper y homogeneous Neumann, the rest are homogeneous Dirichlet.
tests[5] = {}
tests[5]["a"]            = 2.0
tests[5]["b"]            = 5.0
tests[5]["c"]            = {(tests[1]["a"]/12.0-0.5), 0.0}
tests[5]["d"]            = {0.0, (tests[1]["b"]/12.0-0.5)}
tests[5]["mu"]           = {0.0, 1.0}
tests[5]["bcStr"]        = "LxHomoDUxHomoDLyHomoDUyHomoN"
tests[5]["periodicDirs"] = {}
tests[5]["bcLower"]      = { {T="D", V=0.0}, {T="D", V=0.0} }
tests[5]["bcUpper"]      = { {T="D", V=0.0}, {T="N", V=0.0} }
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
-- x homogeneous Dirichlet BCs, y periodic.
tests[7] = {}
tests[7]["a"]            = 2.0
tests[7]["b"]            = {1, 2, 5}
tests[7]["c"]            = {(tests[1]["a"]/12.0-0.5), 0.0} 
tests[7]["d"]            = 0.0
tests[7]["mu"]           = {0.0, 0.0}
tests[7]["bcStr"]        = "LxHomoDUxHomoDyPeriodic"
tests[7]["periodicDirs"] = {2}
tests[7]["bcLower"]      = { {T="D", V=0.0}, {T="P", V=0.0} }
tests[7]["bcUpper"]      = { {T="D", V=0.0}, {T="P", V=0.0} }
-- x periodic, y homogeneous Dirichlet BCs.
tests[8] = {}
tests[8]["a"]            = {1, 2, 5}
tests[8]["b"]            = 2.0
tests[8]["c"]            = 0.0
tests[8]["d"]            = {(tests[1]["b"]/12.0-0.5), 0.0}
tests[8]["mu"]           = {0.0, 0.0}
tests[8]["bcStr"]        = "xPeriodicLyHomoDUyHomoD"
tests[8]["periodicDirs"] = {1}
tests[8]["bcLower"]      = { {T="P", V=0.0}, {T="D", V=0.0} }
tests[8]["bcUpper"]      = { {T="P", V=0.0}, {T="D", V=0.0} }
-- Lower x homogeneous Dirichlet BCs, upper x homogeneous Neumann, y periodic.
tests[9] = {}
tests[9]["a"]            = 5.0
tests[9]["b"]            = {1, 2, 5}
tests[9]["c"]            = {0.0, (tests[1]["a"]/12.0-0.5)} 
tests[9]["d"]            = 0.0
tests[9]["mu"]           = {1.0, 0.0}
tests[9]["bcStr"]        = "LxHomoDUxHomoNyPeriodic"
tests[9]["periodicDirs"] = {2}
tests[9]["bcLower"]      = { {T="D", V=0.0}, {T="P", V=0.0} }
tests[9]["bcUpper"]      = { {T="N", V=0.0}, {T="P", V=0.0} }

-- Function to project onto the DG basis.
local function testFunction(xIn, aIn, bIn, cIn, dIn, muIn, pDirs)
   local x, y = xIn[1], xIn[2]
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
               cosT = amn[m+1][n+1]*math.cos(2.*math.pi*m*x) --*math.cos(2.*math.pi*n*y)
               sinT = bmn[m+1][n+1]*math.sin(2.*math.pi*m*x) --*math.sin(2.*math.pi*n*y)
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

local function createGrid(lo,up,nCells,pDirs)
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

local fromIdx, toIdx = 1, 2

local fineGridIdx = fromIdx
if numCells[toIdx][1] > numCells[fromIdx][1] then fineGridIdx = toIdx end

local poissonSlv = {}
for iT = 1, #tests do

   -- Grids and basis.
   local grid  = {createGrid(lower, upper, numCells[fromIdx], tests[iT]["periodicDirs"]),
                  createGrid(lower, upper, numCells[toIdx], tests[iT]["periodicDirs"])}
   local basis = createBasis(grid[fromIdx]:ndim(), pOrder, basisName)

   -- Fields.
   local phiDG  = {createField(grid[fromIdx],basis), createField(grid[toIdx],basis)}
   local phiFEM = {createField(grid[fromIdx],basis), createField(grid[toIdx],basis)}
   
   -- Create the MG Poisson solver.
   poissonSlv[iT] = Updater.MGpoisson {
      solverType = 'FEM',
      onGrid     = grid[fineGridIdx],
      basis      = basis,
      bcLower    = tests[iT]["bcLower"],
      bcUpper    = tests[iT]["bcUpper"],
   }

   for i = 1,2 do
      fromIdx = i
      toIdx   = 3-i

      for j = 1,#phiDG do phiDG[j]:clear(0.0) end
      for j = 1,#phiFEM do phiFEM[j]:clear(0.0) end
      -- Project the test function onto the first grid.
      local project = createProject(grid[fromIdx], basis)
      project:setFunc(
         function (t, xn)
            local a  = tests[iT]["a"]
            local b  = tests[iT]["b"]
            local c  = tests[iT]["c"]
            local d  = tests[iT]["d"]
            local mu = tests[iT]["mu"]
            local periodicDirs = tests[iT]["periodicDirs"]
            return testFunction(xn, a, b, c, d, mu, periodicDirs)
         end)
      project:advance(0.0, {}, {phiDG[fromIdx]})
      phiDG[fromIdx]:sync()
      phiDG[fromIdx]:write(string.format("phiDGin_%sP%i_Nx%iNy%i-Nx%iNy%i_%s.bp",basisName,pOrder,
                                         numCells[fromIdx][1],numCells[fromIdx][2],
                                         numCells[toIdx][1],numCells[toIdx][2],tests[iT]["bcStr"]),0.0)

      -- Translate the DG field to an FEM field, and output the FEM field.
      poissonSlv[iT]:translateDGtoFEM(phiDG[fromIdx], phiFEM[fromIdx])
      phiFEM[fromIdx]:write(string.format("phiFEMin_%sP%i_Nx%iNy%i-Nx%iNy%i_%s.bp",basisName,pOrder,
                                          numCells[fromIdx][1],numCells[fromIdx][2],
                                          numCells[toIdx][1],numCells[toIdx][2],tests[iT]["bcStr"]),0.0)

      -- Restrict/prolong the FEM field.
      if i==1 then
         poissonSlv[iT]:prolongFEM(phiFEM[fromIdx], phiFEM[toIdx])
      elseif i==2 then
         poissonSlv[iT]:restrictFEM(phiFEM[fromIdx], phiFEM[toIdx])
      end

      -- Output the data on the finer mesh.
      phiFEM[toIdx]:write(string.format("phiFEMout_%sP%i_Nx%iNy%i-Nx%iNy%i_%s.bp",basisName,pOrder,
                                        numCells[fromIdx][1],numCells[fromIdx][2],
                                        numCells[toIdx][1],numCells[toIdx][2],tests[iT]["bcStr"]),0.0)

   end

end
