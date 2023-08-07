-- Gkyl ------------------------------------------------------------------------
--
-- Tests for FEM parallel projection operator.
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
local DecompRegionCalc = require "Lib.CartDecomp"
local Mpi              = require "Comm.Mpi"
local Range      = require "Lib.Range"
local lume       = require "Lib.lume"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats        = Unit.stats

local xsys = require "xsys"
local ffi  = require "ffi"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
   "new, copy, fill, sizeof, typeof, metatype")

function log(msg)
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   if rank == 0 then
      print(msg)
   end
end

function allReduceOneInt(localv)
   local sendbuf, recvbuf = new("int[1]"), new("int[1]")
   sendbuf[0] = localv
   Mpi.Allreduce(sendbuf, recvbuf, 1, Mpi.INT, Mpi.SUM, Mpi.COMM_WORLD)
   return recvbuf[0]
end

local function createGrid(lo, up, nCells, pDirs, decomp)
   pDirs = pDirs or {}
   local gridOut = Grid.RectCart {
      lower = lo,  cells        = nCells,
      upper = up,  periodicDirs = pDirs,
      decomposition = decomp,
   }
   return gridOut
end

local function createBasis(dim, pOrder, bKind)
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

local function createField(grid, basis, nGhost)
   local fld = DataStruct.Field {
      onGrid        = grid,
      numComponents = basis:numBasis(),
      ghost         = {nGhost, nGhost},
      metaData = {polyOrder  = basis:polyOrder(),
                  basisType  = basis:id(),},
   }
   return fld
end

function test_femparproj_1x(nx, polyOrder)
   log("")
   log("Testing 1D parallel FEM projection...")

   local comm = Mpi.COMM_WORLD
   local sz = Mpi.Comm_size(comm)
   if sz < 2 then
      log("Not running test_copyFieldGlobalFromLocal_1D as numProcs less than 2")
      return
   end

   local lower = {-.50}
   local upper = { .50}
   local cells = {nx}
   local numGhost = 1
   log(string.format("nx=%d polyOrder=%d\n", nx, polyOrder))

   local decomp = DecompRegionCalc.CartProd { cuts = {sz}, }
   local grid = createGrid(lower, upper, cells, nil, decomp)

   local noDecomp = DecompRegionCalc.CartProd {
      cuts = { 1 },
      __serTesting = true, -- Hack to create a decomp without cuts.
   }
   local gridGlobal = createGrid(lower, upper, cells, nil, noDecomp)

   local basis = createBasis(grid:ndim(), polyOrder)

   local phi             = createField(grid, basis, numGhost)
   local phiSmooth       = createField(grid, basis, numGhost)
   local phiGlobal       = createField(gridGlobal, basis, numGhost)
   local phiSmoothGlobal = createField(gridGlobal, basis, numGhost)

   local phiIn        = createField(grid, basis, 0)
   local phiGlobalIn  = createField(gridGlobal, basis, 0)
   local globalBuffer = createField(gridGlobal, basis, 0)

   -- Initialize phi.
   local initPhiFunc = function (t,xn)
      local z = xn[1]
      return math.sin(2.*math.pi*z);
   end
   local initPhi = Updater.ProjectOnBasis {
      onGrid   = grid,  basis = basis,
      evaluate = function (t,xn) return 0. end,  -- Set later.
   }
   initPhi:setFunc(initPhiFunc)
   initPhi:advance(0.,{},{phi})

   -- Copy to a local buffer without ghosts (these two are equivalent here).
   phiIn:copyRangeToRange(phi, phiIn:localRange(), phi:localRange())
--   phi:copyRangeToBuffer(phi:localRange(), phiIn:dataPointer())

   -- Gather flat buffers from other processes.
   Mpi.Allgather(phiIn:dataPointer(), phiIn:size(), phiIn:elemCommType(),
                 globalBuffer:dataPointer(), phiIn:size(), globalBuffer:elemCommType(), comm)

   -- Rearrange into a global field without ghosts, i.e. reconstruct the multi-D array
   -- using the collection of flat buffers we gathered from each process.
   local decompRange = grid:decomposedRange()
   local cutsIdxr    = decompRange:cutsIndexer()
   local ones, upper = {}, {}
   for d=1,decompRange:ndim() do ones[d], upper[d] = 1, decompRange:cuts(d) end
   local cutsRange = Range.Range(ones, upper) -- A range over the MPI decomposition.
   for idx in cutsRange:colMajorIter() do  -- MPI processes are column-major ordered.
      local subdomainID  = cutsIdxr(idx)
      local destRange    = decompRange:subDomain(subdomainID)
      local lv, uv       = destRange:lowerAsVec(), destRange:upperAsVec()
      local destRange    = phiGlobalIn:globalExtRange():subRange(lv,uv)
      local bufferOffset = (subdomainID-1)*phiIn:size()
      phiGlobalIn:copyRangeFromBuffer(destRange, globalBuffer:dataPointer()+bufferOffset)
   end

   -- Copy to a global field with ghosts.
   phiGlobal:copyRangeToRange(phiGlobalIn, phiGlobal:localRange(), phiGlobalIn:localRange())

   local smoothUpd = Updater.FemParproj {
      onGrid = gridGlobal,  periodicParallelDir = false,
      basis  = basis,       onField             = phiGlobal,
   }
   local t1 = os.clock()
   smoothUpd:advance(0.,{phiGlobal},{phiSmoothGlobal})
   log(string.format("1D parallel smooth took %g s\n", os.clock()-t1))

   -- Copy portion of the global field belonging to this process.
   local lv = phiSmooth:localRange():lowerAsVec()
   local uv = phiSmooth:localRange():upperAsVec()
   local sourceRange = phiSmoothGlobal:globalExtRange():subRange(lv,uv)
   phiSmooth:copyRangeToRange(phiSmoothGlobal, phiSmooth:localRange(), sourceRange)

   local err = createField(grid, basis, numGhost)
   err:combine(1.0, phi, -1.0, phiSmooth)

--   phi:write("phi-rough-1d.bp", 0.0)
--   phiSmooth:write("phi-smooth-1d.bp", 0.0)
--   err:write("error-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid = grid,   numComponents = 1,
      basis  = basis,  quantity      = "V2",
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   log(string.format("Average RMS error = %g\n", math.sqrt(lv[1])))
   return math.sqrt(lv[1])
end

function test_femparproj_3x(nx, ny, nz, polyOrder)
   log("")
   log("Testing 3D parallel FEM projection...")

   local comm = Mpi.COMM_WORLD
   local sz = Mpi.Comm_size(comm)
   if sz ~= 4 then
      log("Not running test_femparproj_3x as numProcs is not 4")
      return
   end

   local lower = {-2.0, -2.0, -.50}
   local upper = { 2.0,  2.0,  .50}
   local cells = {nx, ny, nz}
   local ncuts = {1, 2, 2}
   local numGhost = 1
   log(string.format("nx=%d ny=%d nz=%d polyOrder=%d\n", nx, ny, nz, polyOrder))

   local decomp = DecompRegionCalc.CartProd { cuts = ncuts, }
   local grid = createGrid(lower, upper, cells, nil, decomp)

   local noDecomp = DecompRegionCalc.CartProd {
      cuts = { 1, 1, 1 },
      __serTesting = true, -- Hack to create a decomp without cuts.
   }
   local gridGlobal = createGrid(lower, upper, cells, nil, noDecomp)

   local basis = createBasis(grid:ndim(), polyOrder)

   local phi             = createField(grid, basis, numGhost)
   local phiSmooth       = createField(grid, basis, numGhost)
   local phiGlobal       = createField(gridGlobal, basis, numGhost)
   local phiSmoothGlobal = createField(gridGlobal, basis, numGhost)

   local phiIn        = createField(grid, basis, 0)
   local phiGlobalIn  = createField(gridGlobal, basis, 0)
   local globalBuffer = createField(gridGlobal, basis, 0)

   -- Initialize phi.
   local initPhiFunc = function (t,xn)
      local x, y, z = xn[1], xn[2], xn[3]
      local mu = {0.2, 0.2}
      local sig = 0.3;
      return math.exp(-(math.pow(x-mu[1],2)+math.pow(y-mu[2],2))/(2.0*sig*sig))*math.sin(2.*math.pi*z);
   end
   local initPhi = Updater.ProjectOnBasis {
      onGrid   = grid,  basis = basis,
      evaluate = function (t,xn) return 0. end,  -- Set later.
   }
   initPhi:setFunc(initPhiFunc)
   initPhi:advance(0.,{},{phi})

   -- Copy to a local buffer without ghosts (these two are equivalent here).
   phiIn:copyRangeToRange(phi, phiIn:localRange(), phi:localRange())
--   phi:copyRangeToBuffer(phi:localRange(), phiIn:dataPointer())

   -- Gather flat buffers from other processes.
   Mpi.Allgather(phiIn:dataPointer(), phiIn:size(), phiIn:elemCommType(),
                 globalBuffer:dataPointer(), phiIn:size(), globalBuffer:elemCommType(), comm)

   -- Rearrange into a global field without ghosts, i.e. reconstruct the multi-D array
   -- using the collection of flat buffers we gathered from each process.
   local decompRange = grid:decomposedRange()
   local cutsIdxr    = decompRange:cutsIndexer()
   local ones, upper = {}, {}
   for d=1,decompRange:ndim() do ones[d], upper[d] = 1, decompRange:cuts(d) end
   local cutsRange = Range.Range(ones, upper) -- A range over the MPI decomposition.
   for idx in cutsRange:colMajorIter() do  -- MPI processes are column-major ordered.
      local subdomainID  = cutsIdxr(idx)
      local destRange    = decompRange:subDomain(subdomainID)
      local lv, uv       = destRange:lowerAsVec(), destRange:upperAsVec()
      local destRange    = phiGlobalIn:globalExtRange():subRange(lv,uv)
      local bufferOffset = (subdomainID-1)*phiIn:size()
      phiGlobalIn:copyRangeFromBuffer(destRange, globalBuffer:dataPointer()+bufferOffset)
   end

   -- Copy to a global field with ghosts.
   phiGlobal:copyRangeToRange(phiGlobalIn, phiGlobal:localRange(), phiGlobalIn:localRange())

   local smoothUpd = Updater.FemParproj {
      onGrid = gridGlobal,   periodicParallelDir = false,
      basis  = basis,        onField             = phiGlobal,
   }
   local t1 = os.clock()
   smoothUpd:advance(0.,{phiGlobal},{phiSmoothGlobal})
   log(string.format("3D parallel smooth took %g s\n", os.clock()-t1))

   -- Copy portion of the global field belonging to this process.
   local lv = phiSmooth:localRange():lowerAsVec()
   local uv = phiSmooth:localRange():upperAsVec()
   local sourceRange = phiSmoothGlobal:globalExtRange():subRange(lv,uv)
   phiSmooth:copyRangeToRange(phiSmoothGlobal, phiSmooth:localRange(), sourceRange)

   local err = createField(grid, basis, numGhost)
   err:combine(1.0, phi, -1.0, phiSmooth)

--   phi:write("phi-rough-3d.bp", 0.0)
--   phiSmooth:write("phi-smooth-3d.bp", 0.0)
--   err:write("error-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid = grid,   numComponents = 1,
      basis  = basis,  quantity      = "V2",
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   log(string.format("Average RMS error = %g\n", math.sqrt(lv[1])))
   return math.sqrt(lv[1])
end

function test_femparproj_3x_zglobalOnly(nx, ny, nz, polyOrder)
   log("")
   log("Testing 3D parallel FEM projection staying local in x-y...")

   local ncuts = {2, 2, 2}
   local comm = Mpi.COMM_WORLD
   local sz = Mpi.Comm_size(comm)
   local numProcs = lume.reduce(ncuts, function(a,b) return a*b end, 1)
   if sz < numProcs then
      log(string.format("Not running test_femparproj_3x_zglobalOnly as numProcs is not %d",numProcs))
      return
   end

   local lower = {-2.0, -2.0, -.50}
   local upper = { 2.0,  2.0,  .50}
   local cells = {nx, ny, nz}
   local numGhost = 1
   log(string.format("nx=%d ny=%d nz=%d polyOrder=%d\n", nx, ny, nz, polyOrder))

   local decomp = DecompRegionCalc.CartProd { cuts = ncuts, }
   local grid = createGrid(lower, upper, cells, nil, decomp)

   -- Get decomp range and other objects needed to creat z and xy communicators.
   local decompRange = grid:decomposedRange()
   local cutsIdxr    = decompRange:cutsIndexer()
   local ones, upper = {}, {}
   for d=1,decompRange:ndim() do ones[d], upper[d] = 1, decompRange:cuts(d) end
   local cutsRange = Range.Range(ones, upper) -- A range over the MPI decomposition.
   local subdomainIdx = {}
   decompRange:cutsInvIndexer()(grid:subGridId(), subdomainIdx) -- Grid ID on this processor.
   local group  = Mpi.Comm_group(comm)
   -- Create a cuts range that only has the subdomains at the same x and y.
   local zCutsRange = cutsRange:deflate({1,1,0}, {subdomainIdx[1],subdomainIdx[2],0})

   -- Create a z communicator.
   local zNumRanks = decompRange:cuts(3)
   local zCommGroupRanks = Lin.IntVec(zNumRanks)
   local j = 0
   for idx in cutsRange:colMajorIter() do  -- MPI processes are column-major ordered.
      if idx[1] == subdomainIdx[1] and idx[2] == subdomainIdx[2] then 
         j = j+1;  zCommGroupRanks[j] = cutsIdxr(idx)-1
      end
   end
   local zGroup = Mpi.Group_incl(group, decompRange:cuts(3), zCommGroupRanks:data());
   local tag    = (subdomainIdx[2]-1)*decompRange:cuts(1)+subdomainIdx[1]-1
   local zComm  = Mpi.Comm_create_group(comm, zGroup, tag);

   -- Create a x-y communicator.
   local xyNumRanks = decompRange:cuts(1)*decompRange:cuts(2) 
   local xyCommGroupRanks = Lin.IntVec(xyNumRanks)
   local j = 0
   for idx in cutsRange:colMajorIter() do  -- MPI processes are column-major ordered.
      if idx[3] == subdomainIdx[3] then 
         j = j+1;  xyCommGroupRanks[j] = cutsIdxr(idx)-1
      end
   end
   local xyGroup = Mpi.Group_incl(group, xyNumRanks, xyCommGroupRanks:data());
   local tag     = subdomainIdx[3]-1
   local xyComm  = Mpi.Comm_create_group(comm, xyGroup, tag);

   Mpi.Group_free(group)
   Mpi.Group_free(zGroup)
   Mpi.Group_free(xyGroup)

   -- Create decomposition object that decomposes in x-y but not z.
   local xyDecomp = DecompRegionCalc.CartProd {
      cuts = { ncuts[1], ncuts[2], 1 },
      comm = xyComm,
   }
   local gridGlobal = createGrid(lower, upper, cells, nil, xyDecomp)

   local basis = createBasis(grid:ndim(), polyOrder)

   -- Allocate fields.
   local phi       = createField(grid, basis, numGhost)
   local phiSmooth = createField(grid, basis, numGhost)
   -- Extra fields needed to solve problem in parallel.
   local phiIn        = createField(grid, basis, 0)
   local phiGlobal    = createField(gridGlobal, basis, numGhost)
   local globalBuffer = createField(gridGlobal, basis, 0)

   -- Initialize phi.
   local initPhiFunc = function (t,xn)
      local x, y, z = xn[1], xn[2], xn[3]
      local mu = {0.2, 0.2}
      local sig = 0.3;
      return math.exp(-(math.pow(x-mu[1],2)+math.pow(y-mu[2],2))/(2.0*sig*sig))*math.sin(2.*math.pi*z);
   end
   local initPhi = Updater.ProjectOnBasis {
      onGrid   = grid,  basis = basis,
      evaluate = function (t,xn) return 0. end,  -- Set later.
   }
   initPhi:setFunc(initPhiFunc)
   initPhi:advance(0.,{},{phi})

   -- Copy to a local buffer without ghosts (these two are equivalent here).
   phiIn:copyRangeToRange(phi, phiIn:localRange(), phi:localRange())
--   phi:copyRangeToBuffer(phi:localRange(), phiIn:dataPointer())

   -- Gather flat buffers from other processes.
   Mpi.Allgather(phiIn:dataPointer(), phiIn:size(), phiIn:elemCommType(),
                 globalBuffer:dataPointer(), phiIn:size(), globalBuffer:elemCommType(), zComm)

   -- Rearrange into a global field, i.e. reconstruct the multi-D array
   -- using the collection of flat buffers we gathered from each process.
   for idx in zCutsRange:colMajorIter() do  -- MPI processes are column-major ordered.
      local subdomainID  = cutsIdxr({subdomainIdx[1],subdomainIdx[2],idx[1]})
      local bufferOffset = (idx[1]-1)*phiIn:size()
      local destRange    = decompRange:subDomain(subdomainID)
      local lv, uv       = destRange:lowerAsVec(), destRange:upperAsVec()
      local destRange    = phiGlobal:localExtRange():subRange(lv,uv)
      phiGlobal:copyRangeFromBuffer(destRange, globalBuffer:dataPointer()+bufferOffset)
   end

   local smoothUpd = Updater.FemParproj {
      onGrid = gridGlobal,   periodicParallelDir = false,
      basis  = basis,        onField             = phiGlobal,
   }
   local t1 = os.clock()
   smoothUpd:advance(0.,{phiGlobal},{phiGlobal})
   log(string.format("3D parallel smooth took %g s\n", os.clock()-t1))

   -- Copy portion of the global field belonging to this process.
   local lv = phiSmooth:localRange():lowerAsVec()
   local uv = phiSmooth:localRange():upperAsVec()
   local sourceRange = phiGlobal:localExtRange():subRange(lv,uv)
   phiSmooth:copyRangeToRange(phiGlobal, phiSmooth:localRange(), sourceRange)

   local err = createField(grid, basis, numGhost)
   err:combine(1.0, phi, -1.0, phiSmooth)

--   phi:write("phi-rough-3d.bp", 0.0)
--   phiSmooth:write("phi-smooth-3d.bp", 0.0)
--   err:write("error-1d.bp", 0.0)

   local calcInt = Updater.CartFieldIntegratedQuantCalc {
      onGrid = grid,   numComponents = 1,
      basis  = basis,  quantity      = "V2",
   }
   local dynVec = DataStruct.DynVector { numComponents = 1 }
   calcInt:advance(0.0, {err}, {dynVec})
   local tm, lv = dynVec:lastData()
   log(string.format("Average RMS error = %g\n", math.sqrt(lv[1])))

   Mpi.Comm_free(zComm)
   Mpi.Comm_free(xyComm)

   return math.sqrt(lv[1])
end

function test_1x_p1()
   log("--- Testing convergence of 1D FEM parallel projection with p=1 ---")
   err1 = test_femparproj_1x(32, 1)
   err2 = test_femparproj_1x(64, 1)
   err3 = test_femparproj_1x(128, 1)
   if err1 and err2 and err3 then
      log(string.format("Order: %g %g\n", err1/err2/2.0, err2/err3/2.0))
      assert_close(4.0, err1/err2/2.0, .05)
      assert_close(4.0, err2/err3/2.0, .05)
   end
   log("")
end

function test_3x_p1()
   log("--- Testing convergence of 3D FEM parallel projection with p=1 ---")
   err1 = test_femparproj_3x(4, 4, 32, 1)
   err2 = test_femparproj_3x(4, 4, 64, 1)
   err3 = test_femparproj_3x(4, 4, 128, 1)
   if err1 and err2 and err3 then
      log(string.format("Order: %g %g\n", err1/err2/2.0, err2/err3/2.0))
      assert_close(4.0, err1/err2/2.0, .05)
      assert_close(4.0, err2/err3/2.0, .05)
   end
   log("")
end

function test_3x_p1_zglobalOnly()
   log("--- Testing convergence of 3D FEM parallel projection with p=1 ---")
   err1 = test_femparproj_3x_zglobalOnly(4, 4, 32, 1)
   err2 = test_femparproj_3x_zglobalOnly(4, 4, 64, 1)
   err3 = test_femparproj_3x_zglobalOnly(4, 4, 128, 1)
   if err1 and err2 and err3 then
      log(string.format("Order: %g %g\n", err1/err2/2.0, err2/err3/2.0))
      assert_close(4.0, err1/err2/2.0, .05)
      assert_close(4.0, err2/err3/2.0, .05)
   end
   log("")
end

test_1x_p1()
test_3x_p1()
test_3x_p1_zglobalOnly()
log("")

totalFail = allReduceOneInt(stats.fail)
totalPass = allReduceOneInt(stats.pass)

if totalFail > 0 then
   log(string.format("\nPASSED %d tests", totalPass))
   log(string.format("**** FAILED %d tests", totalFail))
else
   log(string.format("PASSED ALL %d tests!", totalPass))
end
