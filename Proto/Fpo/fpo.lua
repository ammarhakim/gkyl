-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"
local Updater = require "Updater"
local ffi = require "ffi"
local xsys = require "xsys"

local _ = require "Proto.Fpo.fpoKernelsCdef"

return function(tbl)
   -- Simulation parameters
   local polyOrder = tbl.polyOrder -- polynomial order

   local dragKernelNm = string.format("fpoDragKernelP%d", polyOrder)
   local dragKernelFn = ffi.C[dragKernelNm]
   local diffKernelNm = string.format("fpoDiffKernelP%d", polyOrder)
   local diffKernelFn = ffi.C[diffKernelNm]
   local momsKernelNm = string.format("fpoMomsKernelP%d", polyOrder)
   local momsKernelFn = ffi.C[momsKernelNm]
   local diagKernelNm = string.format("fpoDiagKernelP%d", polyOrder)
   local diagKernelFn = ffi.C[diagKernelNm]

   local cflFrac = tbl.cflFrac and tbl.cflFrac or 1.0
   local cfl = cflFrac*0.5/(2*polyOrder+1) -- CFL number
   local tEnd = tbl.tEnd
   local nFrames = tbl.nFrames
   local updatePotentials = xsys.pickBool(tbl.updatePotentials, true)
   local doughertyPotentials = xsys.pickBool(tbl.doughertyPotentials, false)
   -- Not sure about this logic...
   if doughertyPotentials then updatePotentials = false end
   local writeDiagnostics = xsys.pickBool(tbl.writeDiagnostics, false)
   local cross = xsys.pickBool(tbl.cross, true)
   if cross then
      cross = 1
   else
      cross = 0
   end
   local cells = tbl.cells
   local lower = tbl.lower
   local upper = tbl.upper

   local periodicDirs = tbl.periodicDirs and tbl.periodicDirs or {}
   local periodicX, periodicY = false, false
   for i = 1, #periodicDirs do
      if periodicDirs[i] == 1 then periodicX = true end
      if periodicDirs[i] == 2 then periodicY = true end
   end

   local tmRosen, tmFpo = 0.0, 0.0


   ----------------------
   -- Grids and Fields --
   ----------------------
   local grid = Grid.RectCart {
      lower = {lower[1], lower[2]},
      upper = {upper[1], upper[2]},
      cells = {cells[1], cells[2]},
      periodicDirs = periodicDirs,
   }

   -- basis functions
   local basis = Basis.CartModalSerendipity {
      ndim = grid:ndim(),
      polyOrder = polyOrder
   }

   -- fields
   local function getField()
      return DataStruct.Field {
         onGrid = grid,
         numComponents = basis:numBasis(),
         ghost = {1, 1},
         metaData = {
            polyOrder = basis:polyOrder(),
            basisType = basis:id(),
         },
      }
   end
   local f = getField()
   local f1 = getField()
   local f2 = getField()
   local fe = getField()
   local fs = getField()
   local fNew = getField()
   local h = getField()
   local g = getField()

   local moms = DataStruct.DynVector { numComponents = 4 }
   local diag = DataStruct.DynVector { numComponents = 3 }


   --------------
   -- Updaters --
   --------------
   local function absorbFunc(dir, tm, idxIn, fIn, fOut)
      for i = 1, basis:numBasis() do
         fOut[i] = 0.0
      end
   end
   local function copyFunc(dir, tm, idxIn, fIn, fOut)
      for i = 1, basis:numBasis() do
         fOut[i] = fIn[i]
      end
   end
   local function openFunc(dir, tm, idxIn, fIn, fOut)
      basis:flipSign(dir, fIn, fOut)
   end

   local bcT = Updater.Bc {
      onGrid = grid,
      boundaryConditions = { {copyFunc} },
      dir = 2,
      edge = "upper",
   }
   local bcB = Updater.Bc {
      onGrid = grid,
      boundaryConditions = { {copyFunc} },
      dir = 2,
      edge = "lower",
   }
   local bcL = Updater.Bc {
      onGrid = grid,
      boundaryConditions = { {copyFunc} },
      dir = 1,
      edge = "lower",
   }
   local bcR = Updater.Bc {
      onGrid = grid,
      boundaryConditions = { {copyFunc} },
      dir = 1,
      edge = "upper",
   }

   local function applyBc(fld)
      bcT:advance(0.0, {}, {fld})
      bcB:advance(0.0, {}, {fld})
      bcL:advance(0.0, {}, {fld})
      bcR:advance(0.0, {}, {fld})
      fld:sync()
   end

   -- projection to apply ICs
   local initDist = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = tbl.init,
   }

   local poisson = Updater.FemPerpPoisson {
      onGrid = grid,
      basis = basis,
      bcBottom = { T = "D", V = 0.0 },
      bcTop = { T = "D", V = 0.0 },
      bcLeft = { T = "D", V = 0.0 },
      bcRight = { T = "D", V = 0.0 },
   }

   local function updateRosenbluthDrag(fIn, hOut)
      local tmStart = Time.clock()
      fs:combine(-1.0, fIn)
      if updatePotentials then
         poisson:advance(0.0, {fs}, {hOut})
      end
      tmRosen = tmRosen + Time.clock()-tmStart
   end
   local function updateRosenbluthDiffusion(hIn, gOut)
      local tmStart = Time.clock()
      if updatePotentials then
         poisson:advance(0.0, {hIn}, {gOut})
      end
      tmRosen = tmRosen + Time.clock()-tmStart
   end

   -----------------
   -- Diagnostics --
   -----------------
   local function calcMoms(tCurr, fIn, momVec)
      local dv = Lin.Vec(3)
      dv[1], dv[2] = grid:dx(1), grid:dx(2)
      local vc = Lin.Vec(3)
      local localRange = fIn:localRange()
      local indexer = fIn:genIndexer()
      local out = Lin.Vec(4)
      out[1] = 0.0
      out[2] = 0.0
      out[3] = 0.0
      out[4] = 0.0

      for idxs in localRange:colMajorIter() do
	 grid:setIndex(idxs)
	 grid:cellCenter(vc)
         local fPtr = fIn:get(indexer(idxs))
	 momsKernelFn(dv:data(), vc:data(), fPtr:data(), out:data()) 
      end
      momVec:appendData(tCurr, out)
      return out
   end

   local function calcDiag(tCurr, fIn, hIn, diagVec)
      local dv = Lin.Vec(3)
      dv[1], dv[2] = grid:dx(1), grid:dx(2)
      local vc = Lin.Vec(3)
      local localRange = fIn:localRange()
      local indexer = fIn:genIndexer()
      local out = Lin.Vec(3)
      out[1] = 0.0
      out[2] = 0.0
      out[3] = 0.0

      for idxs in localRange:colMajorIter() do
	 grid:setIndex(idxs)
	 grid:cellCenter(vc)
         local fPtr = fIn:get(indexer(idxs))
         local hPtr = hIn:get(indexer(idxs))
	 diagKernelFn(dv:data(), vc:data(), fPtr:data(), hPtr:data(), out:data()) 
      end
      diagVec:appendData(tCurr, out)
   end

   local function writeData(fr, tm)
      f:write(string.format("f_%d.bp", fr), tm, fr)
      if updatePotentials then
	 h:write(string.format('h_%d.bp', fr), tm, fr)
	 g:write(string.format('g_%d.bp', fr), tm, fr)
      end
      if writeDiagnostics then
	 moms:write(string.format("moms_%d.bp", fr), tm, fr)
         diag:write(string.format("diag_%d.bp", fr), tm, fr)
      end
   end


   --------------------
   -- Initialization --
   --------------------
   initDist:advance(0.0, {}, {f})
   applyBc(f)
   local momVec = calcMoms(0, f, moms)

   -- Check if drag/diff functions are provided
   local initDragFunc = tbl.initDrag and tbl.initDrag or function(t, xn) return 0.0 end
   local initDiffFunc = tbl.initDiff and tbl.initDiff or function(t, xn) return 0.0 end

   -- Overwrite the init functions if the the Dougherty potentials are turned on
   if doughertyPotentials then
      initDragFunc = function (t, z)
	 local ux = momVec[2]/momVec[1]
	 local uy = momVec[3]/momVec[1]
	 return -0.5*((z[1]-ux)^2 + (z[2]-uy)^2)
      end
      initDiffFunc = function (t, z)
	 local ux = momVec[2]/momVec[1]
	 local uy = momVec[3]/momVec[1]
	 local dvth2 = momVec[4]/momVec[1] - (ux^2 + uy^2)
	 return dvth2/2 * (z[1]^2 + z[2]^2) -- /2 is for dimensions
      end
   end

   local initDrag = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = initDragFunc,
      projectOnGhosts = true,
   }
   local initDiff = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = initDiffFunc,
      projectOnGhosts = true,
   }

   -- update Rosenbluth potentials
   if updatePotentials then
      updateRosenbluthDrag(f, h)
      updateRosenbluthDiffusion(h, g)
   else
      initDrag:advance(0.0, {}, {h})
      initDiff:advance(0.0, {}, {g})
   end

   -- write initial conditions
   if writeDiagnostics then
      calcDiag(0, f, h, diag)
   end
   writeData(0, 0.0)
   if updatePotentials == false then
      h:write(string.format('h_%d.bp', 0), 0.0, 0)
      g:write(string.format('g_%d.bp', 0), 0.0, 0)
   end


   local function forwardEuler(dt, fIn, hIn, gIn, fOut)
      local tmStart = Time.clock()

      local dv = Lin.Vec(3)
      dv[1], dv[2] = grid:dx(1), grid:dx(2)
      local localRange = fIn:localRange()
      local indexer = fIn:genIndexer()
      local idxsR, idxsL = {}, {}
      local idxsT, idxsB = {}, {}

      for idxs in localRange:colMajorIter() do
         idxsR[1] = idxs[1]+1
         idxsL[1] = idxs[1]-1
         idxsR[2] = idxs[2]
         idxsL[2] = idxs[2]
         idxsT[1] = idxs[1]
         idxsB[1] = idxs[1]
         idxsT[2] = idxs[2]+1
         idxsB[2] = idxs[2]-1

         local isTopEdge, isBotEdge = 0, 0
         local isLeftEdge, isRightEdge = 0, 0

         if periodicX == false then
            if idxs[1] == 1 then isLeftEdge = 1 end
            if idxs[1] == cells[1] then isRightEdge = 1 end
         end
         if periodicY == false then
            if idxs[2] == 1 then isBotEdge = 1 end
            if idxs[2] == cells[2] then isTopEdge = 1 end
         end

         local fPtr = fIn:get(indexer(idxs))
         local fRPtr = fIn:get(indexer(idxsR))
         local fLPtr = fIn:get(indexer(idxsL))
         local fTPtr = fIn:get(indexer(idxsT))
         local fBPtr = fIn:get(indexer(idxsB))

         local hPtr = hIn:get(indexer(idxs))
         local hRPtr = hIn:get(indexer(idxsR))
         local hLPtr = hIn:get(indexer(idxsL))
         local hTPtr = hIn:get(indexer(idxsT))
         local hBPtr = hIn:get(indexer(idxsB))

         local gPtr = gIn:get(indexer(idxs))
         local gRPtr = gIn:get(indexer(idxsR))
         local gLPtr = gIn:get(indexer(idxsL))
         local gTPtr = gIn:get(indexer(idxsT))
         local gBPtr = gIn:get(indexer(idxsB))

         local fOutPtr= fOut:get(indexer(idxs))

         dragKernelFn(dt, dv:data(),
		      fPtr:data(), fLPtr:data(), fRPtr:data(), fTPtr:data(), fBPtr:data(),
		      hPtr:data(), hLPtr:data(), hRPtr:data(), hTPtr:data(), hBPtr:data(),
		      isTopEdge, isBotEdge, isLeftEdge, isRightEdge,
		      fOutPtr:data())
         diffKernelFn(dt, dv:data(),
		      fPtr:data(), fLPtr:data(), fRPtr:data(), fTPtr:data(), fBPtr:data(),
		      gPtr:data(), gLPtr:data(), gRPtr:data(), gTPtr:data(), gBPtr:data(),
		      isTopEdge, isBotEdge, isLeftEdge, isRightEdge, cross,
		      fOutPtr:data())
      end

      tmFpo = tmFpo + Time.clock()-tmStart
   end

   local function rk3(dt, fIn, fOut)
      -- Stage 1
      updateRosenbluthDrag(fIn, h)
      updateRosenbluthDiffusion(h, g)
      forwardEuler(dt, fIn, h, g, f1)
      applyBc(f1)

      -- Stage 2
      updateRosenbluthDrag(f1, h)
      updateRosenbluthDiffusion(h, g)
      forwardEuler(dt, f1, h, g, fe)
      f2:combine(3.0/4.0, fIn, 1.0/4.0, fe)
      applyBc(f2)

      -- Stage 3
      updateRosenbluthDrag(f2, h)
      updateRosenbluthDiffusion(h, g)
      forwardEuler(dt, f2, h, g, fe)
      fOut:combine(1.0/3.0, fIn, 2.0/3.0, fe)
      applyBc(fOut)
   end

   -- run simulation with RK3
   return function ()
      local tCurr = 0.0
      local step = 1
      local dt = cfl*grid:dx(1)

      local frameInt = tEnd/nFrames
      local nextFrame = 1
      local isDone = false

      local tmStart = Time.clock()
      while not isDone do
         if (tCurr+dt >= tEnd) then
            isDone = true
            dt = tEnd-tCurr
         end
         print(string.format("Step %d at time %g with dt %g ...", step, tCurr, dt))
         rk3(dt, f, fNew)
         f:copy(fNew)

         if writeDiagnostics then
	    calcMoms(tCurr+dt, f, moms)
            updateRosenbluthDrag(f, h)
            calcDiag(tCurr+dt, f, h, diag)
         end

         tCurr = tCurr+dt
         if tCurr >= nextFrame*frameInt or math.abs(tCurr-nextFrame*frameInt) < 1e-10 then
            writeData(nextFrame, tCurr)
            nextFrame = nextFrame+1
         end
         step = step+1
      end
      local tmTotal = Time.clock()-tmStart

      print("")
      print(string.format("Poisson solver took %g secs", tmRosen))
      print(string.format("FPO solver took %g secs", tmFpo))
      print(string.format("Total run-time %g secs", tmTotal))
   end
end
