-- Gkyl ------------------------------------------------------------------------
--
--
--------------------------------------------------------------------------------

local Grid = require "Grid"
local DataStruct = require "DataStruct"
local Time = require "Lib.Time"
local Basis = require "Basis"
local Updater = require "Updater"
local xsys = require "xsys"

local dragStencil = require "Proto.Fpo.dragStencil"
local diffStencil = require "Proto.Fpo.diffStencil"

return function(tbl)

   -- Simulation parameters
   local polyOrder = tbl.polyOrder -- polynomial order
   local cflFrac = tbl.cflFrac and tbl.cflFrac or 1.0
   local cfl = cflFrac*0.5/(2*polyOrder+1) -- CFL number
   local tEnd = tbl.tEnd
   local nFrames = tbl.nFrames
   local updatePotentials = xsys.pickBool(tbl.updatePotentials, true)
   local writeDiagnostics = xsys.pickBool(tbl.writeDiagnostics, false)

   local cells = tbl.cells
   local lower = tbl.lower
   local upper = tbl.upper

   local tmRosen, tmFpo = 0.0, 0.0

   ----------------------
   -- Grids and fields --
   ----------------------
   local grid = Grid.RectCart {
      lower = {lower[1], lower[2]},
      upper = {upper[1], upper[2]},
      cells = {cells[1], cells[2]},
      periodicDirs = {},
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

   local M0 = DataStruct.DynVector { numComponents = 1 }

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
      boundaryConditions = { openFunc },
      dir = 2,
      edge = "upper",
   }
   local bcB = Updater.Bc {
      onGrid = grid,
      boundaryConditions = { openFunc },
      dir = 2,
      edge = "lower",
   }
   local bcL = Updater.Bc {
      onGrid = grid,
      boundaryConditions = { openFunc },
      dir = 1,
      edge = "lower",
   }
   local bcR = Updater.Bc {
      onGrid = grid,
      boundaryConditions = { openFunc },
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
   initDist:advance(0.0, {}, {f})
   applyBc(f)
   f:write("f_0.bp", 0.0, 0)

   local poisson = Updater.FemPerpPoisson {
      onGrid = grid,
      basis = basis,
      bcBottom = { T = "D", V = 0.0 },
      bcTop = { T = "D", V = 0.0 },
      bcLeft = { T = "D", V = 0.0 },
      bcRight = { T = "D", V = 0.0 },
   }
   local initDrag = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = tbl.initDrag,
   }
   local initDiff = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = tbl.initDiff,
   }

   local function updateRosenbluthDrag(fIn, hOut)
      local tmStart = Time.clock()
      fs:combine(-1.0, fIn)
      if updatePotentials then
	 poisson:advance(0.0, {fs}, {hOut})
      else
	 initDrag:advance(0.0, {}, {hOut})
	 applyBc(hOut)
      end
      tmRosen = tmRosen + Time.clock()-tmStart
   end
   local function updateRosenbluthDiffusion(hIn, gOut)
      local tmStart = Time.clock()
      if updatePotentials then
	 poisson:advance(0.0, {hIn}, {gOut})
      else
	 initDiff:advance(0.0, {}, {gOut})
	 applyBc(gOut)
      end
      tmRosen = tmRosen + Time.clock()-tmStart
   end

   -- update Rosenbluth potentials
   updateRosenbluthDrag(f, h)
   updateRosenbluthDiffusion(h, g)
   h:write('h_0.bp', 0.0, 0)
   g:write('g_0.bp', 0.0, 0)

   local M0Calc = Updater.CartFieldIntegratedQuantCalc {
      onGrid = grid,
      basis = basis,
      numComponents = 1,
      quantity = "V"
   }
   if writeDiagnostics then M0:write("f_M0_0.bp", 0.0) end

   local function forwardEuler(dt, fIn, hIn, gIn, fOut)
      local tmStart = Time.clock()

      local dx, dy = grid:dx(1), grid:dx(2)
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

	 local isTopEdge, isBotEdge = false, false
	 local isLeftEdge, isRightEdge = false, false

	 if idxs[1] == 1 then isLeftEdge = true end
	 if idxs[1] == cells[1] then isRightEdge = true end
	 if idxs[2] == 1 then isBotEdge = true end
	 if idxs[2] == cells[2] then isTopEdge = true end

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

	 dragStencil(dt, dx, dy,
	 	     fPtr, fLPtr, fRPtr, fTPtr, fBPtr,
	 	     hPtr, hLPtr, hRPtr, hTPtr, hBPtr,
	 	     fOutPtr)
	 diffStencil(dt, dx, dy,
		     fPtr, fLPtr, fRPtr, fTPtr, fBPtr,
		     gPtr, gLPtr, gRPtr, gTPtr, gBPtr,
		     isTopEdge, isBotEdge, isLeftEdge, isRightEdge,
		     fOutPtr)
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
	 if writeDiagnostics then M0Calc:advance(tCurr+dt, { f }, { M0 }) end
	 if tCurr >= nextFrame*frameInt then
	    f:write(string.format("f_%d.bp", nextFrame), tCurr, nextFrame)
	    if writeDiagnostics then M0:write(string.format("f_M0_%d.bp", nextFrame)) end
	    nextFrame = nextFrame+1
	 end
	 step = step+1
	 tCurr = tCurr+dt
      end
      local tmTotal = Time.clock()-tmStart

      print("")
      print(string.format("Poisson solver took %g secs", tmRosen))
      print(string.format("FPO solver took %g secs", tmFpo))
      print(string.format("Total run-time %g secs", tmTotal))
   end
end
