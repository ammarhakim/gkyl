-- Gkyl ------------------------------------------------------------------------
--
-- Test the Bc updater used to apply boundary conditions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Unit         = require "Unit"
local Grid         = require "Grid"
local DataStruct   = require "DataStruct"
local Basis        = require "Basis"
local Updater      = require "Updater"
local Range        = require "Lib.Range"
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"

local assert_equal = Unit.assert_equal
local assert_close = Unit.assert_close
local stats        = Unit.stats

local m0     = 1.0
local uDrift = {0.75, 0.25, 0.5}
local vt     = 1.0

local function createGrid(lo, up, nCells, pDirs)
   pDirs = pDirs or {}
   local gridOut = Grid.RectCart {
      lower = lo,  cells        = nCells,
      upper = up,  periodicDirs = pDirs,
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
      metaData      = {polyOrder = basis:polyOrder(),
                       basisType = basis:id(),}
   }
   return fld
end

function getGhostRange(edge, dir, global, globalExt)
   local lv, uv = globalExt:lowerAsVec(), globalExt:upperAsVec()
   if edge == "lower" then
      uv[dir] = global:lower(dir)-1   -- For ghost cells on "left".
   else
      lv[dir] = global:upper(dir)+1   -- For ghost cells on "right".
   end
   return Range.Range(lv, uv)
end

function test_bc(polyOrder, bcDir, bcType, lower, upper, numCells, cdim)
   local confLower, confUpper, confNumCells = {}, {}, {}
   for d = 1, cdim do
      confLower[d], confUpper[d], confNumCells[d] = lower[d], upper[d], numCells[d]
   end

   local phaseGrid  = createGrid(lower, upper, numCells)
   local phaseBasis = createBasis(phaseGrid:ndim(), polyOrder)
   local confGrid   = createGrid(confLower, confUpper, confNumCells)
   local confBasis  = createBasis(confGrid:ndim(), polyOrder)

   local pdim = phaseGrid:ndim()
   local vdim = pdim-cdim

   local m0Fld = createField(confGrid, confBasis)
   local distf = createField(phaseGrid, phaseBasis)

   -- Conf space function to project
   local m0Func_1x = function(t, xn)
      local x = xn[1]
      return m0*math.sin(2.*math.pi*x)
   end
   local m0Func_2x = function(t, xn)
      local x, y = xn[1], xn[2]
      return m0*math.sin(2.*math.pi*x)*math.sin(2.*math.pi*y)
   end
   local m0Func_3x = function(t, xn)
      local x, y = xn[1], xn[2], xn[3]
      return m0*math.sin(2.*math.pi*x)*math.sin(2.*math.pi*y)*math.exp(-(z^2)/(2.*(0.15^2)))
   end
   local m0Func
   if cdim==1 then m0Func=m0Func_1x
   elseif cdim==2 then m0Func=m0Func_2x
   elseif cdim==3 then m0Func=m0Func_3x
   end

   -- Velocity space function to project
   local maxwellianFunc_1v = function(t, xn)
      local v = xn[cdim+1]
      local fOut = (1./math.sqrt(2.*math.pi*vt^2))*math.exp(-((v-uDrift[1])^2)/(2*(vt^2)))
      return fOut
   end
   local maxwellianFunc_2v = function(t, xn)
      local vx, vy = xn[cdim+1], xn[cdim+2]
      local fOut = (1./(2.*math.pi*vt^2))*math.exp(-((vx-uDrift[1])^2+(vy-uDrift[2])^2)/(2*(vt^2)))
      return fOut
   end
   local maxwellianFunc_3v = function(t, xn)
      local vx, vy, vz = xn[cdim+1], xn[cdim+2], xn[cdim+3]
      local fOut = (1./((2.*math.pi*vt^2)^(3./2.)))
                  *math.exp(-((vx-uDrift[1])^2+(vy-uDrift[2])^2+(vz-uDrift[3])^2)/(2*(vt^2)))
      return fOut
   end
   local maxwellianFunc
   if vdim==1 then maxwellianFunc=maxwellianFunc_1v
   elseif vdim==2 then maxwellianFunc=maxwellianFunc_2v
   elseif vdim==3 then maxwellianFunc=maxwellianFunc_3v
   end

   local projectConf = Updater.ProjectOnBasis {
      onGrid = confGrid,   evaluate = m0Func,
      basis  = confBasis,  onGhosts = false,
   }
   local projectPhase = Updater.ProjectOnBasis {
      onGrid = phaseGrid,   evaluate = function(t,xn) return m0Func(t,xn)*maxwellianFunc(t,xn) end,
      basis  = phaseBasis,  onGhosts = false,
   }

   -- Do projection.
   projectConf:advance(0.0, {}, {m0Fld})
   projectPhase:advance(0.0, {}, {distf})

   -- Conf BC updaters. 
   local confCopy = function(dir, tm, idxIn, fIn, fOut) -- Use skinLoop = "pointwise".
      for i = 1, confBasis:numBasis() do fOut[i] = fIn[i] end
   end
   local confReflect = function(dir, tm, idxIn, fIn, fOut) -- Use skinLoop = "flip".
      confBasis:flipSign(bcDir, fIn, fOut)
      confBasis:flipSign(bcDir+cdim, fOut, fOut)
   end
   local bcFunc = bcType == "reflect" and confReflect or confCopy
   local confBc = {lower=Updater.Bc{onGrid   = confGrid,   edge               = "lower",
                                    cdim     = cdim,       boundaryConditions = {bcFunc},
                                    dir      = bcDir,      vdir               = cdim+bcDir,
                                    basis    = confBasis,  confBasis          = confBasis,
                                    skinLoop = "pointwise",},
                   upper=Updater.Bc{onGrid   = confGrid,   edge               = "upper",
                                    cdim     = cdim,       boundaryConditions = {bcFunc},
                                    dir      = bcDir,      vdir               = cdim+bcDir,
                                    basis    = confBasis,  confBasis          = confBasis,
                                    skinLoop = "pointwise",}
   }

   -- Phase BC updaters. 
   local phaseCopy = function(dir, tm, idxIn, fIn, fOut) -- Use skinLoop = "pointwise".
      for i = 1, phaseBasis:numBasis() do fOut[i] = fIn[i] end
   end
   local phaseReflect = function(dir, tm, idxIn, fIn, fOut) -- Use skinLoop = "flip".
      phaseBasis:flipSign(bcDir, fIn, fOut)
      phaseBasis:flipSign(bcDir+cdim, fOut, fOut)
   end
   local bcFunc = bcType == "reflect" and phaseReflect or phaseCopy
   local phaseBc = {lower=Updater.Bc{onGrid   = phaseGrid,   edge               = "lower",
                                     cdim     = cdim,        boundaryConditions = {bcFunc},
                                     dir      = bcDir,       vdir               = cdim+bcDir,
                                     basis    = phaseBasis,  confBasis          = confBasis,
                                     skinLoop = "pointwise",},
                    upper=Updater.Bc{onGrid   = phaseGrid,   edge               = "upper",
                                     cdim     = cdim,        boundaryConditions = {bcFunc},
                                     dir      = bcDir,       vdir               = cdim+bcDir,
                                     basis    = phaseBasis,  confBasis          = confBasis,
                                     skinLoop = "pointwise",}
   }


   -- Apply BCs.
   for _, edge in ipairs({"lower","upper"}) do
      confBc[edge]:advance(0., {m0Fld}, {m0Fld})
      phaseBc[edge]:advance(0., {distf}, {distf})
   end

   local idxSkin = Lin.IntVec(pdim)

   -- Check copy BCs.
   local m0Ghost = Lin.Vec(confBasis:numBasis())
   local confGlobal, confGlobalExt = m0Fld:globalRange(), m0Fld:globalExtRange()
   local confLocalExt = m0Fld:localExtRange()
   local m0FldIdxr, m0FldPtr, m0FldSkinPtr = m0Fld:genIndexer(), m0Fld:get(1), m0Fld:get(1)
   for _, edge in ipairs({"lower","upper"}) do
      local confGhost = confLocalExt:intersect(getGhostRange(edge,bcDir,confGlobal,confGlobalExt))   -- Range spanning ghost cells.
      -- Decompose ghost region into threads.
      local confGhostDecomp = LinearDecomp.LinearDecompRange{range = confGhost, numSplit = confGrid:numSharedProcs()}

      local tId = confGrid:subGridSharedId() -- Local thread ID.
      for idx in confGhostDecomp:rowMajorIter(tId) do
         m0Fld:fill(m0FldIdxr(idx), m0FldPtr)

         -- Copy out index into in index
         idx:copyInto(idxSkin)
         idxSkin[bcDir] = edge=="lower" and confGlobal:lower(bcDir) or confGlobal:upper(bcDir)

         m0Fld:fill(m0FldIdxr(idxSkin), m0FldSkinPtr)

         if bcType == "copy" then
            for k = 1, confBasis:numBasis() do
               assert_equal(m0FldPtr[k], m0FldSkinPtr[k], "1x1v_p"..tostring(polyOrder).." "..edge.." edge: Checking m0Fld BC "..bcType)
            end
         elseif bcType == "reflect" then
            confBasis:flipSign(bcDir, m0FldSkinPtr, m0Ghost)
            for k = 1, confBasis:numBasis() do
               assert_equal(m0FldPtr[k], m0Ghost[k], "1x1v_p"..tostring(polyOrder).." "..edge.." edge: Checking m0Fld BC "..bcType)
            end
         end
      end
   end

   local distfGhost = Lin.Vec(phaseBasis:numBasis())
   local phaseGlobal, phaseGlobalExt = distf:globalRange(), distf:globalExtRange()
   local phaseLocalExt = distf:localExtRange()
   local distfIdxr, distfPtr, distfSkinPtr = distf:genIndexer(), distf:get(1), distf:get(1)
   for _, edge in ipairs({"lower","upper"}) do
      local phaseGhost = phaseLocalExt:intersect(getGhostRange(edge,bcDir,phaseGlobal,phaseGlobalExt))   -- Range spanning ghost cells.
      -- Decompose ghost region into threads.
      local phaseGhostDecomp = LinearDecomp.LinearDecompRange{range = phaseGhost, numSplit = phaseGrid:numSharedProcs()}

      local tId = phaseGrid:subGridSharedId() -- Local thread ID.
      for idx in phaseGhostDecomp:rowMajorIter(tId) do
         distf:fill(distfIdxr(idx), distfPtr)

         -- Copy out index into in index
         idx:copyInto(idxSkin)
         idxSkin[bcDir] = edge=="lower" and phaseGlobal:lower(bcDir) or phaseGlobal:upper(bcDir)

         distf:fill(distfIdxr(idxSkin), distfSkinPtr)

         if bcType == "copy" then
            for k = 1, phaseBasis:numBasis() do
               assert_equal(distfPtr[k], distfSkinPtr[k], "1x1v_p"..tostring(polyOrder).." "..edge.." edge: Checking distf BC "..bcType)
            end
         elseif bcType == "reflect" then
            phaseBasis:flipSign(bcDir, distfSkinPtr, distfGhost)
            phaseBasis:flipSign(bcDir+cdim, distfGhost, distfGhost)
            for k = 1, phaseBasis:numBasis() do
               assert_equal(distfPtr[k], distfGhost[k], "1x1v_p"..tostring(polyOrder).." "..edge.." edge: Checking distf BC "..bcType)
            end
         end
      end
   end
end

-- Run tests 1x1v tests.
local lower    = {-0.50, -6.0*vt}
local upper    = { 0.50,  6.0*vt}
local numCells = {8, 8}
test_bc(1, 1, "copy",    lower, upper, numCells, 1)
test_bc(1, 1, "reflect", lower, upper, numCells, 1)
test_bc(2, 1, "copy",    lower, upper, numCells, 1)
test_bc(2, 1, "reflect", lower, upper, numCells, 1)
-- Run tests 1x2v tests.
local lower    = {-0.50, -6.0*vt, -5.0*vt}
local upper    = { 0.50,  6.0*vt,  5.0*vt}
local numCells = {8, 8, 8}
test_bc(1, 1, "copy",    lower, upper, numCells, 1)
test_bc(1, 1, "reflect", lower, upper, numCells, 1)
test_bc(2, 1, "copy",    lower, upper, numCells, 1)
test_bc(2, 1, "reflect", lower, upper, numCells, 1)
-- Run tests 1x3v tests.
local lower    = {-0.50, -6.0*vt, -5.0*vt, -5.50*vt}
local upper    = { 0.50,  6.0*vt,  5.0*vt,  5.50*vt}
local numCells = {8, 8, 8, 8}
test_bc(1, 1, "copy",    lower, upper, numCells, 1)
test_bc(1, 1, "reflect", lower, upper, numCells, 1)
test_bc(2, 1, "copy",    lower, upper, numCells, 1)
test_bc(2, 1, "reflect", lower, upper, numCells, 1)
-- Run tests 2x2v tests.
local lower    = {-0.50, -0.50, -6.0*vt, -5.0*vt}
local upper    = { 0.50,  0.50,  6.0*vt,  5.0*vt}
local numCells = {8, 8, 8, 8}
test_bc(1, 1, "copy",    lower, upper, numCells, 2)
test_bc(1, 1, "reflect", lower, upper, numCells, 2)
test_bc(2, 1, "copy",    lower, upper, numCells, 2)
test_bc(2, 1, "reflect", lower, upper, numCells, 2)
test_bc(1, 2, "copy",    lower, upper, numCells, 2)
test_bc(1, 2, "reflect", lower, upper, numCells, 2)
test_bc(2, 2, "copy",    lower, upper, numCells, 2)
test_bc(2, 2, "reflect", lower, upper, numCells, 2)
-- Run tests 2x3v tests.
local lower    = {-0.50, -0.50, -6.0*vt, -5.0*vt, -5.0*vt}
local upper    = { 0.50,  0.50,  6.0*vt,  5.0*vt,  5.0*vt}
local numCells = {8, 8, 8, 8, 4}
test_bc(1, 1, "copy",    lower, upper, numCells, 2)
test_bc(1, 1, "reflect", lower, upper, numCells, 2)
test_bc(2, 1, "copy",    lower, upper, numCells, 2)
test_bc(2, 1, "reflect", lower, upper, numCells, 2)
test_bc(1, 2, "copy",    lower, upper, numCells, 2)
test_bc(1, 2, "reflect", lower, upper, numCells, 2)
test_bc(2, 2, "copy",    lower, upper, numCells, 2)
test_bc(2, 2, "reflect", lower, upper, numCells, 2)

if stats.fail > 0 then
   print(string.format("\nPASSED %d tests", stats.pass))
   print(string.format("**** FAILED %d tests", stats.fail))
else
   print(string.format("PASSED ALL %d tests!", stats.pass))
end
