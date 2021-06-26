-------------------------------------------------------------------
-- Gkeyll Proto ---------------------------------------------------
-------------------------------------------------------------------
-- Prototype FPO App
-- FPO generally requires more than 2 cells to do an update on a
-- surface which is currently not supported by the HyperDisCont
-- Updater; therefore, we have the prototype App. Additionally, this
-- App is currently hardcoded for 0x3v case.
--
-- We are solving:
-- df/dt = - d/dvi ( dh/dvi f ) + 1/2 d^2/(dvi dvj) ( d^2g/(dvi dvj) f )
-- where
-- nabla^2 h = - 4 pi f; nabla^2 g = h
-------------------------------------------------------------------

local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"
local Updater = require "Updater"
local ffi = require "ffi"
local xsys = require "xsys"

ffi.cdef[[
double erf(double a);
]]

return function(tbl)
   -- Simulation parameters
   local polyOrder = tbl.polyOrder -- polynomial order
   local cflFrac = tbl.cflFrac and tbl.cflFrac or 1.0
   local fixedDt = tbl.fixedDt
   local tEnd = tbl.tEnd
   local nFrames = tbl.nFrames
   local dynamicPotentials = xsys.pickBool(tbl.dynamicPotentials, true)
   local doughertyPotentials = xsys.pickBool(tbl.doughertyPotentials, false)
   local maxwellianPotentials = xsys.pickBool(tbl.maxwellianPotentials, false)
   local writeDiagnostics = xsys.pickBool(tbl.writeDiagnostics, false)
   local doughertyFix = xsys.pickBool(tbl.doughertyFix, false)
   local cells = tbl.cells
   local lower = tbl.lower
   local upper = tbl.upper

   local periodicDirs = tbl.periodicDirs and tbl.periodicDirs or {}
   local hasPeriodicBC = false
   for i = 1, #periodicDirs do
      if periodicDirs[i] == 1 then hasPeriodicBC = true end
      if periodicDirs[i] == 2 then hasPeriodicBC = true end
      if periodicDirs[i] == 3 then hasPeriodicBC = true end
   end

   local tmRosen, tmFpo = 0.0, 0.0


   -------------------------------------------------------------------
   -- Loading C kernels  ---------------------------------------------
   -------------------------------------------------------------------
   local _ = require "Proto.Fpo.vlasov_fpo_cdef"

   local dragKernelNm = string.format("vlasov_fpo_drag_cell_3x_ser_p%d", polyOrder)
   local dragKernelFn = ffi.C[dragKernelNm]

   local schemeType = tbl.schemeType and tbl.schemeType or "I"

   local diffSurfXLSerNm = string.format("vlasov_fpo_diff%s_surfxL_3x_ser_p%d", schemeType, polyOrder)
   local diffSurfXLSerFn = ffi.C[diffSurfXLSerNm]
   local diffSurfXUSerNm = string.format("vlasov_fpo_diff%s_surfxU_3x_ser_p%d", schemeType, polyOrder)
   local diffSurfXUSerFn = ffi.C[diffSurfXUSerNm]

   local diffSurfYLSerNm = string.format("vlasov_fpo_diff%s_surfyL_3x_ser_p%d", schemeType, polyOrder)
   local diffSurfYLSerFn = ffi.C[diffSurfYLSerNm]
   local diffSurfYUSerNm = string.format("vlasov_fpo_diff%s_surfyU_3x_ser_p%d", schemeType, polyOrder)
   local diffSurfYUSerFn = ffi.C[diffSurfYUSerNm]

   local diffSurfZLSerNm = string.format("vlasov_fpo_diff%s_surfzL_3x_ser_p%d", schemeType, polyOrder)
   local diffSurfZLSerFn = ffi.C[diffSurfZLSerNm]
   local diffSurfZUSerNm = string.format("vlasov_fpo_diff%s_surfzU_3x_ser_p%d", schemeType, polyOrder)
   local diffSurfZUSerFn = ffi.C[diffSurfZUSerNm]

   local diffVolSerNm = string.format("vlasov_fpo_diff%s_vol_3x_ser_p%d", schemeType, polyOrder)
   local diffVolSerFn = ffi.C[diffVolSerNm]

   local momsKernelNm = string.format("vlasov_fpo_moms_3x_ser_p%d", polyOrder)
   local momsKernelFn = ffi.C[momsKernelNm]
   local momM0modKernelNm = string.format("vlasov_fpo_momM0mod_3x_ser_p%d", polyOrder)
   local momM0modKernelFn = ffi.C[momM0modKernelNm]
   local momM2KernelNm = string.format("vlasov_fpo_momM2_3x_ser_p%d", polyOrder)
   local momM2KernelFn = ffi.C[momM2KernelNm]
   local diagKernelNm = string.format("vlasov_fpo_diag_3x_ser_p%d", polyOrder)
   local diagKernelFn = ffi.C[diagKernelNm]

   local fStencil7 = ffi.new("stencil7")
   local hStencil7 = ffi.new("stencil7")

   -------------------------------------------------------------------
   -- Grids and Fields -----------------------------------------------
   -------------------------------------------------------------------
   local grid = Grid.RectCart {
      lower = lower,
      upper = upper,
      cells = cells,
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
         syncCorners = hasPeriodicBC,
      }
   end
   local f = getField()
   local f1 = getField()
   local f2 = getField()
   local fe = getField()
   local fNew = getField()
   local h = getField()
   local g = getField()
   local tmp = getField()

   --local test = getField()
   local src = getField()

   local momVec = DataStruct.DynVector { numComponents = 5 }
   local diagVec = DataStruct.DynVector { numComponents = 4 }


   -------------------------------------------------------------------
   -- Boundary conditions --------------------------------------------
   -------------------------------------------------------------------
   local function applyBc(fld)
      fld:sync()
   end


   -------------------------------------------------------------------
   -- Diagnostics and IO ---------------------------------------------
   -------------------------------------------------------------------
   local function calcMoms(fIn)
      local dv = Lin.Vec(3)
      dv[1], dv[2], dv[3] = grid:dx(1), grid:dx(2), grid:dx(3)
      local vc = Lin.Vec(3)
      local localRange = fIn:localRange()
      local indexer = fIn:genIndexer()
      local out = Lin.Vec(5)
      out[1] = 0.0
      out[2] = 0.0
      out[3] = 0.0
      out[4] = 0.0
      out[5] = 0.0

      for idxs in localRange:colMajorIter() do
         grid:setIndex(idxs)
         grid:cellCenter(vc)
         local fPtr = fIn:get(indexer(idxs))
         momsKernelFn(dv:data(), vc:data(), fPtr:data(), out:data())
      end
      return out
   end

   local function calcVth2(fIn)
      local dv = Lin.Vec(3)
      dv[1], dv[2], dv[3] = grid:dx(1), grid:dx(2), grid:dx(3)
      local vc = Lin.Vec(3)
      local localRange = fIn:localRange()
      local indexer = fIn:genIndexer()
      local M0mod, M2 = 0.0, 0.0

      local isXloEdge, isXupEdge = 0, 0
      local isYloEdge, isYupEdge = 0, 0
      local isZloEdge, isZupEdge = 0, 0
      for idxs in localRange:colMajorIter() do
         -- if not hasPeriodicBC then
         --    if idxs[1] == localRange:lower(1) then isXloEdge = 1 end
         --    if idxs[1] == localRange:upper(1) then isXupEdge = 1 end
         --    if idxs[2] == localRange:lower(2) then isYloEdge = 1 end
         --    if idxs[2] == localRange:upper(2) then isYupEdge = 1 end
         --    if idxs[3] == localRange:lower(3) then isZloEdge = 1 end
         --    if idxs[3] == localRange:upper(3) then isZupEdge = 1 end
         -- end
               
         grid:setIndex(idxs)
         grid:cellCenter(vc)
         local fPtr = fIn:get(indexer(idxs))
         M0mod = M0mod + momM0modKernelFn(dv:data(), vc:data(), fPtr:data(),
                                   isXloEdge, isXupEdge, isYloEdge, isYupEdge, isZloEdge, isZupEdge)
         M2 = M2 + momM2KernelFn(dv:data(), vc:data(), fPtr:data())
      end
      return M2/M0mod
   end

   local function calcDiag(tCurr, fIn, hIn, vec)
      local dv = Lin.Vec(3)
      dv[1], dv[2], dv[3] = grid:dx(1), grid:dx(2), grid:dx(3)
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
         local hPtr = hIn:get(indexer(idxs))
         diagKernelFn(dv:data(), vc:data(), fPtr:data(), hPtr:data(), out:data())
      end
      vec:appendData(tCurr, out)
   end

   local function writeData(fr, tm)
      f:write(string.format("f_%d.bp", fr), tm, fr)
      if dynamicPotentials then
         h:write(string.format('h_%d.bp', fr), tm, fr)
         g:write(string.format('g_%d.bp', fr), tm, fr)
      end

      if writeDiagnostics then
         momVec:write("moms.bp", tm, fr)
         diagVec:write("diag.bp", tm, fr)
      end
   end


   -------------------------------------------------------------------
   -- Projection updater ---------------------------------------------
   -------------------------------------------------------------------
   local projectUpd = Updater.ProjectOnBasis {
      onGrid = grid,
      basis = basis,
      evaluate = function(t, xn) return 0. end,
      onGhosts = true,
   }

   projectUpd:setFunc(tbl.init)
   projectUpd:advance(0.0, {}, {f})
   applyBc(f)
   local moms = calcMoms(f)
   momVec:appendData(0.0, moms)

   -------------------------------------------------------------------
   -- Potentials -----------------------------------------------------
   -------------------------------------------------------------------
   local function getDoughertyH(momsIn)
      return function(t, z)
         local ux = momsIn[2]/momsIn[1]
         local uy = momsIn[3]/momsIn[1]
         local uz = momsIn[4]/momsIn[1]
         return -0.5*((z[1]-ux)^2 + (z[2]-uy)^2 + (z[3]-uz)^2)
      end
   end
   local function getMaxwellianH(momsIn)
      return function(t, z)
         local vx, vy, vz = z[1], z[2], z[3]
         local v = math.sqrt(vx*vx + vy*vy + vz*vz)
         local n = momsIn[1]
         local svth = math.sqrt(2)*math.sqrt(momsIn[5]/momsIn[1]/3)
         return n/v*ffi.C.erf(v/svth) * 2
      end
   end

   local function getDoughertyG(momsIn)
      return function(t, z)
         local ux = momsIn[2]/momsIn[1]
         local uy = momsIn[3]/momsIn[1]
         local uz = momsIn[4]/momsIn[1]
         local vth2 = (momsIn[5]/momsIn[1] - ux^2 - uy^2 - uz^2) / 3
         return vth2 * (z[1]^2 + z[2]^2 + z[3]^2)
      end
   end
   local function getDoughertyGmod(vth2)
      return function(t, z)
         return vth2 * (z[1]^2 + z[2]^2 + z[3]^2)
      end
   end
   local function getMaxwellianG(momsIn)
      return function(t, z)
         local vx, vy, vz = z[1], z[2], z[3]
         local v = math.sqrt(vx*vx + vy*vy + vz*vz)
         local n = momsIn[1]
         local svth = math.sqrt(2)*math.sqrt(momsIn[5]/momsIn[1]/3)
         return n*(svth/math.sqrt(math.pi)*math.exp(-v^2/svth^2) +
                      svth*(0.5*svth/v+v/svth)*ffi.C.erf(v/svth))
      end
   end

   local poissonH = Updater.DiscontPoisson {
      onGrid = grid,
      basis = basis,
      bcFunc = getMaxwellianH(moms)
   }
   local poissonG = Updater.DiscontPoisson {
      onGrid = grid,
      basis = basis,
      bcFunc = getMaxwellianG(moms)
   }

   local function updatePotentials(tCurr, fIn, hOut, gOut, isFirst)
      local tmStart = Time.clock()
      if dynamicPotentials or isFirst then
         local moms = calcMoms(fIn)

         if tbl.fixedH then
            projectUpd:setFunc(tbl.fixedH)
            projectUpd:advance(tCurr, {}, {hOut})
         elseif doughertyPotentials then
            projectUpd:setFunc( getDoughertyH(moms) )
            projectUpd:advance(tCurr, {}, {hOut})
         elseif maxwellianPotentials then
            projectUpd:setFunc( getMaxwellianH(moms) )
            projectUpd:advance(tCurr, {}, {hOut})
         else
            tmp:combine(8*math.pi, fIn)
            poissonH:advance(0.0, {tmp}, {hOut})
         end

         if tbl.fixedG then
            projectUpd:setFunc(tbl.fixedG)
            projectUpd:advance(tCurr, {}, {gOut})
         elseif doughertyPotentials then
            if doughertyFix then
               local vth2 = calcVth2(fIn)
               print(moms[5]/moms[1]/3 - vth2)
               projectUpd:setFunc( getDoughertyGmod(vth2) )
            else
               projectUpd:setFunc( getDoughertyG(moms) )
            end
            projectUpd:advance(tCurr, {}, {gOut})
         elseif maxwellianPotentials then
            projectUpd:setFunc( getMaxwellianG(moms) )
            projectUpd:advance(tCurr, {}, {gOut})
         else
            tmp:combine(-1, hOut)
            poissonG:advance(0.0, {tmp}, {gOut})
         end

         -- HACK to add the missing 1/2 for the drag term in scheme II
         -- (i.e., the same drag term can be used for both scheme I and II)
         if schemeType == "II" then
            hOut:scale(0.5)
         end
      end
      tmRosen = tmRosen + Time.clock()-tmStart
   end


   updatePotentials(0.0, f, h, g, true)

   -- write initial conditions
   if writeDiagnostics then
      calcDiag(0, f, h, diagVec)
   end
   writeData(0, 0.0)
   if dynamicPotentials == false then
      h:write(string.format('h.bp'), 0.0, 0)
      g:write(string.format('g.bp'), 0.0, 0)
   end


   -------------------------------------------------------------------
   -- Additional source ----------------------------------------------
   -------------------------------------------------------------------
   local hasConstSource = false
   if tbl.constSource then
      hasConstSource = true
      -- A wrapper to add time dependance for projectOnBasis
      projectUpd:setFunc(function(t,xn) return tbl.constSource(xn) end)
      projectUpd:advance(0.0, {}, {src})
      if writeDiagnostics then
         src:write('src.bp', 0.0, 0)
      end
   end
   --test:accumulate(1, src)
   --test:write('test.bp', 0.0, 0)


   -------------------------------------------------------------------
   -- Forward Euler --------------------------------------------------
   -------------------------------------------------------------------
   local function forwardEuler(dt, fIn, hIn, gIn, fOut)
      local tmStart = Time.clock()

      local dv = Lin.Vec(3)
      dv[1], dv[2], dv[3] = grid:dx(1), grid:dx(2), grid:dx(3)
      local localRange = fIn:localRange()
      local indexer = fIn:genIndexer()
      local idxs_LLC = {}
      local idxs_LCL = {}
      local idxs_LCC = {}
      local idxs_LCU = {}
      local idxs_LUC = {}
      local idxs_CLL = {}
      local idxs_CLC = {}
      local idxs_CLU = {}
      local idxs_CCL = {}
      --local idxs_CCC = {}
      local idxs_CCU = {}
      local idxs_CUL = {}
      local idxs_CUC = {}
      local idxs_CUU = {}
      local idxs_ULC = {}
      local idxs_UCL = {}
      local idxs_UCC = {}
      local idxs_UCU = {}
      local idxs_UUC = {}

      local cflFreq  = 0.0
      local dragFreq, diffFreq
      
      local isXloEdge, isXupEdge
      local isYloEdge, isYupEdge
      local isZloEdge, isZupEdge

      for idxs in localRange:colMajorIter() do
         isXloEdge, isXupEdge = false, false
         isYloEdge, isYupEdge = false, false
         isZloEdge, isZupEdge = false, false

         if not hasPeriodicBC then
            if idxs[1] == localRange:lower(1) then isXloEdge = true end
            if idxs[1] == localRange:upper(1) then isXupEdge = true end
            if idxs[2] == localRange:lower(2) then isYloEdge = true end
            if idxs[2] == localRange:upper(2) then isYupEdge = true end
            if idxs[3] == localRange:lower(3) then isZloEdge = true end
            if idxs[3] == localRange:upper(3) then isZupEdge = true end
         end
         
         idxs_LLC[1], idxs_LLC[2], idxs_LLC[3] = idxs[1]-1, idxs[2]-1, idxs[3]
         idxs_LCL[1], idxs_LCL[2], idxs_LCL[3] = idxs[1]-1, idxs[2], idxs[3]-1
         idxs_LCC[1], idxs_LCC[2], idxs_LCC[3] = idxs[1]-1, idxs[2], idxs[3]
         idxs_LCU[1], idxs_LCU[2], idxs_LCU[3] = idxs[1]-1, idxs[2], idxs[3]+1
         idxs_LUC[1], idxs_LUC[2], idxs_LUC[3] = idxs[1]-1, idxs[2]+1, idxs[3]

         idxs_CLL[1], idxs_CLL[2], idxs_CLL[3] = idxs[1], idxs[2]-1, idxs[3]-1
         idxs_CLC[1], idxs_CLC[2], idxs_CLC[3] = idxs[1], idxs[2]-1, idxs[3]
         idxs_CLU[1], idxs_CLU[2], idxs_CLU[3] = idxs[1], idxs[2]-1, idxs[3]+1
         idxs_CCL[1], idxs_CCL[2], idxs_CCL[3] = idxs[1], idxs[2], idxs[3]-1
         --idxs_CCC[1], idxs_CCC[2], idxs_CCC[3] = idxs[1], idxs[2], idxs[3]
         idxs_CCU[1], idxs_CCU[2], idxs_CCU[3] = idxs[1], idxs[2], idxs[3]+1
         idxs_CUL[1], idxs_CUL[2], idxs_CUL[3] = idxs[1], idxs[2]+1, idxs[3]-1
         idxs_CUC[1], idxs_CUC[2], idxs_CUC[3] = idxs[1], idxs[2]+1, idxs[3]
         idxs_CUU[1], idxs_CUU[2], idxs_CUU[3] = idxs[1], idxs[2]+1, idxs[3]+1

         idxs_ULC[1], idxs_ULC[2], idxs_ULC[3] = idxs[1]+1, idxs[2]-1, idxs[3]
         idxs_UCL[1], idxs_UCL[2], idxs_UCL[3] = idxs[1]+1, idxs[2], idxs[3]-1
         idxs_UCC[1], idxs_UCC[2], idxs_UCC[3] = idxs[1]+1, idxs[2], idxs[3]
         idxs_UCU[1], idxs_UCU[2], idxs_UCU[3] = idxs[1]+1, idxs[2], idxs[3]+1
         idxs_UUC[1], idxs_UUC[2], idxs_UUC[3] = idxs[1]+1, idxs[2]+1, idxs[3]

         fStencil7.C = fIn:get(indexer(idxs)):data()
         fStencil7.xL = fIn:get(indexer(idxs_LCC)):data()
         fStencil7.xU = fIn:get(indexer(idxs_UCC)):data()
         fStencil7.yL = fIn:get(indexer(idxs_CLC)):data()
         fStencil7.yU = fIn:get(indexer(idxs_CUC)):data()
         fStencil7.zL = fIn:get(indexer(idxs_CCL)):data()
         fStencil7.zU = fIn:get(indexer(idxs_CCU)):data()

         hStencil7.C = hIn:get(indexer(idxs)):data()
         hStencil7.xL = hIn:get(indexer(idxs_LCC)):data()
         hStencil7.xU = hIn:get(indexer(idxs_UCC)):data()
         hStencil7.yL = hIn:get(indexer(idxs_CLC)):data()
         hStencil7.yU = hIn:get(indexer(idxs_CUC)):data()
         hStencil7.zL = hIn:get(indexer(idxs_CCL)):data()
         hStencil7.zU = hIn:get(indexer(idxs_CCU)):data()


         local f_LLC = fIn:get(indexer(idxs_LLC)):data()
         local f_LCL = fIn:get(indexer(idxs_LCL)):data()
         local f_LCC = fIn:get(indexer(idxs_LCC)):data()
         local f_LCU = fIn:get(indexer(idxs_LCU)):data()
         local f_LUC = fIn:get(indexer(idxs_LUC)):data()

         local f_CLL = fIn:get(indexer(idxs_CLL)):data()
         local f_CLC = fIn:get(indexer(idxs_CLC)):data()
         local f_CLU = fIn:get(indexer(idxs_CLU)):data()
         local f_CCL = fIn:get(indexer(idxs_CCL)):data()
         local f_CCC = fIn:get(indexer(idxs)):data()
         local f_CCU = fIn:get(indexer(idxs_CCU)):data()
         local f_CUL = fIn:get(indexer(idxs_CUL)):data()
         local f_CUC = fIn:get(indexer(idxs_CUC)):data()
         local f_CUU = fIn:get(indexer(idxs_CUU)):data()

         local f_ULC = fIn:get(indexer(idxs_ULC)):data()
         local f_UCL = fIn:get(indexer(idxs_UCL)):data()
         local f_UCC = fIn:get(indexer(idxs_UCC)):data()
         local f_UCU = fIn:get(indexer(idxs_UCU)):data()
         local f_UUC = fIn:get(indexer(idxs_UUC)):data()


         local g_LLC = gIn:get(indexer(idxs_LLC)):data()
         local g_LCL = gIn:get(indexer(idxs_LCL)):data()
         local g_LCC = gIn:get(indexer(idxs_LCC)):data()
         local g_LCU = gIn:get(indexer(idxs_LCU)):data()
         local g_LUC = gIn:get(indexer(idxs_LUC)):data()

         local g_CLL = gIn:get(indexer(idxs_CLL)):data()
         local g_CLC = gIn:get(indexer(idxs_CLC)):data()
         local g_CLU = gIn:get(indexer(idxs_CLU)):data()
         local g_CCL = gIn:get(indexer(idxs_CCL)):data()
         local g_CCC = gIn:get(indexer(idxs)):data()
         local g_CCU = gIn:get(indexer(idxs_CCU)):data()
         local g_CUL = gIn:get(indexer(idxs_CUL)):data()
         local g_CUC = gIn:get(indexer(idxs_CUC)):data()
         local g_CUU = gIn:get(indexer(idxs_CUU)):data()

         local g_ULC = gIn:get(indexer(idxs_ULC)):data()
         local g_UCL = gIn:get(indexer(idxs_UCL)):data()
         local g_UCC = gIn:get(indexer(idxs_UCC)):data()
         local g_UCU = gIn:get(indexer(idxs_UCU)):data()
         local g_UUC = gIn:get(indexer(idxs_UUC)):data()


         --local srcP = src:get(indexer(idxs))

         local f_out = fOut:get(indexer(idxs)):data()

         dragFreq = dragKernelFn(dt, dv:data(),
                                 fStencil7, hStencil7,
                                 isXloEdge, isXupEdge,
                                 isYloEdge, isYupEdge,
                                 isZloEdge, isZupEdge,
                                 f_out)

         diffSurfXLSerFn(dt, dv:data(),
                         f_LCC, f_LLC, f_LUC, f_LCL, f_LCU,
                         f_CCC, f_CLC, f_CUC, f_CCL, f_CCU,
                         g_LCC, g_LLC, g_LUC, g_LCL, g_LCU,
                         g_CCC, g_CLC, g_CUC, g_CCL, g_CCU,
                         isXloEdge,
                         f_out)
         diffSurfXUSerFn(dt, dv:data(),
                         f_CCC, f_CLC, f_CUC, f_CCL, f_CCU,
                         f_UCC, f_ULC, f_UUC, f_UCL, f_UCU,
                         g_CCC, g_CLC, g_CUC, g_CCL, g_CCU,
                         g_UCC, g_ULC, g_UUC, g_UCL, g_UCU,
                         isXupEdge,
                         f_out)
         diffSurfYLSerFn(dt, dv:data(),
                         f_CLC, f_LLC, f_ULC, f_CLL, f_CLU,
                         f_CCC, f_LCC, f_UCC, f_CCL, f_CCU,
                         g_CLC, g_LLC, g_ULC, g_CLL, g_CLU,
                         g_CCC, g_LCC, g_UCC, g_CCL, g_CCU,
                         isYloEdge,
                         f_out)
         diffSurfYUSerFn(dt, dv:data(),
                         f_CCC, f_LCC, f_UCC, f_CCL, f_CCU,
                         f_CUC, f_LUC, f_UUC, f_CUL, f_CUU,
                         g_CCC, g_LCC, g_UCC, g_CCL, g_CCU,
                         g_CUC, g_LUC, g_UUC, g_CUL, g_CUU,
                         isYupEdge,
                         f_out)
         diffSurfZLSerFn(dt, dv:data(),
                         f_CCL, f_LCL, f_UCL, f_CLL, f_CUL,
                         f_CCC, f_LCC, f_UCC, f_CLC, f_CUC,
                         g_CCL, g_LCL, g_UCL, g_CLL, g_CUL,
                         g_CCC, g_LCC, g_UCC, g_CLC, g_CUC,
                         isZloEdge,
                         f_out)
         diffSurfZUSerFn(dt, dv:data(),
                         f_CCC, f_LCC, f_UCC, f_CLC, f_CUC,
                         f_CCU, f_LCU, f_UCU, f_CLU, f_CUU,
                         g_CCC, g_LCC, g_UCC, g_CLC, g_CUC,
                         g_CCU, g_LCU, g_UCU, g_CLU, g_CUU,
                         isZupEdge,
                         f_out)
         diffFreq = diffVolSerFn(dt, dv:data(),
                                 f_CCC, g_CCC,
                                 f_out)

         --cflFreq = math.max(cflFreq, dragFreq, diffFreq)
         cflFreq = math.max(cflFreq)
      end

      if hasConstSource then
         fOut:accumulate(dt, src)
      end

      tmFpo = tmFpo + Time.clock()-tmStart
      return cflFreq
   end


   -------------------------------------------------------------------
   -- The RK3 time stepping ------------------------------------------
   -------------------------------------------------------------------
   local function rk3(tCurr, dt, fIn, fOut)
      local cflFreq
      local localDt
      -- Stage 1 -----------------------------------------------------
      updatePotentials(tCurr, fIn, h, g)
      cflFreq = forwardEuler(dt, fIn, h, g, f1)
      localDt = cflFrac/cflFreq
      if localDt < 0.9*dt then
         return false, localDt
      end
      applyBc(f1)

      -- Stage 2 -----------------------------------------------------
      updatePotentials(tCurr, f1, h, g)
      cflFreq = forwardEuler(dt, f1, h, g, fe)
      localDt = cflFrac/cflFreq
      if localDt < 0.9*dt then
        return false, localDt
      end
      f2:combine(3.0/4.0, fIn, 1.0/4.0, fe)
      applyBc(f2)

      -- Stage 3 -----------------------------------------------------
      updatePotentials(tCurr, f2, h, g)
      cflFreq = forwardEuler(dt, f2, h, g, fe)
      localDt = cflFrac/cflFreq
      if localDt < 0.9*dt then
        return false, localDt
      end
      fOut:combine(1.0/3.0, fIn, 2.0/3.0, fe)
      applyBc(fOut)

      return true, localDt
   end


   -------------------------------------------------------------------
   -- The Main loop --------------------------------------------------
   -------------------------------------------------------------------
   return function ()
      local tCurr = 0.0
      local step = 1
      local dt = tEnd
      local dynDt, status
      if fixedDt then
         dt = fixedDt
      end

      local frameInt = tEnd/nFrames
      local nextFrame = 1
      local isDone = false

      local tmStart = Time.clock()
      while not isDone do
         if (tCurr+dt >= tEnd) then
            isDone = true
            dt = tEnd-tCurr
         end
         print(string.format("Step %d at time %g with dt %g ...",
                             step, tCurr, dt))
         status, dynDt = rk3(tCurr, dt, f, fNew)
         if fixedDt then
            if dynDt < fixedDt then
               print("'fixedDt' is violating the stability condition, exiting")
               break
            end
         else
            dt = dynDt
         end
         if status then
            f:copy(fNew)

            if writeDiagnostics then
               local moms = calcMoms(f)
               momVec:appendData(tCurr+dt, moms)
               updatePotentials(tCurr+dt, f, h, g)
               calcDiag(tCurr+dt, f, h, diagVec)
            end

            tCurr = tCurr+dt
            if tCurr >= nextFrame*frameInt or math.abs(tCurr-nextFrame*frameInt) < 1e-10 then
               writeData(nextFrame, tCurr)
               nextFrame = nextFrame+1
            end
            step = step+1
         else
            isDone = false
            print("dt too big, retaking")
         end
      end

      local tmTotal = Time.clock()-tmStart

      print("")
      print(string.format("Poisson solver took %g secs", tmRosen))
      print(string.format("FPO solver took %g secs", tmFpo))
      print(string.format("Total run-time %g secs", tmTotal))
   end
end
