local Basis = require "Basis"
local DataStruct = require "DataStruct"
local Grid = require "Grid"
local Updater = require "Updater"

return function(tbl)
   local polyOrder = tbl.polyOrder
   local cflFrac = tbl.cflFrac and tbl.cflFrac or 1.0
   local cfl = cflFrac*0.5/(2*polyOrder+1)
   local tEnd = tbl.tEnd
   local nFrames = tbl.nFrames
   local a = tbl.a
   local cells = tbl.cells
   local lower = tbl.lower
   local upper = tbl.upper
   local periodicDirs = {1}

   -- Select scheme:
   -- 1: Naive (no conection between cells)
   -- 2: Central flux
   -- 3: Upwinding flux
   -- 4: 1-cell recovery (no integration by parts)
   -- 5: 2-cell recovery (standard recovery DG with IBP)
   -- 6: Biased 3-cell recovery
   local schemeId = tbl.schemeId

  
   ----------------------
   -- Grids and Fields --
   ----------------------
   local grid = Grid.RectCart {
      lower = {lower[1]},
      upper = {upper[1]},
      cells = {cells[1]},
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
         ghost = {2, 2},
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
   local fNew = getField()
   
   --------------------
   -- Initialization --
   --------------------
   if tbl.exactInit then 
      local localRange = f:localRange()
      local indexer = f:genIndexer()
      for idxs in localRange:colMajorIter() do
         local fPtr = f:get(indexer(idxs))
         for k = 1, basis:numBasis() do
            fPtr[k] = tbl.exactInit[idxs[1]][k]
         end
      end
   else
      local projectUpd = Updater.ProjectOnBasis {
         onGrid = grid,
         basis = basis,
         evaluate = tbl.init,
         numQuad = 8,
         onGhosts = false,
      }
      projectUpd:advance(0.0, {}, {f})
   end
   f:sync()
   f:write("f_0.bp", 0, 0)
   
   local function forwardEuler(dt, fIn, fOut)
      local dv = grid:dx(1)
      local localRange = fIn:localRange()
      local indexer = fIn:genIndexer()
      local idxsLL, idxsL, idxsR = {}, {}, {}
      
      for idxs in localRange:colMajorIter() do
         idxsLL[1] = idxs[1]-2
         idxsL[1] = idxs[1]-1
         idxsR[1] = idxs[1]+1
         
         local fC = fIn:get(indexer(idxs))
         local fLL = fIn:get(indexer(idxsLL))
         local fL = fIn:get(indexer(idxsL))
         local fR = fIn:get(indexer(idxsR))
         local fO = fOut:get(indexer(idxs))
         
         if schemeId == 1 then  -- Naive
            fO[1] = fC[1] + -(3.464101615137754*fC[2]*a*dt)/dv
            fO[2] = fC[2] + 0.0
         elseif schemeId == 2 then  -- Central
            fO[1] = fC[1] + (0.5*(1.732050807568877*fR[2]+1.732050807568877*fL[2]-3.464101615137754*fC[2]-1.0*fR[1]+fL[1])*a*dt)/dv
            fO[2] = fC[2] + (0.5*(3.0*fR[2]-3.0*fL[2]-1.732050807568877*fR[1]-1.732050807568877*fL[1]+3.464101615137754*fC[1])*a*dt)/dv
         elseif schemeId == 3 then  -- Upwinding
            fO[1] = fC[1] + (1.732050807568877*fL[2]-1.732050807568877*fC[2]+fL[1]-1.0*fC[1])*a*dt/dv
            fO[2] = fC[2] + -(3.0*fL[2]+3.0*fC[2]+1.732050807568877*fL[1]-1.732050807568877*fC[1])*a*dt/dv
         elseif schemeId == 4 then  -- Recovery 1-cell ("volume" update)
            fO[1] = fC[1] + (0.1666666666666667*(3.464101615137754*fR[2]+3.464101615137754*fL[2]-6.928203230275509*fC[2]-3.0*fR[1]+3.0*fL[1])*a*dt)/dv
            fO[2] = fC[2] + (0.2886751345948129*(3.464101615137754*fR[2]-3.464101615137754*fL[2]-3.0*fR[1]-3.0*fL[1]+6.0*fC[1])*a*dt)/dv
         elseif schemeId == 5 then  -- Recovery 2-cell (volume + surface updates)
            if polyOrder == 1 then
               fO[1] = fC[1] + dt*((0.5773502691896258*fR[2]*a)/dv+(0.5773502691896258*fL[2]*a)/dv-(1.154700538379252*fC[2]*a)/dv-(0.5*fR[1]*a)/dv+(0.5*fL[1]*a)/dv)
               fO[2] = fC[2] + dt*((fR[2]*a)/dv-(1.0*fL[2]*a)/dv-(0.8660254037844386*fR[1]*a)/dv-(0.8660254037844386*fL[1]*a)/dv+(1.732050807568877*fC[1]*a)/dv)
            elseif polyOrder == 2 then
               fO[1] = fC[1] + dt*((-(0.489139870078079*fR[3]*a)/dv)+(0.489139870078079*fL[3]*a)/dv+(0.7036456405748563*fR[2]*a)/dv+(0.7036456405748563*fL[2]*a)/dv-(1.407291281149713*fC[2]*a)/dv-(0.5*fR[1]*a)/dv+(0.5*fL[1]*a)/dv)
               fO[2] = fC[2] + dt*((-(0.8472151069828725*fR[3]*a)/dv)-(0.8472151069828725*fL[3]*a)/dv-(1.694430213965745*fC[3]*a)/dv+(1.21875*fR[2]*a)/dv-(1.21875*fL[2]*a)/dv-(0.8660254037844386*fR[1]*a)/dv-(0.8660254037844386*fL[1]*a)/dv+(1.732050807568877*fC[1]*a)/dv)
               fO[3] = fC[3] + dt*((-(1.09375*fR[3]*a)/dv)+(1.09375*fL[3]*a)/dv+(1.573399484396763*fR[2]*a)/dv+(1.573399484396763*fL[2]*a)/dv+(4.599167723621307*fC[2]*a)/dv-(1.118033988749895*fR[1]*a)/dv+(1.118033988749895*fL[1]*a)/dv)
            end
         elseif schemeId == 6 then  -- Recovery 3-cell (volume + surface updates)
            fO[1] = fC[1] + dt*((0.2993668062464727*fR[2]*a)/dv+(0.1336458956457467*fLL[2]*a)/dv+(1.15470053837925*fL[2]*a)/dv-(1.587713240271471*fC[2]*a)/dv-(0.2962962962962963*fR[1]*a)/dv+(0.1203703703703704*fLL[1]*a)/dv+(0.462962962962963*fL[1]*a)/dv-(0.287037037037037*fC[1]*a)/dv)
            fO[2] = fC[2] + dt*((0.5185185185185185*fR[2]*a)/dv-(0.2314814814814815*fLL[2]*a)/dv-(2.462962962962963*fL[2]*a)/dv-(1.712962962962963*fC[2]*a)/dv-(0.5132002392796675*fR[1]*a)/dv-(0.2084875972073649*fLL[1]*a)/dv-(1.21885056828921*fL[1]*a)/dv+(1.940538404776242*fC[1]*a)/dv)
         end
      end
   end
   
   local function rk3(dt, fIn, fOut)
      -- Stage 1
      forwardEuler(dt, fIn, f1)
      f1:sync()
      
      -- Stage 2
      forwardEuler(dt, f1, fe)
      f2:combine(3.0/4.0, fIn, 1.0/4.0, fe)
      f2:sync()
      
      -- Stage 3
      forwardEuler(dt, f2, fe)
      fOut:combine(1.0/3.0, fIn, 2.0/3.0, fe)
      fOut:sync()
   end
   
   -- run simulation with RK3
   return function ()
      local tCurr = 0.0
      local step = 1
      local dt = 1e-5--cfl*grid:dx(1)
      
      local frameInt = tEnd/nFrames
      local nextFrame = 1
      local isDone = false
   
      while not isDone do
         if (tCurr+dt >= tEnd) then
            isDone = true
            dt = tEnd-tCurr
         end
         print(string.format("Step %d at time %g with dt %g ...", step, tCurr, dt))
         rk3(dt, f, fNew)
         f:copy(fNew)
         
         tCurr = tCurr+dt
         if tCurr >= nextFrame*frameInt or math.abs(tCurr-nextFrame*frameInt) < 1e-10 then
            f:write(string.format("f_%d.bp", nextFrame), tCurr, nextFrame)
            nextFrame = nextFrame+1
         end
         step = step+1
      end
   end
end
