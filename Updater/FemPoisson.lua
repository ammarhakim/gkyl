-- Gkyl ------------------------------------------------------------------------
--
-- Wraps FemParPoisson and FemPerpPoisson for common cases in GK
-- To be generalized further...
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto          = require "Lib.Proto"
local UpdaterBase    = require "Updater.Base"
local FemParPoisson  = require "Updater.FemParPoisson"
local FemPerpPoisson = require "Updater.FemPerpPoisson"
local DataStruct     = require "DataStruct"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local CartFieldBinOp = require "Updater.CartFieldBinOp"
local xsys           = require "xsys"
local SkinGhostAvg   = require "Updater.SkinGhostAvg"
local Grid = require("Grid")

-- FEM Poisson solver updater object
local FemPoisson = Proto(UpdaterBase)

function FemPoisson:init(tbl)
   FemPoisson.super.init(self, tbl)

   local function contains(table, element)
      for _, value in pairs(table) do
         if value == element then
            return true
         end
      end
      return false
   end

   self.xLCFS=tbl.xLCFS
   self.grid  = tbl.onGrid
   self.ndim  = self.grid:ndim()
   self.basis = tbl.basis
   local ndim = self.ndim
   self.zContinuous = xsys.pickBool(tbl.zContinuous, true)

   -- Boundary conditions.
   self.bcLower = tbl.bcLower
   self.bcUpper = tbl.bcUpper

   -- Metric coefficients.
   self.gxx = tbl.gxx
   self.gxy = tbl.gxy
   self.gyy = tbl.gyy

   self.smooth = xsys.pickBool(tbl.smooth, false)

   self.slvr = nil
   self._nx = self.grid:numCells(1)
   self._ny = self.grid:numCells(2)
   self._nz = self.grid:numCells(3)
   self.axisymmetric = nil
   self.twistShift=nil
   if self.xLCFS then
      if self._ny==1 then
         self.axisymmetric=true
      else
         self.twistShift=true
         self.onZGhosts=true
      end
   end

   if ndim == 1 then
      self.zSmoother = function(tCurr, src) end
      self.slvr = FemParPoisson {
         onGrid  = self.grid,
         basis   = self.basis,
         bcLower = self.bcLower,
         bcUpper = self.bcUpper,
         smooth  = self.smooth
      }
      -- Set up weak division operator for special case when solve is algebraic.
      self.weakDivide = CartFieldBinOp {
         onGrid    = self.grid,
         weakBasis = self.basis,
         operation = "Divide",
         onGhosts  = true,
      }
   elseif ndim == 2 then
      self.zSmoother = function(tCurr, src) end
      self.slvr = FemPerpPoisson {
        onGrid       = self.grid,
        basis        = self.basis,
        bcLower      = self.bcLower,
        bcUpper      = self.bcUpper,
        constStiff   = self.constStiff,
        gxx          = self.gxx,
        gxy          = self.gxy,
        gyy          = self.gyy,
        smooth       = self.smooth
      }
   elseif ndim == 3 then
      self.slvr = FemPerpPoisson {
        xLCFS        = self.xLCFS,
        onGrid       = self.grid,
        basis        = self.basis,
        bcLower      = self.bcLower,
        bcUpper      = self.bcUpper,
        zContinuous  = self.zContinuous,
        constStiff   = self.constStiff,
        gxx          = self.gxx,
        gxy          = self.gxy,
        gyy          = self.gyy,
        smooth       = self.smooth,
        onZGhosts    = self.onZGhosts
      }
      if self.zContinuous then
         if self.xLCFS then
            self.zDiscontToCont = FemParPoisson {
               onGrid  = self.grid,
               basis   = self.basis,
               bcLower = {{T= "N",V=0.0}},
               bcUpper = {{T= "N",V=0.0}},
               smooth  = true,
            }
            local function createField(grid, basis, ghostCells, vComp, periodicSync)
               vComp = vComp or 1
               local fld = DataStruct.Field {
                  onGrid           = grid,
                  numComponents    = basis:numBasis()*vComp,
                  ghost            = ghostCells,
                  metaData         = {polyOrder = basis:polyOrder(),
                                      basisType = basis:id()},
                  syncPeriodicDirs = periodicSync,
               }
               fld:clear(0.0)
               return fld
            end
            self.srcSOL = createField(self.grid,self.basis,{1,1})
            self:setSplitRanges(self.srcSOL)
            if self.axisymmetric then
               self.zDiscontToContCore = FemParPoisson {
                  onGrid  = self.grid,
                  basis   = self.basis,
                  bcLower = {{T = "P",V=0.0}},
                  bcUpper = {{T = "P",V=0.0}},
                  smooth  = true,
               }
               self.zSmoother = function(tCurr,src) return self:axisymmetricSmoother(tCurr,src) end
            end
            if self.twistShift then
               self.zDiscontToContCore = FemParPoisson {
                  onGrid  = self.grid,
                  basis   = self.basis,
                  bcLower = {{T = "N",V=0.0}},
                  bcUpper = {{T = "N",V=0.0}},
                  smooth  = true,
               }
               self.zSmoother = function(tCurr,src) return self:twistShiftSmoother(tCurr,src) end
               self.skinGhostAvgLower = SkinGhostAvg{
                  grid = self.grid,
                  basis = self.basis,
                  edge='lower',
                  bcDir = 3,
                  advanceArgs = {self.srcSOL},
               }
               self.skinGhostAvgUpper= SkinGhostAvg{
                  grid = self.grid,
                  basis = self.basis,
                  edge='upper',
                  bcDir = 3,
                  advanceArgs = {self.srcSOL},
               }
            end
         else
            self.zSmoother = function(tCurr,src) return self:regularSmoother(tCurr,src) end
            self.zDiscontToCont = FemParPoisson {
               onGrid  = self.grid,
               basis   = self.basis,
               bcLower = {tbl.bcLower[3]},
               bcUpper = {tbl.bcUpper[3]},
               smooth  = true,
            }
         end
      end
   else
      assert(false, "Updater.FemPoisson: Requires ndim<=3")
   end

   -- Option to write development-related diagnostics.
   self.verbose = xsys.pickBool(tbl.verbose, false)

   return self
end


function FemPoisson:regularSmoother(tCurr,src)
   self.zDiscontToCont:advance(tCurr, {src}, {src})
end

function FemPoisson:setSplitRanges(sampleFld)
   local gridRange = self.grid:globalRange()
   local xLCFS=self.xLCFS
   local coordLCFS = {xLCFS-1.e-7}
   idxLCFS    = {-9}
   local xGridIngr = self.grid:childGrid({1})
   local xGrid = Grid.RectCart {
      lower = xGridIngr.lower,  periodicDirs  = xGridIngr.periodicDirs,
      upper = xGridIngr.upper,  decomposition = xGridIngr.decomposition,
      cells = xGridIngr.cells,
   }
   xGrid:findCell(coordLCFS, idxLCFS)
   self.globalCore = gridRange:shorten(1, idxLCFS[1])
   self.globalSOL = gridRange:shortenFromBelow(1, self.grid:numCells(1)-idxLCFS[1]+1)

   local localRange = sampleFld:localRange()
   self.localSOLRange = sampleFld:localRange():intersect(self.globalSOL)
   self.localCoreRange = sampleFld:localRange():intersect(self.globalCore)
end

function FemPoisson:axisymmetricSmoother(tCurr,src)
   self.srcSOL:copy(src)
   self.zDiscontToCont:advance(tCurr, {self.srcSOL}, {self.srcSOL})
   self.zDiscontToContCore:advance(tCurr, {src}, {src})
   src:copyRange(self.localSOLRange, self.srcSOL)
end

function FemPoisson:twistShiftSmoother(tCurr,src)
   self.srcSOL:copy(src)
   self.zDiscontToCont:advance(tCurr, {self.srcSOL}, {self.srcSOL})
   self.skinGhostAvgLower:advance(src)
   self.skinGhostAvgUpper:advance(src)
   self.zDiscontToContCore:advance(tCurr, {src}, {src})
   src:copyRange(self.localSOLRange, self.srcSOL)
end


function FemPoisson:assemble(tCurr, inFld, outFld)
   -- Begin assembling the source vector and, if needed, the stiffness matrix.
   self.zSmoother(tCurr, inFld[1])
   self.slvr:assemble(tCurr, inFld, outFld)
end

function FemPoisson:assemble(tCurr, inFld, outFld)
   -- Begin assembling the source vector and, if needed, the stiffness matrix.
   self.slvr:assemble(tCurr, inFld, outFld)
end

function FemPoisson:solve(tCurr, inFld, outFld)
   -- Assuming the right-side vector (and if needed the stiffness matrix)
   -- has been assembled, this solves the linear problem.
   -- If the assembly initiated an MPI non-blocking reduce, this waits for it.
   self.slvr:solve(tCurr, inFld, outFld)
   self.zSmoother(tCurr, outFld[1])
end

function FemPoisson:_advance(tCurr, inFld, outFld)
   -- Advance method. This assembles the right-side source vector and, if needed,
   -- the stiffness matrix. Then it solves the linear problem.
   self.zSmoother(tCurr, inFld[1])
   if self.ndim == 1 and not self.zContinuous and self.slvr._hasLaplacian == false then
      -- Special case where solve is just algebraic.
      local src = inFld[1]
      local sol = outFld[1]

      self.weakDivide:advance(0, {self.slvr:getModifierWeight(), src}, {sol})
   else
      self.slvr:advance(tCurr, inFld, outFld)
   end
   self.zSmoother(tCurr, outFld[1])
end

function FemPoisson:setLaplacianWeight(weight)
   self.slvr:setLaplacianWeight(weight)
end
function FemPoisson:setModifierWeight(weight)
   self.slvr:setModifierWeight(weight)
end

function FemPoisson:getLaplacianWeight()
   return self.slvr:getLaplacianWeight()
end
function FemPoisson:getModifierWeight()
   return self.slvr:getModifierWeight()
end

function FemPoisson:printDevDiagnostics()
   if self.verbose then self.slvr:printDevDiagnostics() end
end

return FemPoisson
