-- Gkyl ------------------------------------------------------------------------
--
-- Updater to solve Poisson equations on a logically rectangular grid
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require("ffi")
local Proto = require("Lib.Proto")
local UpdaterBase = require("Updater.Base")
local Lin = require("Lib.Linalg")

-- Updater to solve Poisson equations on a logically rectangular grid
local MappedPoisson = Proto(UpdaterBase)

--set up C/Lua interface
ffi.cdef[[
  typedef struct MapPoisson MapPoisson;
  //constructor
  void *new_MapPoisson(int nx_, int ny_, double xl_, double yl_, double xu_, double yu_,
int xbctypel_, int ybctypel_, int xbctypeu_, int ybctypeu_, double xbcu_,
double xbcl_, double ybcu_, double ybcl_);
  //metric function
  void setMetricFuncPointer_MapPoisson(MapPoisson *d, void (*gfunc)(double *xcl, double *g));
  //mapc2p for C++ side
  void setMapcpp_MapPoisson(MapPoisson *d, void (*mapcpp)(double xc, double yc, double *myxp));
  //solver bits
  void wrap_factorize(MapPoisson *d);
  void wrap_phisolve(MapPoisson *d);
  //wrap source and solution solvers
  void wrap_getSrcvalatIJ(MapPoisson *d, int i, int j, double sitrij);
  double wrap_getSolvalatIJ(MapPoisson *d, int i, int j);
]]

--initialization function
function MappedPoisson:init(tbl)

   --intialize table objects with self
   MappedPoisson.super.init(self, tbl)
   self._grid = assert(tbl.onGrid, "Error: must provide grid object with onGrid")
   self._perDirs = tbl.periodicDirs
   self._left = tbl.bcLeft
   self._right = tbl.bcRight
   self._top = tbl.bcTop
   self._bottom = tbl.bcBottom

   --set up metric function
   local myXc, myG = Lin.Vec(2), Lin.Vec(3)
   local function gfunc(xc, g)
      myXc[1], myXc[2] = xc[1], xc[2]
      self._grid:calcMetric(myXc, myG)
      g[1], g[2], g[3] = myG[1], myG[2], myG[3]
   end

   --C++ callable mapc2p function for ruled surface metric calc
   local myxc = Lin.Vec(2)
   local function mapcpp(xc, yc, myxp)
     myxc[1], myxc[2] = xc, yc
     local xp, yp, zp = self._grid:mapc2p(myxc) --works
     --local xp, yp = self._grid:mapc2p(myxc)
     myxp[1], myxp[2], myxp[3] = xp, yp, zp
   end

   --global idx/y so that accessible to solver
   local idx = self._grid:numCells(1)
   local idy = self._grid:numCells(2)

   local xbctu, xbctl, ybctu, ybctl
   --reformat periodic boundary condition inputs
   if (self._perDirs ~= nil) then
      if (self._perDirs[1] == 1 or self._perDirs[2] == 1) then
	       xbctu = 2
	       xbctl = 2
      end
      if (self._perDirs[1] == 2 or self._perDirs[2] == 2) then
	       ybctu = 2
	       ybctl = 2
      end
      if (self._perDirs[1] == 3) then --special periodic bc for calhoun _grid
        xbctu = 2
        xbctl = 2 --long dimension should be x so that x periodic BCs behave normally
        ybctu = 3
        ybctl = 3 --y bcs need to be modified
      end
   end

   local xbcl, xbcu
   --reformat dirichlet boundary condition inputs
   if (self._left ~= nil) then --left side (x)
      if (self._left.T == "D") then
	 xbctl = 0
	 xbcl = self._left.V
      elseif (self._left.T == "N") then
	 xbctl = 1
	 xbcl = self._left.V
      end
   else
      xbcl = 0
   end

   if (self._right ~= nil) then --right side (x)
      if (self._right.T =="D") then
	 xbctu = 0
	 xbcu = self._right.V
      elseif (self._right.T =="N") then
	 xbctu = 1
	 xbcu = self._right.V
      end
   else
      xbcu = 0
   end

   local ybcl, ybcu
   if (self._bottom ~= nil) then --bottom side (y)
      if (self._bottom.T == "D") then
	 ybctl = 0
	 ybcl = self._bottom.V
      elseif (self._bottom.T == "N") then
	 ybctl = 1
	 ybcl = self._bottom.V
      end
   else
      ybcl = 0
   end

   if (self._top ~= nil) then --top side (y)
      if (self._top.T =="D") then
	 ybctu = 0
	 ybcu = self._top.V
      elseif (self._top.T == "N") then
	 ybctu = 1
	 ybcu = self._top.V
      end
   else
      ybcu = 0
   end

   --check to make sure all bcs known after setting values
   assert(xbctl ~= nil, "Error: Left BC not specified")
   assert(xbctu ~= nil, "Error: Right BC not specified")
   assert(ybctl ~= nil, "Error: Bottom BC not specified")
   assert(ybctl ~= nil, "Error: Top BC not specified")

   --initialize function pointer
   self.point = ffi.C.new_MapPoisson(idx, idy, self._grid:lower(1), self._grid:lower(2),
				     self._grid:upper(1), self._grid:upper(2), xbctl, ybctl, xbctu, ybctu, xbcu,
				     xbcl, ybcu, ybcl);
   --set up metric function for c
   ffi.C.setMetricFuncPointer_MapPoisson(self.point, gfunc)
   ffi.C.setMapcpp_MapPoisson(self.point, mapcpp)
   --set up factorization wth pointer
   ffi.C.wrap_factorize(self.point)
end


--solver function
function MappedPoisson:advance(tCurr, dt, inFld, outFld)

   --set up input objects
   local grid = self._grid
   local source = inFld[1]
   local phipot = outFld[1]

   -- set up RHS
   local localRange = source:localRange()
   local indexer = source:indexer()
   for j = localRange:lower(2), localRange:upper(2) do
      for i = localRange:lower(1), localRange:upper(1) do
	 grid:setIndex({i,j})
	 local sitr = source:get(indexer(i,j))
	 ffi.C.wrap_getSrcvalatIJ(self.point, i-1, j-1, sitr[1])
      end
   end

   --solve the problem
   ffi.C.wrap_phisolve(self.point)

   --copy solution into output field
   localRange = phipot:localRange()
   for j = localRange:lower(2), localRange:upper(2) do
      for i = localRange:lower(1), localRange:upper(1) do
	 grid:setIndex({i,j})
	 local pitr = phipot:get(indexer(i,j))
	 pitr[1] = ffi.C.wrap_getSolvalatIJ(self.point, i-1, j-1)
      end
   end

   return true, GKYL_MAX_DOUBLE
end

return MappedPoisson
--------------------------------------------------------------------------------------
