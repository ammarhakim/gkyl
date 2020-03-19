-- Gkyl ------------------------------------------------------------------------
--
-- Uniform Cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries
local DecompRegionCalc = require "Lib.CartDecomp"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"

local cuda = nil
if GKYL_HAVE_CUDA then
   cuda = require "Cuda.RunTime"
end

ffi.cdef [[
  void cellCenter(int ndim, double* lower, int* currIdx, double* dx, double* xc);

  // Uniform Cartesian grid C representation
  typedef struct {
    int32_t ndim;
    int32_t cells[6];
    double lower[6], upper[6];
    double vol, dx[6];
  } RectCart_t;
]]

local rectCartSz = sizeof("RectCart_t")

local function getDevicePointerToGrid(grid)
   if not GKYL_HAVE_CUDA then return nil end

   -- first copy stuff to a C-struct
   local g = ffi.new("RectCart_t")
   g.ndim = grid:ndim()
   g.vol = grid:cellVolume()
   for i = 1, grid:ndim() do
      g.cells[i-1] = grid:numCells(i)
      g.lower[i-1] = grid:lower(i)
      g.upper[i-1] = grid:upper(i)
      g.dx[i-1] = grid:dx(i)
   end
   
   -- copy stuff to device
   local cuGrid, err = cuda.Malloc(rectCartSz)
   cuda.Memcpy(cuGrid, g, rectCartSz, cuda.MemcpyHostToDevice)
   return cuGrid, err
end

-- Determine local domain index. This is complicated by the fact that
-- when using MPI-SHM the processes that live in shmComm all have the
-- same domain index. Hence, we need to use rank from nodeComm and
-- broadcast it to everyone in shmComm.
local function getSubDomIndex(nodeComm, shmComm)
   local idx = ffi.new("int[1]")
   if Mpi.Is_comm_valid(nodeComm) then
      idx[0] = Mpi.Comm_rank(nodeComm)+1 -- sub-domains are indexed from 1
   end
   -- send from rank 0 of shmComm to everyone else: this works as rank
   -- 0 of shmComm is contained in nodeComm
   Mpi.Bcast(idx, 1, Mpi.INT, 0, shmComm)
   return idx[0]
end

-- RectCart --------------------------------------------------------------------
--
-- A uniform Cartesian grid
--------------------------------------------------------------------------------

-- Uniform Cartesian grid
local RectCart = Proto()

function RectCart:init(tbl)
   local cells = tbl.cells
   self._ndim = #cells
   local lo = tbl.lower and tbl.lower or {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
   local up = tbl.upper and tbl.upper or {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}

   self._periodicDirs = tbl.periodicDirs and tbl.periodicDirs or {}
   self._isDirPeriodic = {false, false, false, false, false, false}
   for _, dir in ipairs(self._periodicDirs) do
      self._isDirPeriodic[dir] = true
   end

   self._lower = Lin.Vec(self._ndim)
   self._upper = Lin.Vec(self._ndim)
   self._dx = Lin.Vec(self._ndim)
   self._numCells = Lin.IntVec(self._ndim)
   self._currIdx = Lin.IntVec(self._ndim)

   self._vol = 1.0
   self._gridVol = 1.0
   for d = 1, #cells do
      self._lower[d], self._upper[d] = lo[d], up[d]
      self._numCells[d] = cells[d]
      self._dx[d] = (up[d]-lo[d])/cells[d]
      self._vol = self._vol*self._dx[d]
      self._gridVol = self._gridVol*(up[d]-lo[d])
   end

   -- compute global range
   local l, u = {}, {}
   for d = 1, #cells do
      l[d], u[d] = 1, cells[d]
   end
   self._globalRange = Range.Range(l, u)   
   self._localRange = Range.Range(l, u)
   self._block = 1 -- block number for use in parallel communications
   self._isShared = false
   
   local decomp = tbl.decomposition and tbl.decomposition or nil  -- decomposition
   if decomp then
      assert(decomp:ndim() == self._ndim,
	     "Decomposition dimensions must be same as grid dimensions!")

      self._isShared = decomp:isShared()
      -- in parallel, we need to adjust local range      
      self._commSet = decomp:commSet()
      self._decomposedRange = decomp:decompose(self._globalRange)
      local subDomIdx = getSubDomIndex(self._commSet.nodeComm, self._commSet.sharedComm)
      self._block = subDomIdx
      local localRange = self._decomposedRange:subDomain(subDomIdx)
      self._localRange:copy(localRange)
      self._cuts = {}
      for i = 1, self._ndim do 
	 assert(decomp:cuts(i) <= self._numCells[i],
		"Cannot have more decomposition cuts than cells in any dimension!")
        self._cuts[i] = decomp:cuts(i) 
      end
   else
      -- create a dummy decomp and use it to set the grid
      local cuts = {}
      for i = 1, self._ndim do cuts[i] = 1 end
      local dec1 = DecompRegionCalc.CartProd { cuts = cuts, useShared = true }
      self._commSet = dec1:commSet()
      self._decomposedRange = dec1:decompose(self._globalRange)
      self._block = 1
      self._cuts = cuts
   end
end

-- member functions
function RectCart:id() return "uniform" end
function RectCart:commSet() return self._commSet end 
function RectCart:isShared() return self._isShared end
function RectCart:subGridId() return self._block end
function RectCart:subGridIdByDim(idx) 
   local invIdxr = self._decomposedRange:cutsInvIndexer()
   invIdxr(self._block, idx)
end
function RectCart:numSharedProcs()
   return Mpi.Comm_size(self._commSet.sharedComm)
end
function RectCart:subGridSharedId()
   return Mpi.Comm_rank(self._commSet.sharedComm)+1 -- we are 1-indexed while MPI is 0-indexed
end
function RectCart:decomposedRange() return self._decomposedRange end
function RectCart:ndim() return self._ndim end
function RectCart:lower(dir) return self._lower[dir] end
function RectCart:mid(dir) return self._lower[dir] + (self._upper[dir]-self._lower[dir])/2 end
function RectCart:upper(dir) return self._upper[dir] end
function RectCart:numCells(dir) return self._numCells[dir] end
function RectCart:totalNumCells()
   total = 1
   for d = 1, self._ndim do
      total = total*self._numCells[d] 
   end
   return total
end
function RectCart:localRange() return self._localRange end
function RectCart:globalRange() return self._globalRange end
function RectCart:isDirPeriodic(dir) return self._isDirPeriodic[dir] end
function RectCart:cuts(dir) return self._cuts[dir] end
function RectCart:setIndex(idxIn)
   if type(idxIn) == "cdata" then
     idxIn:copyInto(self._currIdx)
   else 
     for d = 1, self._ndim do
        self._currIdx[d] = idxIn[d]
     end
   end
end
function RectCart:dx(dir) return self._dx[dir] end
function RectCart:getDx(dxOut) 
   if type(dxOut) == "cdata" then
     self._dx:copyInto(dxOut)
   else 
     for d = 1, self._ndim do
        dxOut[d] = self._dx[d]
     end
   end
end
function RectCart:cellCenterInDir(d)
   return self:lower(d) + (self._currIdx[d]-0.5)*self:dx(d)
end
function RectCart:cellLowerInDir(d)
   return self:lower(d) + (self._currIdx[d]-1)*self:dx(d)
end
function RectCart:cellUpperInDir(d)
   return self:lower(d) + (self._currIdx[d])*self:dx(d)
end
function RectCart:cellCenter(xc)
   ffiC.cellCenter(self._ndim, self._lower:data(), self._currIdx:data(), self._dx:data(), xc:data())
end
function RectCart:cellVolume() return self._vol end
function RectCart:gridVolume() return self._gridVol end

function RectCart:write(fName)
   -- nothing to write
end

local _numMetricElems = { 1, 3, 6 }
-- size of vector to store metric elements
function RectCart:numMetricElems()
   return _numMetricElems[self:ndim()]
end

function RectCart:calcMetric(xc, gOut)
   if self._ndim == 1 then
      gOut[1] = 1.0 -- g_11
   elseif self._ndim == 2 then
      -- g_11, g_12, g_22
      gOut[1], gOut[2], gOut[3] = 1.0, 0.0, 1.0
   elseif self._ndim == 3 then
      -- g_11, g_12, g_13, g_22, g_23, g_33
      gOut[1], gOut[2], gOut[3], gOut[4], gOut[5], gOut[6] = 1.0, 0.0, 0.0, 1.0, 0.0, 1.0
   else
      assert(false, "RectCart:calcMetric does not support more than 3 dimensions!")
   end
end

function RectCart:calcContraMetric(xc, gOut)
   if self._ndim == 1 then
      gOut[1] = 1.0 -- g_11
   elseif self._ndim == 2 then
      -- g^11, g^12, g^22
      gOut[1], gOut[2], gOut[3] = 1.0, 0.0, 1.0
   elseif self._ndim == 3 then
      -- g^11, g^12, g^13, g^22, g^23, g^33
      gOut[1], gOut[2], gOut[3], gOut[4], gOut[5], gOut[6] = 1.0, 0.0, 0.0, 1.0, 0.0, 1.0
   else
      assert(false, "RectCart:calcContraMetric does not support more than 3 dimensions!")
   end
end

function RectCart:calcJacobian(xc)
   return 1.0
end

function RectCart:copyHostToDevice()
   return getDevicePointerToGrid(self)
end

return RectCart
