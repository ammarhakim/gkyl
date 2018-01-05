-- Gkyl ------------------------------------------------------------------------
--
-- Dynamically growing 1D vector
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries
local Adios = require "Io.Adios"
local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"

-- Code from Lua wiki to convert table to comma-seperated-values
-- string.
-- Used to escape "'s by toCSV
local function escapeCSV (s)
  if string.find(s, '[,"]') then
    s = '"' .. string.gsub(s, '"', '""') .. '"'
  end
  return s
end
-- Convert from table to CSV string
local function toCSV (tt)
  local s = ""
  -- ChM 23.02.2014: changed pairs to ipairs assumption is that
  -- fromCSV and toCSV maintain data as ordered array
  for _,p in ipairs(tt) do  
    s = s .. "," .. escapeCSV(p)
  end
  return string.sub(s, 2)      -- remove first comma
end

-- Template to copy from table/vector
local copyTempl = xsys.template [[
return function (src, dest)
| for i = 1, NCOMP do
  dest[${i}] = src[${i}]
| end
end 	    
]]

-- DynVector -------------------------------------------------------------------
--
-- Dynamically growing 1D vector. Used to store diagnostics
--------------------------------------------------------------------------------

local DynVector = Proto()

-- Constructor for DynVector
function DynVector:init(tbl)
   self._numComponents = tbl.numComponents and tbl.numComponents or 1 -- default numComponents=1

   -- We store 1 extra element than requested to allow for 1-based indexing of returned values
   local allocator = Alloc.createAllocator(string.format("double[%d]", self._numComponents+1))
   self._timeMesh = Alloc.Double() -- time-mesh
   self._data = allocator() -- data
   self._tmpData = new(string.format("double[%d]", self._numComponents+1)) -- temp storage for single entry

   -- construct various functions from template representations
   self._copyToTempData = loadstring( copyTempl {NCOMP=self._numComponents} )()

   -- write only from rank-0: create sub-communicator and use that for
   -- writing data (perhaps one needs a user-specified write-rank)
   local ranks = Lin.IntVec(1); ranks[1] = 0
   self._ioComm = Mpi.Split_comm(Mpi.COMM_WORLD, ranks)

   -- allocate space for IO buffer
   self._ioBuff = Alloc.Double()
end

function DynVector:numComponents() return self._numComponents end
function DynVector:lastTime() return self._timeMesh:last() end
function DynVector:lastData() return self._timeMesh:last(), self._data:last() end
function DynVector:timeMesh() return self._timeMesh end
function DynVector:data() return self._data end

function DynVector:appendData(t, v)
   self._timeMesh:push(t)
   self._copyToTempData(v, self._tmpData)
   self._data:push(self._tmpData)
end

function DynVector:removeLast()
   local tm = self._timeMesh:popLast()
   local v = self._data:popLast()
   return tm, v
end


function DynVector:clear()
   self._data:clear()
   self._timeMesh:clear()
end

function DynVector:_copy_from_dynvector(buff)
   local c = 1
   for i = 1, self._data:size() do
      local v = self._data[i]
      for n = 1, self._numComponents do
	 buff[c] = v[n]
	 c = c+1
      end
   end
end

function DynVector:write(outNm, tmStamp)
   local comm = self._ioComm
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   if rank ~= 0 then  -- only run on rank 0 ...
      self:clear() -- ... but clear data on all ranks
      return
   end

   -- setup ADIOS for IO
   Adios.init_noxml(comm[0])

   -- create group and set I/O method
   local grpId = Adios.declare_group("DynVector", "", Adios.flag_no)
   Adios.select_method(grpId, "MPI", "", "")
   
   -- ADIOS expects CSV string to specify data shape
   local localTmSz = toCSV( {self._data:size()} )
   local localDatSz = toCSV( {self._data:size(), self._numComponents} )
   
   -- define data to write
   Adios.define_var(
      grpId, "time", "", Adios.double, "", "", "")
   Adios.define_var(
      grpId, "TimeMesh", "", Adios.double, localTmSz, "", "")
   Adios.define_var(
      grpId, "Data", "", Adios.double, localDatSz, "", "")

   local fullNm = GKYL_OUT_PREFIX .. "_" .. outNm -- concatenate prefix

   -- open file to write out group
   local fd = Adios.open("DynVector", fullNm, "w", comm[0])

   -- write data
   local tmStampBuff = new("double[1]"); tmStampBuff[0] = tmStamp
   Adios.write(fd, "time", tmStampBuff)
   Adios.write(fd, "TimeMesh", self._timeMesh:data())
   -- copy data to IO buffer
   self._ioBuff:expand(self._data:size()*self._numComponents)
   self:_copy_from_dynvector(self._ioBuff)
   Adios.write(fd, "Data", self._ioBuff:data())
   
   Adios.close(fd)
   
   Adios.finalize(rank)
   self:clear() -- clear data for next round of IO
end

return { DynVector = DynVector }
