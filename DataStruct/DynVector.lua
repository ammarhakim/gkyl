-- Gkyl ------------------------------------------------------------------------
--
-- Dynamically growing 1D vector.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries.
local Adios       = require "Io.Adios"
local AdiosReader = require "Io.AdiosReader"
local Alloc       = require "Lib.Alloc"
local Lin         = require "Lib.Linalg"
local Mpi         = require "Comm.Mpi"
local Proto       = require "Lib.Proto"

-- Code from Lua wiki to convert table to comma-seperated-values string.
-- Used to escape "'s by toCSV.
local function escapeCSV (s)
  if string.find(s, '[,"]') then
    s = '"' .. string.gsub(s, '"', '""') .. '"'
  end
  return s
end
-- Convert from table to CSV string.
local function toCSV (tt)
  local s = ""
  -- ChM 23.02.2014: changed pairs to ipairs assumption is that
  -- fromCSV and toCSV maintain data as ordered array.
  for _,p in ipairs(tt) do  
    s = s .. "," .. escapeCSV(p)
  end
  return string.sub(s, 2)      -- Remove first comma.
end

-- Template to copy from table/vector.
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

-- Constructor for DynVector.
function DynVector:init(tbl)
   self._numComponents = tbl.numComponents and tbl.numComponents or 1

   -- We store 1 extra element than requested to allow for 1-based
   -- indexing of returned values.
   local allocator = Alloc.createAllocator(
      string.format("double[%d]", self._numComponents+1))
   self._timeMesh = Alloc.Double()
   self._data     = allocator()
   -- Temp storage for single entry.
   self._tmpData  = new(string.format("double[%d]", self._numComponents+1))
   self._tmpTable = {}
   for i = 1, self._numComponents do self._tmpTable[i] = nil end

   -- Construct various functions from template representations.
   self._copyToTempData = loadstring( copyTempl {NCOMP=self._numComponents} )()

   -- Write only from rank-0: create sub-communicator and use that for
   -- writing data (perhaps one needs a user-specified write-rank).
   local ranks  = Lin.IntVec(1); ranks[1] = 0
   self._ioComm = Mpi.Split_comm(Mpi.COMM_WORLD, ranks)

   -- Allocate space for IO buffer.
   self._ioBuff = Alloc.Double()
   self.frNum = -1

   -- If we have meta-data to write out, store it.
   if GKYL_EMBED_INP then
      self._metaData = {
         -- Write out input file contents (encoded as base64 string).
         ["inputfile"] = {
            value = GKYL_INP_FILE_CONTENTS, vType = "string"
         }
      }
   else
      -- Write some dummy text otherwise.
      self._metaData = { ["inputfile"] = { value = "inputfile", vType = "string" } }
   end
   if tbl.metaData then
      -- Store value and its type for each piece of data.
      for k,v in pairs(tbl.metaData) do
         if type(v) == "number" then
            -- Check if this is an integer or float.
            if math.floor(math.abs(v)) == math.abs(v) then
               self._metaData[k] = { value = new("int[1]", v), vType = "integer", }
            else
               self._metaData[k] = { value = new("double[1]", v), vType = "double", }
            end
         elseif type(v) == "string" then
            self._metaData[k] = { value = v, vType = "string" }
         elseif type(v) == "table" then
            assert(type(v[1])=="number", "Io.AdiosCartFieldIo: Metadata table must have elements must be numbers.")
            isInt = (math.floor(math.abs(v[1])) == math.abs(v[1]))
            for _, val in pairs(v) do
               assert(isInt == (math.floor(math.abs(val)) == math.abs(val)), "Io.AdiosCartFieldIo: Metadata table must have elements of the same type (int or double).")
            end
            if isInt then
               self._metaData[k] = { value = new("int[?]", #v), vType = "table", numElements = #v, elementType = Adios.integer }
            else
               self._metaData[k] = { value = new("double[?]", #v), vType = "table", numElements = #v, elementType = Adios.double }
            end
            for i = 1, #v do self._metaData[k].value[i-1] = v[i] end
         end
      end
   end

   -- Used to save the last data, in case we need it later.
   self.tLast, self.dLast = -1., Lin.Vec(self._numComponents)
   self.flushed = false

   self._isFirst = true
end

function DynVector:appendData(t, v)
   self._timeMesh:push(t)
   self._copyToTempData(v, self._tmpData)
   self._data:push(self._tmpData)
end

function DynVector:appendLast(dynV)
   local tm, lv = dynV:lastData()
   self:appendData(tm, lv)
end

function DynVector:assignLastTime(tm)
   self._timeMesh:assignLast(tm)
end

function DynVector:assignLastVal(fact, val)
   for i = 1, self._numComponents do self._tmpData[i] = val[i]*fact end
   self._data:assignLast(self._tmpData)
end

function DynVector:assignLast(fact, dynV)
   local _, lv = dynV:lastData()
   for i = 1, self._numComponents do self._tmpData[i] = lv[i]*fact end
   self._data:assignLast(self._tmpData)
end

function DynVector:accumulateLastOne(fact, dynV)
   local _, lv = dynV:lastData()
   local selflv
   if self:size() == 0 and self.flushed then
      selflv = self.dLast
   else
      selflv = self._data:last()
   end
   for i = 1, self._numComponents do self._tmpData[i] = selflv[i]+lv[i]*fact end
   self._data:assignLast(self._tmpData)
end

function DynVector:accumulateLast(c1, dynV1, ...)
   local args = {...}                  -- Package up rest of args as table.
   self:accumulateLastOne(c1, dynV1)   -- Accumulate first factor*vec:last() pair.
   for i = 1, #args/2 do               -- Accumulate rest of the factor*vec:last() pairs.
      self:accumulateLastOne(args[2*i-1], args[2*i])
   end
end

function DynVector:combineLast(c1, dynV1, ...)
   local args = {...}           -- Package up rest of args as table.
   self:assignLast(c1, dynV1)   -- Assign first factor*vec:last() pair.
   for i = 1, #args/2 do        -- Accumulate rest of the factor*vec:last() pairs.
      self:accumulateLastOne(args[2*i-1], args[2*i])
   end
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

function DynVector:_copy_to_dynvector(buff)
   local c = 1
   for i = 1, self._data:size() do
      local v = self._data[i]
      for n = 1, self._numComponents do
	 v[n] = buff[c]
	 c = c+1
      end
   end
end

function DynVector:clear()
   self._data:clear()
   self._timeMesh:clear()
end

function DynVector:copy(dynV)
   -- Copy over time-mesh and data. Assumes number of components is the same.
   local dataIn = dynV:data()
   local timeMeshIn = dynV:timeMesh()
   local nValIn = dataIn:size()

   self._timeMesh:expand(nValIn-(self:size()-1))
   for i = 1, self._timeMesh:size() do
      self._timeMesh[i] = timeMeshIn[i]
   end

   self._data:expand(nValIn-(self:size()-1))
   for i = 1, self._data:size() do
      for n = 1, self._numComponents do
	 self._data[i][n] = dataIn[i][n]
      end
   end
end

function DynVector:copyLast(dynV)
   local tm, lv = dynV:lastData()
   self._timeMesh:assignLast(tm)
   self._data:assignLast(lv)
end

function DynVector:data() return self._data end
function DynVector:lastTime() return self._timeMesh:last() end

function DynVector:lastTime()
   if self:size() == 0 and self.flushed then 
      return self.tLast
   else
      return self._timeMesh:last()
   end
end

function DynVector:lastData()
   if self:size() == 0 and self.flushed then 
      return self.tLast, self.dLast
   else
      return self._timeMesh:last(), self._data:last() 
   end
end

function DynVector:numComponents() return self._numComponents end

-- Returns time-stamp and frame number.
function DynVector:read(fName)
   local comm = Mpi.COMM_WORLD -- Need to read from all processors.

   local fullNm = GKYL_OUT_PREFIX .. "_" .. fName
   local reader = AdiosReader.Reader(fullNm, comm)

   local timeMesh, data
   if reader:hasVar("TimeMesh") then
      timeMesh = reader:getVar("TimeMesh"):read()
      data     = reader:getVar("Data"):read()
   elseif reader:hasVar("TimeMesh0") then
      timeMesh = reader:getVar("TimeMesh0"):read()
      data     = reader:getVar("Data0"):read()
      varCnt   = 1
      while reader:hasVar("TimeMesh"..varCnt) do
         local timeMeshN = reader:getVar("TimeMesh"..varCnt):read()
         local dataN     = reader:getVar("Data"..varCnt):read()
         for i = 1, timeMeshN:size() do timeMesh:push(timeMeshN[i]) end
         for i = 1, dataN:size() do data:push(dataN[i]) end
         varCnt = varCnt + 1
      end
   end

   -- Copy over time-mesh and data.
   local nVal = data:size()/self._numComponents
   self._timeMesh:expand(nVal)
   for i = 1, self._timeMesh:size() do
      self._timeMesh[i] = timeMesh[i]
   end

   self._data:expand(nVal)   
   self:_copy_to_dynvector(data)

   reader:close()
end

function DynVector:removeLast()
   local tm = self._timeMesh:popLast()
   local v  = self._data:popLast()
   return tm, v
end

function DynVector:size() return self._data:size() end

local function split (inputstr, sep)
   if sep == nil then
      sep = "%s"
   end
   local t={}
   for str in string.gmatch(inputstr, "([^"..sep.."]+)") do
      table.insert(t, str)
   end
   return t
end

function DynVector:timeMesh() return self._timeMesh end

function DynVector:write(outNm, tmStamp, frNum, flushData, appendData)
   local comm = self._ioComm
   local rank = Mpi.Comm_rank(Mpi.COMM_WORLD)
   if rank ~= 0 then  -- Only run on rank 0 ...
      self:clear()    -- ... but clear data on all ranks.
      return
   end

   if not tmStamp then tmStamp = 0.0 end -- Default time-stamp.
   local flushData  = xsys.pickBool(flushData, true)  -- Default flush data on write.
   local appendData = xsys.pickBool(appendData, true) -- Default append data to single file.

   if appendData and (frNum and frNum>=0) then 
      self.frNum = frNum 
   else 
      self.frNum = "" 
      frNum      = frNum or 0
   end

   -- Create group and set I/O method.
   local grpNm = "DynVector"..frNum..outNm
   local grpId = Adios.declare_group(grpNm, "", Adios.flag_no)
   Adios.select_method(grpId, "MPI", "", "")
   
   -- ADIOS expects CSV string to specify data shape
   local localTmSz  = toCSV( {self._data:size()} )
   local localDatSz = toCSV( {self._data:size(), self._numComponents} )
   
   -- Define data to write.
   Adios.define_var(
      grpId, "TimeMesh"..self.frNum, "", Adios.double, localTmSz, "", "")
   Adios.define_var(
      grpId, "Data"..self.frNum, "", Adios.double, localDatSz, "", "")

   if self._isFirst then
      -- Write meta-data for this file.
      for attrNm, v in pairs(self._metaData) do
         if v.vType == "integer" then
            Adios.define_attribute_byvalue(grpId, attrNm, "", Adios.integer, 1, v.value)
         elseif v.vType == "double" then
            Adios.define_attribute_byvalue(grpId, attrNm, "", Adios.double, 1, v.value)
         elseif v.vType == "string" then
            Adios.define_attribute_byvalue(grpId, attrNm, "", Adios.string, 1, v.value)
         elseif v.vType == "table" then
            Adios.define_attribute_byvalue(grpId, attrNm, "", v.elementType, v.numElements, v.value)
         end
      end
      self._isFirst = false
   end

   local fullNm = GKYL_OUT_PREFIX .. "_" .. outNm
   if frNum == 0 or not appendData then
      fd = Adios.open(grpNm, fullNm, "w", comm[0])
   else
      fd = Adios.open(grpNm, fullNm, "u", comm[0])
   end

   Adios.write(fd, "TimeMesh"..self.frNum, self._timeMesh:data())
   -- Copy data to IO buffer.
   self._ioBuff:expand(self._data:size()*self._numComponents)
   self:_copy_from_dynvector(self._ioBuff)
   Adios.write(fd, "Data"..self.frNum, self._ioBuff:data())
   
   Adios.close(fd)

   -- Clear data for next round of IO.
   if flushData then 
      local tLast, dLast = self:lastData() -- Save the last data, in case we need it later.
      self.tLast = tLast
      for i = 1, self._numComponents do self.dLast[i] = dLast[i] end
      self:clear()
      self.flushed = true
   end
end

return { DynVector = DynVector }
