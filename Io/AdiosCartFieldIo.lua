-- Gkyl ------------------------------------------------------------------------
--
-- CartField I/O using ADIOS format.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries.
local Mpi   = require "Comm.Mpi"
local Adios = require "Io.Adios"
local Alloc = require "Lib.Alloc"
local Proto = require "Lib.Proto"
local Lin   = require "Lib.Linalg"

-- Code from Lua wiki to convert table to comma-seperated-values
-- string.
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

-- AdiosCartFieldIo ------------------------------------------------------------
--
-- CartField I/O using ADIOS.
--------------------------------------------------------------------------------

-- CartField I/O using ADIOS.
local AdiosCartFieldIo = Proto()

-- Constructor to make a new uniform grid.
function AdiosCartFieldIo:init(tbl)
   -- By default, write out doubles.
   local elct     = tbl.elemType and tbl.elemType or typeof("double")
   self._elemType = elct -- element type stored in field
   self._method   = tbl.method and tbl.method or "MPI"

   -- Set ADIOS data-types.
   self._elctIoType = Adios.double
   if ffi.istype(new(elct), new("double")) then
      self._elctIoType = Adios.double
   elseif ffi.istype(new(elct), new("float")) then
      self._elctIoType = Adios.real
      isNumberType = true
   elseif ffi.istype(new(elct), new("int")) then
      self._elctIoType = Adios.integer
   else
      self._elctIoType = Adios.byte
   end

   -- Create memory allocator.
   self._allocator = Alloc.Alloc_meta_ctor(elct)
   -- Allocate memory buffer for use in ADIOS I/O.
   self._outBuff   = self._allocator(1) -- This will be resized on an actual write().

   self._writeGhost = xsys.pickBool(tbl.writeGhost, false)

   -- If we have meta-data to write out, store it.
   self._metaData = {
      -- DISABLING OUTPUT FOR INPUT FILE FOR NOW AT THIS BARFS IT THE
      -- STRING IS TOO LONG
      ["inputfile"] = {
      	 value = "inputfile"
      }
      -- -- We always write out input file contents (encoded as base64 string).
      -- ["inputfile"] = {
      -- 	 value = GKYL_INP_FILE_CONTENTS, vType = "string"
      -- }
   }
   if tbl.metaData then
      -- Store value and its type for each piece of data.
      for k,v in pairs(tbl.metaData) do
	 if type(v) == "number" then
	    -- Check if this is an integer or float.
	    if math.floor(math.abs(v)) == math.abs(v) then
	       self._metaData[k] = {
		  value = new("int[1]", v), vType = "integer",
	       }
	    else
	       self._metaData[k] = {
		  value = new("double[1]", v), vType = "double",
	       }	       
	    end
	 elseif type(v) == "string" then
	    self._metaData[k] = {
	       value = v, vType = "string"
	    }
	 end
      end
   end
end

-- Writes field to file.
-- fName: file name
-- tmStamp: time-stamp
-- frNum: frame number
-- writeGhost: Flag to indicate if we should write ghost-cells
-- 
function AdiosCartFieldIo:write(field, fName, tmStamp, frNum, writeGhost)
   local _writeGhost = self._writeGhost
   if writeGhost ~= nil then _writeGhost = writeGhost end
   local comm =  Mpi.getComm(field:grid():commSet().nodeComm)
   -- (the extra getComm() is needed as Lua has no concept of
   -- pointers and hence we don't know before hand if nodeComm is a
   -- pointer or an object)

   -- No need to do anything if communicator is not valid.
   if not Mpi.Is_comm_valid(comm) then return end

   local ndim = field:ndim()
   local localRange, globalRange = field:localRange(), field:globalRange()
   if _writeGhost then 
      localRange  = field:localExtRange() 
      globalRange = field:globalExtRange() 
   end
   
   -- For use in ADIOS output.
   local _adLocalSz, _adGlobalSz, _adOffset = {}, {}, {}
   for d = 1, ndim do
      _adLocalSz[d]  = localRange:shape(d)
      _adGlobalSz[d] = globalRange:shape(d)
      _adOffset[d]   = localRange:lower(d)-1
      if _writeGhost then _adOffset[d] = _adOffset[d] + 1 end
   end
   _adLocalSz[ndim+1]  = field:numComponents()
   _adGlobalSz[ndim+1] = field:numComponents()
   _adOffset[ndim+1]   = 0

   -- Convert tables to comma-seperated-string. For some strange
   -- reasons, this is what ADIOS expects.
   adLocalSz  = toCSV(_adLocalSz)
   adGlobalSz = toCSV(_adGlobalSz)
   adOffset   = toCSV(_adOffset)

   if not frNum then frNum = 5000 end    -- Default frame-number.
   if not tmStamp then tmStamp = 0.0 end -- Default time-stamp.

   -- Resize buffer (only done if needed. Alloc handles this automatically).
   self._outBuff:expand(field:size())

   local rank = Mpi.Comm_rank(comm)

   -- Setup ADIOS for IO.
   Adios.init_noxml(comm)
   --Adios.set_max_buffer_size(16) -- 16 MB chunks	 

   -- Setup group and set I/O method.
   local grpId = Adios.declare_group("CartField", "", Adios.flag_no)
   Adios.select_method(grpId, self._method, "", "")

   -- Global attributes for Gkyl build.
   Adios.define_attribute_byvalue(grpId, "changeset", "", Adios.string, 1, GKYL_GIT_CHANGESET)
   Adios.define_attribute_byvalue(grpId, "builddate", "", Adios.string, 1, GKYL_BUILD_DATE)

   -- Field attributes.
   Adios.define_attribute_byvalue(grpId, "type", "", Adios.string, 1, field:grid():id())
   local gridFullNm = GKYL_OUT_PREFIX .. "_grid.bp"
   Adios.define_attribute_byvalue(grpId, "grid", "", Adios.string, 1, gridFullNm)
   
   local cells = new("int[?]", ndim)
   for d = 1, ndim do cells[d-1] = globalRange:shape(d) end
   Adios.define_attribute_byvalue(grpId, "numCells", "", Adios.integer, ndim, cells)

   local lower = new("double[?]", ndim)
   for d = 1, ndim do 
      lower[d-1] = field:grid():lower(d) 
      if _writeGhost then lower[d-1] = lower[d-1] - field:lowerGhost()*field:grid():dx(d) end
   end
   Adios.define_attribute_byvalue(grpId, "lowerBounds", "", Adios.double, ndim, lower)

   local upper = new("double[?]", ndim)
   for d = 1, ndim do 
      upper[d-1] = field:grid():upper(d) 
      if _writeGhost then upper[d-1] = upper[d-1] + field:upperGhost()*field:grid():dx(d) end
   end
   Adios.define_attribute_byvalue(grpId, "upperBounds", "", Adios.double, ndim, upper)

   -- Write meta-data for this file.
   for attrNm, v in pairs(self._metaData) do
      if v.vType == "integer" then
	 Adios.define_attribute_byvalue(grpId, attrNm, "", Adios.integer, 1, v.value)
      elseif v.vType == "double" then
	 Adios.define_attribute_byvalue(grpId, attrNm, "", Adios.double, 1, v.value)
      elseif v.vType == "string" then
	 Adios.define_attribute_byvalue(grpId, attrNm, "", Adios.string, 1, v.value)
      end
   end

   -- Define data to write.
   Adios.define_var(
      grpId, "frame", "", Adios.integer, "", "", "")
   Adios.define_var(
      grpId, "time", "", Adios.double, "", "", "")
   Adios.define_var(
      grpId, "CartGridField", "", self._elctIoType, adLocalSz, adGlobalSz, adOffset)

   -- Copy field into output buffer (this copy is needed as
   -- field also contains ghost-cell data, and, in addition,
   -- ADIOS expects data to be laid out in row-major order).
   field:_copy_from_field_region(localRange, self._outBuff)

   local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.
   -- Open file to write out group.
   local fd = Adios.open("CartField", fullNm, "w", comm)

   local tmStampBuff = new("double[1]"); tmStampBuff[0] = tmStamp
   Adios.write(fd, "time", tmStampBuff)

   local frNumBuff = new("int[1]"); frNumBuff[0] = frNum
   Adios.write(fd, "frame", frNumBuff)

   Adios.write(fd, "CartGridField", self._outBuff:data())
   Adios.close(fd)
   
   Adios.finalize(rank)
end

-- Read field from file.
-- fName: file name
function AdiosCartFieldIo:read(field, fName, readGhost) --> time-stamp, frame-number
   local _readGhost = self._writeGhost
   if readGhost ~= nil then _readGhost = readGhost end
   local comm    =  Mpi.getComm(field:grid():commSet().nodeComm)
   local shmComm = Mpi.getComm(field:grid():commSet().sharedComm)
   -- (the extra getComm() is needed as Lua has no concept of
   -- pointers and hence we don't know before hand if nodeComm is a
   -- pointer or an object)

   -- Create time stamp and frame number arrays outside of if statement.
   local tmStampBuff = new("double[1]")
   local frNumBuff   = new("int[1]")

   -- Only read in data if communicator is valid.
   if Mpi.Is_comm_valid(comm) then

      local ndim                    = field:ndim()
      local localRange, globalRange = field:localRange(), field:globalRange()

      if readGhost then
         localRange  = field:localExtRange() 
         globalRange = field:globalExtRange() 
      end

      -- For use in ADIOS output.
      local _adLocalSz, _adGlobalSz, _adOffset = {}, {}, {}
      for d = 1, ndim do
	 _adLocalSz[d]  = localRange:shape(d)
	 _adGlobalSz[d] = globalRange:shape(d)
	 _adOffset[d]   = localRange:lower(d)-1
         if _readGhost then _adOffset[d] = _adOffset[d] + 1 end
      end
      _adLocalSz[ndim+1]  = field:numComponents()
      _adGlobalSz[ndim+1] = field:numComponents()
      _adOffset[ndim+1]   = 0

      -- Convert tables to comma-seperated-string. For some strange
      -- reasons, this is what ADIOS expects.
      adLocalSz  = toCSV(_adLocalSz)
      adGlobalSz = toCSV(_adGlobalSz)
      adOffset   = toCSV(_adOffset)

      -- Resize buffer (only done if needed. Alloc handles this automatically).
      self._outBuff:expand(field:size())

      local rank = Mpi.Comm_rank(comm)

      -- Setup ADIOS for IO.
      Adios.init_noxml(comm)
      --Adios.set_max_buffer_size(16) -- 16 MB chunks.

      -- Setup group and set I/O method.
      local grpId = Adios.declare_group("CartField", "", Adios.flag_no)
      Adios.select_method(grpId, self._method, "", "")

      -- Field attributes.
      local cells = new("int[?]", ndim)
      for d = 1, ndim do cells[d-1] = globalRange:shape(d) end
      Adios.define_attribute_byvalue(grpId, "numCells", "", Adios.integer, ndim, cells)

      local lower = new("double[?]", ndim)
      for d = 1, ndim do 
         lower[d-1] = field:grid():lower(d)
         if _readGhost then lower[d-1] = lower[d-1] - field:lowerGhost()*field:grid():dx(d) end
      end
      Adios.define_attribute_byvalue(grpId, "lowerBounds", "", Adios.double, ndim, lower)

      local upper = new("double[?]", ndim)
      for d = 1, ndim do 
         upper[d-1] = field:grid():upper(d) 
         if _readGhost then upper[d-1] = upper[d-1] + field:upperGhost()*field:grid():dx(d) end
      end
      Adios.define_attribute_byvalue(grpId, "upperBounds", "", Adios.double, ndim, upper)

      -- Define data to read.
      Adios.define_var(
	 grpId, "frame", "", Adios.integer, "", "", "")
      Adios.define_var(
	 grpId, "time", "", Adios.double, "", "", "")
      Adios.define_var(
	 grpId, "CartGridField", "", self._elctIoType, adLocalSz, adGlobalSz, adOffset)

      local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.
      -- Open file to read data.
      local fd = Adios.read_open_file(fullNm, comm)

      Adios.schedule_read(fd, Adios.selBoundingBox, "time", 0, 1, tmStampBuff)
      Adios.schedule_read(fd, Adios.selBoundingBox, "frame", 0, 1, frNumBuff)

      -- ADIOS expects input to be const uint64_t* objects, hence vector
      -- types below)
      local start, count = Lin.UInt64Vec(ndim+1), Lin.UInt64Vec(ndim+1)
      for d = 1, ndim do
         start[d] = localRange:lower(d)-1
         count[d] = localRange:shape(d)
         if _readGhost then start[d] = start[d] + 1 end
      end
      count[ndim+1] = field:numComponents()
      start[ndim+1] = 0
      local sel = Adios.selection_boundingbox(ndim+1, start, count)

      Adios.schedule_read(fd, sel, "CartGridField", 0, 1, self._outBuff:data())
      Adios.perform_reads(fd, 1)

      Adios.read_close(fd) -- No reads actually happen unless one closes file!

      -- Copy output buffer into field.
      field:_copy_to_field_region(localRange, self._outBuff)

      Adios.finalize(rank)

   end
   -- If running with shared memory, need to broadcast time stamp and frame number.
   Mpi.Bcast(tmStampBuff, 1, Mpi.DOUBLE, 0, shmComm)
   Mpi.Bcast(frNumBuff, 1, Mpi.INT, 0, shmComm)
   return tmStampBuff[0], frNumBuff[0]
end

return AdiosCartFieldIo
