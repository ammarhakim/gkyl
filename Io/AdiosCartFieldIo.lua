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
   -- Allocate memory buffer for use in ADIOS I/O. Allocate a number 
   -- of them in case user wants to output multiple fields to one file.
   self._outBuff   = {}
   self._maxFields = 50
   for i = 1, self._maxFields do
      self._outBuff[i] = self._allocator(1) -- This will be resized on an actual read/write.
   end

   -- write ghost cells on boundaries of global domain (for BCs)
   self._writeGhost = xsys.pickBool(tbl.writeGhost, false)

   -- If we have meta-data to write out, store it.
   --
   -- WARNING: For now we are not writing out huge input files dues to
   -- a limitation in the ADIOS reader in which the size is limited to
   -- int16_t lenght data. This is an ADIOS problem in the file
   -- bp_utils.c function bp_read_data_from_buffer(). Ammar, 5/16/2020
   if GKYL_EMBED_INP and #GKYL_INP_FILE_CONTENTS < GKYL_MAX_INT16 then
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
	    assert(type(v[1])=="number", "Io.AdiosCartFieldIo: Metadata table elements must be numbers.")
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

   self.grpIds = {}
end

-- Writes field(s) to file.
-- Inputs:
--   fieldsIn: single field or table of key-value pairs of fields
--             to be written out. Fields must live on the same grid.
--   fName:    file name
--   tmStamp:  time-stamp
--   frNum:    frame number
--   writeGhost: Flag to indicate if we should write skin-cells on boundaries of global domain
function AdiosCartFieldIo:write(fieldsIn, fName, tmStamp, frNum, writeGhost)
   local _writeGhost = self._writeGhost
   if writeGhost ~= nil then _writeGhost = writeGhost end

   -- Identify if fieldsIn is a CartField using the self._ndim variable (MF: there ought to be a better way).
   local fieldsTbl = type(fieldsIn._ndim)=="number" and {CartGridField = fieldsIn} or fieldsIn
   local numFields = 0
   for fldNm, _ in pairs(fieldsTbl) do numFields=numFields+1 end
   -- The check below is not actually needed for writes, but it is needed for reads
   -- so we use them for writes too in case this file will be read in (e.g. restart).
   assert(numFields <= self._maxFields, "AdiosCartFieldIo: Cannot read/write more fields than self._maxFields.")

   -- Assume fields are defined on the same grid and grab
   -- grid descriptors from the first field in the table.
   local field
   for _, fld in pairs(fieldsTbl) do
      field = fld
      break
   end

   local comm  = Mpi.getComm(field:grid():commSet().nodeComm)
   -- (the extra getComm() is needed as Lua has no concept of
   -- pointers and hence we don't know before hand if nodeComm is a
   -- pointer or an object)

   -- No need to do anything if communicator is not valid.
   if not Mpi.Is_comm_valid(comm) then return end
   local rank = Mpi.Comm_rank(comm)

   local ndim = field:ndim()
   local localRange, globalRange = field:localRange(), field:globalRange()
   if _writeGhost then 
      -- extend localRange to include ghost cells if on edge of global domain
      for d = 1, ndim do
         if localRange:lower(d) == globalRange:lower(d) then localRange = localRange:extendDir(d, field:lowerGhost(), 0) end
         if localRange:upper(d) == globalRange:upper(d) then localRange = localRange:extendDir(d, 0, field:upperGhost()) end
      end
      -- extend globalRange to include ghost cells
      globalRange = field:globalExtRange() 
   end
   
   -- For use in ADIOS output.
   local _adLocalSz, _adGlobalSz, _adOffset = {}, {}, {}
   for d = 1, ndim do
      _adLocalSz[d]  = localRange:shape(d)
      _adGlobalSz[d] = globalRange:shape(d)
      _adOffset[d]   = localRange:lower(d)-1
      if _writeGhost then _adOffset[d] = _adOffset[d] + field:lowerGhost() end
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
   self._outBuff[1]:expand(localRange:volume()*field:numComponents())

   -- Get group name based on fName with frame and suffix chopped off.
   local grpNm = string.gsub(string.gsub(fName, "_(%d+).bp", ""), ".bp", "")

   -- Setup group and set I/O method. Only need to do once for each grpNm.
   if not self.grpIds[grpNm] then
      self.grpIds[grpNm] = Adios.declare_group(grpNm, "", Adios.flag_no)
      --Adios.set_max_buffer_size(16)
      Adios.select_method(self.grpIds[grpNm], self._method, "", "")
      
      -- Global attributes for Gkyl build.
      Adios.define_attribute_byvalue(self.grpIds[grpNm], "changeset", "", Adios.string, 1, GKYL_GIT_CHANGESET)
      Adios.define_attribute_byvalue(self.grpIds[grpNm], "builddate", "", Adios.string, 1, GKYL_BUILD_DATE)
      
      -- Field attributes.
      Adios.define_attribute_byvalue(self.grpIds[grpNm], "type", "", Adios.string, 1, field:grid():id())
      if self._metaData["grid"] == nil then   -- Otherwise it gets written below with other meta data.
         local gridFullNm = GKYL_OUT_PREFIX .. "_grid.bp"
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "grid", "", Adios.string, 1, gridFullNm)
      end
      
      local cells = new("int[?]", ndim)
      for d = 1, ndim do cells[d-1] = globalRange:shape(d) end
      Adios.define_attribute_byvalue(self.grpIds[grpNm], "numCells", "", Adios.integer, ndim, cells)
      
      local lower = new("double[?]", ndim)
      for d = 1, ndim do 
         lower[d-1] = field:grid():lower(d) 
         if _writeGhost then lower[d-1] = lower[d-1] - field:lowerGhost()*field:grid():dx(d) end
      end
      Adios.define_attribute_byvalue(self.grpIds[grpNm], "lowerBounds", "", Adios.double, ndim, lower)
      
      local upper = new("double[?]", ndim)
      for d = 1, ndim do 
         upper[d-1] = field:grid():upper(d) 
         if _writeGhost then upper[d-1] = upper[d-1] + field:upperGhost()*field:grid():dx(d) end
      end
      Adios.define_attribute_byvalue(self.grpIds[grpNm], "upperBounds", "", Adios.double, ndim, upper)
      
      -- Write meta-data for this file.
      for attrNm, v in pairs(self._metaData) do
         if v.vType == "integer" then
            Adios.define_attribute_byvalue(self.grpIds[grpNm], attrNm, "", Adios.integer, 1, v.value)
         elseif v.vType == "double" then
            Adios.define_attribute_byvalue(self.grpIds[grpNm], attrNm, "", Adios.double, 1, v.value)
         elseif v.vType == "string" then
            Adios.define_attribute_byvalue(self.grpIds[grpNm], attrNm, "", Adios.string, 1, v.value)
         elseif v.vType == "table" then
            Adios.define_attribute_byvalue(self.grpIds[grpNm], attrNm, "", v.elementType, v.numElements, v.value)
         end
      end
      
      -- Define data to write.
      Adios.define_var(
         self.grpIds[grpNm], "frame", "", Adios.integer, "", "", "")
      Adios.define_var(
         self.grpIds[grpNm], "time", "", Adios.double, "", "", "")
      for fldNm, _ in pairs(fieldsTbl) do
         Adios.define_var(
            self.grpIds[grpNm], fldNm, "", self._elctIoType, adLocalSz, adGlobalSz, adOffset)
      end
   end

   local writeRank = field:grid():commSet().writeRank
   -- Don't do anything if the writeRank of the field (which is the rank from which data should be written)
   -- does not match the rank of the global communicator.
   -- This is for cases when the communicator has been split, and the write only happens over
   -- a subset of the domain (and a subset of the global ranks).
   if writeRank ~= Mpi.Comm_rank(Mpi.COMM_WORLD) then return end

   local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.


   -- Open file to write out group.
   local fd = Adios.open(grpNm, fullNm, "w", comm)

   local tmStampBuff = new("double[1]"); tmStampBuff[0] = tmStamp
   Adios.write(fd, "time", tmStampBuff)

   local frNumBuff = new("int[1]"); frNumBuff[0] = frNum
   Adios.write(fd, "frame", frNumBuff)

   for fldNm, fld in pairs(fieldsTbl) do
      -- Copy field into output buffer (this copy is needed as
      -- field also contains ghost-cell data, and, in addition,
      -- ADIOS expects data to be laid out in row-major order).
      fld:_copy_from_field_region(localRange, self._outBuff[1])

      Adios.write(fd, fldNm, self._outBuff[1]:data())
   end
   Adios.close(fd)
end

-- Read field(s) from file.
-- Inputs:
--   fieldsOut: single field or table of key-value pairs of fields to
--              be read in. Fields must live on the same grid.
--   fName:     file name.
--   readGhost:  Flag to indicate if we should read skin-cells on boundaries of global domain.
function AdiosCartFieldIo:read(fieldsOut, fName, readGhost) --> time-stamp, frame-number
   local _readGhost = self._writeGhost
   if readGhost ~= nil then _readGhost = readGhost end

   -- Identify if fieldsIn is a CartField using the self._ndim variable (MF: there ought to be a better way).
   local fieldsTbl = type(fieldsOut._ndim)=="number" and {CartGridField = fieldsOut} or fieldsOut
   local numFields = 0
   for fldNm, _ in pairs(fieldsTbl) do numFields=numFields+1 end
   assert(numFields <= self._maxFields, "AdiosCartFieldIo: Cannot read/write more fields than self._maxFields.")

   -- Assume fields are defined on the same grid and grab
   -- grid descriptors from the first field in the table.
   local field
   for _, fld in pairs(fieldsTbl) do
      field = fld
      break
   end

   local comm    = Mpi.getComm(field:grid():commSet().nodeComm)
   local shmComm = Mpi.getComm(field:grid():commSet().sharedComm)
   -- (the extra getComm() is needed as Lua has no concept of
   -- pointers and hence we don't know before hand if nodeComm is a
   -- pointer or an object)

   -- Create time stamp and frame number arrays outside of if statement.
   local tmStampBuff = new("double[1]")
   local frNumBuff   = new("int[1]")

   -- Only read in data if communicator is valid.
   if Mpi.Is_comm_valid(comm) then

      local ndim = field:ndim()
      local localRange, globalRange = field:localRange(), field:globalRange()
      if _readGhost then 
         -- extend localRange to include ghost cells if on edge of global domain
         for d = 1, ndim do
            if localRange:lower(d) == globalRange:lower(d) then localRange = localRange:extendDir(d, field:lowerGhost(), 0) end
            if localRange:upper(d) == globalRange:upper(d) then localRange = localRange:extendDir(d, 0, field:upperGhost()) end
         end
         -- extend globalRange to include ghost cells
         globalRange = field:globalExtRange() 
      end
      
      -- For use in ADIOS output.
      local _adLocalSz, _adGlobalSz, _adOffset = {}, {}, {}
      for d = 1, ndim do
         _adLocalSz[d]  = localRange:shape(d)
         _adGlobalSz[d] = globalRange:shape(d)
         _adOffset[d]   = localRange:lower(d)-1
         if _readGhost then _adOffset[d] = _adOffset[d] + field:lowerGhost() end
      end
      _adLocalSz[ndim+1] = field:numComponents()
      _adGlobalSz[ndim+1] = field:numComponents()
      _adOffset[ndim+1] = 0

      -- Convert tables to comma-seperated-string. For some strange
      -- reasons, this is what ADIOS expects.
      adLocalSz = toCSV(_adLocalSz)
      adGlobalSz = toCSV(_adGlobalSz)
      adOffset = toCSV(_adOffset)

      -- Resize buffer (only done if needed. Alloc handles this automatically).
      local fldI = 0
      for fldNm, _ in pairs(fieldsTbl) do
         fldI = fldI + 1
         self._outBuff[fldI]:expand(localRange:volume()*field:numComponents())
      end

      local rank = Mpi.Comm_rank(comm)

      -- Get group name based on fName with frame and suffix chopped off.
      local grpNm = string.gsub(string.gsub(fName, "_(%d+).bp", ""), ".bp", "")

      -- Setup group and set I/O method. Only need to do once for each grpNm.
      if not self.grpIds[grpNm] then
         self.grpIds[grpNm] = Adios.declare_group(grpNm, "", Adios.flag_no)
         --Adios.set_max_buffer_size(16)
         Adios.select_method(self.grpIds[grpNm], self._method, "", "")
         
         -- Global attributes for Gkyl build.
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "changeset", "", Adios.string, 1, GKYL_GIT_CHANGESET)
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "builddate", "", Adios.string, 1, GKYL_BUILD_DATE)
         
         -- Field attributes.
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "type", "", Adios.string, 1, field:grid():id())
         local gridFullNm = GKYL_OUT_PREFIX .. "_grid.bp"
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "grid", "", Adios.string, 1, gridFullNm)
         
         local cells = new("int[?]", ndim)
         for d = 1, ndim do cells[d-1] = globalRange:shape(d) end
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "numCells", "", Adios.integer, ndim, cells)
         
         local lower = new("double[?]", ndim)
         for d = 1, ndim do 
            lower[d-1] = field:grid():lower(d) 
            if _writeGhost then lower[d-1] = lower[d-1] - field:lowerGhost()*field:grid():dx(d) end
         end
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "lowerBounds", "", Adios.double, ndim, lower)
         
         local upper = new("double[?]", ndim)
         for d = 1, ndim do 
            upper[d-1] = field:grid():upper(d) 
            if _writeGhost then upper[d-1] = upper[d-1] + field:upperGhost()*field:grid():dx(d) end
         end
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "upperBounds", "", Adios.double, ndim, upper)
         
         -- Write meta-data for this file.
         for attrNm, v in pairs(self._metaData) do
            if v.vType == "integer" then
               Adios.define_attribute_byvalue(self.grpIds[grpNm], attrNm, "", Adios.integer, 1, v.value)
            elseif v.vType == "double" then
               Adios.define_attribute_byvalue(self.grpIds[grpNm], attrNm, "", Adios.double, 1, v.value)
            elseif v.vType == "string" then
               Adios.define_attribute_byvalue(self.grpIds[grpNm], attrNm, "", Adios.string, 1, v.value)
            elseif v.vType == "table" then
               Adios.define_attribute_byvalue(self.grpIds[grpNm], attrNm, "", v.elementType, v.numElements, v.value)
            end
         end
         
         -- Define data to read.
         Adios.define_var(
            self.grpIds[grpNm], "frame", "", Adios.integer, "", "", "")
         Adios.define_var(
            self.grpIds[grpNm], "time", "", Adios.double, "", "", "")
         for fldNm, _ in pairs(fieldsTbl) do
            Adios.define_var(
               self.grpIds[grpNm], fldNm, "", self._elctIoType, adLocalSz, adGlobalSz, adOffset)
         end
      end

      local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.
      -- Open file to read data.
      local fd = Adios.read_open_file(fullNm, comm)

      Adios.schedule_read(fd, Adios.selBoundingBox, "time", 0, 1, tmStampBuff)
      Adios.schedule_read(fd, Adios.selBoundingBox, "frame", 0, 1, frNumBuff)

      -- ADIOS expects input to be const uint64_t* objects, hence vector
      -- types below)
      local start, count = Lin.UInt64Vec(ndim+1), Lin.UInt64Vec(ndim+1)
      for d = 1, ndim do
         local st = localRange:lower(d)-1
         local ct = localRange:shape(d)
         if _readGhost then st = st + field:lowerGhost() end
         start[d] = st
         count[d] = ct
      end
      count[ndim+1] = field:numComponents()
      start[ndim+1] = 0
      local sel = Adios.selection_boundingbox(ndim+1, start, count)

      local fldI = 0
      for fldNm, _ in pairs(fieldsTbl) do
         fldI = fldI + 1
         Adios.schedule_read(fd, sel, fldNm, 0, 1, self._outBuff[fldI]:data())
      end
      Adios.perform_reads(fd, 1)

      Adios.read_close(fd) -- No reads actually happen unless one closes file!

      -- Copy output buffer into field.
      local fldI = 0
      for _, fld in pairs(fieldsTbl) do
         fldI = fldI + 1
         fld:_copy_to_field_region(localRange, self._outBuff[fldI])
      end
   end
   -- If running with shared memory, need to broadcast time stamp and frame number.
   Mpi.Bcast(tmStampBuff, 1, Mpi.DOUBLE, 0, shmComm)
   Mpi.Bcast(frNumBuff, 1, Mpi.INT, 0, shmComm)
   return tmStampBuff[0], frNumBuff[0]
end

return AdiosCartFieldIo
