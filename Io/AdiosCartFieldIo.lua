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

   -- Can specify a table with rank and comm, so only that rank in that comm writes.
   local writeRankCommTbl = tbl.writeRankInComm
   self._writeRank = writeRankCommTbl and writeRankCommTbl[1] or nil
   self._writeComm = writeRankCommTbl and writeRankCommTbl[2] or nil

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
   -- Memory buffer for use in ADIOS I/O.
   self._outBuff = {}

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
--   writeGhost: Flag to indicate if we should write ghost-cells on boundaries of global domain
function AdiosCartFieldIo:write(fieldsIn, fName, tmStamp, frNum, writeGhost)
   local _writeGhost = writeGhost ~= nil and writeGhost or self._writeGhost

   -- Identify if fieldsIn is a CartField using the self._ndim variable (MF: there ought to be a better way).
   local fieldsTbl = type(fieldsIn._ndim)=="number" and {CartGridField = fieldsIn} or fieldsIn
   local numFields = 0
   for fldNm, fld in pairs(fieldsTbl) do 
      numFields=numFields+1 
      if GKYL_USE_GPU and fld:hasCuDev() then
         fld:copyDeviceToHostAsync()
      end
   end

   -- Assume fields are defined on the same grid (and distributed across the same MPI
   -- communicator) and grab grid descriptors from the first field in the table.
   local field
   for _, fld in pairs(fieldsTbl) do
      field = fld
      break
   end

   local dataComm = Mpi.getComm(field:grid():commSet().comm)
   -- (the extra getComm() is needed as Lua has no concept of
   -- pointers and hence we don't know before hand if nodeComm is a
   -- pointer or an object)

   -- No need to do anything if communicator is not valid.
   if not Mpi.Is_comm_valid(dataComm) then return end

   if not frNum then frNum = 5000 end    -- Default frame-number.
   if not tmStamp then tmStamp = 0.0 end -- Default time-stamp.

   local ndim = field:ndim()
   local localRange, globalRange = field:localRange(), field:globalRange()
   if _writeGhost then 
      -- Extend localRange to include ghost cells if on edge of global domain.
      local localExtRange = field:localExtRange()
      for d = 1, ndim do
         if localRange:lower(d)==globalRange:lower(d) then
            local localExtInDirRange = localRange:extendDir(d, field:lowerGhost(), 0)
            localRange = localExtRange:subRange(localExtInDirRange:lowerAsVec(), localExtInDirRange:upperAsVec())
         end
         if localRange:upper(d)==globalRange:upper(d) then
            local localExtInDirRange = localRange:extendDir(d, 0, field:upperGhost())
            localRange = localExtRange:subRange(localExtInDirRange:lowerAsVec(), localExtInDirRange:upperAsVec())
         end
      end
      globalRange = field:globalExtRange()  -- Extend globalRange to include ghost cells.
   end
         
   -- Get group name based on fName with frame and suffix chopped off.
   local grpNm = string.gsub(string.gsub(fName, "_(%d+).bp", ""), ".bp", "")

   -- Setup group and set I/O method. Only need to do once for each grpNm.
   if not self.grpIds[grpNm] then
      self.grpIds[grpNm] = Adios.declare_group(grpNm, "", Adios.flag_no)
      --Adios.set_max_buffer_size(16)
      Adios.select_method(self.grpIds[grpNm], self._method, "", "")
      
      -- Global attributes for Gkyl build.
      if GKYL_GIT_CHANGESET ~= "" then
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "changeset", "", Adios.string, 1, GKYL_GIT_CHANGESET)
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "builddate", "", Adios.string, 1, GKYL_BUILD_DATE)
      end
      
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
      local _adLocalSz, _adGlobalSz, _adOffset = {}, {}, {}
      -- Get local (to this MPI rank) and global shape of dataset. Offset is 0 for now.
      for d = 1, ndim do
         _adLocalSz[d]  = localRange:shape(d)
         _adGlobalSz[d] = globalRange:shape(d)
         _adOffset[d]   = _writeGhost and localRange:lower(d)
            or localRange:lower(d)-field:globalRange():lower(d)
      end
      for fldNm, fld in pairs(fieldsTbl) do
         self._outBuff[fldNm] = self._outBuff[fldNm] or self._allocator(1)
         self._outBuff[fldNm]:expand(localRange:volume()*fld:numComponents())

         _adLocalSz[ndim+1]  = fld:numComponents()
         _adGlobalSz[ndim+1] = fld:numComponents()
         _adOffset[ndim+1]   = 0

         -- Convert tables to comma-separated-string. This is what ADIOS expects.
         local adLocalSz  = toCSV(_adLocalSz)
         local adGlobalSz = toCSV(_adGlobalSz)
         local adOffset   = toCSV(_adOffset)

         Adios.define_var(
            self.grpIds[grpNm], fldNm, "", self._elctIoType, adLocalSz, adGlobalSz, adOffset)
      end
   end

   local writeRank = self._writeRank or field:grid():commSet().writeRank
   local writeComm = self._writeComm or dataComm
   -- Don't do anything if the writeRank of the field (which is the rank from which data should be written)
   -- does not match the rank of the global communicator.
   -- This is for cases when the communicator has been split, and the write only happens over
   -- a subset of the domain (and a subset of the global ranks).
   if writeRank ~= Mpi.Comm_rank(writeComm) then return end

   local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.

   -- Open file to write out group.
   local fd = Adios.open(grpNm, fullNm, "w", dataComm)

   local tmStampBuff = new("double[1]"); tmStampBuff[0] = tmStamp
   Adios.write(fd, "time", tmStampBuff)

   local frNumBuff = new("int[1]"); frNumBuff[0] = frNum
   Adios.write(fd, "frame", frNumBuff)

   for fldNm, fld in pairs(fieldsTbl) do
      -- Copy field into output buffer (this copy is needed as
      -- field also contains ghost-cell data, and, in addition,
      -- ADIOS expects data to be laid out in row-major order).
      fld._zero:copy_to_buffer(self._outBuff[fldNm]:data(), localRange)

      local err = Adios.write(fd, fldNm, self._outBuff[fldNm]:data())
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
   local _readGhost = readGhost ~= nil and readGhost or self._writeGhost

   -- Identify if fieldsIn is a CartField using the self._ndim variable (MF: there ought to be a better way).
   local fieldsTbl = type(fieldsOut._ndim)=="number" and {CartGridField = fieldsOut} or fieldsOut
   local numFields = 0
   for fldNm, _ in pairs(fieldsTbl) do numFields=numFields+1 end

   -- Assume fields are defined on the same grid (and distributed across the same MPI
   -- communicator) and grab grid descriptors from the first field in the table.
   local field
   for _, fld in pairs(fieldsTbl) do
      field = fld
      break
   end

   local comm = Mpi.getComm(field:grid():commSet().comm)
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
         -- Extend localRange to include ghost cells if on edge of global domain.
         local localExtRange = field:localExtRange()
         for d = 1, ndim do
            if localRange:lower(d)==globalRange:lower(d) then
               local localExtInDirRange = localRange:extendDir(d, field:lowerGhost(), 0)
               localRange = localExtRange:subRange(localExtInDirRange:lowerAsVec(), localExtInDirRange:upperAsVec())
            end
            if localRange:upper(d)==globalRange:upper(d) then
               local localExtInDirRange = localRange:extendDir(d, 0, field:upperGhost())
               localRange = localExtRange:subRange(localExtInDirRange:lowerAsVec(), localExtInDirRange:upperAsVec())
            end
         end
         -- extend globalRange to include ghost cells
         globalRange = field:globalExtRange() 
      end
            
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
            if _readGhost then lower[d-1] = lower[d-1] - field:lowerGhost()*field:grid():dx(d) end
         end
         Adios.define_attribute_byvalue(self.grpIds[grpNm], "lowerBounds", "", Adios.double, ndim, lower)
         
         local upper = new("double[?]", ndim)
         for d = 1, ndim do 
            upper[d-1] = field:grid():upper(d) 
            if _readGhost then upper[d-1] = upper[d-1] + field:upperGhost()*field:grid():dx(d) end
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
         local _adLocalSz, _adGlobalSz, _adOffset = {}, {}, {}
         -- Get local (to this MPI rank) and global shape of dataset. Offset is 0 for now.
         for d = 1, ndim do
            _adLocalSz[d]  = localRange:shape(d)
            _adGlobalSz[d] = globalRange:shape(d)
            _adOffset[d]   = _readGhost and localRange:lower(d)
               or localRange:lower(d)-field:globalRange():lower(d)
         end
         for fldNm, fld in pairs(fieldsTbl) do
            self._outBuff[fldNm] = self._outBuff[fldNm] or self._allocator(1)
            self._outBuff[fldNm]:expand(localRange:volume()*fld:numComponents())

            _adLocalSz[ndim+1]  = fld:numComponents()
            _adGlobalSz[ndim+1] = fld:numComponents()
            _adOffset[ndim+1]   = 0

            -- Convert tables to comma-seperated-string. For some strange
            -- reasons, this is what ADIOS expects.
            local adLocalSz  = toCSV(_adLocalSz)
            local adGlobalSz = toCSV(_adGlobalSz)
            local adOffset   = toCSV(_adOffset)
            Adios.define_var(
               self.grpIds[grpNm], fldNm, "", self._elctIoType, adLocalSz, adGlobalSz, adOffset)
         end
      end

      local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.
      -- Open file to read data.
      local fd = Adios.read_open_file(fullNm, comm)

      Adios.schedule_read(fd, Adios.selBoundingBox, "time", 0, 1, tmStampBuff)
      Adios.schedule_read(fd, Adios.selBoundingBox, "frame", 0, 1, frNumBuff)

      -- ADIOS expects input to be const uint64_t* objects, hence vector types below)
      local start, count = Lin.UInt64Vec(ndim+1), Lin.UInt64Vec(ndim+1)

      for fldNm, fld in pairs(fieldsTbl) do
         for d = 1, ndim do
            local ct = localRange:shape(d)
            local st = _readGhost and localRange:lower(d)
               or localRange:lower(d)-field:globalRange():lower(d)
            start[d] = st
            count[d] = ct
         end
         count[ndim+1] = fld:numComponents()
         start[ndim+1] = 0
         local sel = Adios.selection_boundingbox(ndim+1, start, count)

         Adios.schedule_read(fd, sel, fldNm, 0, 1, self._outBuff[fldNm]:data())
      end
      Adios.perform_reads(fd, 1)

      Adios.read_close(fd) -- No reads actually happen unless one closes file!

      -- Copy output buffer into field.
      for fldNm, fld in pairs(fieldsTbl) do
         fld._zero:copy_from_buffer(self._outBuff[fldNm]:data(), localRange)
	 if GKYL_USE_GPU and fld:hasCuDev() then
	    fld:copyHostToDevice()
	 end
      end
   end

   return tmStampBuff[0], frNumBuff[0]
end

return AdiosCartFieldIo
