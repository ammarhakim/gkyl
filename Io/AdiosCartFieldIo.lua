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
local Proto = require "Lib.Proto"
local Lin   = require "Lib.Linalg"
local Adios = require "Io.Adios"
local Alloc = require "Lib.Alloc"
local Mpi   = require "Comm.Mpi"
local lume  = require "Lib.lume"

local AdiosCartFieldIo = Proto()

function AdiosCartFieldIo:init(tbl)
   -- By default, write out doubles.
   local elct     = tbl.elemType and tbl.elemType or typeof("double")
   self._elemType = elct -- element type stored in field

   -- Can specify a table with rank and comm, so only that rank in that comm writes.
   local writeRankCommTbl = tbl.writeRankInComm
   self._writeRank = writeRankCommTbl and writeRankCommTbl[1] or nil
   self._writeComm = writeRankCommTbl and writeRankCommTbl[2] or nil

   -- Set ADIOS data-types.
   self._elctIoType = Adios.type_double
   if ffi.istype(new(elct), new("double")) then
      self._elctIoType = Adios.type_double
   elseif ffi.istype(new(elct), new("float")) then
      self._elctIoType = Adios.type_float
   elseif ffi.istype(new(elct), new("int")) then
      self._elctIoType = Adios.type_int32_t
   else
      assert(false, "AdiosCartFieldIo: data type not supported.")
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
   -- int16_t length data. This is an ADIOS problem in the file
   -- bp_utils.c function bp_read_data_from_buffer(). Ammar, 5/16/2020
   if GKYL_EMBED_INP and #GKYL_INP_FILE_CONTENTS < GKYL_MAX_INT16 then
      self._metaData = {
         -- Write out input file contents (encoded as base64 string).
         inputfile = { value = {data=function(self) return GKYL_INP_FILE_CONTENTS end}, vType = "string",
                       numElements = #GKYL_INP_FILE_CONTENTS, elementType = Adios.type_string, },
      }
   else
      -- Write some dummy text otherwise.
      self._metaData = {
         inputfile = { value = {data=function(self) return "none" end}, vType = "string",
                       numElements = #("none"), elementType = Adios.type_string, },
      }
   end
   if tbl.metaData then
      -- Store value and its type for each piece of data.
      for k,v in pairs(tbl.metaData) do
	 if type(v) == "number" then
	    -- Check if this is an integer or float.
	    if math.floor(math.abs(v)) == math.abs(v) then
	       self._metaData[k] = { value = Lin.IntVec(1), vType = "integer", }
               self._metaData[k].value[1] = v
	    else
	       self._metaData[k] = { value = Lin.Vec(1), vType = "double", }
               self._metaData[k].value[1] = v
	    end
	 elseif type(v) == "string" then
	    self._metaData[k] = { value = {data=function(self) return v end}, vType = "string" }
	 elseif type(v) == "table" then
	    assert(type(v[1])=="number", "Io.AdiosCartFieldIo: Metadata table elements must be numbers.")
            isInt = (math.floor(math.abs(v[1])) == math.abs(v[1]))
            for _, val in pairs(v) do
	       assert(isInt == (math.floor(math.abs(val)) == math.abs(val)), "Io.AdiosCartFieldIo: Metadata table must have elements of the same type (int or double).")
            end
            if isInt then 
               self._metaData[k] = { value = Lin.IntVec(#v), vType = "table", numElements = #v, elementType = Adios.type_int32_t }
            else
               self._metaData[k] = { value = Lin.Vec(#v), vType = "table", numElements = #v, elementType = Adios.type_double }
            end
            for i = 1, #v do self._metaData[k].value[i] = v[i] end
	 end
      end
   end

   -- Time stamp and frame number arrays.
   self.tmStampBuff, self.frNumBuff = Lin.Vec(1), Lin.IntVec(1)

   self.ad_var_fld = {}
   self.localSz, self.globalSz, self.offset = {}, {}, {}
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
      numFields = numFields + 1 
      if GKYL_USE_GPU and fld:hasCuDev() then fld:copyDeviceToHostAsync() end
   end

   -- Assume fields are defined on the same grid (and distributed across the same MPI
   -- communicator) and grab grid descriptors from the first field in the table.
   local field
   for _, fld in pairs(fieldsTbl) do
      field = fld
      break
   end

   local grid     = field:grid()
   local dataComm = grid:commSet().comm

   self.ad_io = self.ad_io or Adios.declare_io(GKYL_ADIOS2_MPI, fName)

   -- No need to do anything if communicator is not valid.
   if not Mpi.Is_comm_valid(dataComm) then return end

   frNum = frNum or -1        -- Default frame-number.
   tmStamp = tmStamp or -1.0  -- Default time-stamp.

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
         
   local attr_names = Adios.available_attributes(self.ad_io)
   local var_names = Adios.available_variables(self.ad_io)

   -- Only need to define attributes and other things once for each adios2_io object.
   if (not lume.any(var_names, function(e) return e=="time" end)) and
      (not lume.any(attr_names, function(e) return e=="numCells" end)) then 
      
      -- Global attributes for Gkyl build.
      if GKYL_GIT_CHANGESET ~= "" then
         Adios.define_attribute(self.ad_io, "changeset", Adios.type_string, GKYL_GIT_CHANGESET)
         Adios.define_attribute(self.ad_io, "builddate", Adios.type_string, GKYL_BUILD_DATE)
      end
      
      -- Field attributes.
      Adios.define_attribute(self.ad_io, "type", Adios.type_string, grid:id())
      if self._metaData["grid"] == nil then   -- Otherwise it gets written below with other meta data.
         local gridFullNm = GKYL_OUT_PREFIX .. "_grid.bp"
         Adios.define_attribute(self.ad_io, "grid", Adios.type_string, gridFullNm)
      end
      
      local cells, lower, upper = globalRange:shape(), grid:lower(), grid:upper()
      if _writeGhost then
         for d = 1, ndim do 
            lower[d] = lower[d] - field:lowerGhost()*grid:dx(d)
            upper[d] = upper[d] + field:upperGhost()*grid:dx(d)
         end
      end
      Adios.define_attribute_array(self.ad_io, "numCells", Adios.type_int32_t, cells, ndim)
      Adios.define_attribute_array(self.ad_io, "lowerBounds", Adios.type_double, lower, ndim)
      Adios.define_attribute_array(self.ad_io, "upperBounds", Adios.type_double, upper, ndim)
      
      -- Write meta-data for this file.
      for attrNm, v in pairs(self._metaData) do
         if v.vType == "integer" then
            Adios.define_attribute(self.ad_io, attrNm, Adios.type_int32_t, v.value:data())
         elseif v.vType == "double" then
            Adios.define_attribute(self.ad_io, attrNm, Adios.type_double, v.value:data())
         elseif v.vType == "string" then
            Adios.define_attribute(self.ad_io, attrNm, Adios.type_string, v.value:data())
         elseif v.vType == "table" then
            Adios.define_attribute_array(self.ad_io, attrNm, v.elementType, v.numElements, v.value:data())
         end
      end
      
      -- Define data to write.
      self.ad_var_time  = Adios.define_variable(self.ad_io, "time", Adios.type_double)
      self.ad_var_frame = Adios.define_variable(self.ad_io, "frame", Adios.type_int32_t)
   else
      self.ad_var_time  = Adios.inquire_variable(self.ad_io, "time")
      self.ad_var_frame = Adios.inquire_variable(self.ad_io, "frame")
   end

   -- Get local (to this MPI rank) and global shape of dataset. Offset is 0 for now.
   for d = 1, ndim do
      self.localSz[d]  = localRange:shape(d)
      self.globalSz[d] = globalRange:shape(d)
      self.offset[d]   = _writeGhost and localRange:lower(d) or localRange:lower(d)-field:globalRange():lower(d)
   end

   local writeRank = self._writeRank or grid:commSet().writeRank
   local writeComm = self._writeComm or dataComm
   -- Don't do anything if the writeRank of the field (which is the rank from which data should be written)
   -- does not match the rank of the global communicator.
   -- This is for cases when the communicator has been split, and the write only happens over
   -- a subset of the domain (and a subset of the global ranks).
   if writeRank ~= Mpi.Comm_rank(writeComm) then return end

   local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.

   -- Open file to write.
   local ad_engine = Adios.open_new_comm(self.ad_io, fullNm, Adios.mode_write, dataComm)

   self.tmStampBuff[1], self.frNumBuff[1] = tmStamp, frNum
   local _ = Adios.put(ad_engine, self.ad_var_time, self.tmStampBuff:data(), Adios.mode_deferred)
   local _ = Adios.put(ad_engine, self.ad_var_frame, self.frNumBuff:data(), Adios.mode_deferred)

   for fldNm, fld in pairs(fieldsTbl) do
      self._outBuff[fldNm] = self._outBuff[fldNm] or self._allocator(1)
      self._outBuff[fldNm]:expand(localRange:volume()*fld:numComponents())

      self.localSz[ndim+1]  = fld:numComponents()
      self.globalSz[ndim+1] = fld:numComponents()
      self.offset[ndim+1]   = 0

      -- Copy field into output buffer (this copy is needed as
      -- field also contains ghost-cell data, and, in addition,
      -- ADIOS expects data to be laid out in row-major order).
      fld._zero:copy_to_buffer(self._outBuff[fldNm]:data(), localRange)

      self.ad_var_fld[fldNm] = self.ad_var_fld[fldNm] or
         (lume.any(var_names, function(e) return e==fldNm end) and
            Adios.inquire_variable(self.ad_io, fldNm) or
            Adios.define_variable(self.ad_io, fldNm, self._elctIoType, ndim+1, self.globalSz, self.offset, self.localSz, true))
      local _ = Adios.put(ad_engine, self.ad_var_fld[fldNm], self._outBuff[fldNm]:data(), Adios.mode_sync)
   end
   local _ = Adios.close(ad_engine)
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

   local grid = field:grid()
   local comm = grid:commSet().comm

   self.ad_io = self.ad_io or Adios.declare_io(GKYL_ADIOS2_MPI, fName)

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
            
      -- Get local (to this MPI rank) and global shape of dataset. Offset is 0 for now.
      for d = 1, ndim do
         self.localSz[d]  = localRange:shape(d)
         self.globalSz[d] = globalRange:shape(d)
         self.offset[d]   = _readGhost and localRange:lower(d) or localRange:lower(d)-field:globalRange():lower(d)
      end

      local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.

      -- Open file to read.
      local ad_engine = Adios.open_new_comm(self.ad_io, fullNm, Adios.mode_readRandomAccess, comm)

      local ad_var_time  = Adios.inquire_variable(self.ad_io, "time")
      local ad_var_frame = Adios.inquire_variable(self.ad_io, "frame")
      local ad_err = Adios.get(ad_engine, ad_var_time, self.tmStampBuff:data(), Adios.mode_deferred)
      local ad_err = Adios.get(ad_engine, ad_var_frame, self.frNumBuff:data(), Adios.mode_deferred)

      for fldNm, fld in pairs(fieldsTbl) do
         self._outBuff[fldNm] = self._outBuff[fldNm] or self._allocator(1)
         self._outBuff[fldNm]:expand(localRange:volume()*fld:numComponents())

         self.localSz[ndim+1]  = fld:numComponents()
         self.globalSz[ndim+1] = fld:numComponents()
         self.offset[ndim+1]   = 0

         local ad_var_fld = Adios.inquire_variable(self.ad_io, fldNm)
         local ad_err = Adios.set_selection(ad_var_fld, ndim+1, self.offset, self.localSz)
         local ad_err = Adios.get(ad_engine, ad_var_fld, self._outBuff[fldNm]:data(), Adios.mode_deferred)
      end

      local _ = Adios.close(ad_engine)

      -- Copy output buffer into field.
      for fldNm, fld in pairs(fieldsTbl) do
         fld._zero:copy_from_buffer(self._outBuff[fldNm]:data(), localRange)
	 if GKYL_USE_GPU and fld:hasCuDev() then fld:copyHostToDevice() end
      end
   end

   return self.tmStampBuff[1], self.frNumBuff[1]
end

return AdiosCartFieldIo
