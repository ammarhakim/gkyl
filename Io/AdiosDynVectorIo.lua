-- Gkyl ------------------------------------------------------------------------
--
-- DynVector I/O using ADIOS format.
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
local Mpi   = require "Comm.Mpi"
local Alloc = require "Lib.Alloc"

local AdiosDynVectorIo = Proto()

function AdiosDynVectorIo:init(tbl)

   -- adios2_io object.
   self.ad_io = tbl.adios_io
   if self.ad_io == nil then
      self.ad = assert(tbl.ioSystem, "AdiosDynVectorIo: must specify the adios object in 'ioSystem', or the adios io object in 'adios_io'.")
   end

   -- Write only from one rank of the specified communicator. Use world comm
   -- and rank 0 if user doesn't specify them.
   self._parentComm = tbl.comm or Mpi.COMM_WORLD
   self.writeRank   = tbl.writeRank or 0
   if self.writeRank ~= Mpi.PROC_NULL then
      local ranks  = Lin.IntVec(1);  ranks[1] = self.writeRank
      self._ioComm = Mpi.Split_comm(self._parentComm, ranks)
   end

   -- Create memory allocator.
   local elct = typeof("double")
   self._allocator = Alloc.Alloc_meta_ctor(elct)
   -- Memory buffer for use in ADIOS I/O.
   self._outBuff = self._allocator(1)

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
            assert(type(v[1])=="number", "Io.AdiosDynVectorIo: Metadata table elements must be numbers.")
            isInt = (math.floor(math.abs(v[1])) == math.abs(v[1]))
            for _, val in pairs(v) do
               assert(isInt == (math.floor(math.abs(val)) == math.abs(val)), "Io.AdiosDynVectorIo: Metadata table must have elements of the same type (int or double).")
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
   self.shape1d, self.start1d = Lin.UInt64Vec(1), Lin.UInt64Vec(1)
   self.shape2d, self.start2d = Lin.UInt64Vec(2), Lin.UInt64Vec(2)

   self.start1d[1] = 0
   self.start2d[1], self.start2d[2] = 0, 0

   self.defineAttrVars = true
end

function AdiosDynVectorIo:write(dynVecIn, fName, tmStamp, frNum, appendData)
   local appendData = xsys.pickBool(appendData, true) -- Default append data to single file.

   frNum = frNum or -1        -- Default frame-number.
   tmStamp = tmStamp or -1.0  -- Default time-stamp.

   -- Create an adios2_object if needed.
   self.ad_io = self.ad_io or Adios.declare_io(self.ad, fName)

   local rank = Mpi.Comm_rank(self._parentComm)
   if rank ~= self.writeRank then  -- Only run on rank that writes ...
      dynVecIn:clear()             -- ... but clear data on all ranks.
      return
   end

   -- Only need to define attributes and other things once for each adios2_io object.
   if self.defineAttrVars then

      -- Global attributes for Gkyl build.
      if GKYL_GIT_CHANGESET ~= "" then
         Adios.define_attribute(self.ad_io, "changeset", Adios.type_string, GKYL_GIT_CHANGESET)
         Adios.define_attribute(self.ad_io, "builddate", Adios.type_string, GKYL_BUILD_DATE)
      end

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

      self.ad_var = {}

      self.defineAttrVars = false
   end

   local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.

   -- Open file to write.
   local ad_engine
   if frNum == 0 or not appendData then
      ad_engine = Adios.open_new_comm(self.ad_io, fullNm, Adios.mode_write, self._ioComm)
   else
      ad_engine = Adios.open_new_comm(self.ad_io, fullNm, Adios.mode_append, self._ioComm)
   end

   local timeVarNm = "TimeMesh" .. frNum
   local dataVarNm = "Data" .. frNum

   self.shape1d[1] = dynVecIn:size()
   self.ad_var[timeVarNm] = self.ad_var[timeVarNm] or
         Adios.define_variable(self.ad_io, timeVarNm, Adios.type_double, 1, self.shape1d, self.start1d, self.shape1d, true)

   self.shape2d[1], self.shape2d[2] = dynVecIn:size(), dynVecIn:numComponents()
   self.ad_var[dataVarNm] = self.ad_var[dataVarNm] or
         Adios.define_variable(self.ad_io, dataVarNm, Adios.type_double, 2, self.shape2d, self.start2d, self.shape2d, true)

   local _ = Adios.put(ad_engine, self.ad_var[timeVarNm], dynVecIn:timeMesh():data(), Adios.mode_deferred)

   -- Copy data to IO buffer. This buffer should be unnecessary, but presently it is needed
   -- because of the extra component in the dynvector used to support base-1 indexing.
   self._outBuff:expand(dynVecIn:size()*dynVecIn:numComponents())
   dynVecIn:_copy_from_dynvector(self._outBuff)
   local _ = Adios.put(ad_engine, self.ad_var[dataVarNm], self._outBuff:data(), Adios.mode_deferred)

   local _ = Adios.close(ad_engine)
end

function AdiosDynVectorIo:read(dynVecIn, fName)

   local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.

   -- Create an adios2_object if needed.
   self.ad_io = self.ad_io or Adios.declare_io(self.ad, fName)

   -- Open file to read.
   local ad_engine = Adios.open_new_comm(self.ad_io, fullNm, Adios.mode_readRandomAccess, self._parentComm)

   local ad_var_time = Adios.inquire_variable(self.ad_io, "TimeMesh")
   local ad_var_data = Adios.inquire_variable(self.ad_io, "Data")
   local varCnt = -1
   while (ad_var_time==nil and ad_var_data==nil) do
      varCnt = varCnt+1
      ad_var_time = Adios.inquire_variable(self.ad_io, "TimeMesh" .. varCnt)
      ad_var_data = Adios.inquire_variable(self.ad_io, "Data" .. varCnt)
   end

   local varShape1d = Adios.variable_shape(ad_var_time)
   self.shape1d[1] = varShape1d[1]
   dynVecIn:timeMesh():expand(varShape1d[1])
   local ad_err = Adios.set_selection(ad_var_time, 1, self.start1d, self.shape1d)
   local ad_err = Adios.get(ad_engine, ad_var_time, dynVecIn:timeMesh():data(), Adios.mode_sync)

   -- Read data into IO buffer. This buffer should be unnecessary, but presently it is needed
   -- because of the extra component in the dynvector used to support base-1 indexing.
   local inBuff = self._allocator(1)
   local varShape2d = Adios.variable_shape(ad_var_data)
   for d=1,2 do self.shape2d[d] = varShape2d[d] end
   local ad_err = Adios.set_selection(ad_var_data, 2, self.start2d, self.shape2d)
   inBuff:expand(varShape2d[1]*varShape2d[2])
   dynVecIn:data():expand(varShape2d[1])
   local ad_err = Adios.get(ad_engine, ad_var_data, inBuff:data(), Adios.mode_sync)

   varCnt = varCnt+1
   local ad_var_time = Adios.inquire_variable(self.ad_io, "TimeMesh" .. varCnt)
   local ad_var_data = Adios.inquire_variable(self.ad_io, "Data" .. varCnt)
   while (ad_var_time~=nil and ad_var_data~=nil) do

      local varShape1d = Adios.variable_shape(ad_var_time)
      self.shape1d[1] = varShape1d[1]
      local currSz = dynVecIn:size()
      dynVecIn:timeMesh():expand(currSz+varShape1d[1])
      local ad_err = Adios.set_selection(ad_var_time, 1, self.start1d, self.shape1d)
      local ad_err = Adios.get(ad_engine, ad_var_time, dynVecIn:timeMesh():data()+currSz, Adios.mode_sync)

      -- Read data into IO buffer. This buffer should be unnecessary, but presently it is needed
      -- because of the extra component in the dynvector used to support base-1 indexing.
      local varShape2d = Adios.variable_shape(ad_var_data)
      for d=1,2 do self.shape2d[d] = varShape2d[d] end
      local ad_err = Adios.set_selection(ad_var_data, 2, self.start2d, self.shape2d)
      inBuff:expand((currSz+varShape2d[1])*varShape2d[2])
      dynVecIn:data():expand(currSz+varShape2d[1])
      local ad_err = Adios.get(ad_engine, ad_var_data, inBuff:data()+currSz*varShape2d[2], Adios.mode_sync)

      varCnt = varCnt+1
      ad_var_time = Adios.inquire_variable(self.ad_io, "TimeMesh" .. varCnt)
      ad_var_data = Adios.inquire_variable(self.ad_io, "Data" .. varCnt)
   end

   dynVecIn:_copy_to_dynvector(inBuff)

   local _ = Adios.close(ad_engine)
end

return AdiosDynVectorIo
