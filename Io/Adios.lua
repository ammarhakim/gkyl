-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for ADIOS
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- don't bother if ADIOS is not built in
assert(GKYL_HAVE_ADIOS, "Gkyl was not built with ADIOS!")

local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

require "Io._AdiosCdef" -- load FFI binding
require "Io._AdiosReadCdef" -- load FFI binding for read API

local _M = {}

-- constants from enum ADIOS_FLAG
_M.flag_yes = ffiC.adios_flag_yes
_M.flag_no = ffiC.adios_flag_no
-- constants from enum ADIOS_DATATYPES
_M.unknown = ffiC.adios_unknown
_M.byte = ffiC.adios_byte
_M.short = ffiC.adios_short
_M.integer = ffiC.adios_integer
_M.long = ffiC.adios_long
_M.unsigned_byte = ffiC.adios_unsigned_byte
_M.unsigned_short = ffiC.adios_unsigned_short
_M.unsigned_integer = ffiC.adios_unsigned_integer
_M.unsigned_long = ffiC.adios_unsigned_long
_M.real = ffiC.adios_real
_M.double = ffiC.adios_double
_M.long_double = ffiC.adios_long_double
_M.string = ffiC.adios_string
_M.complex = ffiC.adios_complex
_M.double_complex = ffiC.adios_double_complex
_M.string_array = ffiC.adios_string_array

-- ADIOS read methods from adios_read
_M.read_method_bp = ffiC.ADIOS_READ_METHOD_BP

-- adios_init_noxml
function _M.init_noxml(comm)
   local err = ffiC.adios_init_noxml(comm)
end
-- adios_set_max_buffer_size
function _M.set_max_buffer_size(szMB)
   local err = ffiC.adios_set_max_buffer_size(szMB)
end

-- adios_declare_group
function _M.declare_group(name, time_index, stats)
   local id = new("int64_t[1]")
   local err = ffiC.adios_declare_group(id, name, time_index, stats)
   return id
end
-- adios_select_method
function _M.select_method(group, method, parameters, base_path)
   local err = ffiC.adios_select_method(group[0], method, parameters, base_path)
end
-- adios_define_var
function _M.define_var(group_id, name, path, typ, dimensions, global_dimensions, local_offsets)
   local varId = ffiC.adios_define_var(
      group_id[0], name, path, typ, dimensions, global_dimensions, local_offsets)
   return varId
end
-- adios_define_attribute
function _M.define_attribute(group, name, path, typ, value, var)
   local err = ffiC.adios_define_attribute(group[0], name, path, typ, value, var)
end
-- adios_define_attribute_byvalue
function _M.define_attribute_byvalue(group, name, path, typ, nelems, values)
   local err = ffiC.adios_define_attribute_byvalue(
      group[0], name, path, typ, nelems, values)
end

-- adios_open
function _M.open(group_name, name, mode, comm)
   local id = new("int64_t[1]")
   local err = ffiC.adios_open(id, group_name, name, mode, comm)
   return id
end
-- adios_group_size
function _M.group_size(fd, data_size)
   local tsz = new("uint64_t[1]")
   ffiC.adios_group_size(fd[0], data_size, tsz)
   return tsz[0]
end

-- adios_write
function _M.write(fd, name, var)
   local err = ffiC.adios_write(fd[0], name, var)
   return err
end

-- adios_read
function _M.read(fd, name, buff, buffSz)
   local err = ffiC.adios_read(fd[0], name, buff, buffSz)
   return err
end

-- adios_close
function _M.close(fd)
   local err = ffiC.adios_close(fd[0])
end
-- adios_finalize
function _M.finalize(rank)
   local err = ffiC.adios_finalize(rank)
end

-- ADIOS v2 read API wrappers

-- adios_read_open_file
function _M.read_open_file(name, comm)
   return ffiC.adios_read_open_file(name, _M.read_method_bp, comm)
end

-- adios_read_close
function _M.read_close(fd)
   return ffiC.adios_read_close(fd)
end

-- adios_inq_var_byid
function _M.inq_var_byid(fd, i)
   return ffiC.adios_inq_var_byid(fd, i)
end

-- adios_type_size
function _M.type_size(vtype, vvalue)
   return ffiC.adios_type_size(vtype, vvalue)
end

-- adios_type_to_string
function _M.type_to_string(vtype)
   return ffi.string(ffiC.adios_type_to_string(vtype))
end

-- adios_selection_boundingbox
function _M.selection_boundingbox(ndim, start, count)
   return ffiC.adios_selection_boundingbox(ndim, start:data(), count:data())
end

-- adios_schedule_read_byid
function _M.schedule_read_byid(fp, sel, varid, from_steps, nsteps, data)
   return ffiC.adios_schedule_read_byid(fp, sel, varid, from_steps, nsteps, data:data())
end

-- adios_perform_reads
function _M.perform_reads(fp, blocking)
   return ffiC.adios_perform_reads(fp, blocking)
end

-- adios_get_attr_byid
function _M.get_attr_byid(fp, attrid) --> type, size, void* to data
   local typePtr, sizePtr = ffi.new("int[1]"), ffi.new("int[1]")
   local voidPtr = ffi.new(typeof("void *[1]"))
   ffiC.adios_get_attr_byid(fp, attrid, typePtr, sizePtr, voidPtr)
   return typePtr[0], sizePtr[0], voidPtr[0]
end

return _M
