-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for ADIOS
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- don't bother if ADIOS is not built in
assert(GKYL_HAVE_ADIOS, "Gkyl was not built with ADIOS!")

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

require "Io._AdiosCdef" -- load FFI binding
require "Io._AdiosReadCdef" -- load FFI binding

local _M = {}

-- constants from enum ADIOS_FLAG
_M.flag_yes = ffi.C.adios_flag_yes
_M.flag_no = ffi.C.adios_flag_no
-- constants from enum ADIOS_DATATYPES
_M.unknown = ffi.C.adios_unknown
_M.byte = ffi.C.adios_byte
_M.short = ffi.C.adios_short
_M.integer = ffi.C.adios_integer
_M.long = ffi.C.adios_long
_M.unsigned_byte = ffi.C.adios_unsigned_byte
_M.unsigned_short = ffi.C.adios_unsigned_short
_M.unsigned_integer = ffi.C.adios_unsigned_integer
_M.unsigned_long = ffi.C.adios_unsigned_long
_M.real = ffi.C.adios_real
_M.double = ffi.C.adios_double
_M.long_double = ffi.C.adios_long_double
_M.string = ffi.C.adios_string
_M.complex = ffi.C.adios_complex
_M.double_complex = ffi.C.adios_double_complex
_M.string_array = ffi.C.adios_string_array

-- ADIOS read methods from adios_read
_M.read_method_bp = ffi.C.ADIOS_READ_METHOD_BP

-- adios_init_noxml
function _M.init_noxml(comm)
   local err = ffi.C.adios_init_noxml(comm)
end
-- adios_set_max_buffer_size
function _M.set_max_buffer_size(szMB)
   local err = ffi.C.adios_set_max_buffer_size(szMB)
end

-- adios_declare_group
function _M.declare_group(name, time_index, stats)
   local id = new("int64_t[1]")
   local err = ffi.C.adios_declare_group(id, name, time_index, stats)
   return id
end
-- adios_select_method
function _M.select_method(group, method, parameters, base_path)
   local err = ffi.C.adios_select_method(group[0], method, parameters, base_path)
end
-- adios_define_var
function _M.define_var(group_id, name, path, typ, dimensions, global_dimensions, local_offsets)
   local varId = ffi.C.adios_define_var(
      group_id[0], name, path, typ, dimensions, global_dimensions, local_offsets)
   return varId
end
-- adios_define_attribute
function _M.define_attribute(group, name, path, typ, value, var)
   local err = ffi.C.adios_define_attribute(group[0], name, path, typ, value, var)
end
-- adios_define_attribute_byvalue
function _M.define_attribute_byvalue(group, name, path, typ, nelems, values)
   local err = ffi.C.adios_define_attribute_byvalue(
      group[0], name, path, typ, nelems, values)
end

-- adios_open
function _M.open(group_name, name, mode, comm)
   local id = new("int64_t[1]")
   local err = ffi.C.adios_open(id, group_name, name, mode, comm)
   return id
end
-- adios_group_size
function _M.group_size(fd, data_size)
   local tsz = new("uint64_t[1]")
   ffi.C.adios_group_size(fd[0], data_size, tsz)
   return tsz[0]
end

-- adios_write
function _M.write(fd, name, var)
   local err = ffi.C.adios_write(fd[0], name, var)
   return err
end

-- adios_read
function _M.read(fd, name, buff, buffSz)
   local err = ffi.C.adios_read(fd[0], name, buff, buffSz)
   return err
end

-- adios_close
function _M.close(fd)
   local err = ffi.C.adios_close(fd[0])
end
-- adios_finalize
function _M.finalize(rank)
   local err = ffi.C.adios_finalize(rank)
end

-- ADIOS v2 read API wrappers
function _M.read_open_file(name, comm)
   return ffi.C.adios_read_open_file(name, _M.read_method_bp, comm)
end

function _M.inq_var_byid(fd, i)
   return ffi.C.adios_inq_var_byid(fd, i)
end

function _M.type_size(vtype, vvalue)
   return ffi.C.adios_type_size(vtype, vvalue)
end

function _M.type_to_string(vtype)
   return ffi.string(ffi.C.adios_type_to_string(vtype))
end

return _M
