-- Gkyl ------------------------------------------------------------------------
--
-- Interface to ADIOS read API
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Mpi = require "Comm.Mpi"
local Adios = require "Io.Adios"
local ffi = require "ffi"
local Proto = require "Lib.Proto"

-- interface hook into the ADIOS read API
local Reader = Proto()

function Reader:init(fName, comm)
   self.fName = fName
   if not comm then comm = Mpi.COMM_WORLD end
   self.fd = Adios.read_open_file(fName, comm)
   assert(fd != nil, string.format(
	     "Unable to open ADIOS file %s for reading", fName))
end

return { Reader = Reader }
