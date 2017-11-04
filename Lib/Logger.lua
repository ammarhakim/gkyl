-- Gkyl ------------------------------------------------------------------------
--
-- Simple logger for use in sims
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Mpi = require "Comm.Mpi"

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- write message to list of output streams
local function writeToFile(outStreams, msg)
   for _, s in ipairs(outStreams) do
      s:write(msg); s:write("\n")
   end
end

-- create logger
local function makeLogger(tbl)
   -- set write rank
   local writeRank = 0
   if tbl.writeRank then writeRank = tbl.writeRank end
   
   -- get hold of communicator if specified or use MPI_COMM_WORLD
   local comm = tbl.comm and tbl.comm or Mpi.COMM_WORLD
   local rank = Mpi.Comm_rank(comm) -- local rank

   local outStreams = {} -- list of streams to which logs are written
   if writeRank == rank then
      outStreams[1] = io.stdout -- always write to console
      -- open log file if needed
      if tbl.logToFile then
	 outStreams[2] = io.open(GKYL_OUT_PREFIX .. "_" .. writeRank .. ".log", "w")
      end
   end

   -- return function to do logging
   return function (msg)
      writeToFile(outStreams, msg)
   end
end

return makeLogger
