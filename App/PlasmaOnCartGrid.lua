-- Gkyl ------------------------------------------------------------------------
--
-- Common entry point into various Gkeyll Apps.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Adios = require "Io.Adios"
local Mpi = require "Comm.Mpi"

return {
   -- Return function to create Moments app objects
   Moments = function()
      local moments = require("App.Moments")
      return moments
   end
}
