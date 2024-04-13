-- Gkyl ------------------------------------------------------------------------
--
-- Plasma solver on a Cartesian grid. Works in arbitrary CDIM/VDIM
-- (VDIM>=CDIM) with either Vlasov, gyrokinetic, fluids and Maxwell,
-- Poisson or specified EM fields.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Infrastructure loads
local DecompRegionCalc = require "Lib.CartDecomp"
local LinearTrigger = require "Lib.LinearTrigger"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local date = require "xsys.date"
local lfs = require "lfs"
local lume = require "Lib.lume"
local xsys = require "xsys"
local ffi = require "ffi"

