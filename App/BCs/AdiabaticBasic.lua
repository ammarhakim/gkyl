-- Gkyl ------------------------------------------------------------------------
--
-- Basic boundary condition for a gyrofluid species, i.e. those that can be
-- applied with Updater/Bc.lua.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local BCsBase    = require "App.BCs.BCsBase"
local DataStruct = require "DataStruct"
local Updater    = require "Updater"
local Mpi        = require "Comm.Mpi"
local Proto      = require "Lib.Proto"
local Time       = require "Lib.Time"

local AdiabaticBasicBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function AdiabaticBasicBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function AdiabaticBasicBC:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.bcKind      = assert(tbl.kind, "AdiabaticBasicBC: must specify the type of BC in 'kind'.")
   self.diagnostics = tbl.diagnostics or {}
   self.saveFlux    = tbl.saveFlux or false
end

function AdiabaticBasicBC:bcAbsorb(dir, tm, idxIn, fIn, fOut)
   -- The idea is that by setting the plasma quantities to zero in the
   -- ghost cell nothing is transported into the domain, and whatever is transported
   -- out is lost. If setting them to exactly zero doesn't work it's likely that the
   -- sound speed & drift velocity diverge, so set them to something small.
   for i = 1, self.basis:numBasis() do fOut[i] = 1.e-10*fIn[i] end   -- Density.
end

function AdiabaticBasicBC:createSolver(mySpecies)
   local bcFunc, skinType
   if self.bcKind == "copy" then
      bcFunc   = function(...) return self:bcCopy(...) end
      skinType = "pointwise"
   elseif self.bcKind == "absorb" then
      bcFunc   = function(...) return self:bcAbsorb(...) end
      skinType = "pointwise"
   end

   self.bcSolver = Updater.Bc {
      onGrid             = self.grid,
      cdim               = self.grid:ndim(),
      dir                = self.bcDir,
      edge               = self.bcEdge,
      boundaryConditions = {bcFunc},
      skinLoop           = skinType,
   }
end

function AdiabaticBasicBC:advance(tCurr, species, inFlds)
   self.bcSolver:advance(tCurr, {}, inFlds)
end

return AdiabaticBasicBC
