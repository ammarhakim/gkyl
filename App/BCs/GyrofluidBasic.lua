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

local GyrofluidBasicBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function GyrofluidBasicBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GyrofluidBasicBC:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.bcKind      = assert(tbl.kind, "GyrofluidBasicBC: must specify the type of BC in 'kind'.")
   self.diagnostics = tbl.diagnostics or {}
   self.saveFlux    = tbl.saveFlux or false 
end

function GyrofluidBasicBC:bcCopy(dir, tm, idxIn, fIn, fOut)
   for i = 1, self.nMoments*self.basis:numBasis() do fOut[i] = fIn[i] end
end

function GyrofluidBasicBC:bcAbsorb(dir, tm, idxIn, fIn, fOut)
   -- The idea is that by setting the plasma quantities to zero in the
   -- ghost cell nothing is transported into the domain, and whatever is transported
   -- out is lost. We can't set them to exactly zero or else the sound speed
   -- and drift velocity would diverge, so we set them to something small.
   local numB = self.basis:numBasis()
   for i = 1, numB do fOut[0*numB+i] = 1.e-10*fIn[0*numB+i] end   -- Mass density.
   for i = 1, numB do fOut[1*numB+i] = 0. end                     -- Momentum density.
   for i = 1, numB do fOut[2*numB+i] = 1.e-10*fIn[2*numB+i] end   -- Energy density.
   for i = 1, numB do fOut[3*numB+i] = 1.e-10*fIn[3*numB+i] end   -- Perpendicular pressure (divided by B).
end

function GyrofluidBasicBC:createSolver(mySpecies)
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

function GyrofluidBasicBC:advance(tCurr, species, inFlds)
   self.bcSolver:advance(tCurr, {}, inFlds)
end

return GyrofluidBasicBC
