local SourceBase     = require "App.Sources.SourceBase"
local DataStruct     = require "DataStruct"
local ffi            = require "ffi"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"

local VmSourceReplenish = Proto(SourceBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmSourceReplenish:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmSourceReplenish:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.
   self.sourceSpecies = assert(tbl.sourceSpecies, "App.VmSource: Must specify names of species to in 'sourceSpecies'.")
   self.sourceLength = assert(tbl.sourceLength, "App.VmSource: Must specify names of species to in 'sourceLength'.")

   self.tmEvalSrc = 0.0
end

function VmSourceReplenish:setConfBasis(basis)
   self.confBasis = basis
end
function VmSourceReplenish:setConfGrid(grid)
   self.confGrid = grid
end

function VmSourceReplenish:advance(tCurr, fIn, species, fRhsOut)
   local tmEvalSourceStart = Time.clock()
   local localEdgeFlux = ffi.new("double[3]")
   localEdgeFlux[0] = 0.0
   localEdgeFlux[1] = 0.0
   localEdgeFlux[2] = 0.0
   for sInd, otherNm in ipairs(self.sourceSpecies) do
      local numConfDims = self.confBasis:ndim()
      assert(numConfDims==1, "VlasovSpecies: The steady state source is available only for 1X.")
      local numConfBasis = self.confBasis:numBasis()
      local lower, upper = Lin.Vec(numConfDims), Lin.Vec(numConfDims)
      lower[1], upper[1] = -1.0, 1.0
      local basisUpper = Lin.Vec(numConfBasis)
      local basisLower = Lin.Vec(numConfBasis)

      self.confBasis:evalBasis(upper, basisUpper)
      self.confBasis:evalBasis(lower, basisLower)

      local flux = species[otherNm]:fluidMoments()[2]
      local fluxIndexer, fluxItr = flux:genIndexer(), flux:get(1)
      for idx in flux:localRangeIter() do
         if idx[1] == self.confGrid:numCells(1) then
            flux:fill(fluxIndexer(idx), fluxItr)
            for k = 1, numConfBasis do
               localEdgeFlux[0] = localEdgeFlux[0] + fluxItr[k]*basisUpper[k]
            end
         elseif idx[1] == 1 then
            flux:fill(fluxIndexer(idx), fluxItr)
            for k = 1, numConfBasis do
               localEdgeFlux[0] = localEdgeFlux[0] - fluxItr[k]*basisLower[k]
            end
         end
      end
   end
   local globalEdgeFlux = ffi.new("double[3]")
   Mpi.Allreduce(localEdgeFlux, globalEdgeFlux, 1, Mpi.DOUBLE, Mpi.MAX, self.confGrid:commSet().comm)
   local densFactor = globalEdgeFlux[0]/self.sourceLength
   fRhsOut:accumulate(densFactor, species[self.speciesName].fSource)
   local tmEvalMomStart = Time.clock()
   self.tmEvalSrc = self.tmEvalSrc + Time.clock() - tmEvalSrcStart
end

function VmSourceReplenish:srcTime()
   return self.tmEvalSource
end

return VmSourceReplenish
