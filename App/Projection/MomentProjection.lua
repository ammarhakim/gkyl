-- Gkyl ------------------------------------------------------------------------
--
-- App support code: MomentProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ProjectionBase   = require "App.Projection.ProjectionBase"
local Proto            = require "Lib.Proto"
local Updater          = require "Updater"
local xsys             = require "xsys"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"

-- Shell class for FV multimoment model projections.
local MomentProjection = Proto(ProjectionBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function MomentProjection:init(tbl)
   self.tbl = tbl
end

function MomentProjection:fullInit(species)
   self.species = species

   self.basis    = species.basis
   self.confGrid = species.grid

   self.cdim = self.confGrid:ndim()
end

----------------------------------------------------------------------
-- Base class for projection of arbitrary function. No re-scaling.
local FunctionProjection = Proto(MomentProjection)

function FunctionProjection:fullInit(species)
   FunctionProjection.super.fullInit(self, species)

   local func = self.tbl.func
   if not func then
      func = self.tbl[1]
   end
   assert(func, "FunctionProjection: Must specify the function")
   assert(type(func) == "function",
	  "The input must be a table containing function")
   self.project = Updater.ProjectOnBasis {
      onGrid   = self.confGrid,
      basis    = self.basis,
      evaluate = func,
      onGhosts = true
   }
   self.initFunc = func
end

function FunctionProjection:advance(time, inFlds, outFlds)
   local distf = outFlds[1]
   self.project:advance(time, {}, {distf})
end
----------------------------------------------------------------------
-- Base class for reading an initial condition.
local ReadInput = Proto(MomentProjection)

function ReadInput:fullInit(species)
   ReadInput.super.fullInit(self, species)

   self.userInputFile = self.tbl.inputFile
end
function ReadInput:run(t, distf)
   self.momIoRead = AdiosCartFieldIo {
      elemType = distf:elemType(),
      method   = "MPI",
      metaData = {polyOrder = self.basis:polyOrder(),
                  basisType = self.basis:id()},
   }

   local readSkin = false
   local tm, fr     = self.momIoRead:read(distf, string.format("%s",self.userInputFile), readSkin)
end
----------------------------------------------------------------------
return {
   FunctionProjection = FunctionProjection,
   ReadInput          = ReadInput,
}
