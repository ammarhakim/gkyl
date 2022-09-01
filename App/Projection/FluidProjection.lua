-- Gkyl ------------------------------------------------------------------------
--
-- App support code: FluidProjection object.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ProjectionBase   = require "App.Projection.ProjectionBase"
local Proto            = require "Lib.Proto"
local Updater          = require "Updater"
local xsys             = require "xsys"
local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"

-- Shell class for fluid projections.
local FluidProjection = Proto(ProjectionBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function FluidProjection:init(tbl)
   self.tbl = tbl
end

function FluidProjection:fullInit(species)
   self.species = species

   self.basis = species.basis
   self.grid  = species.grid

   self.ndim = self.grid:ndim()

   self.mass = species:getMass()

   self.fromFile = self.tbl.fromFile

   self.weakMultiply = Updater.CartFieldBinOp {
      weakBasis = self.basis,  operation = "Multiply",
      onGhosts  = true,
   }
end

----------------------------------------------------------------------
-- Base class for projection of arbitrary function. No re-scaling.
local FunctionProjection = Proto(FluidProjection)

function FunctionProjection:fullInit(species)
   FunctionProjection.super.fullInit(self, species)

   local func = self.tbl.func
   if not func then func = self.tbl[1] end

   assert(func, "FunctionProjection: Must specify the function")
   assert(type(func) == "function", "The input must be a table containing function")

   if self.fromFile then
      self.ioMethod  = "MPI"
      self.writeSkin = true
      self.fieldIo = AdiosCartFieldIo {
         elemType  = species.moments[1]:elemType(),
         method    = self.ioMethod,
         writeSkin = self.writeSkin,
         metaData  = {polyOrder = self.basis:polyOrder(),
                      basisType = self.basis:id()}
      }
   else
      self.project = Updater.ProjectOnBasis {
         onGrid = self.grid,   evaluate = func,
         basis  = self.basis,  onGhosts = true
      }
   end
   self.initFunc = func
end

function FunctionProjection:advance(time, inFlds, outFlds)
   local extField = inFlds[1]
   local momOut   = outFlds[1]
   if self.fromFile then
      local tm, fr = self.fieldIo:read(momOut, self.fromFile)
   else
      self.project:advance(time, {}, {momOut})

      local jacob = extField.geo.jacobGeo
      if jacob then self.weakMultiply:advance(0, {momOut, jacob}, {momOut}) end
   end
end
----------------------------------------------------------------------
-- Base class for reading an initial condition.
local ReadInput = Proto(FluidProjection)

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
   local tm, fr   = self.momIoRead:read(distf, string.format("%s",self.userInputFile), readSkin)
end
----------------------------------------------------------------------
return {
   FluidProjection    = FluidProjection,
   FunctionProjection = FunctionProjection,
   ReadInput          = ReadInput,
}
