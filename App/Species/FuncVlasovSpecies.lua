-- Gkyl ------------------------------------------------------------------------
--
-- A species object with fluid moments specified as functions.
-- 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local AdiosCartFieldIo = require "Io.AdiosCartFieldIo"
local Basis            = require "Basis"
local DataStruct       = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid             = require "Grid"
local LinearTrigger    = require "Lib.LinearTrigger"
local Proto            = require "Lib.Proto"
local SpeciesBase      = require "App.Species.SpeciesBase"
local Time             = require "Lib.Time"
local Updater          = require "Updater"
local xsys             = require "xsys"
local ffi              = require "ffi"

local FuncVlasovSpecies = Proto(SpeciesBase)

function FuncVlasovSpecies:init(tbl)
   FuncVlasovSpecies.super.init(self, tbl)
   self.tbl = tbl
end

function FuncVlasovSpecies:fullInit(appTbl)
   local tbl = self.tbl

   self.evolve    = xsys.pickBool(tbl.evolve, true) -- by default evolve field
   self.confBasis = nil -- Will be set later
   self.charge    = tbl.charge and tbl.charge or 1.0
   self.mass      = tbl.mass and tbl.mass or 1.0   

   self.momDenFunc = tbl.momentumDensity
   self.vdim       = #{ self.momDenFunc(0.0, appTbl.lower) }
end

function FuncVlasovSpecies:getNdim()
   return self.confGrid:ndim()
end
function FuncVlasovSpecies:vdim()
   return self.vdim
end
function FuncVlasovSpecies:setName(nm)
   self.name = nm
end
function FuncVlasovSpecies:setIoMethod(ioMethod)
   self.ioMethod = ioMethod
end
function FuncVlasovSpecies:setConfBasis(basis)
   self.confBasis = basis
end
function FuncVlasovSpecies:setConfGrid(cgrid)
   self.confGrid = cgrid
end

function FuncVlasovSpecies:allocVectorMoment(dim)
   return DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis()*dim,
      ghost         = {1, 1}
   }
end

function FuncVlasovSpecies:alloc(nRkDup)
   self.momDensity = self:allocVectorMoment(self.vdim)
   self.dt          = ffi.new("double[2]")
   self.dtGlobal    = ffi.new("double[2]")
end

function FuncVlasovSpecies:allocMomCouplingFields()
   return { currentDensity = self:allocVectorMoment(self.vdim) }
end

function FuncVlasovSpecies:createSolver()
   self.momDensitySlvr = Updater.ProjectOnBasis {
      onGrid   = self.confGrid,
      basis    = self.confBasis,
      evaluate = self.momDenFunc,
   }
end

function FuncVlasovSpecies:createDiagnostics()
end

function FuncVlasovSpecies:initDist()
   self.momDensitySlvr:advance(0.0, {}, {self.momDensity})
end

function FuncVlasovSpecies:calcCouplingMoments(tCurr, rkIdx)
   if self.evolve then
      self.momDensitySlvr:advance(tCurr, {}, {self.momDensity})
   end
end

function FuncVlasovSpecies:getCharge() return self.charge end
function FuncVlasovSpecies:getMass() return self.mass end
function FuncVlasovSpecies:getMomDensity() return self.momDensity end

function FuncVlasovSpecies:rkStepperFields()
   return { nil, nil, nil }
end

function FuncVlasovSpecies:totalSolverTime()
   return self.momDensitySlvr.totalTime
end

return FuncVlasovSpecies
