-- Gkyl ------------------------------------------------------------------------
--
-- Basic boundary condition for a fluid species, i.e. those that can be
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
--local ConstDiffusionModDecl = require "Eq.constDiffusionData.ConstDiffusionModDecl"

local FluidBasicBC = Proto(BCsBase)

-- Store table passed to it and defer construction to :fullInit().
function FluidBasicBC:init(tbl) self.tbl = tbl end

-- Function initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function FluidBasicBC:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.bcKind      = assert(tbl.kind, "FluidBasicBC: must specify the type of BC in 'kind'.")
   self.diagnostics = tbl.diagnostics or {}
   self.saveFlux    = tbl.saveFlux or false

   if self.bcKind == "dirichlet" or self.bcKind == "neumann" then
      self.bcValue = assert(tbl.value, "FluidBasicBC: must specify the BC value in 'value'.")
   end
end

function FluidBasicBC:setConfBasis(basis) self.basis = basis end
function FluidBasicBC:setConfGrid(grid) self.grid = grid end

function FluidBasicBC:bcDirichlet(dir, tm, idxIn, fIn, fOut)
   -- Impose f=fBC at the boundary.
   if (idxIn[dir] == 1) then
      self.constDiffDirichletBCs[dir][1](self.grid:dx(dir),fIn:data(),self.bcValue,fOut:data())
   else
      self.constDiffDirichletBCs[dir][2](self.grid:dx(dir),fIn:data(),self.bcValue,fOut:data())
   end
end

function FluidBasicBC:bcNeumann(dir, tm, idxIn, fIn, fOut)
   -- Impose f'=fpBC at the boundary.
   if (idxIn[dir] == 1) then
      self.constDiffNeumannBCs[dir][1](self.grid:dx(dir),fIn:data(),self.bcValue,fOut:data())
   else
      self.constDiffNeumannBCs[dir][2](self.grid:dx(dir),fIn:data(),self.bcValue,fOut:data())
   end
end

function FluidBasicBC:createSolver(mySpecies, field, externalField)
   self.nMoments = mySpecies.nMoments

   local bcFunc, skinType
   if self.bcKind == "copy" then
      bcFunc   = function(...) return self:bcCopy(...) end
      skinType = "pointwise"
   elseif self.bcKind == "absorb" then
      bcFunc   = function(...) return self:bcAbsorb(...) end
      skinType = "pointwise"
   elseif self.bcKind == "dirichlet" then
--      self.constDiffDirichletBCs = ConstDiffusionModDecl.selectBCs(self.basis:id(), self.basis:ndim(), 
--                                                                   self.basis:polyOrder(), 2, "Dirichlet")
--      bcFunc   = function(...) return self:bcDirichlet(...) end
--      skinType = "pointwise"
--
      assert(false, "FluidBasicBC: BC kind not available.")
   elseif self.bcKind == "neumann" then
--      self.constDiffNeumannBCs = ConstDiffusionModDecl.selectBCs(self.basis:id(), self.basis:ndim(), 
--                                                                 self.basis:polyOrder(), 2, "Neumann")
--      bcFunc   = function(...) return self:bcNeumann(...) end
--      skinType = "pointwise"
      assert(false, "FluidBasicBC: BC kind not available.")
   else
      assert(false, "FluidBasicBC: BC kind not recognized.")
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

function FluidBasicBC:advance(tCurr, mySpecies, field, externalField, inIdx, outIdx)
   local fIn = mySpecies:rkStepperFields()[outIdx]
   self.bcSolver:advance(tCurr, {}, {fIn})
end

-- ................... Classes meant as aliases to simplify input files ...................... --
local FluidAbsorbBC = Proto(FluidBasicBC)
function FluidAbsorbBC:fullInit(mySpecies)
   self.tbl.kind  = "absorb"
   FluidAbsorbBC.super.fullInit(self, mySpecies)
end

local FluidCopyBC = Proto(FluidBasicBC)
function FluidCopyBC:fullInit(mySpecies)
   self.tbl.kind  = "copy"
   FluidCopyBC.super.fullInit(self, mySpecies)
end

local FluidZeroFluxBC = Proto()
function FluidZeroFluxBC:init(tbl)
   self.tbl      = tbl
   self.tbl.kind = "zeroFlux"
end
-- ................... End of FluidBasicBC alias classes .................... --

return {FluidBasic    = FluidBasicBC,
        FluidAbsorb   = FluidAbsorbBC,
        FluidCopy     = FluidCopyBC,
        FluidZeroFlux = FluidZeroFluxBC,}

