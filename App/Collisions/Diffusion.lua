-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: a (hyper)diffusion term of the form
--   partial_i^n( D_i(x) partial_i^n(f) )
-- where i=1,2,3 depending on cdim, n=diffOrder/2, and D_i(x) is a
-- diffusion coefficient that (can) depend on space and that is different
-- along each direction.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local DataStruct     = require "DataStruct"
local ZeroArray      = require "DataStruct.ZeroArray"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local lume           = require "Lib.lume"
local Updater        = require "Updater"

-- Diffusion ---------------------------------------------------------------
--
-- Add a diffusion term to the right side of an equation.
--------------------------------------------------------------------------------

local Diffusion = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function Diffusion:init(tbl) self.tbl = tbl end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function Diffusion:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.speciesKind = self.speciesKind or tbl.speciesKind -- Should be 'vlasov' or 'gyokinetic'. Other species to be added later.

   -- Diffusion coefficient. This can be a single number, or a table of numbers and/or functions:
   --   1) Single number: applies that diff coefficient along all directions.
   --   2) Table of numbers: must have the same number elements as diffusiveDirs, and applies those diffusion 
   --                        coefficients along the directions indicated in diffusiveDirs.
   --   3) Functions: must be scalar functions with (t,xn) signature that specify the diffusion
   --                 coefficients along the directions indicated in diffusiveDirs.
   self.diffCoeff = assert(tbl.coefficient,
      "App.Diffusion: Must specify the diffusion coefficient (vector) in 'coefficient'.")

   self.diffDirs  = tbl.diffusiveDirs         -- Directions in which to apply diffusion.
   self.diffOrder = tbl.order                 -- Read the diffusion operator order.

   self.timers = {advance = 0.,}
end

function Diffusion:setSpeciesName(nm)   self.speciesName = nm end
function Diffusion:setName(nm)          self.name        = self.speciesName.."_"..nm end
function Diffusion:setConfBasis(basis)  self.confBasis   = basis end
function Diffusion:setConfGrid(grid)    self.confGrid    = grid end
function Diffusion:setPhaseBasis(basis) self.phaseBasis  = basis end
function Diffusion:setPhaseGrid(grid)   self.phaseGrid   = grid end

function Diffusion:createSolver(mySpecies, externalField)

   local cdim = self.confGrid:ndim()
   local diffCoeffType = type(self.diffCoeff)
   if self.diffDirs then
      assert(#self.diffDirs<=cdim, "App.Diffusion: 'diffusiveDirs' cannot have more entries than the simulation's conf dimensions.")
      for d = 1, #self.diffDirs do assert(self.diffDirs[d]<=cdim, "App.Diffusion: cannot apply diffusion in velocity space.") end
   else
      if diffCoeffType == "table" then
         assert(#self.diffCoeff==cdim, "App.Diffusion: 'coefficient' must have same number of entries as there are dimensions if not specifying 'diffusiveDirs'.")
      end
      -- Apply diffusion in all directions.
      self.diffDirs = {}
      for d = 1, cdim do self.diffDirs[d] = d end
   end

   local diffConfRange
   local isVarCoeff = false
   if diffCoeffType == "number" then
      -- Set the diffusion coefficient to the same amplitude in all directions.
      self.coefficient = {
         _zero = ZeroArray.Array(ZeroArray.double, cdim, 1, GKYL_USE_GPU and 1 or 0),
      }
      self.coefficient._zeroDevice = self.coefficient._zero
      self.coefficient._zero:clear(0.)
      for d = 1, #self.diffDirs do self.coefficient._zero:shiftc(self.diffCoeff, self.diffDirs[d]-1) end

      diffConfRange = self.confGrid:localRange() -- Not used.

   elseif diffCoeffType == "table" then

      if self.tbl.diffusiveDirs then
         assert(#self.diffCoeff==#self.tbl.diffusiveDirs, "App.Diffusion: 'coefficient' table must have the same number of entries as 'diffusiveDirs'.")
      else
         assert(#self.diffCoeff==cdim, "App.Diffusion: 'coefficient' table must have the same number of entries as the simulation's conf dimensions if 'diffusiveDirs' is not given.")
      end

      -- Check if any of the diffusion coefficients are functions. If yes, project
      -- the diffusion coefficient vector onto the basis.
      isVarCoeff = lume.any(self.diffCoeff, function(e) return type(e) == "function" end)
      if isVarCoeff == false then
         -- All entries should be numbers, so we'll use constant diffusion.
         for d = 1, #self.diffDirs do assert(type(self.diffCoeff[d]) == "number", "App.Diffusion: 'coefficient' entries can only be numbers or functions.") end

         self.coefficient = {
            _zero = ZeroArray.Array(ZeroArray.double, cdim, 1, GKYL_USE_GPU and 1 or 0),
         }
         self.coefficient._zeroDevice = self.coefficient._zero
         self.coefficient._zero:clear(0.)
         for d = 1, #self.diffDirs do self.coefficient._zero:shiftc(self.diffCoeff[d], self.diffDirs[d]-1) end

         diffConfRange = self.confGrid:localRange() -- Not used.
      else
         -- Some entries must be functions. Turn number entries into functions and use varying diffusion.

         local coeffFuncs = {}
         for d = 1, cdim do coeffFuncs[d] = function(t, xn) return 0. end end
         for d = 1, #self.diffDirs do 
            if type(self.diffCoeff[d]) == "number" then
               coeffFuncs[self.diffDirs[d]] = function(t, zn) return self.diffCoeff[d] end
            else
               coeffFuncs[self.diffDirs[d]] = self.diffCoeff[d]
            end
         end
         local diffCoeffFunc
         if cdim==1 then
            diffCoeffFunc = function(t,xn) return coeffFuncs[1](t,xn) end
         elseif cdim==2 then
            diffCoeffFunc = function(t,xn) return coeffFuncs[1](t,xn), coeffFuncs[2](t,xn) end
         elseif cdim==3 then
            diffCoeffFunc = function(t,xn) return coeffFuncs[1](t,xn), coeffFuncs[2](t,xn), coeffFuncs[3](t,xn) end
         end
         self.coefficient = DataStruct.Field {
            onGrid = self.confGrid,  numComponents = cdim*self.confBasis:numBasis(),
            ghost  = {1, 1},         metaData      = {polyOrder = self.confBasis:polyOrder(),
                                                      basisType = self.confBasis:id(),},
         }
         local projectDiffCoeff = Updater.EvalOnNodes {
            onGrid = self.confGrid,   evaluate = diffCoeffFunc,
            basis  = self.confBasis,  onGhosts = false,
         }
         projectDiffCoeff:advance(0.0, {}, {self.coefficient})

         diffConfRange = self.coefficient:localRange() -- Used to index coefficient.
      end
   end

   -- When using mapped coordinates (default for gyrokinetics), scale the distribution
   -- function by the reciprocal of the configuration-space Jacobian.
   self.fNoJacobGeo = nil
   self.divByJacobGeo = function(tm, fIn, fOut) end
   if externalField.geo then
      if externalField.geo.jacobGeo and externalField.geo.jacobGeoInv then
         assert(isVarCoeff, "App.Diffusion: diffusion coefficient must be a function if using mapped coordinates.") 
         self.fNoJacobGeo = mySpecies:allocDistf()
         self.jacobGeoInv = externalField.geo.jacobGeoInv
         self.confPhaseWeakMultiply = Updater.CartFieldBinOp {
            weakBasis  = self.phaseBasis,  operation = "Multiply",
            fieldBasis = self.confBasis,   onGhosts = true,
         }
         -- Multiply the diffusion coefficient by the jacobian.
         local confWeakMultiply = mySpecies.confWeakMultiply
         confWeakMultiply:advance(0., {self.coefficient, externalField.geo.jacobGeo}, {self.coefficient})
         self.coefficient:sync()

         self.divByJacobGeo = function(tm, fIn, fOut)
            self.confPhaseWeakMultiply:advance(tm, {fIn, self.jacobGeoInv}, {fOut})
         end
      end
   end

   local grid  = self.phaseGrid and self.phaseGrid  or self.confGrid
   local basis = self.phaseGrid and self.phaseBasis or self.confBasis

   -- Diffusion equation solver.
   if self.speciesKind == "vlasov" then
      self.collisionSlvr = Updater.VlasovDiffusion {
         onGrid    = grid,            constDiff    = not isVarCoeff,
         onBasis   = basis,           directions   = self.diffDirs,
         confBasis = self.confBasis,  order        = self.diffOrder,
         confRange = diffConfRange,   zeroFluxDirs = mySpecies.zeroFluxDirections,
      }
   elseif self.speciesKind == "gyrokinetic" then
      self.collisionSlvr = Updater.GyrokineticDiffusion {
         onGrid    = grid,            constDiff    = not isVarCoeff,
         onBasis   = basis,           directions   = self.diffDirs,
         confBasis = self.confBasis,  order        = self.diffOrder,
         confRange = diffConfRange,   zeroFluxDirs = mySpecies.zeroFluxDirections,
      }
   end
end

function Diffusion:advance(tCurr, fIn, population, outFlds)
   local tmStart = Time.clock()

   local fRhsOut, cflRateByCell = outFlds[1], outFlds[2]

   self.fNoJacobGeo = self.fNoJacobGeo or fIn
   self.divByJacobGeo(tCurr, fIn, self.fNoJacobGeo)
   self.collisionSlvr:advance(tCurr, {self.fNoJacobGeo, self.coefficient}, {fRhsOut, cflRateByCell})

   self.timers.advance = self.timers.advance + Time.clock() - tmStart
end

function Diffusion:write(tm, frame) end

-- ................... Classes meant as aliases to simplify input files ...................... --

local DiffusionVlasov = Proto(Diffusion)
function DiffusionVlasov:fullInit(mySpecies)
   self.speciesKind = "vlasov"
   DiffusionVlasov.super.fullInit(self, mySpecies)
end

local DiffusionGyrokinetic = Proto(Diffusion)
function DiffusionGyrokinetic:fullInit(mySpecies)
   self.speciesKind = "gyrokinetic"
   DiffusionGyrokinetic.super.fullInit(self, mySpecies)
end

-- ................... End of aliases ...................... --

return {Diffusion            = Diffusion,
        DiffusionVlasov      = DiffusionVlasov,
        DiffusionGyrokinetic = DiffusionGyrokinetic,}
