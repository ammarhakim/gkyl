-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Gyrokinetic LB Collision operator
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local DataStruct     = require "DataStruct"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local GkLBOconstNuEq = require "Eq.GkLBO"
local xsys           = require "xsys"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"

-- GkLBOCollisions ---------------------------------------------------------------
--
-- Lenard-Bernstein Collision operator
-- Actually dates back to Lord Rayleigh, Philos. Mag. 32, 424 (1891).
--------------------------------------------------------------------------------

local GkLBOCollisions = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkLBOCollisions:init(tbl)
   self.tbl = tbl
end

-- Function to find the index of an element in table.
local function findInd(tbl, el)
   for i, v in ipairs(tbl) do
      if v == el then
         return i
      end
   end
   return #tbl+1    -- If not found return a number larger than the length of the table.
end


-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- previously stored table

   self.cfl = 0.0 -- Will be replaced.

   local collidingSpecies = assert(tbl.collideWith, "App.GkLBOCollisions: Must specify names of species to collide with in 'collideWith'.")

   -- First determine if self-species and/or cross-species collisions take place,
   -- and (if cross-collisions=true) put the names of the other colliding species in a list.
   local selfSpecInd = findInd(collidingSpecies, self.speciesName)
   if selfSpecInd < (#collidingSpecies+1) then
      self.selfCollisions = true                 -- Apply self-species collisions.
      if #collidingSpecies > 1 then
         self.crossCollisions = true             -- Apply cross-species collisions.
         self.crossSpecies    = collidingSpecies
         table.remove(self.crossSpecies, selfSpecInd)
      else
         self.crossCollisions = false            -- Don't apply cross-species collisions.
      end
   else
      self.selfCollisions  = false               -- Don't apply self-species collisions.
      self.crossCollisions = true                -- Apply cross-species collisions.
      self.crossSpecies    = collidingSpecies    -- All species in collidingSpecies must be cross-species.
   end
   
   -- Now establish if user wants constant or spatially varying collisionality.
   -- For constant nu, separate self and cross collision frequencies.
   local collFreqs          = tbl.frequencies -- List of collision frequencies, if using spatially constant nu. 
   if collFreqs then
      self.varNu            = false    -- Not spatially varying nu.
      self.cellConstNu      = true     -- Cell-wise constant nu?
      if self.selfCollisions then
         self.collFreqSelf  = collFreqs[selfSpecInd]
      end
      if self.crossCollisions then
         self.collFreqCross = collFreqs
         table.remove(self.collFreqCross, selfSpecInd)
      end
   else
      --local normNuIn      = tbl.normNu
      local normNuIn      = assert(tbl.normNu, "App.GkLBOCollisions: No constant collision frequencies provided ('frequencies'). For spatially varying nu must specify 'normNu', list of collisionalities normalized by (T_0^(3/2)/n_0) evaluated somewhere in the simulation.") 
      if normNuIn then
         self.varNu       = true    -- Spatially varying nu.
         -- Below normNu is the collisionality normalized by (T_0^(3/2)/n_0) evaluated somewhere in the simulation.
         if self.selfCollisions then
            self.normNuSelf  = normNuIn[selfSpecInd]
         end
         if self.crossCollisions then
            self.normNuCross = normNuIn
            table.remove(self.normNuCross, selfSpecInd)
         end
         -- For now only cell-wise constant nu is implemented.
         -- self.cellConstNu = assert(tbl.cellAvFrequencies, "App.GkLBOCollisions: Must specify 'useCellAverageNu=true/false' for using cellwise constant/expanded spatially varying collisionality.")
         self.cellConstNu = true
      end
   end

   if self.crossCollisions then
      -- IMPORTANT: First implement the option in Eric Shi's 2017 PhD thesis (HeavyIons).
      --            For now assume constant collisionality only.
      -- For now, crossMomOp=HeavyIons is the only option.
      -- self.crossMomOp  = assert(tbl.crossMomOp, "App.GkLBOCollisions: Must specify 'crossMomOp' (HeavyIons, Greene, GreeneSmallAngle, GreeneSmallAngleLimit), formulas used to calculate cross-species primitive moments.")
      self.crossMomOp  = "HeavyIons"
   end

   self.mass = speciesTbl.mass
   self.tmEvalMom = 0.0
end

function GkLBOCollisions:setName(nm)
   self.name = nm
end
function GkLBOCollisions:setSpeciesName(nm)
   self.speciesName = nm
end

function GkLBOCollisions:setCfl(cfl)
   self.cfl = cfl -- what should this be? - AHH
end
function GkLBOCollisions:setConfBasis(basis)
   self.confBasis = basis
end
function GkLBOCollisions:setConfGrid(grid)
   self.confGrid = grid
end
function GkLBOCollisions:setPhaseBasis(basis)
   self.phaseBasis = basis
end
function GkLBOCollisions:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function GkLBOCollisions:createSolver(funcField)
   self.vDim = self.phaseGrid:ndim() - self.confGrid:ndim()

   -- Number of physical velocity dimensions.
   self.vDimPhys = 1.0
   if self.vDim == 2 then
     self.vDimPhys = 3.0
   end

   -- Maximum velocity of the velocity grid (and its square).
   self.vParMax            = self.phaseGrid:upper(self.confGrid:ndim()+1)
   self.vParMaxSq          = self.vParMax^2

   -- Intemediate storage for output of collisions.
   self.collOut = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Inverse of background magnetic field.
   self.bmag    = funcField.geo.bmag
   -- Inverse of background magnetic field.
   self.bmagInv = funcField.geo.bmagInv
      
   -- Zero-flux BCs in the velocity dimensions.
   local zfd = { }
   for d = 1, self.vDim do
      zfd[d] = self.confGrid:ndim() + d
   end

   self.gkLBOconstNuCalcEq = {}
   if self.varNu then
      -- Collisionality, nu.
      self.collFreq = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Updater to compute spatially varying (Spitzer) nu.
      self.spitzerNu = Updater.SpitzerCollisionality {
         onGrid           = self.confGrid,
         confBasis        = self.confBasis,
         useCellAverageNu = self.cellConstNu,
      }
      -- Lenard-Bernstein equation.
      self.gkLBOconstNuCalcEq = GkLBOconstNuEq {
         phaseBasis       = self.phaseBasis,
         confBasis        = self.confBasis,
         vParUpper        = self.vParMax,
         mass             = self.mass,
         useCellAverageNu = self.cellConstNu,
      }
   else
      self.collFreq = 0.0    -- Assigned in advance method.
      -- Lenard-Bernstein equation.
      self.gkLBOconstNuCalcEq = GkLBOconstNuEq {
         phaseBasis = self.phaseBasis,
         confBasis  = self.confBasis,
         vParUpper  = self.vParMax,
         mass       = self.mass,
      }
   end
   self.collisionSlvr = Updater.HyperDisCont {
      onGrid             = self.phaseGrid,
      basis              = self.phaseBasis,
      cfl                = self.cfl,
      equation           = self.gkLBOconstNuCalcEq,
      updateDirections   = zfd, -- only update velocity directions
      zeroFluxDirections = zfd,
   }
   if self.selfCollisions then
      -- Parallel flow velocity.
      self.uPar = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Thermal speed squared, vth=sqrt(T/m).
      self.vthSq = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      self.primMomSelf = Updater.SelfPrimMoments {
         onGrid     = self.confGrid,
         phaseGrid  = self.phaseGrid,
         phaseBasis = self.phaseBasis,
         confBasis  = self.confBasis,
         gkfacs     = {self.mass, self.bmag},
         operator   = "GkLBO",
      }
   end
   if self.crossCollisions then
      -- Flow velocity of the other species (ions if crossMomOp=HeavyIons).
      self.uParOther = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Cross-species flow velocity.
      self.uParCross = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Cross species thermal speed squared.
      self.vthSqCross = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Mixed kinetic energy density.
      self.kinEnergyDens = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Thermal energy density.
      self.thermEnergyDens = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Weak binary operations.
      self.confDiv = Updater.CartFieldBinOp {
         onGrid    = self.confGrid,
         weakBasis = self.confBasis,
         operation = "Divide",
      }
      self.confMultiply = Updater.CartFieldBinOp {
         onGrid    = self.confGrid,
         weakBasis = self.confBasis,
         operation = "Multiply",
      }
      if self.crossMomOp ~= "HeavyIons" then
         -- Thermal speed squared of the other species.
         self.vthSqOther = DataStruct.Field {
            onGrid        = self.confGrid,
            numComponents = self.confBasis:numBasis(),
            ghost         = {1, 1},
         }
         self.primMomCross = Updater.CrossPrimMoments {
            onGrid     = self.confGrid,
            phaseBasis = self.phaseBasis,
            confBasis  = self.confBasis,
            operator   = "GkLBO",
            gkfacs     = {self.mass, self.bmag},
            formulas   = self.crossMomOp,
         }
      end
   end

   -- Number of cells in which number density was negative (somewhere).
   self.primMomLimitCrossingsL = DataStruct.DynVector {
      numComponents = 1,
   }
   self.primMomLimitCrossingsG = DataStruct.DynVector {
      numComponents = 1,
   }
end

function GkLBOCollisions:advance(tCurr, fIn, species, fRhsOut)
   local selfMom = species[self.speciesName]:fluidMoments()

   local tmEvalMomStart = Time.clock()
   if self.selfCollisions then
      -- Compute primitive moments velocity and vthSq=T/m from zeroth,
      -- first and second moments, and distribution function.
      self.primMomSelf:advance(0.0, {selfMom[1], selfMom[2], selfMom[3],fIn},
                                    {self.uPar, self.vthSq})
      self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

      self.gkLBOconstNuCalcEq.primMomCrossLimit = 0.0

      if self.varNu then
         -- Compute the collisionality.
         self.spitzerNu:advance(0.0, {self.mass, self.normNuSelf, selfMom[1], self.vthSq},{self.collFreq})
      else
         self.collFreq = self.collFreqSelf
      end
      -- Compute increment from collisions and accumulate it into output.
      self.collisionSlvr:advance(
         tCurr, {fIn, self.bmagInv, self.uPar, self.vthSq, self.collFreq}, {self.collOut})

      self.primMomLimitCrossingsG:appendData(tCurr, {0.0})
      self.primMomLimitCrossingsL:appendData(tCurr, {self.gkLBOconstNuCalcEq.primMomCrossLimit})

      fRhsOut:accumulate(1.0, self.collOut)
   end

   if self.crossCollisions then
      -- IMPORTANT: For now implement what was done in Eric Shi's PhD Thesis 2017 (HeavyIons).
      --            Assume ions & electrons only, and only include the effect on
      --            the electrons.
      if self.crossMomOp=="HeavyIons" then
         local tmEvalMomStart = Time.clock()
         -- Compute ion flow velocity.
         local ionMom = species[self.crossSpecies[1]]:fluidMoments()
         self.confDiv:advance(0., {ionMom[1], ionMom[2]}, {self.uParOther})
         self.uParCross:copy(self.uParOther)
         -- Compute the cross thermal speed according to:
         -- vDimPhys*M_{0,e} vthSq_{ei} = M_{2,e} - 2u_i dot M_{1,e} + u_i^2 M_{0,e}
         -- vDimPhys=1 if vDim=1, vDimPhys=3 if vDim=2.
         self.confMultiply:advance(0., {self.uParOther, selfMom[2]}, {self.kinEnergyDens})
         self.confMultiply:advance(0., {self.uParOther, self.uParCross}, {self.uParOther})
         self.confMultiply:advance(0., {self.uParOther, selfMom[1]}, {self.uParOther})
         self.thermEnergyDens:combine( 1.0/self.vDimPhys, selfMom[3],
                                      -2.0/self.vDimPhys, self.kinEnergyDens,
                                       1.0/self.vDimPhys, self.uParOther )
         self.confDiv:advance(0., {selfMom[1], self.thermEnergyDens}, {self.vthSqCross}) 
         self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

         if self.varNu then
            -- Compute the collisionality.
            if not self.selfCollisions then
               -- Need to compute the electron temperature.
               self.confDiv:advance(0., {selfMom[1], selfMom[2]}, {self.uPar})
               self.confMultiply:advance(0., {self.uPar, selfMom[2]}, {self.kinEnergyDens})
               self.thermEnergyDens:combine( 1.0/self.vDimPhys, selfMom[3],
                                            -1.0/self.vDimPhys, self.kinEnergyDens )
               self.confDiv:advance(0., {selfMom[1],self.thermEnergyDens}, {self.vthSq})
            end
            self.spitzerNu:advance(0.0, {self.mass, self.normNuCross[1], otherMom[1], self.vthSq}, {self.collFreq})
         else
            self.collFreq = self.collFreqCross[1]
         end
      end
      -- Compute increment from cross-species collisions and accumulate it into output.
      self.collisionSlvr:advance( tCurr, 
         {fIn, self.bmagInv, self.uParCross, self.vthSqCross, self.collFreq}, {self.collOut} )
      fRhsOut:accumulate(1.0, self.collOut)
   end
end

function GkLBOCollisions:write(tm, frame)
   if self.selfCollisions then
      self.uPar:write(string.format("%s_%s_%d.bp", self.speciesName, "uPar", frame), tm, frame)
      self.vthSq:write(string.format("%s_%s_%d.bp", self.speciesName, "vthSq", frame), tm, frame)
      Mpi.Allreduce(self.primMomLimitCrossingsL:data():data(), 
                    self.primMomLimitCrossingsG:data():data(), self.primMomLimitCrossingsG:size(),
                    Mpi.DOUBLE, Mpi.SUM, self.confGrid:commSet().comm)
      self.primMomLimitCrossingsG:write(string.format("%s_%s_%d.bp", self.speciesName, "primMomLimitCrossings", frame), tm, frame, true)
      self.primMomLimitCrossingsL:clear(0.0)
      self.primMomLimitCrossingsG:clear(0.0)
   end
   if self.crossCollisions then
      self.vthSqCross:write(string.format("%s_%s_%d.bp", self.speciesName, "vthSqCross", frame), tm, frame)
   end
end

function GkLBOCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.tmEvalMom
end

return GkLBOCollisions
