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
local function findInd(tblIn, el)
   for i, v in ipairs(tblIn) do
      if v == el then
         return i
      end
   end
   return #tblIn+1    -- If not found return a number larger than the length of the table.
end


-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- previously stored table

   self.cfl = 0.0    -- Will be replaced.

   self.collidingSpecies = assert(tbl.collideWith, "App.GkLBOCollisions: Must specify names of species to collide with in 'collideWith'.")

   -- First determine if self-species and/or cross-species collisions take place,
   -- and (if cross-collisions=true) put the names of the other colliding species in a list.
   local selfSpecInd = findInd(self.collidingSpecies, self.speciesName)
   if selfSpecInd < (#self.collidingSpecies+1) then
      self.selfCollisions = true                 -- Apply self-species collisions.
      if #self.collidingSpecies > 1 then
         self.crossCollisions = true             -- Apply cross-species collisions.
         self.crossSpecies    = self.collidingSpecies
         table.remove(self.crossSpecies, selfSpecInd)
      else
         self.crossCollisions = false            -- Don't apply cross-species collisions.
      end
   else
      self.selfCollisions  = false               -- Don't apply self-species collisions.
      self.crossCollisions = true                -- Apply cross-species collisions.
      self.crossSpecies    = self.collidingSpecies    -- All species in collidingSpecies must be cross-species.
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
      self.varNu       = true    -- Spatially varying nu.
      self.charge      = speciesTbl.charge    -- Charge of this species.
      -- For now only cell-wise constant nu is implemented.
      -- self.cellConstNu = assert(tbl.cellAvFrequencies, "App.GkLBOCollisions: Must specify 'useCellAverageNu=true/false' for using cellwise constant/expanded spatially varying collisionality.")
      self.cellConstNu = true
      -- If no constant collision frequencies provided ('frequencies'), user can specify 'normNu'
      -- list of collisionalities normalized by (T_0^(3/2)/n_0) evaluated somewhere in the
      -- simulation. Otherwise code compute Spitzer collisionality from scratch.
      local normNuIn   = tbl.normNu
      -- normNuSelf, epsilon0 and elemCharge may not used, but are
      -- initialized to avoid if-statements in advance method.
      if normNuIn then
         self.userInputNormNu = true
         if self.selfCollisions then
            self.normNuSelf  = normNuIn[selfSpecInd]
         end
         if self.crossCollisions then
            self.normNuCross = normNuIn
            table.remove(self.normNuCross, selfSpecInd)
         end
         self.epsilon0   = 8.854187817620389850536563031710750260608e-12    -- Farad/meter.
         self.elemCharge = 1.602176487e-19    -- Coulomb
      else
         self.userInputNormNu = false
         if self.selfCollisions then
            self.normNuSelf  = 0.0
         end
         if self.crossCollisions then
            self.normNuCross = self.collidingSpecies
            table.remove(self.normNuCross, selfSpecInd)
            for i, _ in ipairs(self.normNuCross) do self.normNuCross[i] = 0.0 end
         end
         self.epsilon0 = assert(
            tbl.epsilon0, "Updater.GkLBOCollisions: Must specify vacuum permittivity 'epsilon0' ('frequencies' and 'normNu' are not specified, so nu is calculated via Spitzer).")
         self.elemCharge = assert(
            tbl.elemCharge, "Updater.GkLBOCollisions: Must specify elementary charge with 'elemCharge' ('frequencies' and 'normNu' are not specified, so nu is calculated via Spitzer).")
      end
   end

   if self.crossCollisions then
      self.charge        = speciesTbl.charge    -- Charge of this species.
      local betaGreeneIn = tbl.betaGreene       -- Can specify 'betaGreene' free parameter in Grene cross-species collisions.
      if betaGreeneIn then
         self.betaGreene = betaGreeneIn
      else
         self.betaGreene = 0.0   -- Default value.
      end
   end

   self.mass = speciesTbl.mass   -- Mass of this species.

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
   -- Parallel flow velocity times collisionality, summed over species.
   self.nuUParSum = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }
   -- Thermal speed squared times collisionality, summed over species.
   self.nuVtSqSum = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
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
      -- Collisionality, nu, summed over all species pairs.
      self.nuSum = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Updater to compute spatially varying (Spitzer) nu.
      self.spitzerNu = Updater.SpitzerCollisionality {
         onGrid           = self.confGrid,
         confBasis        = self.confBasis,
         useCellAverageNu = self.cellConstNu,
         willInputNormNu  = self.userInputNormNu,
         elemCharge       = self.elemCharge,
         epsilon0         = self.epsilon0,
      }
      -- Weak multiplication to multiply nu(x) with uPar or vtSq.
      self.confMul = Updater.CartFieldBinOp {
         onGrid    = self.confGrid,
         weakBasis = self.confBasis,
         operation = "Multiply",
      }
   else
      self.nuSum = 0.0    -- Assigned in advance method.
   end
   -- Lenard-Bernstein equation.
   self.gkLBOconstNuCalcEq = GkLBOconstNuEq {
      phaseBasis       = self.phaseBasis,
      confBasis        = self.confBasis,
      vParUpper        = self.vParMax,
      varyingNu        = self.varNu,
      useCellAverageNu = self.cellConstNu,
      mass             = self.mass,
   }
   self.collisionSlvr = Updater.HyperDisCont {
      onGrid             = self.phaseGrid,
      basis              = self.phaseBasis,
      cfl                = self.cfl,
      equation           = self.gkLBOconstNuCalcEq,
      updateDirections   = zfd, -- only update velocity directions
      zeroFluxDirections = zfd,
   }
   if self.crossCollisions then
      if self.varNu then
         -- Temporary collisionality field.
         self.collFreq = DataStruct.Field {
            onGrid        = self.confGrid,
            numComponents = self.cNumBasis,
            ghost         = {1, 1},
         }
         -- Cross-species u and vtSq multiplied by collisionality.
         self.nuUParCross = DataStruct.Field {
            onGrid        = self.confGrid,
            numComponents = self.confBasis:numBasis(),
            ghost         = {1, 1},
         }
         self.nuVtSqCross = DataStruct.Field {
            onGrid        = self.confGrid,
            numComponents = self.confBasis:numBasis(),
            ghost         = {1, 1},
         }
      else
         self.collFreq = 0.0
      end
      -- Updater to compute cross-species primitive moments.
      self.primMomCross = Updater.CrossPrimMoments {
         onGrid     = self.confGrid,
         phaseBasis = self.phaseBasis,
         confBasis  = self.confBasis,
         operator   = "GkLBO",
         betaGreene = self.betaGreene,
      }
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

   -- Fetch coupling moments and primitive moments of this species.
   local selfMom     = species[self.speciesName]:fluidMoments()
   local primMomSelf = species[self.speciesName]:selfPrimitiveMoments()

   if self.varNu then
      self.nuSum:clear(0.0)
   else
      self.nuSum = 0.0
   end
   self.nuUParSum:clear(0.0)
   self.nuVtSqSum:clear(0.0)

   local tmEvalMomStart = Time.clock()
   if self.selfCollisions then
      self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

      self.gkLBOconstNuCalcEq.primMomCrossLimit = 0.0

      if self.varNu then
         -- Compute the collisionality.
         self.spitzerNu:advance(tCurr, {self.mass, self.charge, selfMom[1], primMomSelf[2], self.normNuSelf},{self.nuSum})
         self.confMul:advance(tCurr, {self.nuSum, primMomSelf[1]}, {self.nuUParSum})
         self.confMul:advance(tCurr, {self.nuSum, primMomSelf[2]}, {self.nuVtSqSum})
      else
         self.nuSum = self.collFreqSelf
         self.nuUParSum:combine(self.collFreqSelf, primMomSelf[1])
         self.nuVtSqSum:combine(self.collFreqSelf, primMomSelf[2])
      end
   end    -- end if self.selfCollisions.

   if self.crossCollisions then

      local bCorrectionsSelf = species[self.speciesName]:boundaryCorrections()
      local starMomSelf      = species[self.speciesName]:starMoments()

      for sInd, otherNm in ipairs(self.crossSpecies) do

         local mOther            = species[otherNm]:getMass()
         local otherMom          = species[otherNm]:fluidMoments()
         local primMomOther      = species[otherNm]:selfPrimitiveMoments()
         local bCorrectionsOther = species[otherNm]:boundaryCorrections()
         local starMomOther      = species[otherNm]:starMoments()

         local collFreqOther
         if self.varNu then
            -- Compute the collisionality.
            self.spitzerNu:advance(tCurr, {self.mass, self.charge, otherMom[1], primMomSelf[2], self.normNuCross[sInd]}, {self.collFreq})
         else
            self.collFreq = self.collFreqCross[sInd]
            collFreqOther = self.mass*self.collFreq/mOther
         end

         local tmEvalMomStart = Time.clock()
         if (not (species[self.speciesName].momentFlags[5][otherNm] and
                  species[otherNm].momentFlags[5][self.speciesName])) then
            -- Cross-primitive moments for the collision of these two species has not been computed.
            self.primMomCross:advance(tCurr, {self.mass, self.collFreq, selfMom, primMomSelf, bCorrectionsSelf, starMomSelf,
                                              mOther, collFreqOther, otherMom, primMomOther, bCorrectionsOther, starMomOther},
                                             {species[self.speciesName].uParCross[otherNm], species[self.speciesName].vtSqCross[otherNm],
                                              species[otherNm].uParCross[self.speciesName], species[otherNm].vtSqCross[self.speciesName]})

            species[self.speciesName].momentFlags[5][otherNm] = true
            species[otherNm].momentFlags[5][self.speciesName] = true
         end
         self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

         if self.varNu then
            self.confMul:advance(tCurr, {self.collFreq, species[self.speciesName].uParCross[otherNm]}, {self.nuUParCross})
            self.confMul:advance(tCurr, {self.collFreq, species[self.speciesName].vtSqCross[otherNm]}, {self.nuVtSqCross})

            self.nuSum:accumulate(1.0, self.collFreq)
            self.nuUParSum:accumulate(1.0, self.nuUParCross)
            self.nuVtSqSum:accumulate(1.0, self.nuVtSqCross)
         else
            self.nuSum = self.nuSum+self.collFreq
            self.nuUParSum:accumulate(self.collFreq, species[self.speciesName].uParCross[otherNm])
            self.nuVtSqSum:accumulate(self.collFreq, species[self.speciesName].vtSqCross[otherNm])
         end

      end    -- end loop over other species that this species collides with.

   end    -- end if self.crossCollisions.

   -- Compute increment from collisions and accumulate it into output.
   self.collisionSlvr:advance(
      tCurr, {fIn, self.bmagInv, self.nuUParSum, self.nuVtSqSum, self.nuSum}, {self.collOut})

   self.primMomLimitCrossingsG:appendData(tCurr, {0.0})
   self.primMomLimitCrossingsL:appendData(tCurr, {self.gkLBOconstNuCalcEq.primMomCrossLimit})

   fRhsOut:accumulate(1.0, self.collOut)

end

function GkLBOCollisions:write(tm, frame)
--   self.uPar:write(string.format("%s_%s_%d.bp", self.speciesName, "uPar", frame), tm, frame)
--   self.vthSq:write(string.format("%s_%s_%d.bp", self.speciesName, "vthSq", frame), tm, frame)
   if self.selfCollisions then
      Mpi.Allreduce(self.primMomLimitCrossingsL:data():data(), 
                    self.primMomLimitCrossingsG:data():data(), self.primMomLimitCrossingsG:size(),
                    Mpi.DOUBLE, Mpi.SUM, self.confGrid:commSet().comm)
      self.primMomLimitCrossingsG:write(string.format("%s_%s_%d.bp", self.speciesName, "primMomLimitCrossings", frame), tm, frame, true)
      self.primMomLimitCrossingsL:clear(0.0)
      self.primMomLimitCrossingsG:clear(0.0)
   end
--   if self.crossCollisions then
--      if self.crossMomOp ~= "HeavyIons" then
--         self.uParCross:write(string.format("%s_%s_%d.bp", self.speciesName, "uParCross", frame), tm, frame)
--      end
--      self.vthSqCross:write(string.format("%s_%s_%d.bp", self.speciesName, "vthSqCross", frame), tm, frame)
--   end
end

function GkLBOCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.tmEvalMom
end

function GkLBOCollisions:slvrTime()
   return self.collisionSlvr.totalTime
end

function GkLBOCollisions:momTime()
   return self.tmEvalMom
end

return GkLBOCollisions
