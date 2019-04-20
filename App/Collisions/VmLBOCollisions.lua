-- Gkyl ------------------------------------------------------------------------
--
-- PlasmaOnCartGrid support code: Vlasov LB Collision operator
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local CollisionsBase = require "App.Collisions.CollisionsBase"
local DataStruct     = require "DataStruct"
local Proto          = require "Lib.Proto"
local Time           = require "Lib.Time"
local Updater        = require "Updater"
local VmLBOconstNuEq = require "Eq.VmLBO"
local xsys           = require "xsys"
local Lin            = require "Lib.Linalg"
local Mpi            = require "Comm.Mpi"

-- VmLBOCollisions ---------------------------------------------------------------
--
-- Lenard-Bernstein Collision operator
-- Actually dates back to Lord Rayleigh, Philos. Mag. 32, 424 (1891).
--------------------------------------------------------------------------------

local VmLBOCollisions = Proto(CollisionsBase)

-- This ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmLBOCollisions:init(tbl)
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
function VmLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.cfl = 0.0    -- Will be replaced.

   local collidingSpecies = assert(tbl.collideWith, "App.VmLBOCollisions: Must specify names of species to collide with in 'collideWith'.")

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
      self.varNu       = true                 -- Spatially varying nu.
      self.mass        = speciesTbl.mass      -- Mass of this species.
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
         self.elemCharge = 1.602176487e-19    -- Coulomb.
      else
         self.userInputNormNu = false
         if self.selfCollisions then
            self.normNuSelf  = 0.0
         end
         if self.crossCollisions then
            self.normNuCross = collidingSpecies
            table.remove(self.normNuCross, selfSpecInd)
            for i, _ in ipairs(self.normNuCross) do self.normNuCross[i] = 0.0 end 
         end
         self.epsilon0 = assert(
            tbl.epsilon0, "Updater.VmLBOCollisions: Must specify vacuum permittivity 'epsilon0' ('frequencies' and 'normNu' are not specified, so nu is calculated via Spitzer).")
         self.elemCharge = assert(
            tbl.elemCharge, "Updater.VmLBOCollisions: Must specify elementary charge with 'elemCharge' ('frequencies' and 'normNu' are not specified, so nu is calculated via Spitzer).")
      end
   end

   if self.crossCollisions then
      self.mass       = speciesTbl.mass      -- Mass of this species.
      self.charge     = speciesTbl.charge    -- Charge of this species.
      local crossOpIn = tbl.crossOption      -- Can specify 'crossOption' (Greene, GreeneSmallAngle, HeavyIons), formulas used to calculate cross-species primitive moments.
      if crossOpIn then
         self.crossMomOp  = crossOpIn
         if self.crossMomOp=="Greene" then
            local betaGreeneIn = tbl.betaGreene   -- Can specify 'betaGreene' free parameter in Grene cross-species collisions.
            if betaGreeneIn then
               self.betaGreene = betaGreeneIn
            else
               self.betaGreene = 1.0   -- Default value is the heavy ion, quasineutral limit.
            end
         else
            self.betaGreene = 0.0   -- Default value is the heavy ion, quasineutral limit.
         end
      else
         self.crossMomOp  = "Greene"    -- Default to Greene-type formulas.
         self.betaGreene  = 1.0         -- Default value is the heavy ion, quasineutral limit.
      end
   end

   self.tmEvalMom = 0.0
end

function VmLBOCollisions:setName(nm)
   self.name = nm
end
function VmLBOCollisions:setSpeciesName(nm)
   self.speciesName = nm
end

function VmLBOCollisions:setCfl(cfl)
   self.cfl = cfl -- what should this be? - AHH
end
function VmLBOCollisions:setConfBasis(basis)
   self.confBasis = basis
end
function VmLBOCollisions:setConfGrid(grid)
   self.confGrid = grid
end
function VmLBOCollisions:setPhaseBasis(basis)
   self.phaseBasis = basis
end
function VmLBOCollisions:setPhaseGrid(grid)
   self.phaseGrid = grid
end

function VmLBOCollisions:createSolver()
   self.vdim      = self.phaseGrid:ndim() - self.confGrid:ndim()

   self.cNumBasis = self.confBasis:numBasis()

   -- Maximum velocity of the velocity grid (and its square).
   self.vMax   = Lin.Vec(self.vdim)
   for vd = 1,self.vdim do
      self.vMax[vd]   = self.phaseGrid:upper(self.confGrid:ndim()+vd)
   end
   self.vMaxSq = self.vMax[1] 
   for vd = 1,self.vdim do
      if (self.vMaxSq < self.vMax[vd]) then
         self.vMaxSq = self.vMax[vd]
      end
   end
   self.vMaxSq = self.vMaxSq^2

   -- Intemediate storage for output of collisions.
   self.collOut = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Flow velocity in vdim directions.
   self.velocity = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.cNumBasis*self.vdim,
      ghost         = {1, 1},
   }
   -- Thermal speed squared, vth=sqrt(T/m).
   self.vthSq = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.cNumBasis,
      ghost         = {1, 1},
   }

   -- Zero-flux BCs in the velocity dimensions.
   local zfd = { }
   for d = 1, self.vdim do
      zfd[d] = self.confGrid:ndim() + d
   end

   if self.varNu then
      -- Collisionality, nu.
      self.collFreq = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis,
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
   else
      self.collFreq = 0.0    -- Assigned in advance method.
   end
   -- Lenard-Bernstein equation.
   local vmLBOconstNuCalc = VmLBOconstNuEq {
      phaseBasis       = self.phaseBasis,
      confBasis        = self.confBasis,
      vUpper           = self.vMax,
      varyingNu        = self.varNu,
      useCellAverageNu = self.cellConstNu,
   }
   self.collisionSlvr = Updater.HyperDisCont {
      onGrid             = self.phaseGrid,
      basis              = self.phaseBasis,
      cfl                = self.cfl,
      equation           = vmLBOconstNuCalc,
      updateDirections   = zfd, -- only update velocity directions
      zeroFluxDirections = zfd,
   }
   if self.crossCollisions then
      -- Flow velocity in vdim directions of the other species.
      self.uOther = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis*self.vdim,
         ghost         = {1, 1},
      }
      -- Thermal speed squared of the other species..
      self.vthSqOther = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis,
         ghost         = {1, 1},
      }
      -- Cross-species flow velocity in vdim directions.
      self.uCross = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis*self.vdim,
         ghost         = {1, 1},
      }
      -- Cross species thermal speed squared.
      self.vthSqCross = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis,
         ghost         = {1, 1},
      }
      -- Kinetic energy density: u dotted with M_1.
      self.kinEnergyDens = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Thermal energy density: M_2-u dot M_1.
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
      self.confDotProduct = Updater.CartFieldBinOp {
         onGrid    = self.confGrid,
         weakBasis = self.confBasis,
         operation = "DotProduct",
      }
      -- Updater to compute cross-species primitive moments.
      self.primMomCross = Updater.CrossPrimMoments {
         onGrid     = self.confGrid,
         phaseBasis = self.phaseBasis,
         confBasis  = self.confBasis,
         operator   = "VmLBO",
         formulas   = self.crossMomOp,
         betaGreene = self.betaGreene, 
      }
   end

   -- Number of cells in which number density was negative (somewhere).
   self.primMomLimitCrossings = DataStruct.DynVector {
      numComponents = 1,
   }
   self.primMomCrossLimitL = Lin.Vec(1)
   self.primMomCrossLimitG = Lin.Vec(1)
   -- Factor dividing zeroth-coefficient in configuration space cell average.
   self.cellAvFac          = 1.0/math.sqrt(2.0^self.confGrid:ndim())
end

function VmLBOCollisions:advance(tCurr, fIn, species, fRhsOut)
   local selfMom = species[self.speciesName]:fluidMoments()

   local tmEvalMomStart = Time.clock()
   if self.selfCollisions then
      -- Fetch primitive moments velocity and vthSq=T/m.
      local primMomSelf = species[self.speciesName]:selfPrimitiveMoments()
      self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

      -- NOTE: The following code is commented out because Vm users don't seem
      -- to be as worried about limit crossings as Gk users, so counting them
      -- is disabled for now. See the 'write' method as well.
      ---- Determine whether primitive moments cross limits based on
      ---- parallel flow speed and thermal speed squared.
      --self.primMomCrossLimitL[1] = 0
      --self.primMomCrossLimitG[1] = 0
      --local confIndexer          = self.velocity:genIndexer()
      --local uItr                 = self.velocity:get(1)
      --local vthSqItr             = self.vthSq:get(1)
      --for idx in self.velocity:localRangeIter() do
      --   self.velocity:fill(confIndexer(idx), uItr)
      --   self.vthSq:fill(confIndexer(idx), vthSqItr)
      --   local primCrossingFound = false
      --   for vd = 1,self.vdim do
      --      if (math.abs(uItr[(vd-1)*self.cNumBasis+1]*self.cellAvFac)>self.vMax[vd]) then
      --         uCrossingFound = true
      --         break
      --      end
      --   end
      --   local vthSq0 = vthSqItr[1]*self.cellAvFac

      --   if (uCrossingFound or (vthSq0<0) or (vthSq0>self.vMaxSq)) then
      --      self.primMomCrossLimitL[1] = self.primMomCrossLimitL[1]+1
      --   end
      --end
      --Mpi.Allreduce(self.primMomCrossLimitL:data(), self.primMomCrossLimitG:data(), 1,
      --              Mpi.DOUBLE, Mpi.SUM, self.confGrid:commSet().comm)
      --self.primMomLimitCrossings:appendData(tCurr+dt, self.primMomCrossLimitG)

      if self.varNu then
         -- Compute the collisionality.
         self.spitzerNu:advance(0.0, {self.mass, self.charge, selfMom[1], primMomSelf[2], self.normNuSelf},{self.collFreq})
      else
         self.collFreq = self.collFreqSelf
      end
      -- Compute increment from collisions and accumulate it into output.
      self.collisionSlvr:advance(
         tCurr, {fIn, primMomSelf[1], primMomSelf[2], self.collFreq}, {self.collOut})
      -- Barrier over shared communicator before accumulate
      Mpi.Barrier(self.phaseGrid:commSet().sharedComm)

      fRhsOut:accumulate(1.0, self.collOut)
   end

   if self.crossCollisions then
      if not self.selfCollisions then
         -- If applying self-collisions, use the same primitive moments calculated above.
         -- Otherwise compute primitive moments, u and vtSq, of this species.
         self.confDiv:advance(0., {selfMom[1], selfMom[2]}, {self.velocity})
         self.confDotProduct:advance(0., {self.velocity, selfMom[2]}, {self.kinEnergyDens})
         self.thermEnergyDens:combine( 1.0/self.vdim, selfMom[3],
                                      -1.0/self.vdim, self.kinEnergyDens )
         self.confDiv:advance(0., {selfMom[1], self.thermEnergyDens}, {self.vthSq})
      end

      for sInd, otherNm in ipairs(self.crossSpecies) do
         -- Obtain coupling moments of other species.
         local otherMom = species[otherNm]:fluidMoments()
         -- Compute primitive moments, u and vtSq, of other species.
         self.confDiv:advance(0., {otherMom[1], otherMom[2]}, {self.uOther})
         self.confDotProduct:advance(0., {self.uOther, otherMom[2]}, {self.kinEnergyDens})
         self.thermEnergyDens:combine( 1.0/self.vdim, otherMom[3],
                                      -1.0/self.vdim, self.kinEnergyDens )
         self.confDiv:advance(0., {otherMom[1], self.thermEnergyDens}, {self.vthSqOther})

         -- Set the subscripts of the collision term (12 or 21).
         -- For "HeavyIons" species 1 is the lighter one.
         local collSub
         local mRat = species[otherNm]:getMass()/self.mass
         if mRat <= 1 then
            collSub = "12"
         else
            collSub = "21"
         end
         -- Get cross-species primitive moments, u_12 and vtSq_12, or u_21 and vtSq_21.
         -- This updater expects inputs in the order:
         -- collSub, mRat, nSelf, uSelf, vtSqSelf, nOther, uOther, vtSqOther.
         self.primMomCross:advance(0., {collSub, mRat, selfMom[1], self.velocity, self.vthSq,
                                        otherMom[1], self.uOther, self.vthSqOther},
                                       {self.uCross, self.vthSqCross})

         if self.varNu then
            -- Compute the collisionality.
            self.spitzerNu:advance(0., {self.mass, self.charge, otherMom[1], self.vthSq, self.normNuCross[sInd]}, {self.collFreq})
         else
            self.collFreq = self.collFreqCross[sInd]
         end

         -- Compute increment from cross-species collisions and accumulate it into output.
         self.collisionSlvr:advance(
            tCurr, {fIn, self.uCross, self.vthSqCross, self.collFreq}, {self.collOut} )
         -- Barrier over shared communicator before accumulate
         Mpi.Barrier(self.phaseGrid:commSet().sharedComm)

         fRhsOut:accumulate(1.0, self.collOut)
      end
   end
end

function VmLBOCollisions:write(tm, frame)
--   self.velocity:write(string.format("%s_%s_%d.bp", self.speciesName, "u", frame), tm, frame)
--   self.vthSq:write(string.format("%s_%s_%d.bp", self.speciesName, "vthSq", frame), tm, frame)
-- Since this doesn't seem to be as big a problem in Vm as in Gk, we comment this out for now.
--   self.primMomLimitCrossings:write(string.format("%s_%s_%d.bp", self.speciesName, "primMomLimitCrossings", frame), tm, frame)
--   if self.crossCollisions then
--      self.uCross:write(string.format("%s_%s_%d.bp", self.speciesName, "uCross", frame), tm, frame)
--      self.vthSqCross:write(string.format("%s_%s_%d.bp", self.speciesName, "vthSqCross", frame), tm, frame)
--   end
end

function VmLBOCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.tmEvalMom
end

function VmLBOCollisions:slvrTime()
   return self.collisionSlvr.totalTime
end

function VmLBOCollisions:momTime()
   return self.tmEvalMom
end

return VmLBOCollisions
