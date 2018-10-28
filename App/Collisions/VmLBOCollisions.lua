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

-- VmLBOCollisions ---------------------------------------------------------------
--
-- Lenard-Bernstein Collision operator
-- Actually dates back to Lord Rayleigh, Philos. Mag. 32, 424 (1891).
--------------------------------------------------------------------------------

local VmLBOCollisions = Proto(CollisionsBase)

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function VmLBOCollisions:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function VmLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- Previously stored table.

   self.cfl            = 0.0    -- Will be replaced.
   self.selfCollisions = xsys.pickBool(tbl.selfCollisions, true) -- By default, self collisions are on.
   self.crossSpecies   = tbl.crossSpecies

   local constNu       = tbl.collFreq
   if constNu then
      self.varNu       = false    -- Not spatially varying nu.
      self.collFreq    = constNu
      self.cellConstNu = true
   else
      self.varNu       = true    -- Spatially varying nu.
      self.normNu      = assert(tbl.normNu, "App.VmLBOCollisions: Must specify 'normNu', collisionality normalized by (T_0^(3/2)/n_0) evaluated somewhere in the simulation.")
      self.mass        = speciesTbl.mass
      -- For now only cell-wise constant nu is implemented.
      -- self.cellConstNu = assert(tbl.useCellAverageNu, "App.VmLBOCollisions: Must specify 'useCellAverageNu=true/false' for using cellwise constant/expanded spatially varying collisionality.")
      self.cellConstNu = true
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
   self.vdim = self.phaseGrid:ndim() - self.confGrid:ndim()

   -- intemediate storage for output of collisions
   self.collOut = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Flow velocity in vdim directions.
   self.velocity = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis()*self.vdim,
      ghost         = {1, 1},
   }
   -- Thermal speed squared, vth=sqrt(T/m).
   self.vthSq = DataStruct.Field {
      onGrid        = self.confGrid,
      numComponents = self.confBasis:numBasis(),
      ghost         = {1, 1},
   }

   -- Zero-flux BCs in the velocity dimensions.
   local zfd = { }
   for d = 1, self.vdim do
      zfd[d] = self.confGrid:ndim() + d
   end

   if self.varNu then
      -- Collisionality, nu.
      self.nuFld = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.confBasis:numBasis(),
         ghost         = {1, 1},
      }
      -- Updater to compute spatially varying (Spitzer) nu.
      self.spitzerNu = Updater.SpitzerCollisionality {
         onGrid           = self.confGrid,
         confBasis        = self.confBasis,
         normalizedNu     = self.normNu,
         mass             = self.mass,
         useCellAverageNu = self.cellConstNu,
      }
      -- Lenard-Bernstein equation.
      vmLBOconstNuCalc = VmLBOconstNuEq {
         phaseBasis       = self.phaseBasis,
         confBasis        = self.confBasis,
         useCellAverageNu = self.cellConstNu,
      }
   else
      -- Lenard-Bernstein equation.
      vmLBOconstNuCalc = VmLBOconstNuEq {
         nu         = self.collFreq,
         phaseBasis = self.phaseBasis,
         confBasis  = self.confBasis,
      }
   end
   self.collisionSlvr = Updater.HyperDisCont {
      onGrid             = self.phaseGrid,
      basis              = self.phaseBasis,
      cfl                = self.cfl,
      equation           = vmLBOconstNuCalc,
      onlyIncrement      = true,
      updateDirections   = zfd, -- only update velocity directions
      zeroFluxDirections = zfd,
   }
   self.primMomSelf = Updater.SelfPrimMoments {
      onGrid     = self.confGrid,
      phaseGrid  = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confBasis  = self.confBasis,
      operator   = "VmLBO",
   }
end

function VmLBOCollisions:forwardEuler(tCurr, dt, fIn, species, fOut)
   local status, dtSuggested = true, GKYL_MAX_DOUBLE
   local selfMom = species[self.speciesName]:fluidMoments()

   if self.selfCollisions then
      local tmEvalMomStart = Time.clock()
      -- Compute primitive moments velocity and vthSq=T/m from zeroth,
      -- first and second moments, and distribution function.
      self.primMomSelf:advance(0.0, 0.0, {selfMom[1], selfMom[2], selfMom[3],fIn},
                                         {self.velocity,self.vthSq})
      self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

      if self.varNu then
         -- Compute the collisionality.
         self.spitzerNu:advance(0.0, 0.0, {selfMom[1], self.vthSq},{self.nuFld})

         -- Compute increment from collisions and accumulate it into output.
         tmpStatus, tmpDt = self.collisionSlvr:advance(
   	    tCurr, dt, {fIn,self.velocity,self.vthSq,self.nuFld}, {self.collOut})
      else
         -- Compute increment from collisions and accumulate it into output.
         tmpStatus, tmpDt = self.collisionSlvr:advance(
   	    tCurr, dt, {fIn, self.velocity, self.vthSq}, {self.collOut})
      end

      status = status and tmpStatus
      dtSuggested = math.min(dtSuggested, tmpDt)

      fOut:accumulate(dt, self.collOut)
   end
   if self.crossSpecies then
      -- Insert cross collisions here!
   end
   return status, dtSuggested
end

function VmLBOCollisions:write(tm, frame)
   self.velocity:write(string.format("%s_%s_%d.bp", self.speciesName, "u", frame), tm, frame)
   self.vthSq:write(string.format("%s_%s_%d.bp", self.speciesName, "vthSq", frame), tm, frame)
end

function VmLBOCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.tmEvalMom
end

return VmLBOCollisions
