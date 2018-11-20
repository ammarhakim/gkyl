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

   -- intemediate storage for output of collisions
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
      self.nuFld = DataStruct.Field {
         onGrid        = self.confGrid,
         numComponents = self.cNumBasis,
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
         vUpper           = self.vMax,
      }
   else
      -- Lenard-Bernstein equation.
      vmLBOconstNuCalc = VmLBOconstNuEq {
         nu         = self.collFreq,
         phaseBasis = self.phaseBasis,
         confBasis  = self.confBasis,
         vUpper     = self.vMax,
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

   -- Number of cells in which number density was negative (somewhere).
   self.primMomLimitCrossings = DataStruct.DynVector {
      numComponents = 1,
   }
   self.primMomCrossLimitL = Lin.Vec(1)
   self.primMomCrossLimitG = Lin.Vec(1)
   -- Factor dividing zeroth-coefficient in configuration space cell average.
   self.cellAvFac          = 1.0/math.sqrt(2.0^self.confGrid:ndim())
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
-- Since this doesn't seem to be as big a problem in Vm as in Gk, we comment this out for now.
--   self.primMomLimitCrossings:write(string.format("%s_%s_%d.bp", self.speciesName, "primMomLimitCrossings", frame), tm, frame)
end

function VmLBOCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.tmEvalMom
end

return VmLBOCollisions
