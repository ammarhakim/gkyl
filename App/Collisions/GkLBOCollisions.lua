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

-- this ctor simply stores what is passed to it and defers actual
-- construction to the fullInit() method below.
function GkLBOCollisions:init(tbl)
   self.tbl = tbl
end

-- Actual function for initialization. This indirection is needed as
-- we need the app top-level table for proper initialization.
function GkLBOCollisions:fullInit(speciesTbl)
   local tbl = self.tbl -- previously stored table

   self.cfl = 0.0 -- will be replaced
   self.selfCollisions = xsys.pickBool(tbl.selfCollisions, true) -- by default, self collisions are on
   self.crossSpecies   = tbl.crossSpecies

   local constNu       = tbl.collFreq
   if constNu then
      self.varNu       = false    -- Not spatially varying nu.
      self.collFreq    = constNu
      self.cellConstNu = true
   else
      self.varNu       = true    -- Spatially varying nu.
      self.normNu      = assert(tbl.normNu, "App.VmLBOCollisions: Must specify 'normNu', collisionality normalized by (T_0^(3/2)/n_0) evaluated somewhere in the simulation.")
      -- For now only cell-wise constant nu is implemented.
      -- self.cellConstNu = assert(tbl.useCellAverageNu, "App.VmLBOCollisions: Must specify 'useCellAverageNu=true/false' for using cellwise constant/expanded spatially varying collisionality.")
      self.cellConstNu = true
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
   self.vdim = self.phaseGrid:ndim() - self.confGrid:ndim()

   -- Maximum velocity of the velocity grid (and its square).
   self.vParMax            = self.phaseGrid:upper(self.confGrid:ndim()+1)
   self.vParMaxSq          = self.vParMax^2

   -- intemediate storage for output of collisions
   self.collOut = DataStruct.Field {
      onGrid        = self.phaseGrid,
      numComponents = self.phaseBasis:numBasis(),
      ghost         = {1, 1},
   }

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
   -- Inverse of background magnetic field.
   self.bmag    = funcField.geo.bmag
   -- Inverse of background magnetic field.
   self.bmagInv = funcField.geo.bmagInv
      
   -- Zero-flux BCs in the velocity dimensions.
   local zfd = { }
   for d = 1, self.vdim do
      zfd[d] = self.confGrid:ndim() + d
   end

   local gkLBOconstNuCalc
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
      gkLBOconstNuCalc = GkLBOconstNuEq {
         phaseBasis       = self.phaseBasis,
         confBasis        = self.confBasis,
         useCellAverageNu = self.cellConstNu,
         vParUpper        = self.vParMax,
         mass             = self.mass,
      }
   else
      -- Lenard-Bernstein equation.
      gkLBOconstNuCalc = GkLBOconstNuEq {
         nu         = self.collFreq,
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
      equation           = gkLBOconstNuCalc,
      updateDirections   = zfd, -- only update velocity directions
      zeroFluxDirections = zfd,
   }
   self.primMomSelf = Updater.SelfPrimMoments {
      onGrid     = self.confGrid,
      phaseGrid  = self.phaseGrid,
      phaseBasis = self.phaseBasis,
      confBasis  = self.confBasis,
      gkfacs     = {self.mass, self.bmag},
      operator   = "GkLBO",
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

function GkLBOCollisions:advance(tCurr, cflRateByCell, fIn, species, fRhsOut)
   local selfMom = species[self.speciesName]:fluidMoments()

   if self.selfCollisions then
      local tmEvalMomStart = Time.clock()
      -- Compute primitive moments velocity and vthSq=T/m from zeroth,
      -- first and second moments, and distribution function.
      self.primMomSelf:advance(0.0, 0.0, {selfMom[1], selfMom[2], selfMom[3],fIn},
                                         {self.uPar,self.vthSq})
      self.tmEvalMom = self.tmEvalMom + Time.clock() - tmEvalMomStart

      -- Determine whether primitive moments cross limits based on
      -- parallel flow speed and thermal speed squared.
      self.primMomCrossLimitL[1] = 0
      self.primMomCrossLimitG[1] = 0
      local confIndexer       = self.uPar:genIndexer()
      local uParItr           = self.uPar:get(1)
      local vthSqItr          = self.vthSq:get(1)
      for idx in self.uPar:localRangeIter() do
         self.uPar:fill(confIndexer(idx), uParItr)
         self.vthSq:fill(confIndexer(idx), vthSqItr)
         -- Cell average values (mind normalization).
         local uPar0  = uParItr[1]*self.cellAvFac
         local vthSq0 = vthSqItr[1]*self.cellAvFac

         if ((math.abs(uPar0)>self.vParMax) or
             (vthSq0<0) or (vthSq0>self.vParMaxSq)) then
            self.primMomCrossLimitL[1] = self.primMomCrossLimitL[1]+1
         end
      end
      Mpi.Allreduce(self.primMomCrossLimitL:data(), self.primMomCrossLimitG:data(), 1,
                    Mpi.DOUBLE, Mpi.SUM, self.confGrid:commSet().comm)
      self.primMomLimitCrossings:appendData(tCurr, self.primMomCrossLimitG)

      if self.varNu then
         -- Compute the collisionality.
         self.spitzerNu:advance(0.0, nil, {selfMom[1], self.vthSq},{self.nuFld})

         -- Compute increment from collisions and accumulate it into output.
         self.collisionSlvr:advance(
            tCurr, cflRateByCell, {fIn, self.bmagInv, self.uPar, self.vthSq, self.nuFld}, {self.collOut})
      else
         -- Compute increment from collisions and accumulate it into output.
         self.collisionSlvr:advance(
            tCurr, cflRateByCell, {fIn, self.bmagInv, self.uPar, self.vthSq}, {self.collOut})
      end

      fRhsOut:accumulate(1.0, self.collOut)
   end
   if self.crossSpecies then
      -- Insert cross collisions here!
   end
end

function GkLBOCollisions:write(tm, frame)
   self.uPar:write(string.format("%s_%s_%d.bp", self.speciesName, "uPar", frame), tm, frame)
   self.vthSq:write(string.format("%s_%s_%d.bp", self.speciesName, "vthSq", frame), tm, frame)
   self.primMomLimitCrossings:write(string.format("%s_%s_%d.bp", self.speciesName, "primMomLimitCrossings", frame), tm, frame)
end

function GkLBOCollisions:totalTime()
   return self.collisionSlvr.totalTime + self.tmEvalMom
end

return GkLBOCollisions
