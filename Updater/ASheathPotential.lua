-- Gkyl ------------------------------------------------------------------------
--
-- Ambipolar, adiabatic electron Sheath Potential (ASheathPotential):
-- Asumming ambipolar fluxes at the sheath entrance, compute the sheath
-- and the electrostatic potential in the whole domain.
-- The potential is computed from
--    phi = phi_s - (T_e/q_e)*log(n_e/n_{es})
-- which thanks to quasineutrality becomes
--    phi = phi_s - (T_e/q_e)*log(n_i/n_{is})
-- where phi_s and n_{ks} are the potential and density of species k at the
-- sheath entrance, respectively. Since we typically have two sheaths, we
-- compute this in two parts. Continuity is later enforced by an FEM smoothing
-- operator. The sheath potential is
--    phi_s = -(T_e/q_e)*log( sqrt(2*pi)*Gamma_i/(n_e*v_{te}) )
-- which using quasineutrality again gives
--    phi_s = -(T_e/q_e)*( log( sqrt(2*pi)*Gamma_i/(n_i*v_{te}) ) )|_{z=z_s}
-- with Gamma_i being the ion particle flux, and the natural logarithm being
-- evaluated at the sheath entrance.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase = require "Updater.Base"
local DataStruct  = require "DataStruct"
local Lin         = require "Lib.Linalg"
local Proto       = require "Lib.Proto"
local xsys        = require "xsys"
local ModDecl     = require "Updater.aSheathPotentialData.asheath_potential_mod_decl"
local Mpi         = require "Comm.Mpi"
local ffi            = require "ffi"

local ffiC = ffi.C
require "Lib.ZeroUtil"

-- Declaration of gkylzero objects and functions.
ffi.cdef [[
// Object type
typedef struct gkyl_ambi_bolt_potential gkyl_ambi_bolt_potential;

/**
 * Create new updater to compute the electrostatic potential assuming ambipolar
 * sheath particle fluxes and Boltzmann isothermal electrons.
 *
 * @param grid Cartesian grid dynamic field is defined on.
 * @param basis Basis object (configuration space).
 * @param mass_e Electron mass.
 * @param charge_e Electron charge.
 * @param temp_e Electron temperature.
 * @param use_gpu Boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_ambi_bolt_potential* gkyl_ambi_bolt_potential_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, double mass_e, double charge_e, double temp_e,
  bool use_gpu);

/**
 * Compute the ion density and electrostatic potential at the sheath entrance.
 * Below we use jac to refer to the fact that some of these quantities are
 * multiplied by the configuration space Jacobian (e.g. for field aligned
 * coordinates). In reality they are also multiplied by the phase space
 * Jacobian (approx B), but we focus on jac because need to divide by it to
 * get the potential (or multiply by 1/jac).
 *
 * @param Ambipolar, Boltzmann electron sheath potential updater.
 * @param edge Lower (-1) or upper (1) boundary along the field line.
 * @param skin_r Global skin range along the field line.
 * @param ghost_r Corresponding global ghost range along the field line.
 * @param jacob_geo_inv Reciprocal of the configuration space Jacobian.
 * @param gammai Ion particle flux at the sheath entrance times the
 *               conf-space Jacobian.
 * @param m0i Ion number density times the conf-space Jacobian.
 * @param sheath_vals Ion number density and potential at the sheath entrance.
 */
void
gkyl_ambi_bolt_potential_sheath_calc(struct gkyl_ambi_bolt_potential *up, enum gkyl_edge_loc edge,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r, const struct gkyl_array *jacob_geo_inv,
  const struct gkyl_array *gammai, const struct gkyl_array *m0i, struct gkyl_array *sheath_vals);

/**
 * Compute the electrostatic potential in the domain as
 *  phi = phi_s + (T_e/e)*ln(n_i/n_is).
 *
 * @param up Ambipolar, Boltzmann electron sheath potential updater.
 * @param local_r Local range.
 * @param extlocal_r Extended local range.
 * @param jacob_geo_inv Reciprocal of the configuration space Jacobian.
 * @param m0i Ion number density times the conf-space Jacobian.
 * @param sheath_vals Ion number density and potential at the sheath entrance.
 * @param phi electrostatic potential.
 */
void
gkyl_ambi_bolt_potential_phi_calc(struct gkyl_ambi_bolt_potential *up,
  const struct gkyl_range *local_r, const struct gkyl_range *extlocal_r,
  const struct gkyl_array *jacob_geo_inv, const struct gkyl_array *m0i,
  const struct gkyl_array *sheath_vals, struct gkyl_array *phi);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_ambi_bolt_potential_release(gkyl_ambi_bolt_potential *up);
]]

-- Inherit the base Updater from UpdaterBase updater object.
local ASheathPotential = Proto(UpdaterBase)

function ASheathPotential:init(tbl)
   ASheathPotential.super.init(self, tbl) -- setup base object

   local grid  = assert(tbl.onGrid, "Updater.ASheathPotential: Must provide configuration space grid object 'onGrid'.")
   local basis = assert(tbl.basis, "Updater.ASheathPotential: Must provide configuration space basis object 'basis'.")
   local boundaryGrids = assert(tbl.boundaryGrids, "Updater.ASheathPotential: Must provide configuration space boundary grids 'boundaryGrids'.")

   local mElc    = assert(tbl.electronMass, "Updater.ASheathPotential: Must provide electron mass in 'electronMass'.")
   local qElc    = assert(tbl.electronCharge, "Updater.ASheathPotential: Must provide electron charge in 'electronCharge'.")
   local tempElc = assert(tbl.electronTemp, "Updater.ASheathPotential: Must provide electron temperature in 'electronTemp'.")

   local useGPU  = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or 0)

   local dim = basis:ndim()
   self.sheathDir = dim   -- Assume the sheath direction is the last dimension.

   self.boundary = {"lower","upper"}

   self._zero = ffi.gc(ffiC.gkyl_ambi_bolt_potential_new(grid._zero, basis._zero, mElc, qElc, tempElc, useGPU),
                       ffiC.gkyl_ambi_bolt_potential_release)

   self.sheathValues = {lower=DataStruct.Field{onGrid = boundaryGrids["lower"],
                                               numComponents = 2*basis:numBasis(), ghost = {0,0},},
                        upper=DataStruct.Field{onGrid = boundaryGrids["upper"], 
                                               numComponents = 2*basis:numBasis(), ghost = {0,0},}}
   self.bcastSheathQuants = function(bInd) end
   if GKYL_HAVE_MPI then
      -- Need a communicator to broadcast the sheath potential and density along z.
      local commSet   = grid:commSet()
      local worldComm = commSet.comm
      if Mpi.Comm_size(worldComm) > 1 then
         local worldRank = Mpi.Comm_rank(worldComm)
         local zCommRank = 0
         if dim==3 then zCommRank = worldRank % (grid:cuts(1)*grid:cuts(2)) end
         self.zComm = Mpi.Comm_split(worldComm, zCommRank, worldRank)
         local zCommSize = Mpi.Comm_size(self.zComm)
         self.sheathRank = {lower=0, upper=zCommSize-1}
         self.numBoundaryDOFs = {lower=self.sheathValues["lower"]:localRange():volume()*self.sheathValues["lower"]:numComponents(),
                                 upper=self.sheathValues["upper"]:localRange():volume()*self.sheathValues["upper"]:numComponents()}
         self.bcastSheathQuants = function(bInd)
            Mpi.Bcast(self.sheathValues[bInd]:dataPointer(), self.numBoundaryDOFs[bInd], Mpi.DOUBLE, self.sheathRank[bInd], self.zComm)
         end
      end
   end
end

function ASheathPotential:_advance(tCurr, inFlds, outFlds)
   local GammaIon, m0Ion, jacobGeoInv = inFlds[1], inFlds[2], inFlds[3]
   local phi = outFlds[1]

   for _, b in ipairs(self.boundary) do
      -- Loop over boundary grid and compute phi_s in each cell using quadrature.
      local edge = b=="lower" and 0 or 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
      local skinRange = b=="lower" and phi:localGlobalSkinRangeIntersectLower()[self.sheathDir]
                                    or phi:localGlobalSkinRangeIntersectUpper()[self.sheathDir]
      local ghostRange = b=="lower" and phi:localGlobalGhostRangeIntersectLower()[self.sheathDir]
                                     or phi:localGlobalGhostRangeIntersectUpper()[self.sheathDir]
      ffiC.gkyl_ambi_bolt_potential_sheath_calc(self._zero, edge, skinRange, ghostRange, jacobGeoInv._zero, GammaIon[b]._zero, m0Ion._zero, self.sheathValues[b]._zero)

      -- Broadcast the sheath potential and density to other ranks along z.
      self.bcastSheathQuants(b)
   end

   -- Average left and right sheath density and potential.
   self.sheathValues["lower"]:accumulate(1., self.sheathValues["upper"])
   self.sheathValues["lower"]:scale(0.5)
      
   ffiC.gkyl_ambi_bolt_potential_phi_calc(self._zero, phi:localRange(), phi:localExtRange(), jacobGeoInv._zero, m0Ion._zero, self.sheathValues["lower"]._zero, phi._zero)
end

function ASheathPotential:_advanceOnDevice(tCurr, inFlds, outFlds)
   local GammaIon, m0Ion, jacobGeoInv = inFlds[1], inFlds[2], inFlds[3]
   local phi = outFlds[1]

   for _, b in ipairs(self.boundary) do
      -- Loop over boundary grid and compute phi_s in each cell using quadrature.
      local edge = b=="lower" and 0 or 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
      local skinRange = b=="lower" and phi:localGlobalSkinRangeIntersectLower()[self.sheathDir]
                                    or phi:localGlobalSkinRangeIntersectUpper()[self.sheathDir]
      local ghostRange = b=="lower" and phi:localGlobalGhostRangeIntersectLower()[self.sheathDir]
                                     or phi:localGlobalGhostRangeIntersectUpper()[self.sheathDir]
      ffiC.gkyl_ambi_bolt_potential_sheath_calc(self._zero, edge, skinRange, ghostRange, jacobGeoInv._zeroDevice, GammaIon[b]._zeroDevice, m0Ion._zeroDevice, self.sheathValues[b]._zeroDevice)

      -- Broadcast the sheath potential and density to other ranks along z.
      self.bcastSheathQuants(b)
   end

   -- Average left and right sheath density and potential.
   self.sheathValues["lower"]:accumulate(1., self.sheathValues["upper"])
   self.sheathValues["lower"]:scale(0.5)
      
   ffiC.gkyl_ambi_bolt_potential_phi_calc(self._zero, phi:localRange(), phi:localExtRange(), jacobGeoInv._zeroDevice, m0Ion._zeroDevice, self.sheathValues["lower"]._zeroDevice, phi._zeroDevice)
end

return ASheathPotential
