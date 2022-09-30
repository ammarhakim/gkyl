-- Gkyl ------------------------------------------------------------------------
--
-- Modify the electrostatic potential near at the sheath entrance to account
-- for a rarefaction wave that speeds the ions up to the sound speed.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------
-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local Proto          = require "Lib.Proto"
local UpdaterBase    = require "Updater.Base"
local ffi            = require "ffi"

local ffiC = ffi.C
require "Lib.ZeroUtil"

-- Declaration of gkylzero objects and functions.
ffi.cdef [[
// Object type
typedef struct gkyl_sheath_rarefaction_pot gkyl_sheath_rarefaction_pot;

/**
 * Create a new updater to modify the electrostatic potential at the
 * boundary to account for the rarefaction wave.
 *
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param local_range_ext Local extended range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param basis Basis on which coefficients in array are expanded (a device pointer if use_gpu=true).
 * @param grid cartesian grid dynamic field is defined on.
 * @param elem_charge elementary charge (i.e. charge of singly ionized ion).
 * @param mass_e electron mass.
 * @param mass_i ion mass.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_sheath_rarefaction_pot*
gkyl_sheath_rarefaction_pot_new(enum gkyl_edge_loc edge, const struct gkyl_range *local_range_ext,
  const int *num_ghosts, const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid,
  double elem_charge, double mass_e, double mass_i, bool use_gpu);

/**
 * Modify the electrostatic potential at the boundary.
 *
 * @param up updater object.
 * @param moms_e first threee moments of the electron distribution.
 * @param m2par_e v_par^2 moment of the electron distribution.
 * @param moms_i first threee moments of the ion distribution.
 * @param m2par_i v_par^2 moment of the ion distribution.
 * @param phi_wall Wall potential.
 * @param phi Electrostatic potential.
 */
void gkyl_sheath_rarefaction_pot_advance(const struct gkyl_sheath_rarefaction_pot *up,
  const struct gkyl_array *moms_e, const struct gkyl_array *m2par_e,
  const struct gkyl_array *moms_i, const struct gkyl_array *m2par_i,
  const struct gkyl_array *phi_wall, struct gkyl_array *phi);

/**
 * Free memory associated with sheath_rarefaction_pot updater.
 *
 * @param up BC updater.
 */
void gkyl_sheath_rarefaction_pot_release(struct gkyl_sheath_rarefaction_pot *up);
]]

-- Boundary condition updater.
local SheathRarePot = Proto(UpdaterBase)

function SheathRarePot:init(tbl)
   SheathRarePot.super.init(self, tbl) -- Setup base object.

   self._grid    = assert(tbl.onGrid, "Updater.SheathRarePot: Must specify configuration space grid with 'onGrid'.")
   self._edge    = assert(tbl.edge, "Updater.SheathRarePot: Must specify edge to apply BCs with 'edge' (lower', 'upper').")
   local basis   = assert(tbl.basis, "Updater.SheathRarePot: Must specify the conf-space basis in 'basis'.")
   self._phiWall = assert(tbl.phiWall, "Updater.SheathRarePot: Must specify the wall potential in 'phiWall'.")
   local onField = assert(tbl.onField, "Updater.SheathRarePot: Must specify the field we'll apply BCs to in 'onField'.")
   local elem_q  = assert(tbl.elem_charge, "Updater.SheathRarePot: Must specify the elementary charge with 'elem_q'.")
   local mass_e  = assert(tbl.mass_e, "Updater.SheathRarePot: Must specify the electron mass with 'mass_e'.")
   local mass_i  = assert(tbl.mass_i, "Updater.SheathRarePot: Must specify the ion mass with 'mass_i'.")
   local useGPU  = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or 0)

   assert(self._edge == "lower" or self._edge == "upper", "Updater.SheathRarePot: 'edge' must be 'lower' or 'upper'.")

   local edge          = self._edge == 'lower' and 0 or 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
   local localExtRange = onField:localExtRange()
   local numGhostVec   = self._edge == 'lower' and onField._lowerGhostVec or onField._upperGhostVec

   self._zero = ffi.gc(ffiC.gkyl_sheath_rarefaction_pot_new(edge, localExtRange, numGhostVec:data(), basis._zero,
                                                            self._grid._zero, elem_q, mass_e, mass_i, useGPU),
                       ffiC.gkyl_sheath_rarefaction_pot_release)
end

function SheathRarePot:_advance(tCurr, inFld, outFld)
   local momsElc  = assert(inFld[1], "SheathRarePot.advance: Must-specify the electron moments.")
   local m2parElc = assert(inFld[2], "SheathRarePot.advance: Must-specify the electron parallel kinetic energy moment.")
   local momsIon  = assert(inFld[3], "SheathRarePot.advance: Must-specify the ion moments.")
   local m2parIon = assert(inFld[4], "SheathRarePot.advance: Must-specify the ion parallel kinetic energy moment.")
   local phi = assert(outFld[1], "SheathRarePot.advance: Must-specify the output electrostatic potential")
   ffiC.gkyl_sheath_rarefaction_pot_advance(self._zero, momsElc._zero, m2parElc._zero,
      momsIon._zero, m2parIon._zero, self._phiWall._zero, phi._zero)
end

function SheathRarePot:_advanceOnDevice(tCurr, inFld, outFld)
   local momsElc  = assert(inFld[1], "SheathRarePot.advance: Must-specify the electron moments.")
   local m2parElc = assert(inFld[2], "SheathRarePot.advance: Must-specify the electron parallel kinetic energy moment.")
   local momsIon  = assert(inFld[3], "SheathRarePot.advance: Must-specify the ion moments.")
   local m2parIon = assert(inFld[4], "SheathRarePot.advance: Must-specify the ion parallel kinetic energy moment.")
   local phi = assert(outFld[1], "SheathRarePot.advance: Must-specify the output electrostatic potential")
   ffiC.gkyl_sheath_rarefaction_pot_advance(self._zero, momsElc._zeroDevice, m2parElc._zeroDevice,
      momsIon._zeroDevice, m2parIon._zeroDevice, self._phiWall._zeroDevice, phi._zeroDevice)
end

return SheathRarePot
