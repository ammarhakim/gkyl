-- System libraries
local xsys = require "xsys"

-- Gkyl libraries.
local Proto       = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local DataStruct  = require "DataStruct"
local ffi         = require "ffi"
local Lin          = require "Lib.Linalg"

local ffiC = ffi.C
require "Lib.ZeroUtil"

-- Declaration of gkylzero objects and functions.
ffi.cdef [[
// Object type
typedef struct gkyl_geo_gyrokinetic gkyl_geo_gyrokinetic;

struct gkyl_range* gkyl_geo_gyrokinetic_get_nrange(gkyl_geo_gyrokinetic* geo);
struct gkyl_array* gkyl_geo_gyrokinetic_get_mc2p_nodal_fd(gkyl_geo_gyrokinetic* geo);
double* gkyl_geo_gyrokinetic_get_dzc(gkyl_geo_gyrokinetic* geo);


// Type of flux surface
enum gkyl_geo_gyrokinetic_type {
  GKYL_SOL_DN, // SOL of double-null configuration
  GKYL_SOL_SN, // SOL of single-null configuration
  GKYL_PF, // Private flux region
  GKYL_CORE // Core (closed flux-surface)
};  

// Inputs to create a new GK geometry creation object
struct gkyl_geo_gyrokinetic_inp {
  // psiRZ and related inputs  
  const struct gkyl_rect_grid *rzgrid; // RZ grid on which psi(R,Z) is defined
  const struct gkyl_basis *rzbasis; // basis functions for R,Z grid
  const struct gkyl_array *psiRZ; // psi(R,Z) DG representation
  const struct gkyl_range *rzlocal; // local range over which psiRZ is defined
  double B0; // Toroidal Field on axis
  double R0; // Axis

  // Parameters for root finder: leave unset to use defaults
  struct {
    int max_iter; // typically 20
    double eps; // typically 1e-10
  } root_param;

  // Parameters for nmumerical quadrature: leave unset to use default
  struct {
   int max_levels; // typically 6-7    
    double eps; // typically 1e-10
  } quad_param;
};

// Inputs to create geometry for a specific computational grid
struct gkyl_geo_gyrokinetic_geo_inp {
  const struct gkyl_rect_grid *cgrid;
  int* bcs;
  const struct gkyl_basis *cbasis;

  enum gkyl_geo_gyrokinetic_type ftype; // type of geometry
  
  double rclose; // closest R to discrimate
  double rleft; // closest R to discrimate
  double rright; // closest R to discrimate
  double zmin, zmax; // extents of Z for integration

  bool write_node_coord_array; // set to true if nodal coordinates should be written
  const char *node_file_nm; // name of nodal coordinate file
};



// Some cumulative statistics
struct gkyl_geo_gyrokinetic_stat {
  long nquad_cont_calls; // num calls from quadrature
  long nroot_cont_calls; // num calls from root-finder
};  

/**
 * Create new updater to compute the geometry (mapc2p) needed in GK
 * simulations.
 *
 * @param inp Input parameters
 * @param New GK geometry updater
 */
gkyl_geo_gyrokinetic *gkyl_geo_gyrokinetic_new(const struct gkyl_geo_gyrokinetic_inp *inp);

/**
 * Get R(psi,Z) for a specified psi and Z value. Multiple values may
 * be returned (or none). The R(psi,Z) and dR/dZ are stored in the R
 * and dR arrays which be allocated by the caller.
 *
 * @param geo Geometry object
 * @param psi Psi value
 * @param Z Z value
 * @param nmaxroots Maximum number of roots
 * @param R on output, R(psi,Z)
 * @param dR on output, dR/dZ
 */
int gkyl_geo_gyrokinetic_R_psiZ(const gkyl_geo_gyrokinetic *geo, double psi, double Z, int nmaxroots,
  double *R, double *dR);

/**
 * Integrate along a specified psi countour and return its length. The
 * contour must lie completely inside the RZ domain of the psiRZ DG
 * field. The @a rclose parameter is used to select amongst the
 * multiple possible countours with the same psi. Foe example, to
 * select a flux surface on the outboard side of a double-null
 * configuration choose rclose to be Rmax.
 *
 * @param geo Geometry object
 * @param psi Psi value of contour
 * @param zmin Starting z location
 * @param zmax Ending z location
 * @param rclose Value of radial coordinate to discrimate between multiple
 *    contours
 * @return Length of contour
 */
double gkyl_geo_gyrokinetic_integrate_psi_contour(const gkyl_geo_gyrokinetic *geo, double psi,
  double zmin, double zmax, double rclose);

/**
 * Compute physical coordinates (mapc2p)  given computational coordinates
 *
 * @param geo Geometry object
 * @param xn computational coordinates
 * @param ret physical coordinates
 */
void gkyl_geo_gyrokinetic_mapc2p(const gkyl_geo_gyrokinetic *geo, const struct gkyl_geo_gyrokinetic_geo_inp *inp,
    const double *xn, double *ret);

/**
 * Compute geometry (mapc2p) on a specified computational grid. The
 * output array must be pre-allocated by the caller.
 *
 * @param geo Geometry object
 * @param ginp Input structure for creating mapc2p
 * @param mapc2p On output, the DG representation of mapc2p
 */
void gkyl_geo_gyrokinetic_calcgeom(const gkyl_geo_gyrokinetic *geo,
  const struct gkyl_geo_gyrokinetic_geo_inp *ginp, struct gkyl_array *mapc2p, struct gkyl_range *conversion_range);

/**
 * Return cumulative statistics from geometry computations
 *
 * @param geo Geometry object
 * @return Cumulative statistics
 */
struct gkyl_geo_gyrokinetic_stat gkyl_geo_gyrokinetic_get_stat(const gkyl_geo_gyrokinetic *geo);

/**
 * Delete updater.
 *
 * @param geo Geometry object to delete
 */
void gkyl_geo_gyrokinetic_release(gkyl_geo_gyrokinetic *geo);

// Object type
typedef struct gkyl_calc_bmag gkyl_calc_bmag;
typedef struct bmag_ctx bmag_ctx;

/**
 * Create new updater to compute the metric coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_bmag* 
gkyl_calc_bmag_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_rect_grid *cgrid, const struct gkyl_rect_grid *pgrid, const gkyl_geo_gyrokinetic *app, const struct gkyl_geo_gyrokinetic_geo_inp *ginp, bool use_gpu);


/**
 * Advance calc_metric (compute the metric coefficients).
 *
 * @param up calc_metric updater object.
 * @param crange Config-space range.
 * @param XYZ field containing DG rep of cartesian coordinates
 * @param gFld output field where metric coefficients will be placed
 */

void gkyl_calc_bmag_advance(const gkyl_calc_bmag *up, const struct gkyl_range *crange, const struct gkyl_range *crange_ext,
     const struct gkyl_range *prange, const struct gkyl_range *prange_ext, struct gkyl_array *psidg, struct gkyl_array *psibyrdg, struct gkyl_array *psibyr2dg, struct gkyl_array *bphidg, struct gkyl_array* bmag_compdg, struct gkyl_array* mapc2p);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_bmag_release(gkyl_calc_bmag* up);


// Object type
typedef struct gkyl_calc_metric gkyl_calc_metric;

/**
 * Create new updater to compute the metric coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param grid configuration space grid.
 * @param bcs 0 for periodic 1 for non-periodic. Determines whether mapc2p extends to ghost cells
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_metric* gkyl_calc_metric_new(const struct gkyl_basis *cbasis,
  const struct gkyl_rect_grid *grid, const int *bcs, bool use_gpu);

/**
 * Use finite differences to calculate metric coefficients at nodes
 *
 * @param up calc_metric updater object.
 * @param nrange nodal range.
 * @param mc2p_nodal_fd nodal array containing cartesian coordinates at nodes and nearby nodes used for FD
 * @param gFld output field where metric coefficients will be placed
 */
void gkyl_calc_metric_advance(gkyl_calc_metric *up, struct gkyl_range *nrange, struct gkyl_array *mc2p_nodal_fd, double *dzc, struct gkyl_array *gFld, struct gkyl_range *update_range);
//void gkyl_calc_metric_advance(const gkyl_calc_metric *up, const struct gkyl_range *crange, struct gkyl_array *XYZ, struct gkyl_array *gFld);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_metric_release(gkyl_calc_metric* up);
// Object type
typedef struct gkyl_calc_derived_geo gkyl_calc_derived_geo;

/**
 * Create new updater to compute the derived_geo coefficients
 *
 * @param cbasis Basis object (configuration space).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_derived_geo* gkyl_calc_derived_geo_new(const struct gkyl_basis *cbasis,
  const struct gkyl_rect_grid *grid, bool use_gpu);

/**
 * Advance calc_derived_geo (compute the derived_geo coefficients).
 *
 * @param up calc_derived_geo updater object.
 * @param crange Config-space range.
 * @param gFld field containing DG rep of the metric coefficients
 * @param jFld output field where jacobian will be placed
 */

void gkyl_calc_derived_geo_advance(const gkyl_calc_derived_geo *up, const struct gkyl_range *crange,
    struct gkyl_array *gFld, struct gkyl_array *bmagFld, struct gkyl_array *jFld, struct gkyl_array *jinvFld,
    struct gkyl_array *grFld, struct gkyl_array *biFld, struct gkyl_array *cmagFld, struct gkyl_array *jtotFld, struct gkyl_array *jtotinvFld, struct gkyl_array *bmaginvFld, struct gkyl_array *bmaginvsqFld,struct gkyl_array *gxxJFld,  struct gkyl_array *gxyJFld, struct gkyl_array *gyyJFld);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_derived_geo_release(gkyl_calc_derived_geo* up);

]]

-- Boundary condition updater.
local GeometryGyrokinetic = Proto(UpdaterBase)

function GeometryGyrokinetic:init(tbl)
   GeometryGyrokinetic.super.init(self, tbl) -- Setup base object.
   self.grid = tbl.grid
   self.rzGrid = tbl.rzGrid
   self.basis = tbl.basis
   self.rzBasis = tbl.rzBasis
   self.localRange = tbl.localRange
   self.localRangeExt = tbl.localRangeExt
   self.rzLocalRange = tbl.rzLocalRange
   self.rzLocalRangeExt = tbl.rzLocalRangeExt
   self.B0 = tbl.B0
   self.R0 = tbl.R0
   self.calcGeom = tbl.calcGeom
   self.calcBmag = tbl.calcBmag
   self.zmin = tbl.zmin
   self.zmax = tbl.zmax
   self.rclose = tbl.rclose


   


   -- Set up BCs 0=periodic, 1=nonperiodic
   -- BCs are used to set ranges for calculating geometry
   -- and to choose stencils when calculating metrics
   -- Currently we have not thought of a case where we need to calculate geometry with periodic bcs
   self.BCs = Lin.IntVec(self.grid:ndim())
   for i = 1, self.grid:ndim() do
      self.BCs[i] = 1
   end
   self.BCs[2] = 0
   --if self.BCs[1] == 1 then
   --   local lower = self.localRangeExt:lowerAsVec()
   --   local upper = self.localRangeExt:upperAsVec()
   --   lower[1] = lower[1]+1
   --   upper[1] = upper[1]-1
   --   lower[3] = lower[3]+1
   --   upper[3] = upper[3]-1
   --   self.conversionRange = self.localRangeExt:subRange(lower, upper)
   --else 
   --   self.conversionRange = self.localRange
   --end
   self.conversionRange = self.localRange
   

   if self.calcGeom then
      -- Create input struct
      self.inp = ffi.new("struct gkyl_geo_gyrokinetic_inp")
      self.inp.rzgrid = self.rzGrid._zero
      self.inp.rzbasis = self.rzBasis._zero
      self.inp.psiRZ = tbl.psiRZ._zero
      self.inp.rzlocal = self.rzLocalRange
      self.inp.B0 = self.B0
      self.inp.R0 = self.R0

      -- Create mapc2p updater
      self._zero_geom = ffi.gc(ffiC.gkyl_geo_gyrokinetic_new(self.inp), ffiC.gkyl_geo_gyrokinetic_release)

      -- Create geo input struct which will be used by mapc2p updater
      self.ginp = ffi.new("struct gkyl_geo_gyrokinetic_geo_inp")
      self.ginp.cgrid = self.grid._zero
      self.ginp.cbasis = self.basis._zero
      --self.ginp.ftype = GKYL_SOL_DN
      self.ginp.rclose = self.rclose --self.rzGrid:lower(1)
      self.ginp.zmin = self.zmin or self.rzGrid:lower(2)
      self.ginp.zmax = self.zmax or self.rzGrid:upper(2)
      self.ginp.write_node_coord_array = true
      self.ginp.node_file_nm = "grid.gkyl"
      self.ginp.bcs = self.BCs:data()
      self._zero_bmag = ffi.gc(ffiC.gkyl_calc_bmag_new(self.basis._zero, self.rzBasis._zero, self.grid._zero, self.rzGrid._zero, self.geo, self.ginp, false), ffiC.gkyl_calc_bmag_release)
   end


   self._zero_metric = ffi.gc(ffiC.gkyl_calc_metric_new(self.basis._zero, self.grid._zero, self.BCs:data(), false), ffiC.gkyl_calc_metric_release)
   self._zero_derived = ffi.gc(ffiC.gkyl_calc_derived_geo_new(self.basis._zero, self.grid._zero, false), ffiC.gkyl_calc_derived_geo_release)



end

function GeometryGyrokinetic:_advance(tCurr, inFlds, outFlds)

   local psiRZ, psibyrRZ, psibyr2RZ, bphiRZ = inFlds[1], inFlds[2], inFlds[3], inFlds[4]
   local mapc2p_field, gFld, jacobGeo, jacobGeoInv, jacobTot, jacobTotInv, bmagInv, bmagInvSq, gxxJ, gxyJ, gyyJ, grFld, b_i, cmag, bmag = outFlds[1], outFlds[2], outFlds[3], outFlds[4], outFlds[5], outFlds[6], outFlds[7], outFlds[8], outFlds[9], outFlds[10], outFlds[11], outFlds[12], outFlds[13], outFlds[14], outFlds[15]


   -- Fill mapc2p
   if self.calcGeom then 
      ffiC.gkyl_geo_gyrokinetic_calcgeom(self._zero_geom, self.ginp, mapc2p_field._zero, self.conversionRange)
   end

   -- Sync before calculating metric
   mapc2p_field:sync(false)

   -- Fill bmag
   if self.calcBmag then
      ffiC.gkyl_calc_bmag_advance(self._zero_bmag, self.localRange, self.localRangeExt, self.rzLocalRange, self.rzLocalRangeExt, psiRZ._zero, psibyrRZ._zero, psibyr2RZ._zero, bphiRZ._zero,  bmag._zero, mapc2p_field._zero)
   end

   -- Fill metrics
   --ffiC.gkyl_calc_metric_advance(self._zero_metric, self.localRange, mapc2p_field._zero, gFld._zero)
   ffiC.gkyl_calc_metric_advance(self._zero_metric, ffiC.gkyl_geo_gyrokinetic_get_nrange(self._zero_geom), ffiC.gkyl_geo_gyrokinetic_get_mc2p_nodal_fd(self._zero_geom), ffiC.gkyl_geo_gyrokinetic_get_dzc(self._zero_geom), gFld._zero, self.conversionRange)

   -- Fill derived geo. Includes adjustment  of g_zz to preserve cmag = const.
   -- Need to add argument about whether or not to do the adjustment
   ffiC.gkyl_calc_derived_geo_advance(self._zero_derived, self.localRange, gFld._zero, bmag._zero, jacobGeo._zero, jacobGeoInv._zero, grFld._zero, b_i._zero , cmag._zero, jacobTot._zero, jacobTotInv._zero, bmagInv._zero, bmagInvSq._zero, gxxJ._zero, gxyJ._zero, gyyJ._zero)
end

return GeometryGyrokinetic
 
