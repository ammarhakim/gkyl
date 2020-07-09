-- Gkyl ------------------------------------------------------------------------
--
-- Updater to update five-moment source terms. This updater allows
-- both explicit and implicit updates.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

local cuda = nil
if GKYL_HAVE_CUDA then
   cuda = require "Cuda.RunTime"
end

local COL_PIV_HOUSEHOLDER_QR = 0;
local PARTIAL_PIV_LU = 1;

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
"new, copy, fill, sizeof, typeof, metatype")

-- Define C types for storing private data for use in updater
ffi.cdef [[
typedef struct {
  double charge, mass; /* Charge and mass */
  bool evolve;

  double qbym;
} FluidData_t;

typedef struct {
  int8_t nFluids; /* Number of fluids */
  double epsilon0; /* Permittivity of free space */
  double chi_e, chi_m; /* Propagation speed factor for electric field error potential */
  int8_t gravityDir; /* Direction of gravity force */
  double gravity; /* Gravitational acceleration */
  bool hasStatic, hasPressure; /* Flag to indicate if there is: static EB field, pressure */
  bool hasSigma; /* Flag to indicate if there is: sigma */
  bool hasAuxSrc;
  int8_t linSolType; /* Flag to indicate linear solver type for implicit method */
} MomentSrcData_t;

  void gkylFiveMomentSrcRk3(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm, double *sigma, double *auxSrc);
  void gkylFiveMomentSrcTimeCentered(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm, double *sigma, double *auxSrc);
  void gkylFiveMomentSrcTimeCenteredDirect(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm, double *sigma, double *auxSrc);
  void gkylFiveMomentSrcTimeCenteredDirect2(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm, double *sigma, double *auxSrc);
  void gkylFiveMomentSrcExact(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma, double *auxSrc);

  /* CUDA GPU */
  /* Opaque structure holding CUBLAS library context */
  struct cublasContext;
  typedef struct cublasContext *cublasHandle_t;

  typedef struct {
    double *d_lhs;
    double *d_rhs;
    double **d_lhs_ptr;
    double **d_rhs_ptr;
    int *d_info;
    cublasHandle_t handle;
  } GkylMomentSrcDeviceData_t;

  GkylMomentSrcDeviceData_t *cuda_gkylMomentSrcInit(
    const char *scheme, const int nFluids, const int numBlocks, const int numThreads);
  void cuda_gkylMomentSrcDestroy(const GkylMomentSrcDeviceData_t *context);
  void momentSrcAdvanceOnDevice(
      const MomentSrcData_t *sd, const FluidData_t *fd, const double dt,
      GkylCartField_t **fluidFlds, GkylCartField_t *emFld,
      const GkylMomentSrcDeviceData_t *context);
]]

-- Explicit, SSP RK3 scheme
local function updateSrcRk3(self, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
   ffi.C.gkylFiveMomentSrcRk3(self._sd, self._fd, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
end

-- Use an explicit scheme to update momentum and electric field: this
-- is an extension of the standard Boris push algorithm, in which the
-- half time-step electric field update is replaced by an implicit
-- step in which both the velocity and electric are updated.
local function updateSrcModBoris(self, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
   print("updateSrcModBoris")
end

-- Use an implicit scheme to update momentum and electric field
local function updateSrcTimeCentered(self, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
   ffi.C.gkylFiveMomentSrcTimeCentered(self._sd, self._fd, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
end

-- Use an implicit scheme to update momentum and electric field
local function updateSrcTimeCenteredDirect(self, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
   ffi.C.gkylFiveMomentSrcTimeCenteredDirect(self._sd, self._fd, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
end

-- Use an implicit scheme to update momentum and electric field
local function updateSrcTimeCenteredDirect2(self, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
   ffi.C.gkylFiveMomentSrcTimeCenteredDirect2(self._sd, self._fd, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
end

-- Use an implicit scheme to update momentum and electric field
local function updateSrcExact(self, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
   ffi.C.gkylFiveMomentSrcExact(self._sd, self._fd, dt, fPtr, emPtr, staticEmPtr, sigmaPtr, auxSrc)
end

-- Five-moment source updater object
local FiveMomentSrc = Proto(UpdaterBase)

function FiveMomentSrc:init(tbl)
   FiveMomentSrc.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid, "Updater.FiveMomentSrc: Must provide grid object using 'onGrid'")

   self._sd = ffi.new(ffi.typeof("MomentSrcData_t"))   
   -- Read in solver parameters
   self._sd.nFluids = assert(tbl.numFluids, "Updater.FiveMomentSrc: Must specify number of fluid using 'numFluids'")

   assert(#tbl.charge == self._sd.nFluids, "Updater.FiveMomentSrc: Charge table must have " .. self._sd.nFluids .. " elements.")
   assert(#tbl.mass == self._sd.nFluids, "Updater.FiveMomentSrc: Mass table must have " .. self._sd.nFluids .. " elements.")

   self._sd.epsilon0 = assert(tbl.epsilon0, "Updater.FiveMomentSrc: Must specify 'epsilon0'")
   self._sd.chi_e = tbl.elcErrorSpeedFactor and tbl.elcErrorSpeedFactor or 0.0
   self._sd.chi_m = tbl.mgnErrorSpeedFactor and tbl.mgnErrorSpeedFactor or 1.0

   self._sd.gravity = tbl.gravity and tbl.gravity or 0.0
   self._sd.gravityDir = tbl.dir and tbl.dir or 1 -- by default gravity acts in X direction
   
   self._sd.hasStatic = tbl.hasStaticField ~= nil and tbl.hasStaticField or false
   self._sd.hasPressure = tbl.hasPressure ~= nil and tbl.hasPressure or true
   self._sd.hasSigma = tbl.hasSigmaField ~= nil and tbl.hasSigmaField or false

   self._sigmaFld = tbl.sigmaField
   
   self._sd.hasAuxSrc = tbl.hasAuxSourceFunction ~= nil and tbl.hasAuxSourceFunction or false
   self._auxSrcFunc = tbl.auxSourceFunction
   
   local linSolType = tbl.linSolType and tbl.linSolType or "partialPivLu"
   if linSolType == "partialPivLu" then
      self._sd.linSolType = PARTIAL_PIV_LU
   elseif linSolType == "colPivHouseholderQr" then
      self._sd.linSolType = COL_PIV_HOUSEHOLDER_QR
   else
     assert(false, string.format("linSolType %s not supported", linSolType))
   end

   self._qbym = {}
   self._fd = ffi.new("FluidData_t[?]", self._sd.nFluids)
   -- store charge and mass for each fluid
   for n = 1, self._sd.nFluids do
      self._fd[n-1].charge = tbl.charge[n]
      self._fd[n-1].mass = tbl.mass[n]
      self._fd[n-1].qbym = self._fd[n-1].charge / self._fd[n-1].mass
      -- self._fd[n-1].evolve = tbl.evolve ~= nil and tbl.evolve[n] or true
      if (tbl.evolve ~= nil) then
        self._fd[n-1].evolve = tbl.evolve[n]
      else
        self._fd[n-1].evolve = true
      end
      self._qbym[n] = self._fd[n-1].qbym
   end

   -- scheme is one of "ssp-rk3", "modified-boris"  or "time-centered"
   local scheme = tbl.scheme and tbl.scheme or "time-centered"
   self._updateSrc = nil
   if scheme == "ssp-rk3" then
      self._updateSrc = updateSrcRk3
   elseif scheme == "modified-boris" then
      self._updateSrc = updateSrcModBoris
   elseif scheme == "time-centered" then
      self._updateSrc = updateSrcTimeCentered
   elseif scheme == "time-centered-direct" or scheme == "direct" then
      self._updateSrc = updateSrcTimeCenteredDirect
   elseif scheme == "time-centered-direct2" or scheme == "direct2" then
      self._updateSrc = updateSrcTimeCenteredDirect2
   elseif scheme == "exact" then
      self._updateSrc = updateSrcExact
   else
      assert(false, string.format("scheme %s not supported", scheme))
   end
   self.scheme = scheme
end


function FiveMomentSrc:initDevice(tbl)
   local scheme = self.scheme
   self.gpu_scheme = "time-centered"
   if (scheme=="time-centered" or scheme=="time-centered-direct" or
      scheme=="direct") then
      self.gpu_scheme = scheme
   else
      print(string.format("scheme %s is not supported on gpu", scheme))
      print(string.format("falling back to %s", self.gpu_scheme))
   end

   local nFluids = self._sd.nFluids

   local sz_sd = sizeof("MomentSrcData_t")
   self.sd_onDevice, err = cuda.Malloc(sz_sd)
   cuda.Memcpy(self.sd_onDevice, self._sd, sz_sd, cuda.MemcpyHostToDevice)

   local sz_fd = sizeof("FluidData_t") * nFluids
   self.fd_onDevice, err = cuda.Malloc(sz_fd)
   cuda.Memcpy(
      self.fd_onDevice, self._fd, sz_fd , cuda.MemcpyHostToDevice)

   local sz_fluidFlds = sizeof("GkylCartField_t*") * nFluids
   self.d_fluidFlds, err = cuda.Malloc(sz_fluidFlds)

   self.numThreads = tbl.numThreads or GKYL_DEFAULT_NUM_THREADS

   local numCellsLocal = self._onGrid:localRange():volume()
   local numThreads = math.min(self.numThreads, numCellsLocal)
   local numBlocks  = math.ceil(numCellsLocal/numThreads)
   self.device_context = ffi.C.cuda_gkylMomentSrcInit(
      self.gpu_scheme, nFluids, numBlocks, numThreads)
   self.numThreads = numThreads
   self.numBlocks = numBlocks
   self.numCellsLocal = numCellsLocal

   self.first = true

   -- callback into proxy.__gc when the proxy becomes free
   -- FIXME Is this the correct way?
   local prox = newproxy(true)
   self.deviceContextDestroyed = false
   getmetatable(prox).__gc = function()
      if not self.deviceContextDestroyed  then
        ffi.C.cuda_gkylMomentSrcDestroy(self.device_context)
     end
     self.deviceContextDestroyed = true
   end
   self[prox] = true
end


function FiveMomentSrc:_advanceDispatch(tCurr, inFld, outFld, target)
   local grid = self._onGrid
   local dt = self._dt
   local nFluids = self._sd.nFluids

   local ndim = grid:ndim()
   local xc = Lin.Vec(ndim) -- cell center
   
   -- check if correct number of inputs were provided
   assert(#outFld == nFluids+1,
	  "Must have exactly " .. nFluids+1 .. " output fields. Provided " .. #outFld .. " instead.")
   local emFld = outFld[#outFld] -- EM field

   local staticEmFld
   if (self._sd.hasStatic) then
      staticEmFld = inFld[1] -- static EM field
   end

   -- make list of indexer functions
   local fIdxr = {}
   for i = 1, nFluids do
      fIdxr[i] = outFld[i]:genIndexer()
   end
   local emIdxr = emFld:genIndexer()
   local staticEmIdxr
   if (self._sd.hasStatic) then
      staticEmIdxr = staticEmFld:genIndexer()
   end

   local sigmaIdxr
   if (self._sd.hasSigma) then
      sigmaIdxr = self._sigmaFld:genIndexer()
   end

   if target=="gpu" then
      if self.first then
         local fluidFlds = ffi.new("GkylCartField_t*[?]", nFluids)
         for n = 1, nFluids do
            fluidFlds[n-1] = outFld[n]._onDevice
         end
         local sz_fluidFlds = sizeof("GkylCartField_t*") * self._sd.nFluids
         cuda.Memcpy(
            self.d_fluidFlds, fluidFlds, sz_fluidFlds, cuda.MemcpyHostToDevice)
         self.d_emFld = emFld._onDevice
         self.first = false
      end

      local numCellsLocal = emFld:localRange():volume()
      assert(numCellsLocal == self.numCellsLocal)

      ffi.C.momentSrcAdvanceOnDevice(
         self.sd_onDevice, self.fd_onDevice, dt, self.d_fluidFlds, self.d_emFld,
         self.device_context)

      return true, GKYL_MAX_DOUBLE
   elseif target=="cpu" then
      -- allocate stuff to pass to C
      local fDp = ffi.new("double*[?]", nFluids)
      local emDp = ffi.new("double*")
      local staticEmDp = ffi.new("double*")
      local sigmaDp = ffi.new("double*")
      local auxSrcDp = ffi.new("double[?]", 3*(nFluids+1))

      -- loop over local range, updating source in each cell
      for idx in emFld:localRangeIter() do
         grid:setIndex(idx)
         grid:cellCenter(xc)

         -- set pointers to fluids and field
         for i = 1, nFluids do
       fDp[i-1] = outFld[i]:getDataPtrAt(fIdxr[i](idx))
         end
         emDp = emFld:getDataPtrAt(emIdxr(idx))
         if (self._sd.hasStatic) then
            staticEmDp = staticEmFld:getDataPtrAt(staticEmIdxr(idx))
         end
         if (self._sd.hasSigma) then
            sigmaDp = self._sigmaFld:getDataPtrAt(sigmaIdxr(idx))
         end

         if (self._sd.hasAuxSrc) then
            self:_auxSrcFunc(
               xc, tCurr, self._sd.epsilon0, self._qbym, fDp, emDp, auxSrcDp)
         end

         -- update sources
         self._updateSrc(self, dt, fDp, emDp, staticEmDp, sigmaDp, auxSrcDp)
      end

      return true, GKYL_MAX_DOUBLE
   else
      assert(false, string.format("target %s not recognized", target))
   end
end

-- advance method
function FiveMomentSrc:_advance(tCurr, inFld, outFld)
   return self:_advanceDispatch(tCurr, inFld, outFld, "cpu")
end

-- advance method
function FiveMomentSrc:_advanceDevice(tCurr, inFld, outFld)
   return self:_advanceDispatch(tCurr, inFld, outFld, "gpu")
end

return FiveMomentSrc
