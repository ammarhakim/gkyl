-- Gkyl ------------------------------------------------------------------------
--
-- Finite volume wave propagation updater on rectangular, Cartesian
-- grid. Non-uniform grids are supported and an embedded BC can also
-- be included. Algorithm based on LeVeque book and Hakim
-- et. al. papers. Original design from my (AHH) thesis code,
-- Miniwarpx and Warpx which significantly extend ideas present in
-- CLAWPACK. This implementation in LuaJIT is, perhaps, completely
-- novel and unique.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"
local DataStruct = require "DataStruct"

local function printf(...) io.write(string.format(...)) end

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
"new, copy, fill, sizeof, typeof, metatype")

local cuda = nil
if GKYL_HAVE_CUDA then
   cuda = require "Cuda.RunTime"
end

ffi.cdef [[ 
  typedef struct GkylEquationFv_t GkylEquationFv_t;
  typedef struct {
      int updateDirs[6];
      int numUpdateDirs;
      double dt;
      double _cfl;
      double _cflm;
      GkylEquationFv_t *equation;
      GkylCartField_t *dtByCell;
  } GkylWavePropagation_t;
    
  void wavePropagationAdvanceOnDevice(
    const int meqn, const int mwave, const int numBlocks, const int numThreads,
    GkylWavePropagation_t *hyper, GkylCartField_t *qIn, GkylCartField_t *qOut);

  void setDt(GkylWavePropagation_t *hyper, double dt);


/* gkylzero components */
typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

struct gkyl_wave_geom {
  struct gkyl_range range; // range over which geometry is defined
  struct gkyl_array *geom; // geometry in each cell
  struct gkyl_ref_count ref_count;
};

struct gkyl_wave_geom *gkyl_wave_geom_new(const struct gkyl_rect_grid *grid,
  struct gkyl_range *range,
  evalf_t mapc2p, void *ctx);

enum gkyl_wave_limiter {
  GKYL_NO_LIMITER = 1, // to allow default to be 0
  GKYL_MIN_MOD,
  GKYL_SUPERBEE,
  GKYL_VAN_LEER,
  GKYL_MONOTONIZED_CENTERED,
  GKYL_BEAM_WARMING,
  GKYL_ZERO
};

struct gkyl_wave_prop_status {
  int success; // 1 if step worked, 0 otherwise
  double dt_suggested; // suggested time-step
};

typedef struct gkyl_wave_prop gkyl_wave_prop;

struct gkyl_wave_prop_inp {
  const struct gkyl_rect_grid *grid; // grid on which to solve equations
  const struct gkyl_wv_eqn *equation; // equation solver
  enum gkyl_wave_limiter limiter; // limiter to use
  int num_up_dirs; // number of update directions
  int update_dirs[7]; // directions to update FIXME GKYL_MAX_DIM=7
  double cfl; // CFL number to use

  const struct gkyl_wave_geom *geom; // geometry
};

gkyl_wave_prop* gkyl_wave_prop_new(struct gkyl_wave_prop_inp winp);

struct gkyl_wave_prop_status gkyl_wave_prop_advance(const gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *qin, struct gkyl_array *qout);

void gkyl_wave_prop_release(gkyl_wave_prop* up);

void gkyl_wave_geom_release(const struct gkyl_wave_geom* wg);
]]

-- Template for function to compute jump
local calcDeltaTempl = xsys.template([[
return function (ql, qr, delta)
   |for i = 1, MEQN do
   delta[${i}] = qr[${i}] - ql[${i}]
   |end
end
]])

-- Template for function to compute maximum CFL number
local calcCflaTempl = xsys.template([[
return function (cfla, dtdx, s)
   local c = cfla
   |for i = 1, MWAVE do
   c = math.max(c, dtdx*math.abs(s[${i}]))
   |end
   return c
end
]])

-- Template for function to compute first-order Gudonov update
local calcFirstOrderGudTempl = xsys.template([[
return function (dtdx, ql, qr, amdq, apdq)
   |for i = 1, MEQN do
   qr[${i}] = qr[${i}] - dtdx*apdq[${i}]
   |end
   |for i = 1, MEQN do
   ql[${i}] = ql[${i}] - dtdx*amdq[${i}]
   |end
end
]])

-- Template for function to compute dot product of waves
local waveDotProdTempl = xsys.template([[
return function (meqn, waves, waves1, mw)
   local mw1 = mw-1
   return
   |for i = 0, MEQN-2 do
   waves1[${MEQN}*mw1+${i}]*waves[${MEQN}*mw1+${i}]+
   |end
   waves1[${MEQN}*mw1+${MEQN-1}]*waves[${MEQN}*mw1+${MEQN-1}]
end
]])

-- Template for function to rescale waves
local rescaleWaveTempl = xsys.template([[
return function (scale, wave)
   |for i = 0, MEQN-1 do
   wave[${i}] = scale*wave[${i}]
   |end
end
]])

local secondOrderFluxTempl = xsys.template([[
return function (dtdx, s, wave, fs)
   local sfact = 0.5*math.abs(s)*(1-math.abs(s)*dtdx)
   |for i = 0, MEQN-1 do
   fs[${i}] = fs[${i}] + sfact*wave[${i}]
   |end
end
]])

local secondOrderUpdateTempl = xsys.template([[
return function (dtdx, fs, fs1, q)
   | for i = 1, MEQN do
   q[${i}] = q[${i}] - dtdx*(fs1[${i-1}]-fs[${i-1}])
   | end
end
]])

-- limiter functions: See LeVeque book
local limiterFunctions = {}
limiterFunctions["no-limiter"] = function (r)
   return 1
end
limiterFunctions["min-mod"] = function (r)
   return math.max(0, math.min(1, r))
end
limiterFunctions["superbee"] = function (r)
   return math.max(0.0, math.min(1, 2*r), math.min(2.0, r))
end
limiterFunctions["van-leer"] = function (r)
   return (r+math.abs(r))/(1+math.abs(r))
end
limiterFunctions["monotonized-centered"] = function (r)
   local c = (1.0+r)/2
   return math.max(0.0, math.min(c, 2, 2*r))
end
limiterFunctions["beam-warming"] = function (r)
   return r
end
limiterFunctions["zero"] = function (r)
   return 0
end

limiter_id = {}
limiter_id["no-limiter"] = ffi.C.GKYL_NO_LIMITER
limiter_id["min-mod"] = ffi.C.GKYL_MIN_MOD
limiter_id["superbee"] = ffi.C.GKYL_SUPERBEE
limiter_id["van-leer"] = ffi.C.GKYL_VAN_LEER
limiter_id["monotonized-centered"] = ffi.C.GKYL_MONOTONIZED_CENTERED
limiter_id["beam-warming"] = ffi.C.GKYL_BEAM_WARMING
limiter_id["zero"] = ffi.C.GKYL_ZERO

-- Helper object for indexing 1D slice data. The slice spans from
-- [lower, upper] (inclusive) and has `stride` pieces of data stored
-- at each location.
local createSliceData = function (dtype)
   local slice_mt = {
      __new = function (self, lower, upper, stride)
         local n = upper-lower+1
         local v = new(self)
         v._data = Alloc.malloc(typeof(dtype)*n*stride)
         v._sz, v._stride, v._lower = n*stride, stride, lower
         return v
      end,
      __index = function (self, k)
         return self._data+(k-self._lower)*self._stride
      end,
      __gc = function (self)
         Alloc.free(self._data)
      end,
   }
   return metatype(typeof(string.format("struct {int _sz, _stride, _lower; %s *_data; }", dtype)), slice_mt)
end
local SliceData = createSliceData("double")
local SliceDataBool = createSliceData("bool")

-- helper function to zero out contents of SliceData
local function clearSliceData(sd, dtype)
   fill(sd._data, sd._sz*sizeof("double"))
end

-- Wave-propagation updater object
local WavePropagation = Proto(UpdaterBase)

-- constructor
function WavePropagation:init(tbl)
   WavePropagation.super.init(self, tbl) -- setup base object

   self._isFirst = true -- will be reset first time _advance() is called

   -- read data from input table
   self._onGrid = assert(tbl.onGrid, "Updater.WavePropagation: Must provide grid object using 'onGrid'")
   self._equation = assert(tbl.equation, "Updater.WavePropagation: Must provide equation object using 'equation'")
   self._limiterFunc = assert(limiterFunctions[tbl.limiter], "Updater.WavePropagation: Must specify limiter to use")
   self._limiterFuncZero = limiterFunctions["zero"]

   self._ndim = self._onGrid:ndim()
   self._cfl = assert(tbl.cfl, "Updater.WavePropagation: Must specify CFL number using 'cfl'")
   self._cflm = tbl.cflm and tbl.cflm or 1.1*self._cfl
   self._hasSsBnd = xsys.pickBool(tbl.hasSsBnd, false)
   self._inOut = tbl.inOut

   self._updateDirs = {}
   local upDirs = tbl.updateDirections and tbl.updateDirections or {1, 2, 3, 4, 5, 6}
   for d = 1, #upDirs do
      self._updateDirs[d] = upDirs[d] -- update directions
   end

   local meqn, mwave = self._equation:numEquations(), self._equation:numWaves()
   self._localRange = tbl.onGrid:localRange()
   local localRange = self._localRange

   -- allocate space for storing 1D slice data
   self.wavesSlice, self.speedsSlice, self.fsSlice = {}, {}, {}
   for d = 1, self._ndim do
      local l, u = localRange:lower(d)-2, localRange:upper(d)+2
      self.wavesSlice[d] = SliceData(l, u, meqn*mwave)
      self.speedsSlice[d] = SliceData(l, u, mwave)
      self.fsSlice[d] = SliceData(l, u, meqn)
   end

   -- masks for edge/cell locations regarding stari-stepped boundaries
   self.onSsBnd = {}  -- if exactly one of the edge's left/right cells is outside
   self.bothOutSsBnd = {} -- if both of the edge's left/right cells are outside
   self.thisOutSsBnd = {}  -- if this cell is outside
   if (self._hasSsBnd) then
      for d = 1, self._ndim do
         local l, u = localRange:lower(d)-2, localRange:upper(d)+2
         self.onSsBnd[d] = SliceDataBool(l, u, 1)
         self.bothOutSsBnd[d] = SliceDataBool(l, u, 1)
         self.thisOutSsBnd[d] = SliceDataBool(l, u, 1)
      end
   end

   -- construct various functions from template representations
   self._calcDelta = loadstring( calcDeltaTempl {MEQN = meqn} )()
   self._calcCfla = loadstring( calcCflaTempl {MWAVE = mwave} )()
   self._calcFirstOrderGud = loadstring( calcFirstOrderGudTempl {MEQN = meqn} )()
   self._waveDotProd = loadstring( waveDotProdTempl {MEQN = meqn} )()
   self._rescaleWave = loadstring( rescaleWaveTempl {MEQN = meqn} )()
   self._secondOrderFlux = loadstring( secondOrderFluxTempl {MEQN = meqn} )()
   self._secondOrderUpdate = loadstring( secondOrderUpdateTempl {MEQN = meqn} )()

   self._limiter_id = limiter_id[tbl.limiter]
end

function WavePropagation:initDevice(tbl)
   self.dtByCell = DataStruct.Field {
      onGrid = self._onGrid,
      numComponents = 1,
      ghost = {1, 1},
      createDeviceCopy = true,
   }
   self.dtPtr = cuAlloc.Double(1)  -- FIXME pass in address of self.dt?
   local hyper = ffi.new("GkylWavePropagation_t")
   hyper.updateDirs = ffi.new("int[6]", self._updateDirs)
   hyper.numUpdateDirs = #self._updateDirs
   hyper.equation = self._equation._onDevice
   hyper._cfl = self._cfl
   hyper._cflm = self._cflm
   hyper.dtByCell = self.dtByCell._onDevice
   self._onHost = hyper
   local sz = sizeof("GkylWavePropagation_t")
   self._onDevice, err = cuda.Malloc(sz)
   cuda.Memcpy(self._onDevice, hyper, sz, cuda.MemcpyHostToDevice)

   self.numThreads = tbl.numThreads or GKYL_DEFAULT_NUM_THREADS
   self._useSharedDevice = xsys.pickBool(tbl.useSharedDevice, false)
end

-- Limit waves: this code closely follows the example of CLAWPACK and
-- my (AHH) thesis code Miniwarpx.
function WavePropagation:limitWaves(lower, upper, wavesSlice, speedsSlice, onSsBnd)
   local meqn, mwave = self._equation:numEquations(), self._equation:numWaves()
   for mw = 1, mwave do
      local dotr = self._waveDotProd(meqn, wavesSlice[lower-1], wavesSlice[lower], mw)
      for i = lower, upper do
         local dotl = dotr
         local wnorm2 = self._waveDotProd(meqn, wavesSlice[i], wavesSlice[i], mw)
         dotr = self._waveDotProd(meqn, wavesSlice[i], wavesSlice[i+1], mw)
         if wnorm2 > 0 then
            local r = speedsSlice[i][mw-1] > 0 and dotl/wnorm2 or dotr/wnorm2
            local wlimitr
            if (self._hasSsBnd and onSsBnd[i][0]) then
               wlimitr = self._limiterFuncZero(r)
            else
               wlimitr = self._limiterFunc(r)
            end
            self._rescaleWave(wlimitr, wavesSlice[i]+(mw-1)*meqn)
         end
      end
   end
end

local isOutside = function (inOutPtr)
   return inOutPtr[1] < 0.
end

-- advance method
function WavePropagation:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local dt = self._dt
   local qIn = assert(inFld[1], "WavePropagation.advance: Must-specify an input field")
   local qOut = assert(outFld[1], "WavePropagation.advance: Must-specify an output field")

   if (self._equation._zero_wv and true) then
      -- TODO move this to init FIXME which grid/range to use?
      if self._isFirst then
         local winp = ffi.new("struct gkyl_wave_prop_inp")
         -- winp.grid = self._onGrid._zero
         winp.grid = qOut._grid._zero
         winp.equation = self._equation._zero_wv
         winp.limiter = self._limiter_id
         winp.num_up_dirs = #self._updateDirs
         for d = 1, #self._updateDirs do
            winp.update_dirs[d-1] = self._updateDirs[d]-1
         end
         winp.cfl = self._cfl
         local mapc2p, ctx = nil, nil
         winp.geom = ffi.C.gkyl_wave_geom_new(
            winp.grid, qOut._localExtRange, mapc2p, ctx);

         self._zero = ffi.gc(
            ffi.C.gkyl_wave_prop_new(winp), ffi.C.gkyl_wave_prop_release)

         ffi.C.gkyl_wave_geom_release(winp.geom)

         self._isFirst = false
      end

      local g0_status = ffi.C.gkyl_wave_prop_advance(
         self._zero, 0, dt, qOut._localRange, qIn._zero, qOut._zero)
      local success = g0_status.success == 1 and true or false
      return self:reduceStatusDt(success, g0_status.dt_suggested)
   end

   local equation = self._equation -- equation to solve
   local meqn, mwave = equation:numEquations(), equation:numWaves()
   local localRange = self._localRange

   local qInIdxr, qOutIdxr = qIn:genIndexer(), qOut:genIndexer() -- indexer functions into fields

   local cfl, cflm = self._cfl, self._cflm
   local cfla = 0.0 -- actual CFL number used

   local delta = Lin.Vec(meqn)
   local waves = Lin.Mat(mwave, meqn)
   local s = Lin.Vec(mwave)
   local amdq, apdq = Lin.Vec(meqn), Lin.Vec(meqn)

   local idxp, idxm = Lin.IntVec(grid:ndim()), Lin.IntVec(grid:ndim())

   -- pointers to (re)use inside inner loop
   local qInL, qInR = qIn:get(1), qIn:get(1)
   local qOutL, qOutR = qOut:get(1), qOut:get(1)
   local q1 = qOut:get(1)

   local hasSsBnd = self._hasSsBnd
   local inOut
   local inOutL, inOutR
   local inOutIdxr
   if (hasSsBnd) then
      inOut = self._inOut
      inOutL, inOutR = inOut:get(1), inOut:get(1)
      inOutIdxr = self._inOut:genIndexer()
   end

   local calcDelta = self._calcDelta
   local calcFirstOrderGud = self._calcFirstOrderGud
   local calcCfla = self._calcCfla
   local secondOrderFlux = self._secondOrderFlux
   local secondOrderUpdate = self._secondOrderUpdate

   qOut:copy(qIn) -- update only adds increments, so set qOut = qIn
   -- update specified directions
   for d = 1, #self._updateDirs do
      local dir = self._updateDirs[d]
      local dtdx = dt/grid:dx(dir)

      local wavesSlice, speedsSlice, fsSlice = self.wavesSlice[dir], self.speedsSlice[dir], self.fsSlice[dir]
      local onSsBnd = self.onSsBnd[dir]
      local bothOutSsBnd = self.bothOutSsBnd[dir]
      local thisOutSsBnd = self.thisOutSsBnd[dir]
      -- lower/upper bounds in direction 'dir': these are edge indices
      local dirLoIdx, dirUpIdx = localRange:lower(dir)-1, localRange:upper(dir)+2

      if self._isFirst then
         self._perpRange[dir] = localRange:shorten(dir) -- range orthogonal to 'dir'
      end
      local perpRange = self._perpRange[dir]

      -- outer loop is over directions orthogonal to 'dir' and inner
      -- loop is over 1D slice in `dir`.
      for idx in perpRange:rowMajorIter() do
         idx:copyInto(idxp); idx:copyInto(idxm)

         -- fill masks along this direction
         if hasSsBnd then
            for i = localRange:lower(dir)-1, localRange:upper(dir)+2 do -- this loop is over edges
               idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'
               inOut:fill(inOutIdxr(idxm), inOutL)
               inOut:fill(inOutIdxr(idxp), inOutR)
               local isOutsideL = isOutside(inOutL)
               local isOutsideR = isOutside(inOutR)
               onSsBnd[i][0] =  (isOutsideL and (not isOutsideR)) or
                  ((not isOutsideL) and isOutsideR)
               bothOutSsBnd[i][0] =  isOutsideL and isOutsideR
               thisOutSsBnd[i][0] = isOutsideR
            end
         end

         for i = dirLoIdx, dirUpIdx do -- this loop is over edges
            idxm[dir], idxp[dir]  = i-1, i -- cell left/right of edge 'i'

            if not (hasSsBnd and bothOutSsBnd[i][0]) then
               qIn:fill(qInIdxr(idxm), qInL); qIn:fill(qInIdxr(idxp), qInR)
               calcDelta(qInL, qInR, delta) -- jump across interface

               equation:rp(dir, delta, qInL, qInR, waves, s) -- compute waves and speeds
               equation:qFluctuations(dir, qInL, qInR, waves, s, amdq, apdq) -- compute fluctuations

               qOut:fill(qOutIdxr(idxm), qOutL); qOut:fill(qOutIdxr(idxp), qOutR)
               calcFirstOrderGud(dtdx, qOutL, qOutR, amdq, apdq) -- first-order Gudonov updates
               cfla = calcCfla(cfla, dtdx, s) -- actual CFL value

               -- copy waves data for use in limiters
               copy(wavesSlice[i], waves:data(), sizeof("double")*meqn*mwave)
               copy(speedsSlice[i], s:data(), sizeof("double")*mwave)
            end
         end

         -- return if time-step was too large
         if cfla > cflm then return self:reduceStatusDt(false, dt*cfl/cfla) end

         -- limit waves before computing second-order updates
         self:limitWaves(localRange:lower(dir), localRange:upper(dir)+1, wavesSlice, speedsSlice, onSsBnd)

         local dirLoIdx2, dirUpIdx2 = localRange:lower(dir), localRange:upper(dir)+1 -- one more edge than cells
         -- compute second order correction fluxes
         clearSliceData(fsSlice)
         for i = dirLoIdx2, dirUpIdx2 do -- this loop is over edges
            if not (hasSsBnd and bothOutSsBnd[i][0]) then
               for mw = 0, mwave-1 do
                  secondOrderFlux(dtdx, speedsSlice[i][mw], wavesSlice[i]+meqn*mw, fsSlice[i])
               end
            end
         end

         local dirLoIdx3, dirUpIdx3 = localRange:lower(dir), localRange:upper(dir)
         -- add them to solution
         for i = dirLoIdx3, dirUpIdx3 do -- this loop is over cells
            idxm[dir] = i -- cell index
            if not (hasSsBnd and thisOutSsBnd[i][0]) then
               qOut:fill(qOutIdxr(idxm), q1)
               secondOrderUpdate(dtdx, fsSlice[i], fsSlice[i+1], q1)
            end
         end
      end
   end

   self._isFirst = false -- no longer first time
   return self:reduceStatusDt(true, dt*cfl/cfla)
end

function WavePropagation:_advanceOnDevice(tCurr, inFld, outFld)
   local qIn = assert(
      inFld[1], "WavePropagation.advanceOnDevice: Must specify an input field")
   local qOut = assert(
      outFld[1],
      "WavePropagation.advanceOnDevice: Must specify an output field")

   local numCellsLocal = qOut:localRange():volume()
   local numEdgesLocal = qOut:localEdgeRange():volume()
   local numThreads = math.min(self.numThreads, numCellsLocal)
   local numBlocks  = math.ceil(numCellsLocal/numThreads)
   local meqn, mwave = self._equation:numEquations(), self._equation:numWaves()

   if self._useSharedDevice then
      -- TODO implement
      -- ffi.C.advanceOnDevice_shared(
      --    numBlocks, numThreads, qIn:numComponents(), self._onDevice,
      --    qIn._onDevice, qOut._onDevice)
   else
      ffi.C.wavePropagationAdvanceOnDevice(
         meqn, mwave, numBlocks, numThreads, self._onDevice, qIn._onDevice,
         qOut._onDevice)
   end

   -- self.dtByCell:deviceReduce('min', self.dtPtr)
end

-- FIXME no cflRate needed here
function WavePropagation:setDtAndCflRate(dt)
   WavePropagation.super.setDtAndCflRate(self, dt, nil)

   if self._onDevice then
      ffi.C.setDt(self._onDevice, dt)
   end
end

return WavePropagation
