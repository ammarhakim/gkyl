-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute diagnostics from a field: diagnostics can be
-- computed using functions that return a single scalar.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Base = require "Updater.Base"
local Mpi = require "Comm.Mpi"
local ffi = require "ffi"

-- operators 
local function op_sum(a, b) return a+b end
local function op_max(a, b) return math.max(a, b) end
local function op_min(a, b) return math.min(a, b) end

-- initial value based on operator
local function initVal(op)
   if op == Mpi.MAX then
      return -GKYL_MAX_DOUBLE
   elseif op == Mpi.MIN then
      return GKYL_MAX_DOUBLE
   end
   return 0.0
end

-- Time-step restriction updater object
local CalcDiagnostic = {}

-- Constructor: read data from tbl
function CalcDiagnostic:new(tbl)
   local self = setmetatable({}, CalcDiagnostic)
   Base.setup(self, tbl) -- setup base object

   -- read operator and operator function
   self._op, self._opFunc = Mpi.SUM, op_sum
   if tbl.operator then
      if tbl.operator == "MAX" then
	 self._op, self._opFunc = Mpi.MAX, op_max
      elseif tbl.op_max == "MIN" then
	 self._op, self._opFunc = Mpi.MIN, op_min
      end
   end

   -- diagnostic function
   self._diagnostic = assert(tbl.diagnostic, "Updater.CalcDiagnostic: Must specify diagnostic using 'diagnostic'")

   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(CalcDiagnostic, { __call = function (self, o) return self.new(self, o) end })

-- advance method computes diagniostic across MPI ranks
local function advance(self, tCurr, dt, inFld, outFld)
   local qIn = assert(inFld[1], "CalcDiagnostic.advance: Must specify an output field")
   local dynOut = assert(outFld[1], "CalcDiagnostic.advance: Must specify an output field")

   local v = initVal(self._op) -- initial value
   
   local indexer = qIn:genIndexer()
   -- loop, computing diagnostics in each cell
   for idx in qIn:localRangeIter() do
      local qItr = qIn:get(indexer(idx))
      v = self._opFunc(v, self._diagnostic(tCurr+dt, qItr))
   end

   if Mpi.Is_comm_valid(self._nodeComm) then
      -- all-reduce across nodes
      local valLocal, val = ffi.new("double[1]"), ffi.new("double[1]")
      valLocal[0] = v
      Mpi.Allreduce(valLocal, val, 1, Mpi.DOUBLE, self._op, self._nodeComm)
      dynOut:appendData(tCurr+dt, { val[0] })
   end

   return true, GKYL_MAX_DOUBLE
end

-- Methods in updater
CalcDiagnostic.__index = { advance = Base.advanceFuncWrap(advance) }

return CalcDiagnostic
