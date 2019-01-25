
--wall-clock time, monotonic time and sleeping for Windows, Linux and OSX.
--Written by Cosmin Apreutesei. Public Domain.

local ffi = require'ffi'

local M = {}
local ffiC = ffi.C

if ffi.os == 'Windows' then

   ffi.cdef[[
	 void time_GetSystemTimeAsFileTime(uint64_t*) asm("GetSystemTimeAsFileTime");
	 int  time_QueryPerformanceCounter(int64_t*) asm("QueryPerformanceCounter");
	 int  time_QueryPerformanceFrequency(int64_t*) asm("QueryPerformanceFrequency");
	 void time_Sleep(uint32_t ms) asm("Sleep");
	   ]]

   local t = ffi.new'uint64_t[1]'
   local DELTA_EPOCH_IN_100NS = 116444736000000000ULL

   function M.time()
      ffiC.time_GetSystemTimeAsFileTime(t)
      return tonumber(t[0] - DELTA_EPOCH_IN_100NS) / 10^7
   end

   assert(ffiC.time_QueryPerformanceFrequency(t) ~= 0)
   local qpf = tonumber(t[0])

   function M.clock()
      assert(ffiC.time_QueryPerformanceCounter(t) ~= 0)
      return tonumber(t[0]) / qpf
   end

   function M.sleep(s)
      ffiC.time_Sleep(s * 1000)
   end

elseif ffi.os == 'Linux' or ffi.os == 'OSX' then

   ffi.cdef[[
	 typedef struct {
	    long s;
	    long ns;
			} time_timespec;
	 
	 int time_nanosleep(time_timespec*, time_timespec *) asm("nanosleep");
	   ]]

   local EINTR = 4
   
   local t = ffi.new'time_timespec'
   
   function M.sleep(s)
      local int, frac = math.modf(s)
      t.s = int
      t.ns = frac * 10^9
      local ret = ffiC.time_nanosleep(t, t)
      while ret == -1 and ffi.errno() == EINTR do --interrupted
	 ret = ffiC.time_nanosleep(t, t)
      end
      assert(ret == 0)
   end

   if ffi.os == 'Linux' then

      ffi.cdef[[
		int time_clock_gettime(int clock_id, time_timespec *tp) asm("clock_gettime");
		]]

      local CLOCK_REALTIME = 0
      local CLOCK_MONOTONIC = 1

      local clock_gettime = ffi.load'rt'.time_clock_gettime

      local function tos(t)
	 return tonumber(t.s) + tonumber(t.ns) / 10^9
      end

      function M.time()
	 assert(clock_gettime(CLOCK_REALTIME, t) == 0)
	 return tos(t)
      end

      function M.clock()
	 assert(clock_gettime(CLOCK_MONOTONIC, t) == 0)
	 return tos(t)
      end

   elseif ffi.os == 'OSX' then

      ffi.cdef[[
	    typedef struct {
	       long    s;
	       int32_t us;
			   } time_timeval;

	    typedef struct {
	       uint32_t numer;
	       uint32_t denom;
			   } time_mach_timebase_info_data_t;

	    int      time_gettimeofday(time_timeval*, void*) asm("gettimeofday");
	    int      time_mach_timebase_info(time_mach_timebase_info_data_t* info) asm("mach_timebase_info");
	    uint64_t time_mach_absolute_time(void) asm("mach_absolute_time");
	      ]]

      local t = ffi.new'time_timeval'

      function M.time()
	 assert(ffiC.time_gettimeofday(t, nil) == 0)
	 return tonumber(t.s) + tonumber(t.us) / 10^6
      end

      --NOTE: this appears to be pointless on Intel Macs. The timebase fraction
      --is always 1/1 and mach_absolute_time() does dynamic scaling internally.
      local timebase = ffi.new'time_mach_timebase_info_data_t'
      assert(ffiC.time_mach_timebase_info(timebase) == 0)
      local scale = tonumber(timebase.numer) / tonumber(timebase.denom) / 10^9
      function M.clock()
	 return tonumber(ffiC.time_mach_absolute_time()) * scale
      end

   end --OSX

end --Linux or OSX

return M
