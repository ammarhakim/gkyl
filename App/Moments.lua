-- Gkyl ------------------------------------------------------------------------
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Infrastructure loads
local DecompRegionCalc = require "Lib.CartDecomp"
local Grid = require "Grid"
local LinearTrigger = require "Lib.LinearTrigger"
local Logger = require "Lib.Logger"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Time = require "Lib.Time"
local date = require "xsys.date"
local xsys = require "xsys"

-- Euler equation initialization
local Euler = function(tbl)
   return G0.Moments.Eq.Euler.new(tbl)
end

-- IsoEuler equation initialization
local IsoEuler = function(tbl)
   return G0.Moments.Eq.IsoEuler.new(tbl)
end

-- TenMoment equation initialization
local TenMoment = function(tbl)
   return G0.Moments.Eq.TenMoment.new(tbl)
end

-- Moment species initialization
local Species = function(tbl)
   return G0.Moments.Species.new(tbl)
end

-- Moment field initialization
local Field = function(tbl)
   return G0.Moments.Field.new(tbl)
end

-- top-level moments App
local App = Proto()

-- Initialize the App object
function App:init(tbl)
   self.label = "Moments"
   self.tStart = tbl.tStart and tbl.tStart or 0.0
   self.tEnd  = tbl.tEnd
   self.nFrame = tbl.nFrame

   self.comm = Mpi.COMM_WORLD
   local commSz = Mpi.Comm_size(self.comm)

   if not tbl.cuts then
      if commSz > 1 then
	 -- build cuts table if not specified explicitly
	 local cuts = DecompRegionCalc.makeCuts(#tbl.cells, commSz, tbl.cells)
	 tbl.decompCuts = cuts
      end
   end

   -- construct G0 Moments App object
   self.g0App = G0.Moments.App.new(tbl)

   -- grid object
   self.grid = Grid.RectCart {
      lower = tbl.lower,
      upper = tbl.upper,
      cells = tbl.cells,
      periodicDirs = tbl.periodicDirs
   }

   -- ghost cell layout
   self.nghost = self.g0App:nghost()

end

-- Run simulation
function App:run()

   local maxSteps = GKYL_MAX_INT

   local log = Logger { logToFile = false }
   log(date(false):fmt()); log("\n") -- Time-stamp for sim start.

   if GKYL_GIT_CHANGESET then
      log(string.format("Gkyl built with %s\n", GKYL_GIT_CHANGESET))
   end
   if GKYL_BUILD_DATE then
      log(string.format("Gkyl built on %s\n", GKYL_BUILD_DATE))
   end

   local tStart, tEnd = self.tStart, self.tEnd
   local nFrame = self.nFrame
   local tCurr = tStart

   local timers = {
      total = 0.0,
      init = 0.0,
      timeloop = 0.0,
   }
   timers.total = Time.clock()
   timers.init  = timers.total
   log(string.format("Initializing %s simulation ...\n", self.label))

   -- Triggers
   local ioTrigger = LinearTrigger(0, tEnd, nFrame)
   local logTrigger = LinearTrigger(0, tEnd, 10)
   local logTrigger1p = LinearTrigger(0, tEnd, 100)

   local tenth = 0
   if tStart > 0 then
      tenth = math.floor(tStart/tEnd*10)
   end

   local p1c = 0
   if tStart > 0 then
      p1c = math.floor(tStart/tEnd*100) % 10
   end

   local ioFrame = 0
   -- function to write data to file
   local function writeData(tCurr)
      if ioTrigger(tCurr) then
	 self.g0App:write(tCurr, ioFrame)
	 self.g0App:write_integrated_mom()
	 self.g0App:write_field_energy()

	 ioFrame = ioFrame + 1
      end
   end

   local logCount = 0 -- This is needed to avoid initial log message.
   -- For writing out step information
   local function writeStepMessage(step, tCurr, dt_next)
      if logTrigger(tCurr) then
	 if logCount > 0 then
	    log (string.format(
		    " Step %6d at time  %#11.8g.  Time step  %.6e.  Completed %g%s\n", 
		    step, tCurr, dt_next, tenth*10, "%"))
	 else
	    logCount = logCount+1
	 end
	 tenth = tenth+1
      end
      if logTrigger1p(tCurr) then
	 log(string.format("%d", p1c))
	 p1c = (p1c+1)%10
      end
   end

   -- Apply initial conditions
   self.g0App:apply_ic(tCurr)
   self.g0App:calc_integrated_mom(tCurr)
   self.g0App:calc_field_energy(tCurr)

   writeData(tCurr)

   timers.init = Time.clock()-timers.init
   log(string.format("Initialization completed in %g sec\n\n", timers.init))

   local first = true
   local dt = tEnd - tStart
   local step = 1
   -- Main simulation loop
   timers.timeloop = Time.clock()
   while true do
      local upStatus = self.g0App:update(dt)

      if not upStatus.success then
	 log(string.format("*** Simulation failed at step %d, time %g! Aborting\n", step, tCurr))
	 break
      end

      if first then 
	 log(string.format(" Step 0 at time %g. Time step %g. Completed 0%%\n", tCurr, upStatus.dt_actual))
	 first = false
      end
      writeStepMessage(step, tCurr, upStatus.dt_suggested)

      tCurr = tCurr + upStatus.dt_actual
      dt = upStatus.dt_suggested
      
      self.g0App:calc_integrated_mom(tCurr)
      self.g0App:calc_field_energy(tCurr)

      if (tCurr >= tEnd) then break end

      writeData(tCurr)
      step = step + 1
   end
   local tmEnd = Time.clock()
   timers.timeloop = tmEnd - timers.timeloop
   timers.total = tmEnd - timers.total

   self.g0App:stat_write()

   log(string.format("\n\nTotal number of time-steps %s\n", step))
   log(string.format("Main-loop took %g secs to run\n", timers.timeloop))
   log(string.format("Simulation took a total of %g secs to run\n", timers.total))
   log(string.format("See log file for details\n"))
end

-- run using G0 run method
function App:g0run(numSteps)
   self.g0App:run(numSteps)
end

return {
   App = App,
   Species = Species,
   Field = Field,

   -- various boundary conditions for species
   SpeciesBc = {
      bcCopy = G0.SpeciesBc.bcCopy,
      bcWall = G0.SpeciesBc.bcWall,
      bcReflect = G0.SpeciesBc.bcReflect,
      bcAbsorb = G0.SpeciesBc.bcAbsorb,
      bcNoSlip = G0.SpeciesBc.bcNoSlip,
      bcWedge = G0.SpeciesBc.bcWedge,
      bcFunc = G0.SpeciesBc.bcFunc,
      bcFixedFunc = G0.SpeciesBc.bcFixedFunc,
      bcZeroFlux = G0.SpeciesBc.bcZeroFlux,
   },

   -- various boundary conditions for fields
   FieldBc = {
      bcCopy = G0.FieldBc.bcCopy,
      bcWall = G0.FieldBc.bcWall,
      bcPECWall = G0.FieldBc.bcPECWall,
      bcSymWall = G0.FieldBc.bcSymWall,
      bcReservoir = G0.FieldBc.bcReservoir,
      bcWedge = G0.FieldBc.bcWedge,
      bcFunc = G0.FieldBc.bcFunc,
   },

   -- supported equation systems
   Eq = {
      Euler = Euler,
      IsoEuler = IsoEuler,
      TenMoment = TenMoment
   }
}
