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
local argparse = require "Lib.argparse"

-- argument parser for run command
local runParser = argparse()
   :name("run")
   :description("Run command")

runParser:option("-s --steps", "Maximum number of steps to use", GKYL_MAX_INT, tonumber)

local function runParserEval(cline, paramList)
   local opts = runParser:parse(cline)
   paramList.maxSteps = opts.steps
end

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

-- SR Euler equation initialization
local SrEuler = function(tbl)
   return G0.Moments.Eq.SrEuler.new(tbl)
end

-- Cold-fluid equation initialization
local ColdFluid = function(tbl)
   return G0.Moments.Eq.ColdFluid.new(tbl)
end

-- Ideal MHD equation initialization
local Mhd = function(tbl)
   return G0.Moments.Eq.Mhd.new(tbl)
end

-- Moment species initialization
local Species = function(tbl)
   return G0.Moments.Species.new(tbl)
end

-- Moment field initialization
local Field = function(tbl)
   return G0.Moments.Field.new(tbl)
end

-- Run simulation: "app" is a pre-conructed Moment App object
local function runSimulation(app, paramList)
   local log = app.log
   local maxSteps = paramList.maxSteps

   local tStart, tEnd = app.tStart, app.tEnd
   local nFrame = app.nFrame
   local tCurr = tStart

   local timers = {
      total = 0.0,
      init = 0.0,
      timeloop = 0.0,
   }
   timers.total = Time.clock()
   timers.init  = timers.total
   log(string.format("Initializing %s simulation ...\n", app.label))

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
	 app.g0App:write(tCurr, ioFrame)
	 app:writeFieldDiagnostics(tCurr)
	 app:writeSpeciesDiagnostics(tcurr)

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
   app.g0App:apply_ic(tCurr)
   app:calcFieldDiagnostics(tCurr)
   app:calcSpeciesDiagnostics(tCurr)

   writeData(tCurr)

   timers.init = Time.clock()-timers.init
   log(string.format("Initialization completed in %g sec\n\n", timers.init))

   local first = true
   local dt = tEnd - tStart
   local step = 0
   -- Main simulation loop
   timers.timeloop = Time.clock()
   while (tCurr < tEnd) and (step < maxSteps) do
      local upStatus = app.g0App:update(dt)

      if not upStatus.success then
	 log(string.format("*** Simulation failed at step %d, time %g! Aborting\n", step, tCurr))
	 break
      end

      if first then 
	 log(string.format(" Step 0 at time %g. Time step %g. Completed 0%%\n", tCurr, upStatus.dt_actual))
	 first = false
      end

      tCurr = tCurr + upStatus.dt_actual
      dt = upStatus.dt_suggested

      writeStepMessage(step, tCurr, upStatus.dt_suggested)
      
      app:calcFieldDiagnostics(tCurr)
      app:calcSpeciesDiagnostics(tCurr)
      writeData(tCurr)

      step = step + 1
   end
   local tmEnd = Time.clock()
   timers.timeloop = tmEnd - timers.timeloop
   timers.total = tmEnd - timers.total

   app.g0App:stat_write()

   log(string.format("\n\nTotal number of time-steps %s\n", step))
   log(string.format("Main-loop took %g secs to run\n", timers.timeloop))
   log(string.format("Simulation took a total of %g secs to run\n", timers.total))
   log(string.format("See log file for details\n"))
end

-- Does nothing
local function nullSimulation(app, paramList)
   local log = Logger { logToFile = false }
   log(string.format("*** Incorrect simulation command specified!\n"))
end

-- top-level moments App
local App = Proto()

-- Initialize the App object
function App:init(tbl)
   self.log = Logger { logToFile = false }

   local log = self.log
   log(date(false):fmt()); log("\n") -- Time-stamp for sim start.

   if GKYL_GIT_CHANGESET then
      log(string.format("Gkyl built with %s\n", GKYL_GIT_CHANGESET))
   end
   if GKYL_BUILD_DATE then
      log(string.format("Gkyl built on %s\n", GKYL_BUILD_DATE))
   end

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

   local cmdParsers = { 
      run = runParser
   }

 -- list of parsers
   local parsers = {
      run = runParserEval,
   }

   -- by default, run the whole simulation
   self.runMethod = runSimulation
   self.paramList = { maxSteps = GKYL_MAX_INT }

   local cline = { }
   -- constuct run method based on commands
   if #GKYL_COMMANDS > 0 then
      local cmd = GKYL_COMMANDS[1]
      if #GKYL_COMMANDS > 1 then
	 for i = 2, #GKYL_COMMANDS do
	    cline[i-1] = GKYL_COMMANDS[i]
	 end
      end

      -- fetch parser and run it
      local parserFunc = parsers[cmd]
      if parserFunc then
	 parserFunc(cline, self.paramList)
      else
	 self.runMethod = nullSimulation
      end
   end
   
end

-- Compute field diagnostics
function App:calcFieldDiagnostics(tm)
   self.g0App:calc_field_energy(tm)
end

-- Write field diagnostics
function App:writeFieldDiagnostics(tm)
   self.g0App:write_field_energy()
end

-- Compute field diagnostics
function App:calcSpeciesDiagnostics(tm)
   self.g0App:calc_integrated_mom(tm)
end

-- Write field diagnostics
function App:writeSpeciesDiagnostics(tm)
   self.g0App:write_integrated_mom()
end


-- Run simulation
function App:run()
   self.runMethod(self, self.paramList)
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
      TenMoment = TenMoment,
      SrEuler = SrEuler,
      ColdFluid = ColdFluid,
      Mhd = Mhd
   }
}
