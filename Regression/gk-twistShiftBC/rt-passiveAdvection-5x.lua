-- Gkyl ---------------------------------------------------------------------------
--
-- Passively advect a structure in the z direction which then gets
-- twist-shifted at the z-boundaries.
--
--
local Plasma = require ("App.PlasmaOnCartGrid").PassiveAdvection()

-- Velocities.
local ux = 0.
local uy = 0.
local uz = 1.

-- Domain size.
local Lx = 1.
local Ly = 1.
local Lz = 1.

local function gaussian(t, xn)
   local r2 = (xn[1]-0.5)^2 + (xn[2]-0.5)^2 + (xn[3]-0.5)^2
   return math.exp(-50*r2)
end
local function sphere(t, xn)
   local r2 = (xn[1]-0.5)^2 + (xn[2]-0.5)^2 + (xn[3]-0.5)^2
   if r2 < 0.25^2 then
      return 1.0
   end
   return 1.0e-5
end
local function cube(t, xn)
   local rx2, ry2, rz2 = (xn[1]-Lx/2)^2, (xn[2]-Ly/2)^2, (xn[3]-Lz/2)^2
   if rx2 < (Lx/4)^2 and ry2 < (Ly/4)^2 and rz2 < (Lz/4)^2 then
      return 1.0
   end
   return 1.0e-10
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd        = 3,                -- End time.
   nFrame      = 1,                -- Number of output frames.
   lower       = { 0,  0,  0, -Lx/2., 0.},        -- Configuration space lower left.
   upper       = {Lx, Ly, Lz,  Lx/2., Ly},     -- Configuration space upper right.
   cells       = {12, 12, 12, 2, 2},     -- Configuration space cells.
   basis       = "serendipity",    -- One of "serendipity" or "maximal-order".
   polyOrder   = 1,                -- Polynomial order.
   timeStepper = "rk3",            -- One of "rk2" or "rk3".
   --maximumDt = 0.01,
   
   -- Decomposition for configuration space.
   decompCuts = {1, 1, 1, 1, 1},    -- Cuts in each configuration direction.
   useShared  = false,    -- If to use shared memory.

   -- Boundary conditions for configuration space.
   periodicDirs = {1, 2}, -- Periodic directions

   -- Fluid species.
   fluidUpLin = Plasma.Species {
      -- Initial conditions.
      init = function(t, xn)
        local n = cube(t,xn)
--        local vx, vy, vz = ux, uy, uz
        return n, 0., 0., xn[4], 0., 0.
      end,
        
      evolve = true, -- Evolve species?
      diagnostics = {"intMom"},
      bcz = {Plasma.TwistShiftBC{shiftFunction=function(t,xn) return 1*(xn[1]-.5) end},
             Plasma.TwistShiftBC{shiftFunction=function(t,xn) return 1*(xn[1]-.5) end}},
   },

--   fluidDownLin = Plasma.Species {
--      -- Initial conditions.
--      init = function(t, xn)
--        local n = cube(t,xn)
--        local vx, vy, vz = ux, uy, -uz
--        return n, vx, vy, vz
--      end,
--        
--      evolve = true, -- Evolve species?
--      diagnostics = {"intMom"},
--      bcz = {Plasma.TwistShiftBC{shiftFunction=function(t,xn) return 1*(xn[1]-.5) end},
--             Plasma.TwistShiftBC{shiftFunction=function(t,xn) return 1*(xn[1]-.5) end}},
--   },

--   fluidUpQuad = Plasma.Species {
--      -- Initial conditions.
--      init = function(t, xn)
--        local n = cube(t,xn)
--        local vx, vy, vz = ux, uy, uz
--        return n, vx, vy, vz
--      end,
--        
--      evolve = true, -- Evolve species?
--      diagnostics = {"intMom"},
--      bcz = {Plasma.TwistShiftBC{shiftFunction=function(t,xn) return 1*(xn[1]-.5)^2 end},
--             Plasma.TwistShiftBC{shiftFunction=function(t,xn) return 1*(xn[1]-.5)^2 end}},
--   },
}
-- Run application.
plasmaApp:run()
