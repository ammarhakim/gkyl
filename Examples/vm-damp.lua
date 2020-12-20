-- Gkyl ------------------------------------------------------------------------
-- 1x1v simulation of collisionless damping of an electron Langmuir wave,
-- using kinetic ions and electrons. We use normalize units.
--------------------------------------------------------------------------------
local Plasma = require("App.PlasmaOnCartGrid").VlasovMaxwell()

permitt  = 1.0   -- Permittivity of free space.
permeab  = 1.0   -- Permeability of free space.
eV       = 1.0   -- Elementary charge, or Joule-eV conversion factor.
elcMass  = 1.0   -- Electron mass.
ionMass  = 1.0   -- Ion mass.

nElc = 1.0    -- Electron number density.
nIon = nElc   -- Ion number density.
Te   = 1.0    -- Electron temperature.
Ti   = Te     -- Ion temperature.

vtElc   = math.sqrt(eV*Te/elcMass)                   -- Electron thermal speed.
vtIon   = math.sqrt(eV*Ti/ionMass)                   -- Ion thermal speed.
wpe     = math.sqrt((eV^2)*nElc/(permitt*elcMass))   -- Plasma frequency.
lambdaD = vtElc/wpe                                  -- Debye length.

-- Amplitude and wavenumber of sinusoidal perturbation.
pertA = 1.0e-3
pertK = .750/lambdaD

-- Maxwellian in (x,vx)-space, given the density (denU), bulk flow
-- velocity (flowU), mass and temperature (temp).
local function maxwellian1D(x, vx, den, flowU, mass, temp)
   local v2   = (vx - flowU)^2
   local vtSq = temp/mass
   return (den/math.sqrt(2*math.pi*vtSq))*math.exp(-v2/(2*vtSq))
end

plasmaApp = Plasma.App {
   tEnd         = 20.0/wpe,           -- End time.
   nFrame       = 20,                 -- Number of output frames.
   lower        = {-math.pi/pertK},   -- Lower boundary of configuration space.
   upper        = { math.pi/pertK},   -- Upper boundary of configuration space.
   cells        = {64},               -- Configuration space cells.
   polyOrder    = 1,                  -- Polynomial order.
   periodicDirs = {1},                -- Periodic directions.
   
   elc = Plasma.Species {
      charge = -eV, mass = elcMass,
      lower = {-6.0*vtElc},      -- Velocity space lower boundary.
      upper = { 6.0*vtElc},      -- Velocity space upper boundary.
      cells = {64},              -- Number of cells in velocity space.
      init = function (t, xn)    -- Initial conditions.
	 local x, v = xn[1], xn[2]
	 return (1+pertA*math.cos(pertK*x))*maxwellian1D(x, v, nElc, 0.0, elcMass, Te) 
      end,
      evolve = true, -- Evolve species?
   },

   ion = Plasma.Species {
      charge = eV, mass = ionMass,
      lower = {-6.0*vtIon},      -- Velocity space lower boundary.
      upper = { 6.0*vtIon},      -- Velocity space upper boundary.
      cells = {64},              -- Number of cells in velocity space.
      init  = function (t, xn)   -- Initial conditions.
	 local x, v = xn[1], xn[2]
	 return maxwellian1D(x, v, nIon, 0.0, ionMass, Ti) 
      end,
      evolve = true, -- Evolve species?
   },

   field = Plasma.Field {
      epsilon0 = permitt, mu0 = permeab,
      init = function (t, xn)   -- Initial conditions.
         local Ex, Ey, Ez = -pertA*math.sin(pertK*xn[1])/pertK, 0.0, 0.0
         local Bx, By, Bz = 0.0, 0.0, 0.0
         return Ex, Ey, Ez, Bx, By, Bz
      end,
      evolve = true, -- Evolve field?
   },
}
-- Run application.
plasmaApp:run()
