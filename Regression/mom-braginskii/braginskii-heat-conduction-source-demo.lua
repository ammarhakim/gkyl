-- BraginskiiHeatConductionSource demo in 1d with external field.
-- Note that we cannot do 2d for now because the corner cells are not properly
-- synchronized by Gkeyll for now.

local Moments = require("App.PlasmaOnCartGrid").Moments()
local Euler = require "Eq.Euler"

local gasGamma  = 2
local epsilon0  = 1.
local elcMass   = 1.
local ionMass   = 100.
local ionCharge = 1.
local elcCharge = -ionCharge
local kappaPara = {0.01, 0}
local kappaPerp = {0.0, 0.01}

local lower        = {-0.5}
local upper        = {0.5}
local cells        = {64}
local periodicDirs = {1}

local cfl = 0.9
local tStart = 0
local tEnd = 9
local nFrame = 6

local function calcT (t, xn)
   local x, y = xn[1], xn[2]

   local Pi   = math.pi
   local Lx   = {upper[1] - lower[1]}
   local kNum = {2*Pi/Lx[1]}

   return 1 + 0.1 * math.sin(kNum[1]*x)
end

local function initElc (t, xn)
   local x, y = xn[1], xn[2]
   local n = 1
   return n * elcMass, 0, 0, 0, n * calcT(t, xn) / (gasGamma-1)
end

local function initIon (t, xn)
   local x, y = xn[1], xn[2]
   local n = 1
   return n * ionMass, 0, 0, 0, n * calcT(t, xn) / (gasGamma-1)
end

local function initEmf0 (t, xn)
   local x, y = xn[1], xn[2]
   return 0, 0, 0, 1, 1, 0, 0, 0
end

local function initEmf (t, xn)
   local x, y = xn[1], xn[2]
   return 0, 0, 0, 0, 0, 0, 0, 0
end

momentApp = Moments.App {
   logToFile = true,

   tEnd = tEnd,
   nFrame = nFrame,
   lower = lower,
   upper = upper,
   cells = cells,
   timeStepper = "fvDimSplit",

   periodicDirs = periodicDirs,

   elc = Moments.Species {
      charge = elcCharge, mass = elcMass,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      init = initElc,
      evolve = false,
      forceWrite = true,
      forceApplyBc = true,
      forceApplyBc = true,
   },

   ion = Moments.Species {
      charge = ionCharge, mass = ionMass,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      forceInv = false,
      init = initIon,
      evolve = false,
      forceWrite = true,
      forceApplyBc = true,
   },

   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = epsilon0,
      init = initEmf,
      evolve = false,
   },

   heatConduction = Moments.BraginskiiHeatConductionSource {
      species = {"elc", "ion"},
      gasGamma = gasGamma,
      kappaMode = "constant",
      kappaPara = kappaPara,
      kappaPerp = kappaPerp,
      cfl = cfl,
   },

   -- Specified, time-dependent electromagnetic fields.
   externalField = Moments.ExternalField {
      emFunc = initEmf0,
      evolve = false,
   },
}

-- run application
momentApp:run()
