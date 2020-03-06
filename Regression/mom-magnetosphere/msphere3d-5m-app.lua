local Constants = require "Constants"
local buildMsphereApp = require("App.Special.Msphere")

local lightSpeed = Constants.SPEED_OF_LIGHT / 50
local massRatio = 25
local pressureRatio = 5
local ionInertialLengthSW = 400e3
local rhoSW = 5e6 * Constants.PROTON_MASS

local mu0 = Constants.MU0
local epsilon0 = 1/mu0/(lightSpeed^2)

local ionMass = Constants.PROTON_MASS
local ionCharge = ionMass / ionInertialLengthSW /
                  math.sqrt(mu0 * rhoSW)

local app = buildMsphereApp {
   -- Earth parameters
   planetRadius = 6400e3,
   planetBx0 = 0, planetBy0 = 0, planetBz0 = -3.12e-5,

   -- solar wind parameters
   rhoIn = rhoSW, pIn = 3.45e-12,
   vxIn = 400e3, vyIn = 0, vzIn = 0,
   BxIn = 0, ByIn = 0, BzIn = -5e-9,

   -- multifluid parameters
   mass = {ionMass, ionMass / massRatio},
   charge = {ionCharge, -ionCharge},
   moments = {5, 5},
   pressureFractions = {1 / (1 + 1 / pressureRatio),
                        1 / (1 + pressureRatio)},

   mu0 = mu0, epsilon0 = epsilon0,

  -- domain and grid parameters
   xlo = -5, xup = 5, nx = 80,
   ylo = -5, yup = 5, ny = 80,
   zlo = -5, zup = 5, nz = 80,

   tEnd = 2700,
   tFrame = 30,
}

app:run()
