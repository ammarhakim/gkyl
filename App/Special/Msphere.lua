-- Gkyl ------------------------------------------------------------------------

local Proto = require("Proto")
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"
local TenMoment = require "Eq.TenMoment"
local Constants = require "Constants"
local Logger = require "Logger"
local common = require "App.Special.Msphere_common"

local logger = Logger {
   logToFile = true
}

local log = function(...)
   logger(string.format(...))
   logger("\n")
end

local function setdefault(tbl, key, val)
   if tbl[key] == nil then
      tbl[key] = val
   end
   return tbl[key]
end

local function buildApp(tbl)
   local moments = tbl.moments
   local nSpecies = #moments

   local mass = tbl.mass
   local totalMass = 0
   tbl.massFractions = {}
   for s = 1, nSpecies do
      totalMass = totalMass + mass[s]
   end
   for s = 1, nSpecies do
      tbl.massFractions[s] = mass[s] / totalMass
   end

   local R = tbl.planetRadius
   local xlo, xup, nx = tbl.xlo * R, tbl.xup * R, tbl.nx
   local ylo, yup, ny = tbl.ylo * R, tbl.yup * R, tbl.ny
   local zlo, zup, nz = tbl.zlo * R, tbl.zup * R, tbl.nz
   local lower = {xlo, ylo, zlo}
   local upper = {xup, yup, zup}
   local cells = {nx, ny, nz}

   local coordinateMap = nil

   local useNonUniformGrid = setdefault(tbl, "useNonUniformGrid", false)
   if useNonUniformGrid then
      local wx, fwx, rx = tbl.wx, tbl.fwx, tbl.rx
      local wy, fwy, ry = tbl.wy, tbl.fwy, tbl.ry
      local wz, fwz, rz = tbl.wz, tbl.fwz, tbl.rz

      coordinateMap = {
         common.buildGridFlatTop1d (xlo, xup, nx, 0, wx, fwx, rx),
         common.buildGridFlatTop1d (ylo, yup, ny, 0, wy, fwy, ry),
         common.buildGridFlatTop1d (zlo, zup, nz, 0, wz, fwz, rz),
      }

      lower = {0, 0, 0}
      upper = {1, 1, 1}
   end

   local tEnd = tbl.tEnd
   local tFrame = tbl.tFrame
   local nFrame = math.floor(tEnd / tFrame)

   local valuesIn = common.calcValuesIn(tbl)

   local initMirdip = common.buildInitMirdip(tbl)

   local species = {}
   local speciesNames = {}
   for s = 1, nSpecies do
      local moment = tbl.moments[s]

      local model
      if moment == 5 then
         model = Euler
      elseif moment == 10 then
         model = TenMoment
      end

      species[s] = Moments.Species {
         charge = tbl.charge[s],
         mass = tbl.mass[s],
         equation = model {gasGamma = tbl.gasGamma},
         equationInv = model {gasGamma = tbl.gasGamma, numericalFlux = "lax"},
         forceInv = false,
         evolve = true,

         init = function(t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            return initMirdip.initFluid(x, y, z, s)
         end,

         bcx = {
            model.bcConst(unpack(valuesIn[s])),
            model.bcCopy
         },
         bcy = {model.bcCopy, model.bcCopy},
         bcz = {model.bcCopy, model.bcCopy},

         hasSsBnd = true,
         inOutFunc = common.buildInOutFunc(tbl),
         ssBc = {Moments.Species.bcCopy}
      }
      speciesNames[s] = "fluid"..s
   end

   local momentApp = Moments.App {
      logToFile = true,

      tEnd = tEnd,
      nFrame = nFrame,

      lower = lower,
      upper = upper,
      cells = cells,
      coordinateMap = coordinateMap,

      timeStepper = "fvDimSplit",

      -- TODO any number of species
      ion = species[1],
      elc = species[2],

      field = Moments.Field {
         epsilon0 = tbl.epsilon0,
         mu0 = tbl.mu0,
         evolve = true,
         init = function(t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            return initMirdip.initField(x, y, z, nSpecies + 1)
         end,
         bcx = {
            Moments.Field.bcConst(unpack(valuesIn[nSpecies + 1])),
            Moments.Field.bcCopy
         },
         bcy = {Moments.Field.bcCopy, Moments.Field.bcCopy},
         bcz = {Moments.Field.bcCopy, Moments.Field.bcCopy},
         hasSsBnd = true,
         inOutFunc = common.buildInOutFunc(tbl),
         ssBc = {Moments.Species.bcReflect}
      },

      emSource = Moments.CollisionlessEmSource {
         species = {"ion", "elc"},
         timeStepper = "direct",
         hasStaticField = true,
         staticEmFunction = function(t, xn)
            local x, y, z = xn[1], xn[2], xn[3]
            return initMirdip.calcStaticEB(x, y, z)
         end
      }
   }
   return momentApp
end

return buildApp

