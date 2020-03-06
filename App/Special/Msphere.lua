-- Gkyl ------------------------------------------------------------------------
local Moments = require("App.PlasmaOnCartGrid").Moments
local Euler = require "Eq.Euler"
local TenMoment = require "Eq.TenMoment"
local Logger = require "Logger"
local common = require "App.Special.Msphere_common"

local function buildApp(tbl)
   common.setdefaultObject(tbl)

   local R = tbl.planetRadius
   local xlo, xup, nx = tbl.xlo * R, tbl.xup * R, tbl.nx
   local ylo, yup, ny = tbl.ylo * R, tbl.yup * R, tbl.ny
   local zlo, zup, nz = tbl.zlo * R, tbl.zup * R, tbl.nz
   local lower = {xlo, ylo, zlo}
   local upper = {xup, yup, zup}
   local cells = {nx, ny, nz}

   local coordinateMap = nil
   local useNonUniformGrid = tbl.useNonUniformGrid
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

   local initMirdip = common.buildInitMirdip(tbl)
   local valuesIn = common.calcValuesIn(tbl)

   local nSpecies = #tbl.moments
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

         bcx = {model.bcConst(unpack(valuesIn[s])), model.bcCopy},
         bcy = {model.bcCopy, model.bcCopy},
         bcz = {model.bcCopy, model.bcCopy},

         hasSsBnd = true,
         inOutFunc = common.buildInOutFunc(tbl),
         ssBc = {Moments.Species.bcCopy}
      }
      speciesNames[s] = "fluid"..s
   end

   local field = Moments.Field {
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
   }

   local emSource = Moments.CollisionlessEmSource {
      species = species,
      timeStepper = "direct",
      hasStaticField = true,
      staticEmFunction = function(t, xn)
         local x, y, z = xn[1], xn[2], xn[3]
         return initMirdip.calcStaticEB(x, y, z)
      end
   }

   local mt = {
      logToFile = true,

      lower = lower,
      upper = upper,
      cells = cells,
      coordinateMap = coordinateMap,

      tEnd = tEnd,
      nFrame = nFrame,

      timeStepper = "fvDimSplit",

      field = field,
      emSource = emSource,
   }
   for s = 1, nSpecies do
      mt[speciesNames[s]] = species[s]
   end

   return Moments.App (mt)
end

return buildApp
