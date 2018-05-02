-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the updater modules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl modules
local Bc = require "Updater.Bc"
local BgkCollisions = require "Updater.BgkCollisions"
local CalcDiagnostic = require "Updater.CalcDiagnostic"
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"
local DistFuncIntegratedMomentCalc = require "Updater.DistFuncIntegratedMomentCalc"
local DistFuncMomentCalc = require "Updater.DistFuncMomentCalc"
local DistFuncPrimMomentCalc = require "Updater.DistFuncPrimMomentCalc"
local FiveMomentSrc = require "Updater.FiveMomentSrc"
local HyperDisCont = require "Updater.HyperDisCont"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local WavePropagation = require "Updater.WavePropagation"
local FemPoisson = require "Updater.FemPoisson"
local FemParPoisson = require "Updater.FemParPoisson"
local FemPerpPoisson = require "Updater.FemPerpPoisson"
local FemParPoisson = require "Updater.FemParPoisson"
local ConfToPhase = require "Updater.ConfToPhase"
local SolidSurface = require "Updater.SolidSurface"

return {
   Bc = Bc,
   BgkCollisions = BgkCollisions,
   CalcDiagnostic = CalcDiagnostic,
   CartFieldIntegratedQuantCalc = CartFieldIntegratedQuantCalc,
   ConfToPhase = ConfToPhase,
   DistFuncIntegratedMomentCalc = DistFuncIntegratedMomentCalc,
   DistFuncMomentCalc = DistFuncMomentCalc,
   DistFuncPrimMomentCalc = DistFuncPrimMomentCalc,
   FemPoisson = FemPoisson,
   FemParPoisson = FemParPoisson,
   FemPerpPoisson = FemPerpPoisson,
   FiveMomentSrc = FiveMomentSrc,
   HyperDisCont = HyperDisCont,
   ProjectOnBasis = ProjectOnBasis,
   SolidSurface = SolidSurface,
   WavePropagation = WavePropagation,
}
