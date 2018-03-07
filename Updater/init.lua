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
local FiveMomentSrc = require "Updater.FiveMomentSrc"
local HyperDisCont = require "Updater.HyperDisCont"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local VlasovDisCont = require "Updater.VlasovDisCont"
local WavePropagation = require "Updater.WavePropagation"
local FemPerpPoisson = require "Updater.FemPerpPoisson"

return {
   Bc = Bc,
   BgkCollisions = BgkCollisions,
   CalcDiagnostic = CalcDiagnostic,
   CartFieldIntegratedQuantCalc = CartFieldIntegratedQuantCalc,
   DistFuncIntegratedMomentCalc = DistFuncIntegratedMomentCalc,
   DistFuncMomentCalc = DistFuncMomentCalc,
   FiveMomentSrc = FiveMomentSrc,
   HyperDisCont = HyperDisCont,
   ProjectOnBasis = ProjectOnBasis,
   VlasovDisCont = VlasovDisCont,   
   WavePropagation = WavePropagation,
   FemPerpPoisson = FemPerpPoisson,
}
