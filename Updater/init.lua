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
local CartFieldBinOp = require "Updater.CartFieldBinOp"
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"
local ConfToPhase = require "Updater.ConfToPhase"
local DistFuncIntegratedMomentCalc = require "Updater.DistFuncIntegratedMomentCalc"
local DistFuncMomentCalc = require "Updater.DistFuncMomentCalc"
local FemGyroaverage = require "Updater.FemGyroaverage"
local FemParPoisson = require "Updater.FemParPoisson"
local FemPerpPoisson = require "Updater.FemPerpPoisson"
local FemPoisson = require "Updater.FemPoisson"
local FiveMomentSrc = require "Updater.FiveMomentSrc"
local TenMomentSrc = require "Updater.TenMomentSrc"
local TenMomentRelax = require "Updater.TenMomentRelax"
local HyperDisCont = require "Updater.HyperDisCont"
local LagrangeFix = require "Updater.LagrangeFix"
local MappedPoisson = require "Updater.MappedPoisson"
local MaxwellianOnBasis = require "Updater.MaxwellianOnBasis"
local PositivityRescale = require "Updater.PositivityRescale"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local SelfPrimMoments = require "Updater.SelfPrimMoments"
local SolidSurface = require "Updater.SolidSurface"
local VoronovIonization = require "Updater.VoronovIonization"
local WavePropagation = require "Updater.WavePropagation"

return {
   Bc = Bc,
   BgkCollisions = BgkCollisions,
   CalcDiagnostic = CalcDiagnostic,
   CartFieldBinOp = CartFieldBinOp,
   CartFieldIntegratedQuantCalc = CartFieldIntegratedQuantCalc,
   ConfToPhase = ConfToPhase,
   DistFuncIntegratedMomentCalc = DistFuncIntegratedMomentCalc,
   DistFuncMomentCalc = DistFuncMomentCalc,
   FemGyroaverage = FemGyroaverage,
   FemParPoisson = FemParPoisson,
   FemPerpPoisson = FemPerpPoisson,
   FemPoisson = FemPoisson,
   FiveMomentSrc = FiveMomentSrc,
   TenMomentSrc = TenMomentSrc,
   TenMomentRelax = TenMomentRelax,
   HyperDisCont = HyperDisCont,
   LagrangeFix = LagrangeFix,
   MappedPoisson = MappedPoisson,
   MaxwellianOnBasis = MaxwellianOnBasis,
   PositivityRescale = PositivityRescale,
   ProjectOnBasis = ProjectOnBasis,
   SelfPrimMoments = SelfPrimMoments,
   SolidSurface = SolidSurface,
   VoronovIonization = VoronovIonization,
   WavePropagation = WavePropagation,
}
