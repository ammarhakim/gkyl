-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the updater modules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl modules
local Bc = require "Updater.Bc"
local BgkCollisions = require "Updater.BgkCollisions"
local CartFieldBinOp = require "Updater.CartFieldBinOp"
local CartFieldIntegratedQuantCalc = require "Updater.CartFieldIntegratedQuantCalc"
local CartFieldInterpolate = require "Updater.CartFieldInterpolate"
local CellAveMaxwellian = require "Updater.CellAveMaxwellian"
local ChargeExchange = require "Updater.ChargeExchange"
local ConfToPhase = require "Updater.ConfToPhase"
local CrossPrimMoments = require "Updater.CrossPrimMoments"
local DiscontGenPoisson = require "Updater.DiscontGenPoisson"
local DiscontPoisson = require "Updater.DiscontPoisson"
local DistFuncMomentCalc = require "Updater.DistFuncMomentCalc"
local EvalOnNodes = require "Updater.EvalOnNodes"
local EvaluateBronoldFehskeBC = require "Updater.EvaluateBronoldFehskeBC"
local FemGyroaverage = require "Updater.FemGyroaverage"
local FemParPoisson = require "Updater.FemParPoisson"
local FemPerpPoisson = require "Updater.FemPerpPoisson"
local FemPoisson = require "Updater.FemPoisson"
local FiveMomentSrc = require "Updater.FiveMomentSrc"
local HyperDisCont = require "Updater.HyperDisCont"
local HyperDisContCellBased = require "Updater.HyperDisContCellBased"
local Ionization = require "Updater.Ionization"
local IntegratedDGMoment = require "Updater.IntegratedDGMoment"
local Ionization = require "Updater.Ionization"
local IterPoisson = require "Updater.IterPoisson"
local LagrangeFix = require "Updater.LagrangeFix"
local MGpoisson = require "Updater.MGpoisson"
local MappedPoisson = require "Updater.MappedPoisson"
local MaxwellianOnBasis = require "Updater.MaxwellianOnBasis"
local PositivityCheck = require "Updater.PositivityCheck"
local PositivityRescale = require "Updater.PositivityRescale"
local ProjectFluxFunc = require "Updater.ProjectFluxFunc"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local SelfPrimMoments = require "Updater.SelfPrimMoments"
local SeparateVectorComponents = require "Updater.SeparateVectorComponents"
local SolidSurface = require "Updater.SolidSurface"
local SpitzerCollisionality = require "Updater.SpitzerCollisionality"
local StairSteppedBc = require "Updater.StairSteppedBc"
local TenMomentRelax = require "Updater.TenMomentRelax"
local TenMomentSrc = require "Updater.TenMomentSrc"
local WavePropagation = require "Updater.WavePropagation"
local AxisymmetricFiveMomentSrc = require "Updater.AxisymmetricFiveMomentSrc"
local AxisymmetricPhMaxwellSrc = require "Updater.AxisymmetricPhMaxwellSrc"
local FiveMomentFrictionSrc = require "Updater.FiveMomentFrictionSrc"
local BraginskiiHeatConduction = require "Updater.BraginskiiHeatConduction"
local BraginskiiViscosityDiffusion = require "Updater.BraginskiiViscosityDiffusion"
local AnisotropicDiffusion = require "Updater.AnisotropicDiffusion"

return {
   Bc = Bc,
   BgkCollisions = BgkCollisions,
   CartFieldBinOp = CartFieldBinOp,
   CartFieldIntegratedQuantCalc = CartFieldIntegratedQuantCalc,
   CartFieldInterpolate = CartFieldInterpolate,
   ChargeExchange = ChargeExchange,
   CellAveMaxwellian = CellAveMaxwellian,
   ConfToPhase = ConfToPhase,
   CrossPrimMoments = CrossPrimMoments,
   DiscontGenPoisson = DiscontGenPoisson,
   DiscontPoisson = DiscontPoisson,
   DistFuncMomentCalc = DistFuncMomentCalc,
   EvalOnNodes = EvalOnNodes,
   EvaluateBronoldFehskeBC = EvaluateBronoldFehskeBC,
   FemGyroaverage = FemGyroaverage,
   FemParPoisson = FemParPoisson,
   FemPerpPoisson = FemPerpPoisson,
   FemPoisson = FemPoisson,
   FiveMomentSrc = FiveMomentSrc,
   HyperDisCont = HyperDisCont,
   HyperDisContCellBased = HyperDisContCellBased,
   Ionization = Ionization,
   IntegratedDGMoment = IntegratedDGMoment,
   Ionization = Ionization,
   IterPoisson = IterPoisson,
   LagrangeFix = LagrangeFix,
   MGpoisson = MGpoisson,
   MappedPoisson = MappedPoisson,
   MaxwellianOnBasis = MaxwellianOnBasis,
   PositivityCheck = PositivityCheck,
   PositivityRescale = PositivityRescale,
   ProjectFluxFunc = ProjectFluxFunc,
   ProjectOnBasis = ProjectOnBasis,
   SelfPrimMoments = SelfPrimMoments,
   SeparateVectorComponents = SeparateVectorComponents,
   SolidSurface = SolidSurface,
   SpitzerCollisionality = SpitzerCollisionality,
   StairSteppedBc = StairSteppedBc,
   TenMomentRelax = TenMomentRelax,
   TenMomentSrc = TenMomentSrc,
   WavePropagation = WavePropagation,
   AxisymmetricFiveMomentSrc = AxisymmetricFiveMomentSrc,
   AxisymmetricPhMaxwellSrc = AxisymmetricPhMaxwellSrc,
   FiveMomentFrictionSrc = FiveMomentFrictionSrc,
   BraginskiiHeatConduction = BraginskiiHeatConduction,
   BraginskiiViscosityDiffusion = BraginskiiViscosityDiffusion,
   AnisotropicDiffusion = AnisotropicDiffusion,
}
