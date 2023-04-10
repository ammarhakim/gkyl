-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the updater modules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl modules
local ASheathPotential = require "Updater.ASheathPotential"
local AnisotropicDiffusion = require "Updater.AnisotropicDiffusion"
local AxisymmetricFiveMomentSrc = require "Updater.AxisymmetricFiveMomentSrc"
local AxisymmetricPhMaxwellSrc = require "Updater.AxisymmetricPhMaxwellSrc"
local Bc = require "Updater.Bc"
local BgkCollisions = require "Updater.BgkCollisions"
local BiMaxwellianOnBasis = require "Updater.BiMaxwellianOnBasis"
local BraginskiiHeatConduction = require "Updater.BraginskiiHeatConduction"
local BraginskiiViscosityDiffusion = require "Updater.BraginskiiViscosityDiffusion"
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
local FiveMomentFrictionSrc = require "Updater.FiveMomentFrictionSrc"
local FiveMomentSrc = require "Updater.FiveMomentSrc"
local HyperDisCont = require "Updater.HyperDisCont"
local HyperDisContCellBased = require "Updater.HyperDisContCellBased"
local IntegratedDGMoment = require "Updater.IntegratedDGMoment"
local Ionization = require "Updater.Ionization"
local Ionization = require "Updater.Ionization"
local IterMaxwellianFix = require "Updater.IterMaxwellianFix"
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
local SqrtOnBasis = require "Updater.SqrtOnBasis"
local StairSteppedBc = require "Updater.StairSteppedBc"
local TenMomentRelax = require "Updater.TenMomentRelax"
local TenMomentSrc = require "Updater.TenMomentSrc"
local TwistShiftBC = require "Updater.TwistShiftBC"
local WavePropagation = require "Updater.WavePropagation"

return {
   ASheathPotential = ASheathPotential,
   AnisotropicDiffusion = AnisotropicDiffusion,
   AxisymmetricFiveMomentSrc = AxisymmetricFiveMomentSrc,
   AxisymmetricPhMaxwellSrc = AxisymmetricPhMaxwellSrc,
   Bc = Bc,
   BgkCollisions = BgkCollisions,
   BiMaxwellianOnBasis = BiMaxwellianOnBasis,
   BraginskiiHeatConduction = BraginskiiHeatConduction,
   BraginskiiViscosityDiffusion = BraginskiiViscosityDiffusion,
   CartFieldBinOp = CartFieldBinOp,
   CartFieldIntegratedQuantCalc = CartFieldIntegratedQuantCalc,
   CartFieldInterpolate = CartFieldInterpolate,
   CellAveMaxwellian = CellAveMaxwellian,
   ChargeExchange = ChargeExchange,
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
   FiveMomentFrictionSrc = FiveMomentFrictionSrc,
   FiveMomentSrc = FiveMomentSrc,
   HyperDisCont = HyperDisCont,
   HyperDisContCellBased = HyperDisContCellBased,
   IntegratedDGMoment = IntegratedDGMoment,
   Ionization = Ionization,
   Ionization = Ionization,
   IterMaxwellianFix = IterMaxwellianFix,
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
   SqrtOnBasis = SqrtOnBasis,
   StairSteppedBc = StairSteppedBc,
   TenMomentRelax = TenMomentRelax,
   TenMomentSrc = TenMomentSrc,
   TwistShiftBC = TwistShiftBC,
   WavePropagation = WavePropagation,
}
