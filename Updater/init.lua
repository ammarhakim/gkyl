-- Gkyl ------------------------------------------------------------------------
--
-- Dispatch into the updater modules
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl modules
local Bc = require "Updater.Bc"
local CalcDiagnostic = require "Updater.CalcDiagnostic"
local DistFuncMomentCalc = require "Updater.DistFuncMomentCalc"
local FiveMomentSrc = require "Updater.FiveMomentSrc"
local LinearHyperbolicDisCont = require "Updater.LinearHyperbolicDisCont"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local WavePropagation = require "Updater.WavePropagation"

-- system modules
local xsys = require "xsys"

return {
   Bc = Bc,
   CalcDiagnostic = CalcDiagnostic,
   DistFuncMomentCalc = DistFuncMomentCalc,
   FiveMomentSrc = FiveMomentSrc,
   LinearHyperbolicDisCont = LinearHyperbolicDisCont,
   ProjectOnBasis = ProjectOnBasis,
   WavePropagation = WavePropagation,
}
