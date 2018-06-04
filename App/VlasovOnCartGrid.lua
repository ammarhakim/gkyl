-- Gkyl ------------------------------------------------------------------------
--
-- LEGACY: Vlasov solver on a Cartesian grid. Works in arbitrary CDIM/VDIM
-- (VDIM>=CDIM) with either Maxwell, Poisson or specified EM fields.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- use PlasmaOnCartGrid infrastructure now
local Plasma = require "App.PlasmaOnCartGrid"

return {
   App             = Plasma.App,
   Species         = Plasma.VlasovSpecies,
   EmField         = Plasma.MaxwellField,
   FuncField       = Plasma.FuncMaxwellField,
   NoField         = Plasma.NoField,
   BgkCollisions   = Plasma.BgkCollisions,   
   VmLBOCollisions = Plasma.VmLBOCollisions,   
}
