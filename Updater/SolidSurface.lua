-- Gkyl ------------------------------------------------------------------------
--
-- Functions for plasma-surface interactions
--
--------------------------------------------------------------------------------

local _M = {}

-- 1 - probability for quantum-mechanical reflection by the surface potential
-- Taken from Bronold, 2018, Electron kinetics at the plasma interface
function _M.StickingProb(E, nu, elcAffinity, condBandElcMass)
   -- direction cosine inside due to the conservation of energy and momentum
   local eta = math.sqrt(1 - (E-elcAffinity)/(condBandElcMass*E) * (1-nu*nu))
   if E > elcAffinity then
      local k = math.sqrt(E-elcAffinity)*nu
      local p = math.sqrt(condBandElcMass*E)*eta

      return (4*condBandElcMass*k*p) / (condBandElcMass*k + p)^2
   else
      return 0.0
   end
end

return _M
