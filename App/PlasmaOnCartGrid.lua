-- Gkyl ------------------------------------------------------------------------
--
-- Common entry point into various Gkeyll Apps.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------


return {
   -- Return function to create Moments app objects
   Moments = function()
      local moments = require("App.Moments")
      return moments
   end
}
