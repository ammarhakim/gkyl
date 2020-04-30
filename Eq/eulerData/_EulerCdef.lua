-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in the Euler equations
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

typedef struct {
   double _gasGamma; /* Gas constant */
   double _rpTime; /* Time spent in RP */
   int _numWaves; /* Number of waves in system */
   int _rpType; /* Type of RP to use */
   double _fl[6], _fr[6]; /* Storage for left/right fluxes ([6] as we want to index from 1) */
} EulerEq_t;
 
]]
