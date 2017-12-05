-- Gkyl ------------------------------------------------------------------------
--
-- Lua file-system library
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- As the library is preloaded, this allows user code to require
-- "lfs". In reality, one never needs to require "lfs" but external
-- code may expect this.
return lfs
