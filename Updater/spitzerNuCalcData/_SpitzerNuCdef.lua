-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in calculating Spitzer collisionality.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

double gkyl_ipow(double base, int exp);

void SpitzerNuCellAvScale1xSer_P1(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvScale2xSer_P1(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvScale3xSer_P1(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvScale1xSer_P2(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvScale2xSer_P2(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvScale3xSer_P2(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvScale1xSer_P3(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvScale2xSer_P3(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvBuild1xSer_P1(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvBuild2xSer_P1(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvBuild3xSer_P1(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvBuild1xSer_P2(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvBuild2xSer_P2(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvBuild3xSer_P2(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvBuild1xSer_P3(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvBuild2xSer_P3(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvScale1xMax_P1(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvScale2xMax_P1(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvScale3xMax_P1(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvScale1xMax_P2(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvScale2xMax_P2(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvScale3xMax_P2(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvScale1xMax_P3(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvScale2xMax_P3(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvBuild1xMax_P1(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvBuild2xMax_P1(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvBuild3xMax_P1(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvBuild1xMax_P2(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvBuild2xMax_P2(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvBuild3xMax_P2(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 


void SpitzerNuCellAvBuild1xMax_P3(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu); 

void SpitzerNuCellAvBuild2xMax_P3(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu);



]] 
