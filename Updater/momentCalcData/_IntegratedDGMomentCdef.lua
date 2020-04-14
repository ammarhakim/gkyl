-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for various C functions used in calculating integrated moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"

ffi.cdef [[

  void IntDGMoment1xSer_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x4_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x4_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x5_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x5Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x4_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x5_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x5Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x6_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x6Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 


  void IntDGMoment1xSer_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x4_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x4_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x5_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x5Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x4_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x5_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x5Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x6_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x6Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 


  void IntDGMoment1xSer_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xSer_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xSer_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xSer_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x4_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xSer_x4Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x4_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x4Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x5_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xSer_x5Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x4_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x4Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x5_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x5Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x6_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xSer_x6Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vSer_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vSer_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vSer_v3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vSer_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vSer_v3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vSer_v3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 


  void IntDGMoment1xMax_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x4_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x4_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x5_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x5Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x4_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x5_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x5Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x6_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x6Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 


  void IntDGMoment1xMax_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x4_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x4_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x5_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x5Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x4_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x5_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x5Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x6_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x6Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 


  void IntDGMoment1xMax_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xMax_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xMax_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xMax_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x4_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xMax_x4Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x4_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x4Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x5_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xMax_x5Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x4_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x4Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x5_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x5Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x6_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xMax_x6Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vMax_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vMax_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vMax_v3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vMax_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vMax_v3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vMax_v3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 


  void IntDGMoment1xTensor_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x4_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x4_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x5_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x5Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_one_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_xSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_xi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x4_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x4Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x5_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x5Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x6_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x6Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_vi_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_vSq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_intM_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v1_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v1Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v2_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v2Sq_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v3_P1(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v3Sq_P1(const double *w, const double *dx, const double *fld, double *out); 


  void IntDGMoment1xTensor_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x4_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x4_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x5_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x5Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_one_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_xSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_xi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x4_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x4Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x5_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x5Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x6_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x6Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_vi_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_vSq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_intM_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v1_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v1Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v2_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v2Sq_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v3_P2(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v3Sq_P2(const double *w, const double *dx, const double *fld, double *out); 


  void IntDGMoment1xTensor_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1xTensor_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2xTensor_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3xTensor_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x4_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment4xTensor_x4Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x4_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x4Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x5_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment5xTensor_x5Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_one_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_xSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_xi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x4_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x4Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x5_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x5Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x6_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment6xTensor_x6Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x1vTensor_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x2vTensor_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment1x3vTensor_v3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x2vTensor_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment2x3vTensor_v3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_vi_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_vSq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_intM_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v1_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v1Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v2_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v2Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v3_P3(const double *w, const double *dx, const double *fld, double *out); 

  void IntDGMoment3x3vTensor_v3Sq_P3(const double *w, const double *dx, const double *fld, double *out); 

]]
