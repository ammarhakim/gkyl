#ifndef CART_FIELD_INTEGRATED_QUANT_CALC_IMPL_H
#define CART_FIELD_INTEGRATED_QUANT_CALC_IMPL_H

extern "C" {
    void gkylCartFieldIntQuantV(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
    void gkylCartFieldIntQuantAbsV(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
    void gkylCartFieldIntQuantV2(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
    void gkylCartFieldIntQuantGradPerpV2(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
}

#endif // CART_FIELD_INTEGRATED_QUANT_CALC_IMPL_H
