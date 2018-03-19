#ifndef VLASOV_INCR_DECL_H
#define VLASOV_INCR_DECL_H

extern "C" {
    // aOut = a*aIn
    void vlasovIncr(unsigned n, const double *aIn, double a, double *aOut);
}

#endif // VLASOV_INCR_DECL_H
