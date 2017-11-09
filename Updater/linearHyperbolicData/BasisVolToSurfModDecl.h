#ifndef BASIS_VOL_TO_SURF_MOD_DECL_H
#define BASIS_VOL_TO_SURF_MOD_DECL_H

extern "C"
{
    void ModalMaxOrderBasisSurf1DP1_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf1DP1_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf1DP2_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf1DP2_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf1DP3_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf1DP3_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf1DP4_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf1DP4_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 

    void ModalMaxOrderBasisSurf2DP1_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP1_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP1_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP1_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP2_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP2_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP2_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP2_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP3_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP3_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP3_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP3_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP4_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP4_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP4_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf2DP4_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 

    void ModalMaxOrderBasisSurf3DP1_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP1_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP1_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP1_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP1_Upper3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP1_Lower3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP2_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP2_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP2_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP2_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP2_Upper3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP2_Lower3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP3_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP3_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP3_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP3_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP3_Upper3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP3_Lower3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP4_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP4_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP4_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP4_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP4_Upper3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalMaxOrderBasisSurf3DP4_Lower3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 


 
    void ModalSerendipBasisSurf1DP1_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf1DP1_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf1DP2_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf1DP2_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf1DP3_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf1DP3_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf1DP4_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf1DP4_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 

    void ModalSerendipBasisSurf2DP1_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP1_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP1_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP1_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP2_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP2_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP2_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP2_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP3_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP3_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP3_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP3_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP4_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP4_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP4_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf2DP4_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 

    void ModalSerendipBasisSurf3DP1_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP1_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP1_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP1_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP1_Upper3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP1_Lower3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP2_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP2_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP2_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP2_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP2_Upper3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP2_Lower3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP3_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP3_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP3_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP3_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP3_Upper3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP3_Lower3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP4_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP4_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP4_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP4_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP4_Upper3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut); 
    void ModalSerendipBasisSurf3DP4_Lower3(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut);     
}

#endif // BASIS_VOL_TO_SURF_MOD_DECL_H
