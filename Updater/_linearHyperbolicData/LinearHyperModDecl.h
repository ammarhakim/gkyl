#ifndef LINEAR_HYPER_MOD_DELC_H
#define LINEAR_HYPER_MOD_DELC_H

extern "C"
{
    void LinearHyperMax1DP1_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax1DP2_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax1DP3_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax1DP4_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);

    void LinearHyperMax2DP1_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax2DP1_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax2DP2_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax2DP2_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax2DP3_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax2DP3_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax2DP4_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax2DP4_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);

    void LinearHyperMax3DP1_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP1_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP1_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP2_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP2_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP2_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP3_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP3_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP3_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP4_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP4_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperMax3DP4_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);

    void LinearHyperSer1DP1_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer1DP2_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer1DP3_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer1DP4_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);

    void LinearHyperSer2DP1_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer2DP1_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer2DP2_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer2DP2_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer2DP3_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer2DP3_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer2DP4_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer2DP4_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);

    void LinearHyperSer3DP1_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP1_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP1_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP2_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP2_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP2_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP3_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP3_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP3_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP4_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP4_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);
    void LinearHyperSer3DP4_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut);    
}    
    
#endif // LINEAR_HYPER_MOD_DELC_H
