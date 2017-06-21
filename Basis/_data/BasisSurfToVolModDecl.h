#ifndef BASIS_SURF_TO_VOL_MOD_DECL_H
#define BASIS_SURF_TO_VOL_MOD_DECL_H

extern "C"
{
    void ModalSer1DP1_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer1DP1_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer1DP2_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer1DP2_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer1DP3_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer1DP3_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer1DP4_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer1DP4_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 

    void ModalSer2DP1_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP1_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP1_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP1_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP2_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP2_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP2_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP2_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP3_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP3_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP3_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP3_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP4_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP4_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP4_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer2DP4_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 

    void ModalSer3DP1_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP1_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP1_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP1_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP1_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP1_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP2_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP2_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP2_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP2_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP2_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP2_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP3_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP3_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP3_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP3_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP3_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP3_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP4_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP4_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP4_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP4_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP4_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalSer3DP4_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 


 
    void ModalMax1DP1_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax1DP1_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax1DP2_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax1DP2_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax1DP3_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax1DP3_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax1DP4_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax1DP4_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 

    void ModalMax2DP1_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP1_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP1_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP1_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP2_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP2_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP2_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP2_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP3_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP3_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP3_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP3_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP4_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP4_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP4_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax2DP4_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 

    void ModalMax3DP1_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP1_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP1_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP1_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP1_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP1_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP2_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP2_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP2_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP2_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP2_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP2_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP3_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP3_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP3_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP3_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP3_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP3_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP4_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP4_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP4_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP4_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP4_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
    void ModalMax3DP4_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut); 
}

#endif // BASIS_SURF_TO_VOL_MOD_DECL_H
