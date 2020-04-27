extern "C" {
    struct _ProblemState {
        double dl, ul, pl;
        double dr, ur, pr;
        double lower, upper;
        int ncell;
        double tEnd;
        double gas_gamma;
        double disLoc;
    };

    void solveRiemannProblem(_ProblemState _ps,
      double *density, double *velocity, double *pressure, double *internalenergy);
}
