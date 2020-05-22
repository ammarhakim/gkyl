/** Compute exact solution to Euler Riemann problem. */

#include <GkExactEulerRpImpl.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

/** Problem state */
struct ProblemState {
    /** Density, velocity, pressure, sound-speed on left */
    double dl, ul, pl, cl;
    /** Density, velocity, pressure, sound-speed on right */
    double dr, ur, pr, cr;
    /** Lower and upper bounds */
    double lower, upper;
    /** Domain length, time at which solution is needed */
    double domlen, tEnd;
    /** Ratio of specific heat */
    double gas_gamma;
    /** Location of discontinuity */
    double disLoc;
};

/** Solution of Riemann problem */
class Solution {
  public:
    /**
     * @param ncells Number of cells in domain.
     */
    Solution(unsigned ncell)
      : ncell(ncell),
        density(ncell),
        velocity(ncell),
        pressure(ncell),
        internalEnergy(ncell)
    {
    }
    
    /** Number of cells */
    unsigned ncell;
    /** Density, velocity, pressure and internal-energy */
    std::vector<double> density, velocity, pressure, internalEnergy;
};

void
prefun(const ProblemState& ps, double& F, double& FD,
  double P, double DK, double PK, double CK) {

  double PRATIO, QRT, AK, BK;
  double gas_gamma = ps.gas_gamma;
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

  if (P<=PK)
  {
// rarefaction wave
    PRATIO = P/PK;
    F = G4*CK*(pow(PRATIO,G1) - 1.0);
    FD = (1.0/(DK*CK))*pow(PRATIO,-G2);
  }
  else
  {
// shock wave
    AK  = G5/DK;
    BK  = G6*PK;
    QRT = sqrt(AK/(BK + P));
    F   = (P - PK)*QRT;
    FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT;
  }
}

/**
 * @param ps Problem state
 * @return Guess for pressure in star region
 */
double
guessp(const ProblemState& ps) {
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  double gas_gamma = ps.gas_gamma;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

  double CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ,
    PTL, PTR, QMAX, QUSER, UM;

  QUSER = 2.0;

  CUP  = 0.25*(DL + DR)*(CL + CR);
  PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP;
  PPV  = std::max(0.0, PPV);
  PMIN = std::min(PL,  PR);
  PMAX = std::max(PL,  PR);
  QMAX = PMAX/PMIN;

  if ((QMAX <= QUSER) && ((PMIN<=PPV) && (PPV<=PMAX)))
  {
    //std::cout << "First if" << std::endl;
    PM = PPV;
  }
  else 
  {
    if (PPV<PMIN) 
    {
      //std::cout << "Second if" << std::endl;
      PQ  = pow(PL/PR,G1);
      UM  = (PQ*UL/CL + UR/CR + G4*(PQ - 1.0))/(PQ/CL + 1.0/CR);
      PTL = 1.0 + G7*(UL - UM)/CL;
      PTR = 1.0 + G7*(UM - UR)/CR;
      PM  = 0.5*(PL*pow(PTL,G3) + PR*pow(PTR,G3));
    }
    else 
    {
      //std::cout << "Second else" << std::endl;
      GEL = sqrt((G5/DL)/(G6*PL + PPV));
      GER = sqrt((G5/DR)/(G6*PR + PPV));
      PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER);
    }
  }
  return PM;
}

void
starpu(const ProblemState& ps, double &pm, double& um) {
  double CHANGE, FL, FLD, FR, FRD, P, POLD, 
    PSTART, TOLPRE, U, UDIFF, PSCALE;

  TOLPRE = 1.0e-6;
  unsigned NRITER = 20;

  double gas_gamma = ps.gas_gamma;
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;
  PSCALE = 1.0;

// compute initial guess
  PSTART = guessp(ps);
  //- std::cout << "Initial guess " << PSTART << std::endl;

  POLD = PSTART;
  UDIFF = UR-UL;

  bool converged = false;
  for (unsigned i=1; i<NRITER; ++i) {
    prefun(ps, FL, FLD, POLD, DL, PL, CL);
    prefun(ps, FR, FRD, POLD, DR, PR, CR);
    P = POLD - (FL + FR + UDIFF)/(FLD + FRD);
    CHANGE = 2.0*fabs((P - POLD)/(P + POLD));
    //- std::cout << "Iteration " << i << " Change " << CHANGE << std::endl;
    if (CHANGE<=TOLPRE) {
      converged = true;
      break;
    }
    if (P<0.0) P = TOLPRE;
    POLD = P;
  }
  if (converged) {
    //- std::cout << "Newton iteration converged ..." << std::endl;
  }
  else {
    std::cout << "Newton iteration did not converge" << std::endl;
    exit(1);
  }
// compute velocity in star region
  U = 0.5*(UL + UR + FR - FL);
  //- std::cout << "Converged pressure " << P << std::endl;
  //- std::cout << "Converged velocity " << U << std::endl;
  pm = P;
  um = U;
}

void
sample(const ProblemState& ps, double PM, double UM, double S,
  double& D, double& U, double& P) {
  double gas_gamma = ps.gas_gamma;
// compute constants related to gamma
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

  double C, CML, CMR, PML, PMR, SHL, SHR, SL, SR, STL, STR;

  //std::cout << "S = " << S << std::endl;
  //std::cout << "UM = " << UM << std::endl;
  if (S<=UM)
  {
// left of contact discontinuity
    //std::cout << "Left of contact .." << std::endl;

    if (PM<=PL)
    {
      //std::cout << "Left rarefaction ...." << std::endl;

// left rarefaction
      SHL = UL-CL;
      if (S<=SHL)
      {
        //std::cout << "Left data state!" << std::endl;

// left data state
        D = DL;
        U = UL;
        P = PL;
      }
      else 
      {
        CML = CL*pow(PM/PL, G1);
        STL = UM-CML;

        if (S>STL)
        {
          //std::cout << "Star left state!" << std::endl;

// Star left state
          D = DL*pow(PM/PL, 1/GAMMA);
          U = UM;
          P = PM;
        } 
        else 
        {
          //std::cout << "Inside left fan!" << std::endl;
// inside left fan
          U = G5*(CL + G7*UL + S);
          C = G5*(CL + G7*(UL - S));
          D = DL*pow(C/CL, G4);
          P = PL*pow(C/CL, G3);
        }
      }
    }
    else
    {
      //std::cout << "Left shock ..." << std::endl;
// left shock
      PML = PM/PL;
      SL = UL - CL*sqrt(G2*PML + G1);

      if (S<=SL)
      {
        //std::cout << "Left data state!" << std::endl;
// point is left data state
        D = DL;
        U = UL;
        P = PL;
      }
      else
      {
        //std::cout << "Star left state!" << std::endl;
// point is star left state
        D = DL*(PML + G6)/(PML*G6 + 1.0);
        U = UM;
        P = PM;
      }
    }
  }
  else
  {
// right of contact discontinuity
    //std::cout << "Right of contact discontinuity ...." << std::endl;
    if (PM>PR)
    {
// right shock
      //std::cout << "Right shock ..." << std::endl;

      PMR = PM/PR;
      SR  = UR + CR*sqrt(G2*PMR + G1);

      if (S>=SR)
      {
        //std::cout << "Right data state!" << std::endl;
// right data state
        D = DR;
        U = UR;
        P = PR;
      }
      else
      {
        //std::cout << "Star right state!" << std::endl;
// right star state
        D = DR*(PMR + G6)/(PMR*G6 + 1.0);
        U = UM;
        P = PM;
      }
    }
    else
    {
      //std::cout << "Right rarefaction ..." << std::endl;
// right rarefaction
      SHR = UR+CR;

      if (S>=SHR)
      {
        //std::cout << "Right data state!" << std::endl;
// right data state
        D = DR;
        U = UR;
        P = PR;
      }
      else
      {
        CMR = CR*pow(PM/PR, G1);
        STR = UM + CMR;

        if (S<=STR)
        {
          //std::cout << "Star right state!" << std::endl;
// star right state
          D = DR*pow(PM/PR, 1.0/GAMMA);
          U = UM;
          P = PM;
        }
        else
        {
          //std::cout << "Inside left fan!" << std::endl;
// inside left fan
          U = G5*(-CR + G7*UR + S);
          C = G5*(CR - G7*(UR - S));
          D = DR*pow(C/CR, G4);
          P = PR*pow(C/CR, G3);
        }
      }
    }
  }
}

void
sampleWithVacuum(const ProblemState& ps, double S, double& D, double& U, double& P) {
  double gas_gamma = ps.gas_gamma;
// compute constants related to gamma
  double GAMMA, G1, G2, G3, G4, G5, G6, G7, G8;
  GAMMA = gas_gamma;
  G1 = (gas_gamma - 1)/(2*gas_gamma);
  G2 = (gas_gamma + 1)/(2*gas_gamma);
  G3 = 2*gas_gamma/(gas_gamma - 1);
  G4 = 2/(gas_gamma - 1);
  G5 = 2/(gas_gamma + 1);
  G6 = (gas_gamma - 1)/(gas_gamma + 1);
  G7 = (gas_gamma - 1)/2;
  G8 = gas_gamma - 1;
// left, right states
  double DL, UL, PL, CL, DR, UR, PR, CR;
  DL = ps.dl;
  UL = ps.ul;
  PL = ps.pl;
  CL = ps.cl;
  DR = ps.dr;
  UR = ps.ur;
  PR = ps.pr;
  CR = ps.cr;

// assume vacuum is to the right (JUST FOR NOW)
  double sstar = UL + 2*CL/G8;
  double relVel = UL-CL;
  if (S<=relVel)
  {
// left of rarefaction
    D = DL;
    U = UL;
    P = PL;
  }
  else
  {
    if ((relVel<S) && (S<=sstar))
    {
// inside rarefaction fan
      D = DL*pow(G5 + G6/CL*(UL-S), G4);
      U = G5*(CL + G7*UL + S);
      P = PL*pow(G5 + G6/CL*(UL-S), G3);
    }
    else
    {
// inside vacuum
      D = 0.0;
      U = 0.0;
      P = 0.0;
    }
  }
}

void
exactEulerRp(const ProblemState& ps, Solution& solution) {
  double g1, g2, g3, g4, g5, g6, g7, g8;

  double gas_gamma = ps.gas_gamma;
// compute constants related to gamma
  g1 = (gas_gamma - 1)/(2*gas_gamma);
  g2 = (gas_gamma + 1)/(2*gas_gamma);
  g3 = 2*gas_gamma/(gas_gamma - 1);
  g4 = 2/(gas_gamma - 1);
  g5 = 2/(gas_gamma + 1);
  g6 = (gas_gamma - 1)/(gas_gamma + 1);
  g7 = (gas_gamma - 1)/2;
  g8 = gas_gamma - 1;

// check if a vacuum is generated
  if (g4*(ps.cl+ps.cr) <= (ps.ur-ps.ul)) {
    std::cout << "Initial conditions will lead to vacuum. Aborting ..." << std::endl;
    exit(1);
  }

// compute pressure and velocity in "star" region
  double pm, um;
  starpu(ps, pm, um);

// cell spacing
  double dx = ps.domlen/solution.ncell;
// compute solution at each grid point
  for (unsigned i=0; i<solution.ncell; ++i) {
    double xpos = ps.lower + (i+0.5)*dx;
    double s = (xpos-ps.disLoc)/ps.tEnd;
// compute solution at (x,t) = (xpos-ps.disLoc, ps.tEnd)
    double dsol, usol, psol;
    sample(ps, pm, um, s, dsol, usol, psol);
    //std::cout << xpos << " " << dsol << std::endl;
// copy solution into array
    solution.density[i] = dsol;
    solution.velocity[i] = usol;
    solution.pressure[i] = psol;
    solution.internalEnergy[i] = psol/dsol/g8;
  }
}

void
exactEulerRpWithVacuum(const ProblemState& ps, Solution& solution) {

  double gas_gamma = ps.gas_gamma;
  double g8 = gas_gamma - 1;

// cell spacing
  double dx = ps.domlen/solution.ncell;
// compute solution at each grid point
  for (unsigned i=0; i<solution.ncell; ++i) {
    double xpos = ps.lower + (i+0.5)*dx;
    double s = (xpos-ps.disLoc)/ps.tEnd;
// compute solution at (x,t) = (xpos-ps.disLoc, ps.tEnd)
    double dsol, usol, psol;
    sampleWithVacuum(ps, s, dsol, usol, psol);
    //std::cout << xpos << " " << dsol << std::endl;
// copy solution into array
    solution.density[i] = dsol;
    solution.velocity[i] = usol;
    solution.pressure[i] = psol;
    if (dsol < std::numeric_limits<double>::epsilon())
      solution.internalEnergy[i] = 0.0;
    else
      solution.internalEnergy[i] = psol/dsol/g8;
  }
}

void
showInput(const ProblemState& ps) {
  std::cout << "gamma " << ps.gas_gamma << std::endl;

  std::cout << "dl " << ps.dl << std::endl;
  std::cout << "ul " << ps.ul << std::endl;
  std::cout << "pl " << ps.pl << std::endl;

  std::cout << "dr " << ps.dr << std::endl;
  std::cout << "ur " << ps.ur << std::endl;
  std::cout << "pr " << ps.pr << std::endl;

  std::cout << "lower " << ps.lower << std::endl;
  std::cout << "upper " << ps.upper << std::endl;
  std::cout << "domelen " << ps.domlen << std::endl;

  std::cout << "disLoc " << ps.disLoc << std::endl;
}

void solveRiemannProblem(const _ProblemState _ps,
  double *density, double *velocity, double *pressure, double *internalenergy) {

  ProblemState ps;
  
  ps.dl = _ps.dl; ps.ul = _ps.ul; ps.pl = _ps.pl;
  ps.dr = _ps.dr; ps.ur = _ps.ur; ps.pr = _ps.pr;
  ps.lower = _ps.lower; ps.upper = _ps.upper;
  ps.gas_gamma = _ps.gas_gamma;
  ps.tEnd = _ps.tEnd;
  ps.disLoc = _ps.disLoc;
  ps.domlen = ps.upper-ps.lower;

  //showInput(ps);
  
  unsigned ncell = _ps.ncell;

  // compute sound speeds in each region
  if (ps.dl != 0)
    ps.cl = std::sqrt(ps.gas_gamma*ps.pl/ps.dl);
  else
    ps.cl = 0.0;

  if (ps.dr != 0)
    ps.cr = std::sqrt(ps.gas_gamma*ps.pr/ps.dr);
  else
    ps.cr = 0.0;

// allocate arrays for solution
  Solution solution(ncell);

// compute solution
  if ((ps.dl != 0.0) && (ps.dr != 0.0))
    exactEulerRp(ps, solution);
  else if (ps.dr == 0)
    exactEulerRpWithVacuum(ps, solution);
  else
    exactEulerRpWithVacuum(ps, solution);

// copy solution to output arrays
  for (auto i=0; i<ncell; ++i) {
    density[i] = solution.density[i];
    velocity[i] = solution.velocity[i];
    pressure[i] = solution.pressure[i];
    internalenergy[i] = solution.internalEnergy[i];
  }
}
