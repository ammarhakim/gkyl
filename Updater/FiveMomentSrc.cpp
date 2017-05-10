// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for five-moment source terms
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <FiveMomentSrc.h>
#include <iostream>
#include <vector>
#include <string>

// Makes indexing cleaner
static const unsigned MX = 1;
static const unsigned MY = 2;
static const unsigned MZ = 3;

static const unsigned EX = 0;
static const unsigned EY = 1;
static const unsigned EZ = 2;
static const unsigned BX = 3;
static const unsigned BY = 4;
static const unsigned BZ = 5;

void
gkylFiveMomentSrcRk3(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em)
{
  unsigned nFluids = sd->nFluids;
  std::vector<double> f1(5), em1(8), curr(3);
  std::vector<double> f2(5), em2(8);

  // B fields don't change
  em1[BX] = em[BX]; em1[BY] = em[BY]; em1[BZ] = em[BZ];
  em2[BX] = em[BX]; em2[BY] = em[BY]; em2[BZ] = em[BZ];

  //------------> RK Stage 1
  curr[0] = curr[1] = curr[2] = 0.0;
  for (unsigned n=0; n<nFluids; ++n)
  { // update fluids
    double *f = ff[n];
    double qmdt = dt*fd[n].charge/fd[n].mass;
    double qmdte = qmdt/sd->epsilon0;

    f1[MX] = f[MX] + qmdt*(em[EX] + f[MY]*em[BZ] - f[MZ]*em[BY]);
    f1[MY] = f[MY] + qmdt*(em[EY] + f[MZ]*em[BX] - f[MX]*em[BZ]);
    f1[MZ] = f[MZ] + qmdt*(em[EZ] + f[MX]*em[BY] - f[MY]*em[BX]);

    curr[EX] += qmdte*f[MX];
    curr[EY] += qmdte*f[MY];
    curr[EZ] += qmdte*f[MZ];
  }
  // update field
  em1[EX] = em[EX] - curr[EX];
  em1[EY] = em[EY] - curr[EY];
  em1[EZ] = em[EZ] - curr[EZ];

  //------------> RK Stage 2
  curr[0] = curr[1] = curr[2] = 0.0;
  for (unsigned n=0; n<nFluids; ++n)
  { // update fluids
    double *f = ff[n];
    double qmdt = dt*fd[n].charge/fd[n].mass;
    double qmdte = qmdt/sd->epsilon0;

    f2[MX] = 0.75*f[MX] + 0.25*(f1[MX] + qmdt*(em1[EX] + f1[MY]*em1[BZ] - f1[MZ]*em1[BY]));
    f2[MY] = 0.75*f[MY] + 0.25*(f1[MY] + qmdt*(em1[EY] + f1[MZ]*em1[BX] - f1[MX]*em1[BZ]));
    f2[MZ] = 0.75*f[MY] + 0.25*(f1[MZ] + qmdt*(em1[EZ] + f1[MX]*em1[BY] - f1[MY]*em1[BX]));

    curr[EX] += qmdte*f1[MX];
    curr[EY] += qmdte*f1[MY];
    curr[EZ] += qmdte*f1[MZ];
  }
  // update field
  em2[EX] = 0.75*em[EX] + 0.25*(em1[EX] - curr[EX]);
  em2[EY] = 0.75*em[EY] + 0.25*(em1[EY] - curr[EY]);
  em2[EZ] = 0.75*em[EZ] + 0.25*(em1[EZ] - curr[EZ]);

  double one3 = 1.0/3.0, two3 = 2.0/3.0;
  //------------> RK Stage 3
  curr[0] = curr[1] = curr[2] = 0.0;
  for (unsigned n=0; n<nFluids; ++n)
  { // update fluids
    double *f = ff[n];
    double qmdt = dt*fd[n].charge/fd[n].mass;
    double qmdte = qmdt/sd->epsilon0;

    f[MX] = one3*f[MX] + two3*(f2[MX] + qmdt*(em2[EX] + f2[MY]*em2[BZ] - f2[MZ]*em2[BY]));
    f[MY] = one3*f[MY] + two3*(f2[MY] + qmdt*(em2[EY] + f2[MZ]*em2[BX] - f2[MX]*em2[BZ]));
    f[MZ] = one3*f[MY] + two3*(f2[MZ] + qmdt*(em2[EZ] + f2[MX]*em2[BY] - f2[MY]*em2[BX]));

    curr[EX] += qmdte*f2[MX];
    curr[EY] += qmdte*f2[MY];
    curr[EZ] += qmdte*f2[MZ];
  }
  // update field
  em[EX] = one3*em[EX] + two3*(em1[EX] - curr[EX]);
  em[EY] = one3*em[EY] + two3*(em1[EY] - curr[EY]);
  em[EZ] = one3*em[EZ] + two3*(em1[EZ] - curr[EZ]);
  
}
