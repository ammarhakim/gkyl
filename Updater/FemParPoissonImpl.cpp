// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for FEM Poisson solver with solve only in the last configuration space direction
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <FemParPoissonImpl.h>
#include <FemMatrices.h>

// std includes
#include <string>
#include <vector>

using namespace Eigen;
static const int DIRICHLET_BC = 0;
static const int DIRICHLET_VARIABLE_BC = 2;
static const int LO = 0;
static const int HI = 1;

double take_lastPar(const double &in,const double &b) { return b; }
void vectorSumPar(double *in, double *inout, int *len, MPI_Datatype *dptr)
{
  int i;
  for(i=0; i< *len; ++i) {
    inout[i] += in[i];
  }
}

FemParPoisson::FemParPoisson(int nz_, int ndim_, int polyOrder_, 
                       double dz_, bool z_periodic_, 
                       bcdataPar_t bc_[2], bool writeMatrix_)
  : nz(nz_), ndim(ndim_), 
    polyOrder(polyOrder_), dz(dz_), z_periodic(z_periodic_),
    writeMatrix(writeMatrix_)
{
//#define DMSG(s) std::cout << s << std::endl;
#define DMSG(s) ;

  DMSG("Inside FEM Poisson initialize");

// copy to input to class structures
  for(int i=0; i<2; i++) {
    bc[i] = bc_[i];
    if(!bc[i].isSet) {
      bc[i].type = -1;
    }
    bc1d[i] = bc[i];
  }

// prepare solver
  int nglobal = getNumParGlobalNodes(nz, ndim, polyOrder, z_periodic);
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  int nlocal_xy = 1;
  if(ndim>1) nlocal_xy = getNumLocalNodes(ndim-1, polyOrder);

  analyzed_ = false; // flag so that stiffness matrix only analyzed once

  // setup boundary condition indexing
  setupBoundaryIndices(bc, ndim, polyOrder);

// initialize global source vector
  globalSrc = VectorXd(nglobal);

// initialize modal-nodal transformation matrices
 
  MatrixXd localMass = MatrixXd::Zero(nlocal, nlocal);
  localMassModToNod = MatrixXd::Zero(nlocal, nlocal);
  localNodToMod = MatrixXd::Zero(nlocal, nlocal);
  localModToNod = MatrixXd::Zero(nlocal, nlocal);

  getParNodToModMatrix(localNodToMod, ndim, polyOrder);
  getMassMatrix(localMass, ndim, polyOrder);
  
  localModToNod = localNodToMod.transpose().inverse().transpose();
  localMassModToNod = localMass*localModToNod;

  if (writeMatrix)
  {
    std::string outName = "poisson-nodtomod"; 
    //outName += std::to_string(ndim) + "d";
    saveMarket(localNodToMod, outName);
    outName = "poisson-modtonod"; // wrong
    //outName += std::to_string(ndim) + "d";
    saveMarket(localModToNod, outName);
    outName = "poisson-massmodtonod";
    //outName += std::to_string(ndim) + "d";
    saveMarket(localMassModToNod, outName);
    outName = "poisson-mass";   
    //outName += std::to_string(ndim) + "d";
    saveMarket(localMass, outName);
    outName = "poisson-check";   
    //outName += std::to_string(ndim) + "d";
    saveMarket(localModToNod*localNodToMod, outName);
  }
  localMass.resize(0,0);

  // create MPI type for eigen triplet vector
  MPI_Type_contiguous(sizeof(Triplet<double>), MPI_BYTE, &MPI_triplet_t);
  MPI_Type_commit(&MPI_triplet_t);

  MPI_Op_create((MPI_User_function *) vectorSumPar, true, &MPI_vectorSum_op);
}

FemParPoisson::~FemParPoisson() 
{
  stiffMat.resize(0,0);
  globalSrc.resize(0);
  sourceModVec.resize(0);
  x.resize(0);
}

// get start and end indices for each boundary
void FemParPoisson::setupBoundaryIndices(bcdataPar_t bc[2], int ndim, int polyOrder)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  int nlocal_xy = 1;
  if(ndim>1) nlocal_xy = getNumLocalNodes(ndim-1, polyOrder);
  std::vector<int> lgMap(nlocal);
  // back boundary
  getParLocalToGlobalInteriorBoundary(lgMap,0,nz,ndim,polyOrder,z_periodic);
  bc[LO].istart = lgMap[0];
  bc[LO].iend = bc[LO].istart + nlocal_xy - 1;

  // front boundary
  getParLocalToGlobalInteriorBoundary(lgMap,nz-1,nz,ndim,polyOrder,z_periodic);
  bc[HI].istart = lgMap[polyOrder];
  bc[HI].iend = bc[HI].istart + nlocal_xy - 1;
}

void FemParPoisson::makeGlobalParStiffnessMatrix(double *laplacianWeight, double *modifierWeight, int idz)
{
  int nglobal = getNumParGlobalNodes(nz, ndim, polyOrder, z_periodic);
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  int nlocal_xy = 1;
  if(ndim>1) nlocal_xy = getNumLocalNodes(ndim-1, polyOrder);

  int nonzeros = nlocal*(std::pow(2.0, 1.0)+1); // estimate number of nonzeros per row

  MatrixXd localStiff = MatrixXd::Zero(nlocal, nlocal);
  MatrixXd localMass = MatrixXd::Zero(nlocal, nlocal);

  std::vector<int> lgMap(nlocal);
  stiffTripletList.reserve(nonzeros*nglobal); // estimate number of nonzero elements

  getParStiffnessMatrix(localStiff, laplacianWeight, ndim, polyOrder, dz);
  getMassMatrix(localMass, modifierWeight, ndim, polyOrder);

  getParLocalToGlobalInteriorBoundary(lgMap,idz,nz,ndim,polyOrder,z_periodic);
  
//ake triplets for constructing Eigen SparseMatrix
  for (unsigned k=0; k<nlocal; ++k)
  {
    for (unsigned m=0; m<nlocal; ++m) {
      double val = -localStiff(k,m) + localMass(k,m);
      unsigned globalIdx_k = lgMap[k];
      unsigned globalIdx_m = lgMap[m];

      stiffTripletList.push_back(Triplet<double>(globalIdx_k, globalIdx_m, val));
    }
  }
}
  
void FemParPoisson::finishGlobalParStiffnessMatrix()
{
  int nglobal = getNumParGlobalNodes(nz, ndim, polyOrder, z_periodic);
  int nlocal_z = 1;
  if(ndim==3) nlocal_z += polyOrder;

// construct SparseMatrix from triplets
  SparseMatrix<double,RowMajor> stiffMatRowMajor(nglobal, nglobal);
  stiffMatRowMajor.setFromTriplets(stiffTripletList.begin(), stiffTripletList.end());
  stiffTripletList.resize(0);

// handle Dirichlet BCs
// note: we do not need to do anything for Neumann BCs
  int nDirichlet = 0;
// set rows corresponding to Dirichlet BCs to zero
  for (int side=0; side<2; ++side)
  {
    if ((bc[side].isSet && bc[side].type == DIRICHLET_BC) ||
        (bc[side].isSet && bc[side].type == DIRICHLET_VARIABLE_BC))
    {
      int start = bc[side].istart;
      int end = bc[side].iend;
      int nRows = end-start+1;
      stiffMatRowMajor.middleRows(start, nRows) = SparseMatrix<double,RowMajor>(nRows, stiffMatRowMajor.cols());   

      nDirichlet += nRows;
    }
  }

// create column major copy of stiffMat so that we can zero columns
  stiffMat = SparseMatrix<double,ColMajor>(stiffMatRowMajor);
// de-allocate row major copy
  stiffMatRowMajor.resize(0,0);

// create matrix for modifying source when there are Dirichlet BCs
// this matrix will have nonzero entries only for the blocks consisting
// of Dirichlet cols and non-Dirichlet rows
  SparseMatrix<double,ColMajor> sourceModMat(nglobal, nglobal);

// sparse matrix for identity blocks in rows/cols corresponding to Dirichlet BCs
  SparseMatrix<double,ColMajor> dirichletIdentity(nglobal, nglobal);
  std::vector<Triplet<double> > identityTripletList;
  identityTripletList.reserve(nDirichlet); // estimate number of nonzero elements

  for (int side=0; side<2; ++side)
  {
    if ((bc[side].isSet && bc[side].type == DIRICHLET_BC) ||
        (bc[side].isSet && bc[side].type == DIRICHLET_VARIABLE_BC))
    {
      int start = bc[side].istart;
      int end = bc[side].iend;
      int nCols = end-start+1;
  // copy cols corresponding to Dirichlet BCs to sourceModMat
      sourceModMat.middleCols(start, nCols) = stiffMat.middleCols(start, nCols);
  
  // in stiffMat, set cols corresponding to Dirichlet BCs to zero
      stiffMat.middleCols(start, nCols) = SparseMatrix<double,ColMajor>(stiffMat.rows(), nCols);
  // set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
      for(int i=start; i<start+nCols; i++) {
        identityTripletList.push_back(Triplet<double>(i, i, 1.)); 
      } 
    }
  }

  // this ensures no duplicates, so that entries for these rows are at most equal to 1
  dirichletIdentity.setFromTriplets(identityTripletList.begin(), identityTripletList.end(), take_lastPar);

  stiffMat+=dirichletIdentity;

// create vector of Dirichlet values
  SparseVector<double> dirichletVec(nglobal);
  dirichletVec.reserve(nDirichlet);
  for (int side=0; side<2; ++side)
  {
    if (bc[side].isSet && bc[side].type == DIRICHLET_BC)
    {
      for(unsigned i=bc[side].istart; i<=bc[side].iend; i++) {
        dirichletVec.coeffRef(i) = bc[side].value;
      }
    }
  }

// calculate vector to subtract from source
  sourceModVec = sourceModMat*dirichletVec;

  DMSG("Finished initializing stiffness matrices");

// compute step: reorder and factorize stiffMat to prepare for solve 
  if(!analyzed_) {
    // only do analyzePattern once, assuming structure of matrix doesn't change
    solver.analyzePattern(stiffMat);
    analyzed_ = true;
  }
  solver.factorize(stiffMat);
  
  DMSG("Solver compute step complete");

  if (writeMatrix)
  {
    std::string outName = "poisson-stiffnessMatrix";
    //outName += std::to_string(ndim) + "d";
    saveMarket(stiffMat, outName);
  }
}

// clear out existing stuff in source vector: this is required
// otherwise successive calls to advance() will accumulate into source
// from prevous calls, which is of course not what we want.
void FemParPoisson::zeroGlobalSrc() {
  int nglobal = getNumParGlobalNodes(nz, ndim, polyOrder, z_periodic);
  globalSrc.setZero(nglobal);
}

// called within an indexer loop over idz
void FemParPoisson::createGlobalSrc(double* localSrcPtr, int idz, double intSrcVol)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);

  std::vector<int> lgMap(nlocal);
  std::vector<double> localSrc(nlocal), localMassSrc(nlocal);

  // copy data from input
  for (unsigned k=0; k<nlocal; ++k)
    localSrc[k] = localSrcPtr[k];

  // adjust src cell-average 
  localSrc[0] = localSrc[0] - intSrcVol;

  getParLocalToGlobalInteriorBoundary(lgMap,idz,nz,ndim,polyOrder,z_periodic);

// evaluate local mass matrix times local source
  for (unsigned k=0; k<nlocal; ++k)
  {
    localMassSrc[k] = 0.0;
    for (unsigned m=0; m<nlocal; ++m)
    {
      localMassSrc[k] += localMassModToNod(k,m)*localSrc[m];
    }
  }
// map into global source vector
  for(unsigned i=0; i<nlocal; ++i)
  {
    globalSrc.coeffRef(lgMap[i]) += localMassSrc[i];
  }
  if (writeMatrix)
  {
    std::string outName = "poisson-src-beforeBCs-";
    //outName += std::to_string(ndim) + "d";
    saveMarket(globalSrc, outName);
  }
}

void FemParPoisson::allreduceGlobalSrc(MPI_Comm comm)
{
  int nglobal = getNumParGlobalNodes(nz, ndim, polyOrder, z_periodic);
  // all reduce (sum) globalSrc
  MPI_Allreduce(MPI_IN_PLACE, globalSrc.data(), nglobal, MPI_DOUBLE, MPI_vectorSum_op, comm);
}
void FemParPoisson::allgatherGlobalStiff(MPI_Comm comm)
{
  // all gather (concatenate) stiffTripletList
  std::vector<Triplet<double> > stiffTripletListGathered;
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  stiffTripletListGathered.resize(stiffTripletList.size()*nprocs); // this resize is required.. just need to make sure to reserve enough space
  MPI_Allgather(stiffTripletList.data(), stiffTripletList.size(), MPI_triplet_t, stiffTripletListGathered.data(), stiffTripletList.size(), MPI_triplet_t, comm);
  stiffTripletList = stiffTripletListGathered;
}

void FemParPoisson::solve()
{
  int nglobal = getNumParGlobalNodes(nz, ndim, polyOrder, z_periodic);
  int nlocal_xy = 1;
  if(ndim>1) nlocal_xy = getNumLocalNodes(ndim-1, polyOrder);
// replace dirichlet nodes of global source with dirichlet values
  for (unsigned side=0; side<2; ++side)
  {
    if (bc[side].isSet && bc[side].type == DIRICHLET_BC)
    {
      for(unsigned i=bc[side].istart; i<=bc[side].iend; ++i) {
        globalSrc.coeffRef(i) = bc[side].value;
      }
    }
  }

  if (writeMatrix)
  {
    std::string outName = "poisson-src";
    //outName += std::to_string(ndim) + "d";
    saveMarket(globalSrc, outName);
  }
 
  x = Eigen::VectorXd(nglobal);
  
// solve linear system(s)
// modify non-dirichlet rows of source by subtracting sourceModVec
  //if(ndim>1 && polyOrder==1) {
  //  // for 2d or 3d p=1, we solve a 1d system for each local x-y node and concatenate the results
  //  int nglobal1d = getNumParGlobalNodes(nz, 1, polyOrder, z_periodic);
  //  for(int i=0; i<nlocal_xy; i++) {
  //    // solve in 1d at x-y node component
  //    x.segment(i*nlocal_xy, nglobal1d) = solver.solve(globalSrc.segment(i*nlocal_xy, nglobal1d)-sourceModVec); 
  //  }
  //}
  //else {
    x = solver.solve(globalSrc-sourceModVec);
  //}

  if (writeMatrix)
  {
    std::string outName = "poisson-sol";
    //outName += std::to_string(ndim) + "d";
    saveMarket(x, outName);
  }
}

void FemParPoisson::getSolution(double* solPtr, int idz)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  std::vector<int> lgMap(nlocal);
  getParLocalToGlobalInteriorBoundary(lgMap, idz,nz,ndim,polyOrder, z_periodic);

  // transform global nodal x to local modal sol
  for (unsigned k=0; k<nlocal; ++k) {
    solPtr[k] = 0.0;
    for (unsigned m=0; m<nlocal; ++m) {
      solPtr[k] += localNodToMod(k,m)*x.coeffRef(lgMap[m]);
    }
  }
}

void FemParPoisson::getNodalSolution(double* solPtr, int idz)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  std::vector<int> lgMap(nlocal);
  getParLocalToGlobalInteriorBoundary(lgMap, idz,nz,ndim,polyOrder, z_periodic);

  // transform global nodal x to local nodal sol
  for (unsigned k=0; k<nlocal; ++k) {
    solPtr[k] = x.coeffRef(lgMap[k]);
  }
}

int FemParPoisson::getNumParGlobalNodes(int nz, int ndim, int polyOrder, bool z_periodic)
{
  int numGlobalNodes = -1;
  int nlocal_xy = 1;
  if(ndim>1) nlocal_xy = getNumLocalNodes(ndim-1, polyOrder);

  if(ndim==1) {
    if (polyOrder == 1) {
      numGlobalNodes = nz + 1;
    }
    else if (polyOrder == 2) {
      numGlobalNodes = 2*nz + 1;
    }
  }
  else if(ndim==2) {
    if (polyOrder == 1) {
      numGlobalNodes = 2*(nz+1);
    }
    else if (polyOrder == 2) {
      numGlobalNodes = 2*(2*nz+1)+nz+1;
    }
  }
  else if(ndim==3) {
    if (polyOrder == 1) {
      numGlobalNodes = 4*(nz+1);
    }
    else if (polyOrder == 2) {
      numGlobalNodes = 4*(2*nz+1)+4*(nz+1);
    }
  }
  // if z-periodic, subtract number of front boundary nodes
  if(z_periodic) numGlobalNodes -= nlocal_xy;

  return numGlobalNodes;
}

int FemParPoisson::getNumLocalNodes(int ndim, int p) 
{
  int numNodes = -1;
  if(ndim==1) {
    numNodes = polyOrder + 1;
  }
  if(ndim==2) {
    if(p==1) numNodes = 4;
    else if(p==2) numNodes = 8;
  }
  else if(ndim==3) {
    if(p==1) numNodes = 8;
    else if(p==2) numNodes = 20;
  }
  return numNodes;
}

void FemParPoisson::getParLocalToGlobalInteriorBoundary(std::vector<int>& lgMap, int iz, int nz, int ndim, int polyOrder, bool z_periodic)
{
// note: this creates a global mapping in z, but not in x and y.

  int ninterior1d = polyOrder*nz - 1;

  if(ndim==1) {
    if (polyOrder == 1)
    {
      lgMap[0] = iz - 1; // 0
      lgMap[1] = iz;     // 1
    }
    else if (polyOrder == 2)
    {
      lgMap[0] = 2*iz - 1; // 0
      lgMap[1] = 2*iz;     // 1
      lgMap[2] = 2*iz + 1; // 2
    }
    // remap boundaries to end of mapping
    if(iz==0) {
      lgMap[0] = ninterior1d;
    }
    else if (iz==nz-1) {
      if(z_periodic) lgMap[polyOrder] = ninterior1d;
      else lgMap[polyOrder] = ninterior1d+1;
    }
  }
  else if (ndim==2) {
    if (polyOrder == 1)
    {
      lgMap[0] = iz - 1; // 0
      lgMap[1] = iz;     // 1
      int offset = nz - 1; // == ninterior1d
      lgMap[2] = offset + iz - 1; // 2
      lgMap[3] = offset + iz;     // 3
 
      // remap boundaries
      int ninterior2d = 2*(nz-1);
      if(iz==0) {
        lgMap[0] = ninterior2d;
        lgMap[2] = ninterior2d+1;
      }
      else if(iz==nz-1) {
        if(z_periodic) {
          lgMap[1] = ninterior2d;
          lgMap[3] = ninterior2d+1;
        } else {
          lgMap[1] = ninterior2d+2;
          lgMap[3] = ninterior2d+3;
        } 
      }
    }
    else if (polyOrder == 2)
    {
      lgMap[0] = 2*iz - 1; // 0
      lgMap[1] = 2*iz;     // 1
      lgMap[2] = 2*iz + 1; // 2
      int offset = 2*nz - 1; // == ninterior1d
      lgMap[3] = offset + iz - 1; // 3
      lgMap[4] = offset + iz;     // 4
      offset += nz - 1; // == ninterior1d for p=1
      lgMap[5] = offset + 2*iz - 1; // 5
      lgMap[6] = offset + 2*iz;     // 6
      lgMap[7] = offset + 2*iz + 1; // 7

      // remap boundaries
      int ninterior2d = 2*(2*nz+1)+nz+1 - 6;
      if(iz==0) {
        lgMap[0] = ninterior2d;
        lgMap[3] = ninterior2d+1;
        lgMap[5] = ninterior2d+2;
      }
      else if(iz==nz-1) {
        if(z_periodic) {
          lgMap[2] = ninterior2d;
          lgMap[4] = ninterior2d+1;
          lgMap[7] = ninterior2d+2;
        } else {
          lgMap[2] = ninterior2d+3;
          lgMap[4] = ninterior2d+4;
          lgMap[7] = ninterior2d+5;
        } 
      }
    }
  }
  else if (ndim==3) {
    if (polyOrder == 1)
    {
      lgMap[0] = iz - 1; // 0
      lgMap[1] = iz;     // 1
      int offset = nz - 1; // == ninterior1d
      lgMap[2] = offset + iz - 1; // 2
      lgMap[3] = offset + iz;     // 3
      lgMap[4] = 2*offset + iz - 1; // 4
      lgMap[5] = 2*offset + iz;     // 5
      lgMap[6] = 3*offset + iz - 1; // 6
      lgMap[7] = 3*offset + iz;     // 7

      // remap boundaries
      int ninterior3d = 4*(nz+1) - 8;
      if(iz==0) {
        lgMap[0] = ninterior3d;
        lgMap[2] = ninterior3d+1;
        lgMap[4] = ninterior3d+2;
        lgMap[6] = ninterior3d+3;
      }
      else if(iz==nz-1) {
        if(z_periodic) {
          lgMap[1] = ninterior3d;
          lgMap[3] = ninterior3d+1;
          lgMap[5] = ninterior3d+2;
          lgMap[7] = ninterior3d+3;
        } else {
          lgMap[1] = ninterior3d+4;
          lgMap[3] = ninterior3d+5;
          lgMap[5] = ninterior3d+6;
          lgMap[7] = ninterior3d+7;
        } 
      }
    }
    else if (polyOrder == 2)
    {
      lgMap[0] = 2*iz - 1; // 0
      lgMap[1] = 2*iz;     // 1
      lgMap[2] = 2*iz + 1; // 2
      int offset = 2*nz - 1; // == ninterior1d
      lgMap[3] = offset + iz - 1; // 3
      lgMap[4] = offset + iz;     // 4
      offset += nz - 1; // == ninterior1d for p=1
      lgMap[5] = offset + 2*iz - 1; // 5
      lgMap[6] = offset + 2*iz;     // 6
      lgMap[7] = offset + 2*iz + 1; // 7

      offset += 2*nz - 1; // == ninterior1d

      lgMap[8] = offset + iz - 1; // 8
      lgMap[9] = offset + iz;     // 9
      offset += nz - 1; // == ninterior1d for p=1
      lgMap[10] = offset + iz - 1; // 10
      lgMap[11] = offset + iz;     // 11

      offset += nz - 1; // == ninterior1d for p=1

      lgMap[12] = offset + 2*iz - 1; // 12
      lgMap[13] = offset + 2*iz;     // 13
      lgMap[14] = offset + 2*iz + 1; // 14
      offset += 2*nz - 1; // == ninterior1d
      lgMap[15] = offset + iz - 1; // 15
      lgMap[16] = offset + iz;     // 16
      offset += nz - 1; // == ninterior1d for p=1
      lgMap[17] = offset + 2*iz - 1; // 17
      lgMap[18] = offset + 2*iz;     // 18
      lgMap[19] = offset + 2*iz + 1; // 19

      // remap boundaries
      int ninterior3d = 4*(2*nz+1)+4*(nz+1) - 16;
      if(iz==0) {
        lgMap[0] = ninterior3d;
        lgMap[3] = ninterior3d+1;
        lgMap[5] = ninterior3d+2;
        lgMap[8] = ninterior3d+3;
        lgMap[10] = ninterior3d+4;
        lgMap[12] = ninterior3d+5;
        lgMap[15] = ninterior3d+6;
        lgMap[17] = ninterior3d+7;
      }
      else if(iz==nz-1) {
        if(z_periodic) {
          lgMap[2] = ninterior3d;
          lgMap[4] = ninterior3d+1;
          lgMap[7] = ninterior3d+2;
          lgMap[9] = ninterior3d+3;
          lgMap[11] = ninterior3d+4;
          lgMap[14] = ninterior3d+5;
          lgMap[16] = ninterior3d+6;
          lgMap[19] = ninterior3d+7;
        } else {
          lgMap[2] = 8 + ninterior3d;
          lgMap[4] = 8 + ninterior3d+1;
          lgMap[7] = 8 + ninterior3d+2;
          lgMap[9] = 8 + ninterior3d+3;
          lgMap[11] = 8 + ninterior3d+4;
          lgMap[14] = 8 + ninterior3d+5;
          lgMap[16] = 8 + ninterior3d+6;
          lgMap[19] = 8 + ninterior3d+7;
        } 
      }
    }
  }
}

// C wrappers for interfacing with FemParPoisson class
extern "C" void* new_FemParPoisson(int nz, int ndim, int polyOrder, double dz, bool z_periodic, bcdataPar_t bc[2], bool writeMatrix)
{
  FemParPoisson *f = new FemParPoisson(nz, ndim, polyOrder, dz, z_periodic, bc, writeMatrix);
  return reinterpret_cast<void*>(f);
}

extern "C" void delete_FemParPoisson(FemParPoisson* f)
{
  delete f;
}

extern "C" void makeParGlobalStiff(FemParPoisson* f, double *laplacianWeight, double *modifierWeight, int idz)
{
  f->makeGlobalParStiffnessMatrix(laplacianWeight, modifierWeight, idz);
}

extern "C" void finishParGlobalStiff(FemParPoisson* f)
{
  f->finishGlobalParStiffnessMatrix();
}

extern "C" void createParGlobalSrc(FemParPoisson* f, double* localSrcPtr, int idz, double intSrcVol)
{
  f->createGlobalSrc(localSrcPtr, idz, intSrcVol);
}

extern "C" void zeroParGlobalSrc(FemParPoisson* f)
{
  f->zeroGlobalSrc();
} 

extern "C" void allreduceParGlobalSrc(FemParPoisson* f, MPI_Comm comm)
{
  f->allreduceGlobalSrc(comm);
}

extern "C" void allgatherParGlobalStiff(FemParPoisson* f, MPI_Comm comm)
{
  f->allgatherGlobalStiff(comm);
}

extern "C" void solvePar(FemParPoisson* f)
{
  f->solve();
}

extern "C" void getSolutionPar(FemParPoisson* f, double* ptr, int idz)
{
  f->getSolution(ptr, idz);
}

extern "C" void getNodalSolutionPar(FemParPoisson* f, double* ptr, int idz)
{
  f->getNodalSolution(ptr, idz);
}

