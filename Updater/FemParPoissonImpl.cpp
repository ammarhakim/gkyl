// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for FEM Poisson solver with solve only in the last configuration space direction
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <FemParPoissonImpl.h>

// std includes
#include <string>
#include <vector>

using namespace Eigen;
static const int DIRICHLET_BC = 0;
static const int NEUMANN_BC = 1;
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
                       bcdataPar_t bc_[2], bool writeMatrix_,
                       double laplacianWeight_, double modifierConstant_) 
  : nz(nz_), ndim(ndim_), 
    polyOrder(polyOrder_), dz(dz_), z_periodic(z_periodic_),
    writeMatrix(writeMatrix_),
    laplacianWeight(laplacianWeight_),
    modifierConstant(modifierConstant_)
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

  adjustSource = false;
  if(z_periodic && modifierConstant==0.0) adjustSource = true;

// prepare solver
  int nglobal = getNumParGlobalNodes(nz, ndim, polyOrder, z_periodic);
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  int nlocal_xy = 1;
  if(ndim>1) nlocal_xy = getNumLocalNodes(ndim-1, polyOrder);

  // setup boundary indices and make global stiffness matrix(es)
  setupBoundaryIndices(bc, ndim, polyOrder);
  //if(ndim>1 && polyOrder==1) {
  //  // trick: use a 2d or 3d p=1 global par stiffness matrix that has a diagonal 
  //  // block structure with blocks for each of the x-y nodes.
  //  // each block is just a 1d p=1 global par stiffness matrix.
  //  // we make a single 1d p=1 global par stiffness matrix, and we will use it multiple times.
  //  // the 2d or 3d p=1 global mass matrix will have the same structure.
  //  // this trick is effectively a pre-factorization for the solver.
  //  setupBoundaryIndices(bc1d, 1, polyOrder);
  //  makeGlobalParStiffnessMatrix(stiffMat, sourceModVec, 1, polyOrder, bc1d);
  //} 
  //else {
    // make full global stiffness matrix
    makeGlobalParStiffnessMatrix(stiffMat, sourceModVec, ndim, polyOrder, bc);
  //}
  DMSG("Finished initializing stiffness matrices");
  stiffMat.makeCompressed();

// compute step: reorder and factorize stiffMat to prepare for solve 
  solver.compute(stiffMat);
  
  DMSG("Solver compute step complete");


  if (writeMatrix)
  {
    std::string outName = "poisson-stiffnessMatrix";
    outName += std::to_string(ndim) + "d";
    saveMarket(stiffMat, outName);
  }

// initialize global source vector
  globalSrc = VectorXd(nglobal);

// initialize modal-nodal transformation matrices
 
  MatrixXd localMass = MatrixXd::Zero(nlocal, nlocal);
  localMassModToNod = MatrixXd::Zero(nlocal, nlocal);
  localNodToMod = MatrixXd::Zero(nlocal, nlocal);
  localModToNod = MatrixXd::Zero(nlocal, nlocal);

  getNodToModMatrix(localNodToMod, ndim, polyOrder);
  getMassMatrix(localMass, ndim, polyOrder);
  
  localModToNod = localNodToMod.transpose().inverse().transpose();
  localMassModToNod = localMass*localModToNod;

  if (writeMatrix)
  {
    std::string outName = "poisson-nodtomod"; 
    outName += std::to_string(ndim) + "d";
    saveMarket(localNodToMod, outName);
    outName = "poisson-modtonod"; // wrong
    outName += std::to_string(ndim) + "d";
    saveMarket(localModToNod, outName);
    outName = "poisson-massmodtonod";
    outName += std::to_string(ndim) + "d";
    saveMarket(localMassModToNod, outName);
    outName = "poisson-mass";   
    outName += std::to_string(ndim) + "d";
    saveMarket(localMass, outName);
    outName = "poisson-check";   
    outName += std::to_string(ndim) + "d";
    saveMarket(localModToNod*localNodToMod, outName);
  }
  localMass.resize(0,0);

  // create MPI type for globalSrc vector
  MPI_Type_contiguous(nglobal, MPI_DOUBLE, &MPI_vector_t);
  MPI_Type_commit(&MPI_vector_t);

  MPI_Op_create((MPI_User_function *) vectorSumPar, true, &MPI_vectorSum_op);
}


void FemParPoisson::allreduceGlobalSrc(MPI_Comm comm)
{
  int nglobal = getNumParGlobalNodes(nz, ndim, polyOrder, z_periodic);
  MPI_Allreduce(MPI_IN_PLACE, globalSrc.data(), nglobal, MPI_DOUBLE, MPI_vectorSum_op, comm);
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

void FemParPoisson::makeGlobalParStiffnessMatrix(
     Eigen::SparseMatrix<double,Eigen::ColMajor>& stiffMat, 
     Eigen::VectorXd& sourceModVec_,
     int ndim, int polyOrder, bcdataPar_t bc[2])
{
  int nglobal = getNumParGlobalNodes(nz, ndim, polyOrder, z_periodic);
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  int nlocal_xy = 1;
  if(ndim>1) nlocal_xy = getNumLocalNodes(ndim-1, polyOrder);

  int nonzeros = nlocal*(std::pow(2.0, 1.0*ndim)+1); // estimate number of nonzeros per row

  MatrixXd localStiff = MatrixXd::Zero(nlocal, nlocal);
  MatrixXd localMass = MatrixXd::Zero(nlocal, nlocal);

  std::vector<int> lgMap(nlocal);
  std::vector<Triplet<double> > tripletList, identityTripletList;
  tripletList.reserve(nonzeros*nglobal); // estimate number of nonzero elements

  getParStiffnessMatrix(localStiff, ndim, polyOrder, dz);
  getMassMatrix(localMass, ndim, polyOrder);

  // loop over global region
  for(int idz=0; idz<nz; idz++) {
    getParLocalToGlobalInteriorBoundary(lgMap,idz,nz,ndim,polyOrder,z_periodic);
    
// make triplets for constructing Eigen SparseMatrix
    for (unsigned k=0; k<nlocal; ++k)
    {
      for (unsigned m=0; m<nlocal; ++m) {
        double val = -laplacianWeight*localStiff(k,m) + modifierConstant*localMass(k,m);
        unsigned globalIdx_k = lgMap[k];
        unsigned globalIdx_m = lgMap[m];

        tripletList.push_back(Triplet<double>(globalIdx_k, globalIdx_m, val));
      }
    }
  }

// construct SparseMatrix from triplets
  SparseMatrix<double,RowMajor> stiffMatRowMajor(nglobal, nglobal);
  stiffMatRowMajor.setFromTriplets(tripletList.begin(), tripletList.end());

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

  if (modifierConstant!=0.0 && laplacianWeight==0.0) {
    stiffMat+=modifierConstant*dirichletIdentity;
  } else {
    stiffMat+=dirichletIdentity;
  }

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
  sourceModVec_ = sourceModMat*dirichletVec;
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
    outName += std::to_string(ndim) + "d";
    saveMarket(globalSrc, outName);
  }
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
    outName += std::to_string(ndim) + "d";
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
    outName += std::to_string(ndim) + "d";
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

void FemParPoisson::getParStiffnessMatrix(Eigen::MatrixXd& stiff, int ndim, int polyOrder, double dz)
{
  double dfacz2 = 4.0/(dz*dz);

  if(ndim==1) {
    if(polyOrder == 1) 
    {
      stiff << 0.5*dfacz2,-0.5*dfacz2,-0.5*dfacz2,0.5*dfacz2;   
    }
    else if (polyOrder == 2) 
    {
      stiff << 1.166666666666667*dfacz2,-1.333333333333333*dfacz2,0.1666666666666667*dfacz2,-1.333333333333333*dfacz2,2.666666666666667*dfacz2,-1.333333333333333*dfacz2,0.1666666666666667*dfacz2,-1.333333333333333*dfacz2,1.166666666666667*dfacz2;
    }
  }
  else if(ndim==2) {
    if(polyOrder == 1) 
    {
      //// for 2d p=1 assume the par stiffness matrix has diagonal block structure.
      //// each diagonal block corresponds to the 1d stiffness matrix 
      //// in the plane of one of the two x nodes
      //int nlocal_1d = getNumLocalNodes(1, polyOrder);
      //MatrixXd stiff1d_p1(nlocal_1d, nlocal_1d);
      //MatrixXd zeros = MatrixXd::Zero(nlocal_1d, nlocal_1d);
      //getParStiffnessMatrix(stiff1d_p1, 1, polyOrder, dz);
      //// initialize stiffness matrix with block structure
      //stiff << stiff1d_p1, zeros,
      //         zeros,      stiff1d_p1;
     
      stiff << 0.3333333333333333*dfacz2,-0.3333333333333333*dfacz2,0.1666666666666667*dfacz2,-0.1666666666666667*dfacz2,-0.3333333333333333*dfacz2,0.3333333333333333*dfacz2,-0.1666666666666667*dfacz2,0.1666666666666667*dfacz2,0.1666666666666667*dfacz2,-0.1666666666666667*dfacz2,0.3333333333333333*dfacz2,-0.3333333333333333*dfacz2,-0.1666666666666667*dfacz2,0.1666666666666667*dfacz2,-0.3333333333333333*dfacz2,0.3333333333333333*dfacz2;
    }
    else if (polyOrder == 2)
    {
      stiff << 0.5777777777777777*dfacz2,-0.8888888888888888*dfacz2,0.3111111111111111*dfacz2,0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.1888888888888889*dfacz2,-0.4444444444444444*dfacz2,0.2555555555555555*dfacz2,-0.8888888888888888*dfacz2,1.777777777777778*dfacz2,-0.8888888888888888*dfacz2,0.0,0.0,-0.4444444444444444*dfacz2,0.8888888888888888*dfacz2,-0.4444444444444444*dfacz2,0.3111111111111111*dfacz2,-0.8888888888888888*dfacz2,0.5777777777777777*dfacz2,-0.06666666666666667*dfacz2,0.06666666666666667*dfacz2,0.2555555555555555*dfacz2,-0.4444444444444444*dfacz2,0.1888888888888889*dfacz2,0.06666666666666667*dfacz2,0.0,-0.06666666666666667*dfacz2,0.5333333333333333*dfacz2,-0.5333333333333333*dfacz2,0.06666666666666667*dfacz2,0.0,-0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.0,0.06666666666666667*dfacz2,-0.5333333333333333*dfacz2,0.5333333333333333*dfacz2,-0.06666666666666667*dfacz2,0.0,0.06666666666666667*dfacz2,0.1888888888888889*dfacz2,-0.4444444444444444*dfacz2,0.2555555555555555*dfacz2,0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.5777777777777777*dfacz2,-0.8888888888888888*dfacz2,0.3111111111111111*dfacz2,-0.4444444444444444*dfacz2,0.8888888888888888*dfacz2,-0.4444444444444444*dfacz2,0.0,0.0,-0.8888888888888888*dfacz2,1.777777777777778*dfacz2,-0.8888888888888888*dfacz2,0.2555555555555555*dfacz2,-0.4444444444444444*dfacz2,0.1888888888888889*dfacz2,-0.06666666666666667*dfacz2,0.06666666666666667*dfacz2,0.3111111111111111*dfacz2,-0.8888888888888888*dfacz2,0.5777777777777777*dfacz2;
    }
  }
  else if(ndim==3) {
    if(polyOrder==1) {
      //// for 3d p=1 assume the par stiffness matrix has diagonal block structure.
      //// each diagonal block corresponds to the 1d stiffness matrix 
      //// in the plane of one of the four x-y nodes
      //int nlocal_1d = getNumLocalNodes(1, polyOrder);
      //MatrixXd stiff1d_p1(nlocal_1d, nlocal_1d);
      //MatrixXd zeros = MatrixXd::Zero(nlocal_1d, nlocal_1d);
      //getParStiffnessMatrix(stiff1d_p1, 1, polyOrder, dz);
      //// initialize stiffness matrix with block structure
      //stiff << stiff1d_p1,	zeros,		zeros,		zeros, 
      //         zeros,		stiff1d_p1,  	zeros,		zeros,
      //         zeros, 		zeros, 		stiff1d_p1,	zeros,
      //         zeros,		zeros,		zeros,		stiff1d_p1;
      
      stiff << 0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,0.05555555555555555*dfacz2,-0.05555555555555555*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,-0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,-0.05555555555555555*dfacz2,0.05555555555555555*dfacz2,0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,0.05555555555555555*dfacz2,-0.05555555555555555*dfacz2,0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,-0.05555555555555555*dfacz2,0.05555555555555555*dfacz2,-0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,0.05555555555555555*dfacz2,-0.05555555555555555*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,-0.05555555555555555*dfacz2,0.05555555555555555*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,-0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,0.05555555555555555*dfacz2,-0.05555555555555555*dfacz2,0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,-0.05555555555555555*dfacz2,0.05555555555555555*dfacz2,-0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,-0.1111111111111111*dfacz2,0.1111111111111111*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2;
    }
    else if(polyOrder==2) {
      stiff << 0.362962962962963*dfacz2,-0.5925925925925926*dfacz2,0.2296296296296296*dfacz2,-0.06666666666666667*dfacz2,0.06666666666666667*dfacz2,0.1703703703703704*dfacz2,-0.2962962962962963*dfacz2,0.1259259259259259*dfacz2,-0.06666666666666667*dfacz2,0.06666666666666667*dfacz2,-0.08888888888888889*dfacz2,0.08888888888888889*dfacz2,0.1703703703703704*dfacz2,-0.2962962962962963*dfacz2,0.1259259259259259*dfacz2,-0.08888888888888889*dfacz2,0.08888888888888889*dfacz2,0.1074074074074074*dfacz2,-0.1481481481481481*dfacz2,0.04074074074074074*dfacz2,-0.5925925925925926*dfacz2,1.185185185185185*dfacz2,-0.5925925925925926*dfacz2,0.0,0.0,-0.2962962962962963*dfacz2,0.5925925925925926*dfacz2,-0.2962962962962963*dfacz2,0.0,0.0,0.0,0.0,-0.2962962962962963*dfacz2,0.5925925925925926*dfacz2,-0.2962962962962963*dfacz2,0.0,0.0,-0.1481481481481481*dfacz2,0.2962962962962963*dfacz2,-0.1481481481481481*dfacz2,0.2296296296296296*dfacz2,-0.5925925925925926*dfacz2,0.362962962962963*dfacz2,0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.1259259259259259*dfacz2,-0.2962962962962963*dfacz2,0.1703703703703704*dfacz2,0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.08888888888888889*dfacz2,-0.08888888888888889*dfacz2,0.1259259259259259*dfacz2,-0.2962962962962963*dfacz2,0.1703703703703704*dfacz2,0.08888888888888889*dfacz2,-0.08888888888888889*dfacz2,0.04074074074074074*dfacz2,-0.1481481481481481*dfacz2,0.1074074074074074*dfacz2,-0.06666666666666667*dfacz2,0.0,0.06666666666666667*dfacz2,0.3555555555555556*dfacz2,-0.3555555555555556*dfacz2,-0.06666666666666667*dfacz2,0.0,0.06666666666666667*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,-0.08888888888888889*dfacz2,0.0,0.08888888888888889*dfacz2,0.1777777777777778*dfacz2,-0.1777777777777778*dfacz2,-0.08888888888888889*dfacz2,0.0,0.08888888888888889*dfacz2,0.06666666666666667*dfacz2,0.0,-0.06666666666666667*dfacz2,-0.3555555555555556*dfacz2,0.3555555555555556*dfacz2,0.06666666666666667*dfacz2,0.0,-0.06666666666666667*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,0.08888888888888889*dfacz2,0.0,-0.08888888888888889*dfacz2,-0.1777777777777778*dfacz2,0.1777777777777778*dfacz2,0.08888888888888889*dfacz2,0.0,-0.08888888888888889*dfacz2,0.1703703703703704*dfacz2,-0.2962962962962963*dfacz2,0.1259259259259259*dfacz2,-0.06666666666666667*dfacz2,0.06666666666666667*dfacz2,0.362962962962963*dfacz2,-0.5925925925925926*dfacz2,0.2296296296296296*dfacz2,-0.08888888888888889*dfacz2,0.08888888888888889*dfacz2,-0.06666666666666667*dfacz2,0.06666666666666667*dfacz2,0.1074074074074074*dfacz2,-0.1481481481481481*dfacz2,0.04074074074074074*dfacz2,-0.08888888888888889*dfacz2,0.08888888888888889*dfacz2,0.1703703703703704*dfacz2,-0.2962962962962963*dfacz2,0.1259259259259259*dfacz2,-0.2962962962962963*dfacz2,0.5925925925925926*dfacz2,-0.2962962962962963*dfacz2,0.0,0.0,-0.5925925925925926*dfacz2,1.185185185185185*dfacz2,-0.5925925925925926*dfacz2,0.0,0.0,0.0,0.0,-0.1481481481481481*dfacz2,0.2962962962962963*dfacz2,-0.1481481481481481*dfacz2,0.0,0.0,-0.2962962962962963*dfacz2,0.5925925925925926*dfacz2,-0.2962962962962963*dfacz2,0.1259259259259259*dfacz2,-0.2962962962962963*dfacz2,0.1703703703703704*dfacz2,0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.2296296296296296*dfacz2,-0.5925925925925926*dfacz2,0.362962962962963*dfacz2,0.08888888888888889*dfacz2,-0.08888888888888889*dfacz2,0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.04074074074074074*dfacz2,-0.1481481481481481*dfacz2,0.1074074074074074*dfacz2,0.08888888888888889*dfacz2,-0.08888888888888889*dfacz2,0.1259259259259259*dfacz2,-0.2962962962962963*dfacz2,0.1703703703703704*dfacz2,-0.06666666666666667*dfacz2,0.0,0.06666666666666667*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,-0.08888888888888889*dfacz2,0.0,0.08888888888888889*dfacz2,0.3555555555555556*dfacz2,-0.3555555555555556*dfacz2,0.1777777777777778*dfacz2,-0.1777777777777778*dfacz2,-0.06666666666666667*dfacz2,0.0,0.06666666666666667*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,-0.08888888888888889*dfacz2,0.0,0.08888888888888889*dfacz2,0.06666666666666667*dfacz2,0.0,-0.06666666666666667*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,0.08888888888888889*dfacz2,0.0,-0.08888888888888889*dfacz2,-0.3555555555555556*dfacz2,0.3555555555555556*dfacz2,-0.1777777777777778*dfacz2,0.1777777777777778*dfacz2,0.06666666666666667*dfacz2,0.0,-0.06666666666666667*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,0.08888888888888889*dfacz2,0.0,-0.08888888888888889*dfacz2,-0.08888888888888889*dfacz2,0.0,0.08888888888888889*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,-0.06666666666666667*dfacz2,0.0,0.06666666666666667*dfacz2,0.1777777777777778*dfacz2,-0.1777777777777778*dfacz2,0.3555555555555556*dfacz2,-0.3555555555555556*dfacz2,-0.08888888888888889*dfacz2,0.0,0.08888888888888889*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,-0.06666666666666667*dfacz2,0.0,0.06666666666666667*dfacz2,0.08888888888888889*dfacz2,0.0,-0.08888888888888889*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,0.06666666666666667*dfacz2,0.0,-0.06666666666666667*dfacz2,-0.1777777777777778*dfacz2,0.1777777777777778*dfacz2,-0.3555555555555556*dfacz2,0.3555555555555556*dfacz2,0.08888888888888889*dfacz2,0.0,-0.08888888888888889*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,0.06666666666666667*dfacz2,0.0,-0.06666666666666667*dfacz2,0.1703703703703704*dfacz2,-0.2962962962962963*dfacz2,0.1259259259259259*dfacz2,-0.08888888888888889*dfacz2,0.08888888888888889*dfacz2,0.1074074074074074*dfacz2,-0.1481481481481481*dfacz2,0.04074074074074074*dfacz2,-0.06666666666666667*dfacz2,0.06666666666666667*dfacz2,-0.08888888888888889*dfacz2,0.08888888888888889*dfacz2,0.362962962962963*dfacz2,-0.5925925925925926*dfacz2,0.2296296296296296*dfacz2,-0.06666666666666667*dfacz2,0.06666666666666667*dfacz2,0.1703703703703704*dfacz2,-0.2962962962962963*dfacz2,0.1259259259259259*dfacz2,-0.2962962962962963*dfacz2,0.5925925925925926*dfacz2,-0.2962962962962963*dfacz2,0.0,0.0,-0.1481481481481481*dfacz2,0.2962962962962963*dfacz2,-0.1481481481481481*dfacz2,0.0,0.0,0.0,0.0,-0.5925925925925926*dfacz2,1.185185185185185*dfacz2,-0.5925925925925926*dfacz2,0.0,0.0,-0.2962962962962963*dfacz2,0.5925925925925926*dfacz2,-0.2962962962962963*dfacz2,0.1259259259259259*dfacz2,-0.2962962962962963*dfacz2,0.1703703703703704*dfacz2,0.08888888888888889*dfacz2,-0.08888888888888889*dfacz2,0.04074074074074074*dfacz2,-0.1481481481481481*dfacz2,0.1074074074074074*dfacz2,0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.08888888888888889*dfacz2,-0.08888888888888889*dfacz2,0.2296296296296296*dfacz2,-0.5925925925925926*dfacz2,0.362962962962963*dfacz2,0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.1259259259259259*dfacz2,-0.2962962962962963*dfacz2,0.1703703703703704*dfacz2,-0.08888888888888889*dfacz2,0.0,0.08888888888888889*dfacz2,0.1777777777777778*dfacz2,-0.1777777777777778*dfacz2,-0.08888888888888889*dfacz2,0.0,0.08888888888888889*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,-0.06666666666666667*dfacz2,0.0,0.06666666666666667*dfacz2,0.3555555555555556*dfacz2,-0.3555555555555556*dfacz2,-0.06666666666666667*dfacz2,0.0,0.06666666666666667*dfacz2,0.08888888888888889*dfacz2,0.0,-0.08888888888888889*dfacz2,-0.1777777777777778*dfacz2,0.1777777777777778*dfacz2,0.08888888888888889*dfacz2,0.0,-0.08888888888888889*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,-0.2222222222222222*dfacz2,0.2222222222222222*dfacz2,0.06666666666666667*dfacz2,0.0,-0.06666666666666667*dfacz2,-0.3555555555555556*dfacz2,0.3555555555555556*dfacz2,0.06666666666666667*dfacz2,0.0,-0.06666666666666667*dfacz2,0.1074074074074074*dfacz2,-0.1481481481481481*dfacz2,0.04074074074074074*dfacz2,-0.08888888888888889*dfacz2,0.08888888888888889*dfacz2,0.1703703703703704*dfacz2,-0.2962962962962963*dfacz2,0.1259259259259259*dfacz2,-0.08888888888888889*dfacz2,0.08888888888888889*dfacz2,-0.06666666666666667*dfacz2,0.06666666666666667*dfacz2,0.1703703703703704*dfacz2,-0.2962962962962963*dfacz2,0.1259259259259259*dfacz2,-0.06666666666666667*dfacz2,0.06666666666666667*dfacz2,0.362962962962963*dfacz2,-0.5925925925925926*dfacz2,0.2296296296296296*dfacz2,-0.1481481481481481*dfacz2,0.2962962962962963*dfacz2,-0.1481481481481481*dfacz2,0.0,0.0,-0.2962962962962963*dfacz2,0.5925925925925926*dfacz2,-0.2962962962962963*dfacz2,0.0,0.0,0.0,0.0,-0.2962962962962963*dfacz2,0.5925925925925926*dfacz2,-0.2962962962962963*dfacz2,0.0,0.0,-0.5925925925925926*dfacz2,1.185185185185185*dfacz2,-0.5925925925925926*dfacz2,0.04074074074074074*dfacz2,-0.1481481481481481*dfacz2,0.1074074074074074*dfacz2,0.08888888888888889*dfacz2,-0.08888888888888889*dfacz2,0.1259259259259259*dfacz2,-0.2962962962962963*dfacz2,0.1703703703703704*dfacz2,0.08888888888888889*dfacz2,-0.08888888888888889*dfacz2,0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.1259259259259259*dfacz2,-0.2962962962962963*dfacz2,0.1703703703703704*dfacz2,0.06666666666666667*dfacz2,-0.06666666666666667*dfacz2,0.2296296296296296*dfacz2,-0.5925925925925926*dfacz2,0.362962962962963*dfacz2;
    }
  }
}

void FemParPoisson::getMassMatrix(Eigen::MatrixXd& mass, int ndim, int polyOrder)
{
  if (ndim==1) {
    if(polyOrder == 1) {
      mass << 0.6666666666666666,0.3333333333333333,0.3333333333333333,0.6666666666666666;
    }
    else if(polyOrder == 2) {
      mass << 0.2666666666666667,0.1333333333333333,-0.06666666666666667,0.1333333333333333,1.066666666666667,0.1333333333333333,-0.06666666666666667,0.1333333333333333,0.2666666666666667;
    }
  }
  else if(ndim==2) {
    if(polyOrder==1) {
      //// for 2d p=1 assume the mass matrix has diagonal block structure.
      //// each diagonal block corresponds to the 1d mass matrix 
      //// in the plane of one of the two x nodes
      //int nlocal_1d = getNumLocalNodes(1, polyOrder);
      //MatrixXd mass1d_p1(nlocal_1d, nlocal_1d);
      //MatrixXd zeros = MatrixXd::Zero(nlocal_1d, nlocal_1d);
      //getMassMatrix(mass1d_p1, 1, polyOrder);
      //// initialize massness matrix with block structure
      //mass << mass1d_p1, zeros,
      //         zeros,      mass1d_p1;

      mass << 0.4444444444444444,0.2222222222222222,0.2222222222222222,0.1111111111111111,0.2222222222222222,0.4444444444444444,0.1111111111111111,0.2222222222222222,0.2222222222222222,0.1111111111111111,0.4444444444444444,0.2222222222222222,0.1111111111111111,0.2222222222222222,0.2222222222222222,0.4444444444444444;
    }
    else if(polyOrder==2) {
      mass << 0.1333333333333333,-0.1333333333333333,0.04444444444444445,-0.1333333333333333,-0.1777777777777778,0.04444444444444445,-0.1777777777777778,0.06666666666666667,-0.1333333333333333,0.7111111111111111,-0.1333333333333333,0.4444444444444444,0.4444444444444444,-0.1777777777777778,0.3555555555555556,-0.1777777777777778,0.04444444444444445,-0.1333333333333333,0.1333333333333333,-0.1777777777777778,-0.1333333333333333,0.06666666666666667,-0.1777777777777778,0.04444444444444445,-0.1333333333333333,0.4444444444444444,-0.1777777777777778,0.7111111111111111,0.3555555555555556,-0.1333333333333333,0.4444444444444444,-0.1777777777777778,-0.1777777777777778,0.4444444444444444,-0.1333333333333333,0.3555555555555556,0.7111111111111111,-0.1777777777777778,0.4444444444444444,-0.1333333333333333,0.04444444444444445,-0.1777777777777778,0.06666666666666667,-0.1333333333333333,-0.1777777777777778,0.1333333333333333,-0.1333333333333333,0.04444444444444445,-0.1777777777777778,0.3555555555555556,-0.1777777777777778,0.4444444444444444,0.4444444444444444,-0.1333333333333333,0.7111111111111111,-0.1333333333333333,0.06666666666666667,-0.1777777777777778,0.04444444444444445,-0.1777777777777778,-0.1333333333333333,0.04444444444444445,-0.1333333333333333,0.1333333333333333;
    }
  }
  else if(ndim==3) {
    if(polyOrder==1) {
      //// for 3d p=1 assume the mass matrix has diagonal block structure.
      //// each diagonal block corresponds to the 1d mass matrix 
      //// in the plane of one of the four x-y nodes
      //int nlocal_1d = getNumLocalNodes(1, polyOrder);
      //MatrixXd mass1d_p1(nlocal_1d, nlocal_1d);
      //MatrixXd zeros = MatrixXd::Zero(nlocal_1d, nlocal_1d);
      //getMassMatrix(mass1d_p1, 1, polyOrder);
      //// initialize massness matrix with block structure
      //mass << mass1d_p1,	zeros,		zeros,		zeros, 
      //         zeros,		mass1d_p1,  	zeros,		zeros,
      //         zeros, 		zeros, 		mass1d_p1,	zeros,
      //         zeros,		zeros,		zeros,		mass1d_p1;

      mass << 0.2962962962962963,0.1481481481481481,0.1481481481481481,0.07407407407407407,0.1481481481481481,0.07407407407407407,0.07407407407407407,0.03703703703703703,0.1481481481481481,0.2962962962962963,0.07407407407407407,0.1481481481481481,0.07407407407407407,0.1481481481481481,0.03703703703703703,0.07407407407407407,0.1481481481481481,0.07407407407407407,0.2962962962962963,0.1481481481481481,0.07407407407407407,0.03703703703703703,0.1481481481481481,0.07407407407407407,0.07407407407407407,0.1481481481481481,0.1481481481481481,0.2962962962962963,0.03703703703703703,0.07407407407407407,0.07407407407407407,0.1481481481481481,0.1481481481481481,0.07407407407407407,0.07407407407407407,0.03703703703703703,0.2962962962962963,0.1481481481481481,0.1481481481481481,0.07407407407407407,0.07407407407407407,0.1481481481481481,0.03703703703703703,0.07407407407407407,0.1481481481481481,0.2962962962962963,0.07407407407407407,0.1481481481481481,0.07407407407407407,0.03703703703703703,0.1481481481481481,0.07407407407407407,0.1481481481481481,0.07407407407407407,0.2962962962962963,0.1481481481481481,0.03703703703703703,0.07407407407407407,0.07407407407407407,0.1481481481481481,0.07407407407407407,0.1481481481481481,0.1481481481481481,0.2962962962962963;
    }
    else if(polyOrder==2) {
      mass << 0.2074074074074074,-0.237037037037037,0.162962962962963,-0.237037037037037,-0.1925925925925926,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.237037037037037,-0.1925925925925926,-0.1925925925925926,-0.1333333333333333,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,-0.1333333333333333,0.1481481481481481,-0.1333333333333333,0.1259259259259259,-0.237037037037037,0.4740740740740741,-0.237037037037037,0.2962962962962963,0.2962962962962963,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.2962962962962963,0.2962962962962963,0.1481481481481481,0.1481481481481481,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.1481481481481481,0.1481481481481481,-0.1333333333333333,0.1185185185185185,-0.1333333333333333,0.162962962962963,-0.237037037037037,0.2074074074074074,-0.1925925925925926,-0.237037037037037,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1925925925925926,-0.237037037037037,-0.1333333333333333,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1333333333333333,-0.1925925925925926,0.1259259259259259,-0.1333333333333333,0.1481481481481481,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.4740740740740741,0.237037037037037,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.2962962962962963,0.1481481481481481,0.2962962962962963,0.1481481481481481,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.237037037037037,0.1185185185185185,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.237037037037037,0.4740740740740741,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.1481481481481481,0.2962962962962963,0.1481481481481481,0.2962962962962963,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.1185185185185185,0.237037037037037,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.237037037037037,-0.1925925925925926,0.2074074074074074,-0.237037037037037,0.162962962962963,-0.1925925925925926,-0.1333333333333333,-0.237037037037037,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.1259259259259259,-0.1925925925925926,-0.1333333333333333,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.2962962962962963,0.2962962962962963,-0.237037037037037,0.4740740740740741,-0.237037037037037,0.1481481481481481,0.1481481481481481,0.2962962962962963,0.2962962962962963,-0.1333333333333333,0.1185185185185185,-0.1333333333333333,0.1481481481481481,0.1481481481481481,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1925925925925926,-0.237037037037037,0.162962962962963,-0.237037037037037,0.2074074074074074,-0.1333333333333333,-0.1925925925925926,-0.1925925925925926,-0.237037037037037,0.1259259259259259,-0.1333333333333333,0.1481481481481481,-0.1333333333333333,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.2962962962962963,0.1481481481481481,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.4740740740740741,0.237037037037037,0.237037037037037,0.1185185185185185,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.2962962962962963,0.1481481481481481,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.1481481481481481,0.2962962962962963,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.237037037037037,0.4740740740740741,0.1185185185185185,0.237037037037037,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.1481481481481481,0.2962962962962963,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.2962962962962963,0.1481481481481481,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.237037037037037,0.1185185185185185,0.4740740740740741,0.237037037037037,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.2962962962962963,0.1481481481481481,-0.237037037037037,0.2962962962962963,-0.1925925925925926,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.1481481481481481,0.2962962962962963,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.1185185185185185,0.237037037037037,0.237037037037037,0.4740740740740741,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.1481481481481481,0.2962962962962963,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,-0.1333333333333333,0.1481481481481481,-0.1333333333333333,0.1259259259259259,-0.237037037037037,-0.1925925925925926,-0.1925925925925926,-0.1333333333333333,0.2074074074074074,-0.237037037037037,0.162962962962963,-0.237037037037037,-0.1925925925925926,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.1481481481481481,0.1481481481481481,-0.1333333333333333,0.1185185185185185,-0.1333333333333333,0.2962962962962963,0.2962962962962963,0.1481481481481481,0.1481481481481481,-0.237037037037037,0.4740740740740741,-0.237037037037037,0.2962962962962963,0.2962962962962963,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1333333333333333,-0.1925925925925926,0.1259259259259259,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,-0.237037037037037,-0.1333333333333333,-0.1925925925925926,0.162962962962963,-0.237037037037037,0.2074074074074074,-0.1925925925925926,-0.237037037037037,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.237037037037037,0.1185185185185185,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.2962962962962963,0.1481481481481481,0.2962962962962963,0.1481481481481481,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.4740740740740741,0.237037037037037,-0.237037037037037,0.2962962962962963,-0.1925925925925926,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.1185185185185185,0.237037037037037,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.1481481481481481,0.2962962962962963,0.1481481481481481,0.2962962962962963,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.237037037037037,0.4740740740740741,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.1481481481481481,-0.1333333333333333,0.1259259259259259,-0.1925925925925926,-0.1333333333333333,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,-0.1333333333333333,-0.237037037037037,-0.1925925925925926,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.237037037037037,-0.1925925925925926,0.2074074074074074,-0.237037037037037,0.162962962962963,-0.1333333333333333,0.1185185185185185,-0.1333333333333333,0.1481481481481481,0.1481481481481481,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.1481481481481481,0.1481481481481481,0.2962962962962963,0.2962962962962963,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.2962962962962963,0.2962962962962963,-0.237037037037037,0.4740740740740741,-0.237037037037037,0.1259259259259259,-0.1333333333333333,0.1481481481481481,-0.1333333333333333,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1333333333333333,-0.1925925925925926,-0.1925925925925926,-0.237037037037037,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1925925925925926,-0.237037037037037,0.162962962962963,-0.237037037037037,0.2074074074074074;
    }
  }
}

void FemParPoisson::getNodToModMatrix(Eigen::MatrixXd& mat, int ndim, int polyOrder)
{
  if(ndim==1) {
    if(polyOrder==1) {
      mat << 0.7071067811865475,0.7071067811865475,-0.408248290463863,0.408248290463863;
    }
    else if (polyOrder==2) {
      mat << 0.2357022603955159,0.9428090415820636,0.2357022603955159,-0.4082482904638631,0.0,0.4082482904638631,0.210818510677892,-0.421637021355784,0.210818510677892;
    }
  }
  else if(ndim==2) {
    if(polyOrder==1) {
      mat << 0.5,0.5,0.5,0.5,-0.2886751345948129,-0.2886751345948129,0.2886751345948129,0.2886751345948129,-0.2886751345948129,0.2886751345948129,-0.2886751345948129,0.2886751345948129,0.1666666666666667,-0.1666666666666667,-0.1666666666666667,0.1666666666666667;
    }
    else if (polyOrder==2) {
      mat << -0.1666666666666667,0.6666666666666666,-0.1666666666666666,0.6666666666666666,0.6666666666666666,-0.1666666666666667,0.6666666666666666,-0.1666666666666667,-0.09622504486493762,-0.3849001794597505,-0.09622504486493762,4.782436703424177e-17,3.326912489338558e-17,0.09622504486493759,0.3849001794597505,0.0962250448649376,-0.09622504486493762,4.158640611673197e-18,0.09622504486493762,-0.3849001794597504,0.3849001794597505,-0.09622504486493762,-4.158640611673197e-18,0.09622504486493762,0.1666666666666667,0.0,-0.1666666666666667,8.317281223346394e-18,-8.317281223346394e-18,-0.1666666666666667,0.0,0.1666666666666667,0.149071198499986,2.079320305836599e-17,0.149071198499986,-0.2981423969999719,-0.298142396999972,0.149071198499986,-3.425906660311779e-17,0.149071198499986,0.149071198499986,-0.2981423969999719,0.149071198499986,1.975354290544769e-17,-1.975354290544769e-17,0.149071198499986,-0.298142396999972,0.149071198499986,-0.08606629658238704,-2.079320305836599e-17,0.08606629658238706,0.1721325931647741,-0.1721325931647741,-0.08606629658238704,0.0,0.08606629658238706,-0.08606629658238704,0.1721325931647741,-0.08606629658238706,-2.079320305836599e-17,6.237960917509796e-18,0.08606629658238706,-0.1721325931647741,0.08606629658238706;
    }
  }
  else if (ndim==3) {
    if(polyOrder==1) {
      mat << 0.3535533905932737,0.3535533905932737,0.3535533905932737,0.3535533905932737,0.3535533905932738,0.3535533905932738,0.3535533905932738,0.3535533905932738,-0.2041241452319315,-0.2041241452319315,0.2041241452319315,0.2041241452319315,-0.2041241452319316,-0.2041241452319316,0.2041241452319316,0.2041241452319316,-0.2041241452319315,-0.2041241452319315,-0.2041241452319315,-0.2041241452319315,0.2041241452319316,0.2041241452319316,0.2041241452319316,0.2041241452319315,-0.2041241452319315,0.2041241452319315,-0.2041241452319315,0.2041241452319315,-0.2041241452319316,0.2041241452319316,-0.2041241452319316,0.2041241452319316,0.1178511301977579,0.1178511301977579,-0.1178511301977579,-0.1178511301977579,-0.117851130197758,-0.117851130197758,0.117851130197758,0.1178511301977579,0.1178511301977579,-0.1178511301977579,-0.1178511301977579,0.1178511301977579,0.117851130197758,-0.117851130197758,-0.117851130197758,0.117851130197758,0.1178511301977579,-0.1178511301977579,0.1178511301977579,-0.1178511301977579,-0.117851130197758,0.117851130197758,-0.117851130197758,0.117851130197758,-0.06804138174397718,0.06804138174397718,0.06804138174397718,-0.06804138174397717,0.06804138174397718,-0.0680413817439772,-0.0680413817439772,0.0680413817439772;
    }
    else if (polyOrder==2) {
      mat << -0.3535533905932737,0.4714045207910317,-0.3535533905932737,0.4714045207910317,0.4714045207910317,-0.3535533905932737,0.4714045207910317,-0.3535533905932737,0.4714045207910317,0.4714045207910317,0.4714045207910317,0.4714045207910317,-0.3535533905932737,0.4714045207910317,-0.3535533905932737,0.4714045207910317,0.4714045207910317,-0.3535533905932737,0.4714045207910317,-0.3535533905932737,0.06804138174397717,-0.2721655269759087,0.06804138174397717,0.0,0.0,-0.06804138174397717,0.2721655269759087,-0.06804138174397717,-0.2721655269759087,-0.2721655269759087,0.2721655269759087,0.2721655269759087,0.06804138174397717,-0.2721655269759087,0.06804138174397717,0.0,0.0,-0.06804138174397717,0.2721655269759087,-0.06804138174397717,0.06804138174397717,-0.2721655269759087,0.06804138174397717,-0.2721655269759087,-0.2721655269759087,0.06804138174397717,-0.2721655269759087,0.06804138174397717,0.0,0.0,0.0,0.0,-0.06804138174397717,0.2721655269759087,-0.06804138174397717,0.2721655269759087,0.2721655269759087,-0.06804138174397717,0.2721655269759087,-0.06804138174397717,0.06804138174397717,0.0,-0.06804138174397717,-0.2721655269759087,0.2721655269759087,0.06804138174397717,0.0,-0.06804138174397717,-0.2721655269759087,0.2721655269759087,-0.2721655269759087,0.2721655269759087,0.06804138174397717,0.0,-0.06804138174397717,-0.2721655269759087,0.2721655269759087,0.06804138174397717,0.0,-0.06804138174397717,0.0392837100659193,0.1571348402636772,0.0392837100659193,0.0,0.0,-0.0392837100659193,-0.1571348402636772,-0.0392837100659193,0.0,0.0,0.0,0.0,-0.0392837100659193,-0.1571348402636772,-0.0392837100659193,0.0,0.0,0.0392837100659193,0.1571348402636772,0.0392837100659193,0.0392837100659193,0.0,-0.0392837100659193,0.0,0.0,-0.0392837100659193,0.0,0.0392837100659193,0.1571348402636772,-0.1571348402636772,-0.1571348402636772,0.1571348402636772,0.0392837100659193,0.0,-0.0392837100659193,0.0,0.0,-0.0392837100659193,0.0,0.0392837100659193,0.0392837100659193,0.0,-0.0392837100659193,0.1571348402636772,-0.1571348402636772,0.0392837100659193,0.0,-0.0392837100659193,0.0,0.0,0.0,0.0,-0.0392837100659193,0.0,0.0392837100659193,-0.1571348402636772,0.1571348402636772,-0.0392837100659193,0.0,0.0392837100659193,0.105409255338946,0.0,0.105409255338946,-0.210818510677892,-0.210818510677892,0.105409255338946,0.0,0.105409255338946,0.0,0.0,0.0,0.0,0.105409255338946,0.0,0.105409255338946,-0.210818510677892,-0.210818510677892,0.105409255338946,0.0,0.105409255338946,0.105409255338946,0.0,0.105409255338946,0.0,0.0,0.105409255338946,0.0,0.105409255338946,-0.210818510677892,-0.210818510677892,-0.210818510677892,-0.210818510677892,0.105409255338946,0.0,0.105409255338946,0.0,0.0,0.105409255338946,0.0,0.105409255338946,0.105409255338946,-0.210818510677892,0.105409255338946,0.0,0.0,0.105409255338946,-0.210818510677892,0.105409255338946,0.0,0.0,0.0,0.0,0.105409255338946,-0.210818510677892,0.105409255338946,0.0,0.0,0.105409255338946,-0.210818510677892,0.105409255338946,-0.06804138174397717,0.0,0.06804138174397717,0.0,0.0,0.06804138174397717,0.0,-0.06804138174397717,0.0,0.0,0.0,0.0,0.06804138174397717,0.0,-0.06804138174397717,0.0,0.0,-0.06804138174397717,0.0,0.06804138174397717,-0.06085806194501844,0.0,-0.06085806194501844,0.1217161238900369,0.1217161238900369,-0.06085806194501844,0.0,-0.06085806194501844,0.0,0.0,0.0,0.0,0.06085806194501844,0.0,0.06085806194501844,-0.1217161238900369,-0.1217161238900369,0.06085806194501844,0.0,0.06085806194501844,-0.06085806194501844,0.0,-0.06085806194501844,0.0,0.0,0.06085806194501844,0.0,0.06085806194501844,0.1217161238900369,0.1217161238900369,-0.1217161238900369,-0.1217161238900369,-0.06085806194501844,0.0,-0.06085806194501844,0.0,0.0,0.06085806194501844,0.0,0.06085806194501844,-0.06085806194501844,0.0,0.06085806194501844,0.1217161238900369,-0.1217161238900369,-0.06085806194501844,0.0,0.06085806194501844,0.0,0.0,0.0,0.0,-0.06085806194501844,0.0,0.06085806194501844,0.1217161238900369,-0.1217161238900369,-0.06085806194501844,0.0,0.06085806194501844,-0.06085806194501844,0.0,0.06085806194501844,0.0,0.0,-0.06085806194501844,0.0,0.06085806194501844,0.1217161238900369,-0.1217161238900369,0.1217161238900369,-0.1217161238900369,-0.06085806194501844,0.0,0.06085806194501844,0.0,0.0,-0.06085806194501844,0.0,0.06085806194501844,-0.06085806194501844,0.1217161238900369,-0.06085806194501844,0.0,0.0,0.06085806194501844,-0.1217161238900369,0.06085806194501844,0.0,0.0,0.0,0.0,-0.06085806194501844,0.1217161238900369,-0.06085806194501844,0.0,0.0,0.06085806194501844,-0.1217161238900369,0.06085806194501844,-0.06085806194501844,0.1217161238900369,-0.06085806194501844,0.0,0.0,-0.06085806194501844,0.1217161238900369,-0.06085806194501844,0.0,0.0,0.0,0.0,0.06085806194501844,-0.1217161238900369,0.06085806194501844,0.0,0.0,0.06085806194501844,-0.1217161238900369,0.06085806194501844,0.03513641844631532,0.0,-0.03513641844631532,-0.07027283689263064,0.07027283689263064,0.03513641844631532,0.0,-0.03513641844631532,0.0,0.0,0.0,0.0,-0.03513641844631532,0.0,0.03513641844631532,0.07027283689263064,-0.07027283689263064,-0.03513641844631532,0.0,0.03513641844631532,0.03513641844631532,0.0,-0.03513641844631532,0.0,0.0,-0.03513641844631532,0.0,0.03513641844631532,-0.07027283689263064,0.07027283689263064,0.07027283689263064,-0.07027283689263064,0.03513641844631532,0.0,-0.03513641844631532,0.0,0.0,-0.03513641844631532,0.0,0.03513641844631532,0.03513641844631532,-0.07027283689263064,0.03513641844631532,0.0,0.0,-0.03513641844631532,0.07027283689263064,-0.03513641844631532,0.0,0.0,0.0,0.0,-0.03513641844631532,0.07027283689263064,-0.03513641844631532,0.0,0.0,0.03513641844631532,-0.07027283689263064,0.03513641844631532;
    }
  }
}

// C wrappers for interfacing with FemParPoisson class
extern "C" void* new_FemParPoisson(int nz, int ndim, int polyOrder, double dz, bool z_periodic, bcdataPar_t bc[2], bool writeMatrix, double laplacianWeight, double modifierConstant)
{
  FemParPoisson *f = new FemParPoisson(nz, ndim, polyOrder, dz, z_periodic, bc, writeMatrix, laplacianWeight, modifierConstant);
  return reinterpret_cast<void*>(f);
}

extern "C" void delete_FemParPoisson(FemParPoisson* f)
{
  delete f;
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

