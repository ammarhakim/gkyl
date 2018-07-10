// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for 2D FEM Poisson solver
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <FemPerpPoissonImpl.h>

// std includes
#include <string>
#include <vector>

using namespace Eigen;
static const int DIRICHLET_BC = 0;
static const int NEUMANN_BC = 1;
static const int DIRICHLET_VARIABLE_BC = 2;
static const int DX = 0;
static const int DY = 1;
static const int LO = 0;
static const int HI = 1;

double take_last(const double &in,const double &b) { return b; }
void vectorSum(double *in, double *inout, int *len, MPI_Datatype *dptr)
{
  int i;
  for(i=0; i< *len; ++i) {
    inout[i] += in[i];
  }
}

FemPerpPoisson::FemPerpPoisson(int nx_, int ny_, int ndim_, int polyOrder_, 
                       double dx_, double dy_, bool periodicFlgs_[2], 
                       bcdata_t bc_[2][2], bool writeMatrix_,
                       double laplacianWeight_, double modifierConstant_) 
  : nx(nx_), ny(ny_), ndim(ndim_), 
    polyOrder(polyOrder_), dx(dx_), dy(dy_), 
    writeMatrix(writeMatrix_),
    laplacianWeight(laplacianWeight_),
    modifierConstant(modifierConstant_)
{
//#define DMSG(s) std::cout << s << std::endl;
#define DMSG(s) ;

  DMSG("Inside FEM Poisson initialize");

// copy to input to class structures
  for(int i=0; i<2; i++) {
    periodicFlgs[i] = periodicFlgs_[i];
    for(int j=0; j<2; j++) {
      bc[i][j] = bc_[i][j];
      if(!bc[i][j].isSet) {
        bc[i][j].type = -1;
      }
      bc2d[i][j] = bc[i][j];
      bc2d_z0[i][j] = bc[i][j];
    }
  }

  allPeriodic = false;
  if(periodicFlgs[0] && periodicFlgs[1]) allPeriodic = true;
  adjustSource = false;
  if(allPeriodic && modifierConstant==0.0) adjustSource = true;
  cornerval = 0.0;

// prepare solver
  int nglobal = getNumPerpGlobalNodes(nx, ny, ndim, polyOrder, periodicFlgs);
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  int nlocal_z = 1;
  if(ndim==3) nlocal_z += polyOrder;

  // setup boundary indices and make global stiffness matrix(es)
  setupBoundaryIndices(bc, ndim, polyOrder);
  if(ndim==3 && polyOrder==1) {
    // trick: use a 3d p=1 global perp stiffness matrix that has a diagonal 
    // block structure with two blocks (one for each of z=-1 and z=1).
    // each block is just a 2d p=1 global perp stiffness matrix.
    // so we make a single 2d p=1 global perp stiffness matrix, and we will use it twice.
    // the 3d p=1 global mass matrix will have the same structure.
    // this trick is effectively a pre-factorization for the solver.
    setupBoundaryIndices(bc2d, 2, polyOrder);
    makeGlobalPerpStiffnessMatrix(stiffMat, sourceModVec, 2, polyOrder, bc2d);
  } 
  else {
    // make global stiffness matrix
    makeGlobalPerpStiffnessMatrix(stiffMat, sourceModVec, ndim, polyOrder, bc);
  }
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
    outName = "poisson-modtonod"; 
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

  MPI_Op_create((MPI_User_function *) vectorSum, true, &MPI_vectorSum_op);
}


void FemPerpPoisson::allreduceGlobalSrc(MPI_Comm comm)
{
  int nglobal = getNumPerpGlobalNodes(nx, ny, ndim, polyOrder, periodicFlgs);
  MPI_Allreduce(MPI_IN_PLACE, globalSrc.data(), nglobal, MPI_DOUBLE, MPI_vectorSum_op, comm);
}

FemPerpPoisson::~FemPerpPoisson() 
{
  stiffMat.resize(0,0);
  stiffMat_z0.resize(0,0);
  globalSrc.resize(0);
  sourceModVec.resize(0);
  sourceModVec_z0.resize(0);
  x.resize(0);
}

// get start and end indices for each boundary
// note: corners not included!
void FemPerpPoisson::setupBoundaryIndices(bcdata_t bc[2][2], int ndim, int polyOrder)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  std::vector<int> lgMap(nlocal);
  // left boundary
    // get index of start and end of left boundary in global mapping
    // start
    getPerpLocalToGlobalInteriorBLRT(lgMap,0,0,nx,ny,ndim,polyOrder,periodicFlgs);
    if(polyOrder==1) {
      bc[DX][LO].cornerstart[0] = lgMap[0];
      bc[DX][LO].istart[0] = lgMap[2];
      if(ndim==3) {  
        bc[DX][LO].cornerstart[1] = lgMap[4];
        bc[DX][LO].istart[1] = lgMap[6];
      }
    }
    else if(polyOrder==2) {
      bc[DX][LO].cornerstart[0] = lgMap[0];
      bc[DX][LO].istart[0] = lgMap[3];
      if(ndim==3) {
        bc[DX][LO].cornerstart[1] = lgMap[8];
        bc[DX][LO].cornerstart[2] = lgMap[12];
        bc[DX][LO].istart[1] = lgMap[10];
        bc[DX][LO].istart[2] = lgMap[15];
      }
    }
    // end
    getPerpLocalToGlobalInteriorBLRT(lgMap,0,ny-1,nx,ny,ndim,polyOrder,periodicFlgs);
    if(polyOrder==1) {
      bc[DX][LO].cornerend[0] = lgMap[2];
      bc[DX][LO].iend[0] = lgMap[0];
      if(ndim==3) {  
        bc[DX][LO].cornerend[1] = lgMap[6];
        bc[DX][LO].iend[1] = lgMap[4];
      }
    }
    else if(polyOrder==2) {
      bc[DX][LO].cornerend[0] = lgMap[5];
      bc[DX][LO].iend[0] = lgMap[3];
      if(ndim==3) {
        bc[DX][LO].cornerend[1] = lgMap[10];
        bc[DX][LO].cornerend[2] = lgMap[17];
        bc[DX][LO].iend[1] = lgMap[8];
        bc[DX][LO].iend[2] = lgMap[15];
      }
    }
  // right boundary
    // get index of start and end of right boundary in global mapping
    // start
    getPerpLocalToGlobalInteriorBLRT(lgMap,nx,0,nx,ny,ndim,polyOrder,periodicFlgs);
    if(polyOrder==1) {
      bc[DX][HI].cornerstart[0] = lgMap[0];
      bc[DX][HI].istart[0] = lgMap[2];
      if(ndim==3) {
        bc[DX][HI].cornerstart[1] = lgMap[4];
        bc[DX][HI].istart[1] = lgMap[6];
      }
    }
    else if(polyOrder==2) {
      bc[DX][HI].cornerstart[0] = lgMap[0];
      bc[DX][HI].istart[0] = lgMap[3];
      if(ndim==3) {
        bc[DX][HI].cornerstart[1] = lgMap[8];
        bc[DX][HI].cornerstart[2] = lgMap[12];
        bc[DX][HI].istart[1] = lgMap[10];
        bc[DX][HI].istart[2] = lgMap[15];
      }
    }
    // end
    getPerpLocalToGlobalInteriorBLRT(lgMap,nx,ny-1,nx,ny,ndim,polyOrder,periodicFlgs);
    if(polyOrder==1) {
      bc[DX][HI].cornerend[0] = lgMap[2];
      bc[DX][HI].iend[0] = lgMap[0];
      if(ndim==3) {
        bc[DX][HI].cornerend[1] = lgMap[6];
        bc[DX][HI].iend[1] = lgMap[4];
      }
    }
    else if(polyOrder==2) {
      bc[DX][HI].cornerend[0] = lgMap[5];
      bc[DX][HI].iend[0] = lgMap[3];
      if(ndim==3) {
        bc[DX][HI].cornerend[1] = lgMap[10];
        bc[DX][HI].cornerend[2] = lgMap[17];
        bc[DX][HI].iend[1] = lgMap[8];
        bc[DX][HI].iend[2] = lgMap[15];
      }
    }
  // bottom boundary
    // get index of start and end of bottom boundary in global mapping
    // start
    getPerpLocalToGlobalInteriorBLRT(lgMap,0,0,nx,ny,ndim,polyOrder,periodicFlgs);
    if(polyOrder==1) {
      bc[DY][LO].cornerstart[0] = lgMap[0];
      bc[DY][LO].istart[0] = lgMap[1];
      if(ndim==3) { 
        bc[DY][LO].cornerstart[1] = lgMap[4];
        bc[DY][LO].istart[1] = lgMap[5];
      }
    }
    else if(polyOrder==2) {
      bc[DY][LO].cornerstart[0] = lgMap[0];
      bc[DY][LO].istart[0] = lgMap[1];
      if(ndim==3) {
        bc[DY][LO].cornerstart[1] = lgMap[8];
        bc[DY][LO].cornerstart[2] = lgMap[12];
        bc[DY][LO].istart[1] = lgMap[9];
        bc[DY][LO].istart[2] = lgMap[13];
      }
    }
    // end
    getPerpLocalToGlobalInteriorBLRT(lgMap,nx-1,0,nx,ny,ndim,polyOrder,periodicFlgs);
    if(polyOrder==1) {
      bc[DY][LO].cornerend[0] = lgMap[1];
      bc[DY][LO].iend[0] = lgMap[0];
      if(ndim==3) {
        bc[DY][LO].cornerend[1] = lgMap[5];
        bc[DY][LO].iend[1] = lgMap[4];
      }
    }
    else if(polyOrder==2) {
      bc[DY][LO].cornerend[0] = lgMap[2];
      bc[DY][LO].iend[0] = lgMap[1];
      if(ndim==3) {
        bc[DY][LO].cornerend[1] = lgMap[9];
        bc[DY][LO].cornerend[2] = lgMap[14];
        bc[DY][LO].iend[1] = lgMap[8];
        bc[DY][LO].iend[2] = lgMap[13];
      }
    }
  // top boundary
    // get index of start and end of top boundary in global mapping
    // start
    getPerpLocalToGlobalInteriorBLRT(lgMap,0,ny,nx,ny,ndim,polyOrder,periodicFlgs);
    if(polyOrder==1) {
      bc[DY][HI].cornerstart[0] = lgMap[0];
      bc[DY][HI].istart[0] = lgMap[1];
      if(ndim==3) {  
        bc[DY][HI].cornerstart[1] = lgMap[4];
        bc[DY][HI].istart[1] = lgMap[5];
      }
    }
    else if(polyOrder==2) {
      bc[DY][HI].cornerstart[0] = lgMap[0];
      bc[DY][HI].istart[0] = lgMap[1];
      if(ndim==3) {
        bc[DY][HI].cornerstart[1] = lgMap[8];
        bc[DY][HI].cornerstart[2] = lgMap[12];
        bc[DY][HI].istart[1] = lgMap[9];
        bc[DY][HI].istart[2] = lgMap[13];
      }
    }
    // end
    getPerpLocalToGlobalInteriorBLRT(lgMap,nx-1,ny,nx,ny,ndim,polyOrder,periodicFlgs);
    if(polyOrder==1) {
      bc[DY][HI].cornerend[0] = lgMap[1];
      bc[DY][HI].iend[0] = lgMap[0];
      if(ndim==3) {
        bc[DY][HI].cornerend[1] = lgMap[5];
        bc[DY][HI].iend[1] = lgMap[4];
      }
    }
    else if(polyOrder==2) {
      bc[DY][HI].cornerend[0] = lgMap[2];
      bc[DY][HI].iend[0] = lgMap[1];
      if(ndim==3) {
        bc[DY][HI].cornerend[1] = lgMap[9];
        bc[DY][HI].cornerend[2] = lgMap[14];
        bc[DY][HI].iend[1] = lgMap[8];
        bc[DY][HI].iend[2] = lgMap[13];
      }
    }

  // check
//  for (int d=0; d<2; ++d)
//  {
//    for (int side=0; side<2; ++side)
//    {
//      printf("Direction %d, side %d = %d : %d", d, side, bc[d][side].istart[0], bc[d][side].iend[0]);
//      if(ndim==3) {
//        for(int i=0; i<polyOrder; i++) {
//          printf(", %d : %d", bc[d][side].istart[i+1], bc[d][side].iend[i+1]);
//        }
//      }
//      printf("\n");
//    }
//  }
//  printf("Corners: \n");
//  printf("DX: Bottom left: %d    ", bc[DX][LO].cornerstart[0]);
//  printf("Bottom right: %d    ", bc[DX][HI].cornerstart[0]);
//  printf("Top left: %d    ", bc[DX][LO].cornerend[0]);
//  printf("Top right %d\n", bc[DX][HI].cornerend[0]);
//  printf("DY: Bottom left: %d    ", bc[DY][LO].cornerstart[0]);
//  printf("Bottom right: %d    ", bc[DY][LO].cornerend[0]);
//  printf("Top left: %d    ", bc[DY][HI].cornerstart[0]);
//  printf("Top right %d\n", bc[DY][HI].cornerend[0]);
    
}

void FemPerpPoisson::makeGlobalPerpStiffnessMatrix(
     Eigen::SparseMatrix<double,Eigen::ColMajor>& stiffMat, 
     Eigen::VectorXd& sourceModVec_,
     int ndim_in, int polyOrder_in, bcdata_t bc_in[2][2])
{
  int nglobal = getNumPerpGlobalNodes(nx, ny, ndim_in, polyOrder_in, periodicFlgs);
  int nlocal = getNumLocalNodes(ndim_in, polyOrder_in);
  int nlocal_z = 1;
  if(ndim_in==3) nlocal_z += polyOrder_in;

  int nonzeros = nlocal*(std::pow(2.0, 1.0*2)+1); // estimate number of nonzeros per row

  MatrixXd localStiff = MatrixXd::Zero(nlocal, nlocal);
  MatrixXd localMass = MatrixXd::Zero(nlocal, nlocal);

  std::vector<int> lgMap(nlocal);
  std::vector<Triplet<double> > tripletList, identityTripletList;
  tripletList.reserve(nonzeros*nglobal); // estimate number of nonzero elements

  getPerpStiffnessMatrix(localStiff, ndim_in, polyOrder_in, dx, dy);
  getMassMatrix(localMass, ndim_in, polyOrder_in);

  // loop over global region
  for(int idy=0; idy<ny; idy++) {
    for(int idx=0; idx<nx; idx++) {
      getPerpLocalToGlobalInteriorBLRT(lgMap,idx,idy,nx,ny,ndim_in,polyOrder_in,periodicFlgs);
      
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
  }

// construct SparseMatrix from triplets
  SparseMatrix<double,RowMajor> stiffMatRowMajor(nglobal, nglobal);
  stiffMatRowMajor.setFromTriplets(tripletList.begin(), tripletList.end());

// handle Dirichlet BCs
// note: we do not need to do anything for Neumann BCs
  int nDirichlet = 0;
// set rows corresponding to Dirichlet BCs to zero
  for (int d=0; d<2; ++d)
  {
    for (int side=0; side<2; ++side)
    {
      if ((bc_in[d][side].isSet && bc_in[d][side].type == DIRICHLET_BC) ||
          (bc_in[d][side].isSet && bc_in[d][side].type == DIRICHLET_VARIABLE_BC))
      {
        for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
          int start = bc_in[d][side].istart[zblock];
          int end = bc_in[d][side].iend[zblock];
          int nRows = end-start+1;
          stiffMatRowMajor.middleRows(start, nRows) = SparseMatrix<double,RowMajor>(nRows, stiffMatRowMajor.cols());   

          nDirichlet += nRows;
        }
      }
    }
  }

// handle Dirichlet corners
  // bottom left
  // also make this corner effectively Dirichlet if all periodic BCs
  if((bc_in[DX][LO].type==DIRICHLET_BC || bc_in[DY][LO].type==DIRICHLET_BC) || adjustSource) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][LO].cornerstart[zblock];
      stiffMatRowMajor.middleRows(corner, 1) = SparseMatrix<double,RowMajor>(1, stiffMatRowMajor.cols());
      nDirichlet += 1;
    }
  }
  // bottom right
  if((bc_in[DX][HI].type==DIRICHLET_BC || bc_in[DY][LO].type==DIRICHLET_BC)) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][LO].cornerend[zblock];
      stiffMatRowMajor.middleRows(corner, 1) = SparseMatrix<double,RowMajor>(1, stiffMatRowMajor.cols());
      nDirichlet += 1;
    }
  }
  // top left
  if(bc_in[DX][LO].type==DIRICHLET_BC || bc_in[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][HI].cornerstart[zblock];
      stiffMatRowMajor.middleRows(corner, 1) = SparseMatrix<double,RowMajor>(1, stiffMatRowMajor.cols());
      nDirichlet += 1;
    }
  }
  // top right
  if(bc_in[DX][HI].type==DIRICHLET_BC || bc_in[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][HI].cornerend[zblock];
      stiffMatRowMajor.middleRows(corner, 1) = SparseMatrix<double,RowMajor>(1, stiffMatRowMajor.cols());
      nDirichlet += 1;
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

  for (int d=0; d<2; ++d)
  {
    for (int side=0; side<2; ++side)
    {
      if ((bc_in[d][side].isSet && bc_in[d][side].type == DIRICHLET_BC) ||
          (bc_in[d][side].isSet && bc_in[d][side].type == DIRICHLET_VARIABLE_BC))
      {
        for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
          int start = bc_in[d][side].istart[zblock];
          int end = bc_in[d][side].iend[zblock];
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
    }
  }

// handle Dirichlet corners
  // bottom left
  // also make this corner effectively Dirichlet if all periodic BCs
  if((bc_in[DX][LO].type==DIRICHLET_BC || bc_in[DY][LO].type==DIRICHLET_BC) || adjustSource) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][LO].cornerstart[zblock];
// copy cols corresponding to Dirichlet BCs to sourceModMat
      sourceModMat.middleCols(corner, 1) = stiffMat.middleCols(corner, 1);
// in stiffMat, set cols corresponding to Dirichlet BCs to zero
      stiffMat.middleCols(corner, 1) = SparseMatrix<double,ColMajor>(stiffMat.rows(), 1);
// set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
      identityTripletList.push_back(Triplet<double>(corner, corner, 1.));
    }
  }
  // bottom right
  if((bc_in[DX][HI].type==DIRICHLET_BC || bc_in[DY][LO].type==DIRICHLET_BC)) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][LO].cornerend[zblock];
// copy cols corresponding to Dirichlet BCs to sourceModMat
      sourceModMat.middleCols(corner, 1) = stiffMat.middleCols(corner, 1);
// in stiffMat, set cols corresponding to Dirichlet BCs to zero
      stiffMat.middleCols(corner, 1) = SparseMatrix<double,ColMajor>(stiffMat.rows(), 1);
// set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
      identityTripletList.push_back(Triplet<double>(corner, corner, 1.));
    }
  }
  // top left
  if(bc_in[DX][LO].type==DIRICHLET_BC || bc_in[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][HI].cornerstart[zblock];
// copy cols corresponding to Dirichlet BCs to sourceModMat
      sourceModMat.middleCols(corner, 1) = stiffMat.middleCols(corner, 1);
// in stiffMat, set cols corresponding to Dirichlet BCs to zero
      stiffMat.middleCols(corner, 1) = SparseMatrix<double,ColMajor>(stiffMat.rows(), 1);
// set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
      identityTripletList.push_back(Triplet<double>(corner, corner, 1.));
    }
  }
  // top right
  if(bc_in[DX][HI].type==DIRICHLET_BC || bc_in[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][HI].cornerend[zblock];
// copy cols corresponding to Dirichlet BCs to sourceModMat
      sourceModMat.middleCols(corner, 1) = stiffMat.middleCols(corner, 1);
// in stiffMat, set cols corresponding to Dirichlet BCs to zero
      stiffMat.middleCols(corner, 1) = SparseMatrix<double,ColMajor>(stiffMat.rows(), 1);
// set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
      identityTripletList.push_back(Triplet<double>(corner, corner, 1.));
    }
  }

  // this ensures no duplicates, so that entries for these rows are at most equal to 1
  dirichletIdentity.setFromTriplets(identityTripletList.begin(), identityTripletList.end(), take_last);

  if (modifierConstant!=0.0 && laplacianWeight==0.0) {
    stiffMat+=modifierConstant*dirichletIdentity;
  } else {
    stiffMat+=dirichletIdentity;
  }

// create vector of Dirichlet values
  SparseVector<double> dirichletVec(nglobal);
  dirichletVec.reserve(nDirichlet);
  for (int d=0; d<2; ++d)
  {
    for (int side=0; side<2; ++side)
    {
      if (bc_in[d][side].isSet && bc_in[d][side].type == DIRICHLET_BC)
      {
        for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
          for(unsigned i=bc_in[d][side].istart[zblock]; i<=bc_in[d][side].iend[zblock]; i++) {
            dirichletVec.coeffRef(i) = bc_in[d][side].value;
          }
        }
      }
    }
  }

// handle Dirichlet corners
  // bottom left
  // also make this corner effectively Dirichlet if all periodic BCs
  if((bc_in[DX][LO].type==DIRICHLET_BC || bc_in[DY][LO].type==DIRICHLET_BC) || adjustSource) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][LO].cornerstart[zblock];
      dirichletVec.coeffRef(corner) = bc_in[DX][LO].type==DIRICHLET_BC ? bc_in[DX][LO].value : bc_in[DY][LO].value;
      if(adjustSource) dirichletVec.coeffRef(corner) = cornerval;
    }
  }
  // bottom right
  if((bc_in[DX][HI].type==DIRICHLET_BC || bc_in[DY][LO].type==DIRICHLET_BC)) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][LO].cornerend[zblock];
      dirichletVec.coeffRef(corner) = bc_in[DX][HI].type==DIRICHLET_BC ? bc_in[DX][HI].value : bc_in[DY][LO].value;
    }
  }
  // top left
  if(bc_in[DX][LO].type==DIRICHLET_BC || bc_in[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][HI].cornerstart[zblock];
      dirichletVec.coeffRef(corner) = bc_in[DX][LO].type==DIRICHLET_BC ? bc_in[DX][LO].value : bc_in[DY][HI].value;
    }
  }
  // top right
  if(bc_in[DX][HI].type==DIRICHLET_BC || bc_in[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc_in[DY][HI].cornerend[zblock];
      dirichletVec.coeffRef(corner) = bc_in[DX][HI].type==DIRICHLET_BC ? bc_in[DX][HI].value : bc_in[DY][HI].value;
    }
  }

// calculate vector to subtract from source
  sourceModVec_ = sourceModMat*dirichletVec;
}


// clear out existing stuff in source vector: this is required
// otherwise successive calls to advance() will accumulate into source
// from prevous calls, which is of course not what we want.
void FemPerpPoisson::zeroGlobalSrc() {
  int nglobal = getNumPerpGlobalNodes(nx, ny, ndim, polyOrder, periodicFlgs);
  globalSrc.setZero(nglobal);
}

// called within an indexer loop over idx and idy
void FemPerpPoisson::createGlobalSrc(double* localSrcPtr, int idx, int idy, double intSrcVol)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);

  std::vector<int> lgMap(nlocal);
  std::vector<double> localSrc(nlocal), localMassSrc(nlocal);

  // copy data from input
  for (unsigned k=0; k<nlocal; ++k)
    localSrc[k] = localSrcPtr[k];

  // adjust src cell-average 
  localSrc[0] = localSrc[0] - intSrcVol;
  
  getPerpLocalToGlobalInteriorBLRT(lgMap,idx,idy,nx,ny,ndim,polyOrder,periodicFlgs);

  //printf("(%d,%d): ", idx,idy);
  //for(unsigned k=0; k<nlocal; ++k) 
  //{
  //  if(k<8)
  //    printf("%d  ", lgMap[k]);
  //  else if (k<12)
  //    printf("%d  ", lgMap[k]-getNumPerpGlobalNodes(nx,ny,2,polyOrder,periodicFlgs));
  //  else
  //    printf("%d  ", lgMap[k]-getNumPerpGlobalNodes(nx,ny,2,polyOrder,periodicFlgs)-getNumPerpGlobalNodes(nx,ny,2,1,periodicFlgs));
  //}
  //printf("\n");

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


void FemPerpPoisson::solve()
{
  int nglobal = getNumPerpGlobalNodes(nx, ny, ndim, polyOrder, periodicFlgs);
  int nlocal_z = 1;
  if(ndim==3) nlocal_z += polyOrder;
// replace dirichlet nodes of global source with dirichlet values
  for (unsigned d=0; d<2; ++d)
  {
    for (unsigned side=0; side<2; ++side)
    {
      if (bc[d][side].isSet && bc[d][side].type == DIRICHLET_BC)
      {
        for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
          //printf("dim %d side %d is dirichlet: start = %d, end = %d, value = %f\n", d, side, bc[d][side].istart[0], bc[d][side].iend[0], bc[d][side].value);
          for(unsigned i=bc[d][side].istart[zblock]; i<=bc[d][side].iend[zblock]; ++i) {
            globalSrc.coeffRef(i) = bc[d][side].value;
          }
        }
      }
    }
  }

// handle Dirichlet corners
// also make corners effectively Dirichlet if all periodic BCs
  // bottom left
  if((bc[DX][LO].type==DIRICHLET_BC || bc[DY][LO].type==DIRICHLET_BC) || adjustSource) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][LO].cornerstart[zblock];
      globalSrc.coeffRef(corner) = bc[DX][LO].type==DIRICHLET_BC ? bc[DX][LO].value : bc[DY][LO].value;
      if(adjustSource) globalSrc.coeffRef(corner) = cornerval;
    }
  }
  // bottom right
  if((bc[DX][HI].type==DIRICHLET_BC || bc[DY][LO].type==DIRICHLET_BC)) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][LO].cornerend[zblock];
      globalSrc.coeffRef(corner) = bc[DX][HI].type==DIRICHLET_BC ? bc[DX][HI].value : bc[DY][LO].value;
    }
  }
  // top left
  if(bc[DX][LO].type==DIRICHLET_BC || bc[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][HI].cornerstart[zblock];
      globalSrc.coeffRef(corner) = bc[DX][LO].type==DIRICHLET_BC ? bc[DX][LO].value : bc[DY][HI].value;
    }
  }
  // top right
  if(bc[DX][HI].type==DIRICHLET_BC || bc[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][HI].cornerend[zblock];
      globalSrc.coeffRef(corner) = bc[DX][HI].type==DIRICHLET_BC ? bc[DX][HI].value : bc[DY][HI].value;
    }
  }

  if (writeMatrix)
  {
    std::string outName = "poisson-src";
    outName += std::to_string(ndim) + "d";
    saveMarket(globalSrc, outName);
  }
 
  x = Eigen::VectorXd(nglobal);
  Eigen::VectorXd x_zm1, x_z0;
  
// solve linear system(s)
// modify non-dirichlet rows of source by subtracting sourceModVec
  if(ndim==3 && polyOrder==1) {
    // for 3d p=1, we solve a 2d system for each local z node and concatenate the results
    int nglobal2d = getNumPerpGlobalNodes(nx, ny, 2, polyOrder, periodicFlgs);
    // solve for z=-1 component
    x_zm1 = solver.solve(globalSrc.head(nglobal2d)-sourceModVec); 
    // solve for z=1 component
    x << x_zm1, solver.solve(globalSrc.tail(nglobal2d)-sourceModVec);
  }
  else {
    x = solver.solve(globalSrc-sourceModVec);
  }

  if (writeMatrix)
  {
    std::string outName = "poisson-sol";
    outName += std::to_string(ndim) + "d";
    saveMarket(x, outName);
  }
}

void FemPerpPoisson::getSolution(double* solPtr, int idx, int idy)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  std::vector<int> lgMap(nlocal);
  getPerpLocalToGlobalInteriorBLRT(lgMap, idx,idy,nx,ny,ndim,polyOrder, periodicFlgs);

  // transform global nodal x to local modal sol
  for (unsigned k=0; k<nlocal; ++k) {
    solPtr[k] = 0.0;
    for (unsigned m=0; m<nlocal; ++m) {
      solPtr[k] += localNodToMod(k,m)*x.coeffRef(lgMap[m]);
    }
  }
}

void FemPerpPoisson::getNodalSolution(double* solPtr, int idx, int idy)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  std::vector<int> lgMap(nlocal);
  getPerpLocalToGlobalInteriorBLRT(lgMap, idx,idy,nx,ny,ndim,polyOrder, periodicFlgs);

  // transform global nodal x to local nodal sol
  for (unsigned k=0; k<nlocal; ++k) {
    solPtr[k] = x.coeffRef(lgMap[k]);
  }
}

int FemPerpPoisson::getNumPerpGlobalNodes(int nx, int ny, int ndim, int polyOrder, bool periodicFlgs[2])
{
  int numGlobalNodes = -1;
  bool x_periodic = periodicFlgs[0];
  bool y_periodic = periodicFlgs[1];

  if(ndim==2) {
    if (polyOrder == 1) {
      numGlobalNodes = (nx+1)*(ny+1);
      // if x-periodic, subtract number of left boundary nodes
      if(x_periodic) numGlobalNodes -= (ny - 1);
      // if y-periodic, subtract number of top boundary nodes
      if(y_periodic) numGlobalNodes -= (nx + 1);
      // if x- and y-periodic, subtract an additional node for bottom right corner
      if(x_periodic && y_periodic) numGlobalNodes -= 1;
    }
    else if (polyOrder == 2)
    {
      // there are 3 owned nodes per cell, 2*nx and 2*ny edge nodes along
      // the top and right edges and the 1 accounts for the node on the
      // top-right corner.
      numGlobalNodes = 3*nx*ny + 2*nx + 2*ny + 1;
      // if x-periodic, subtract number of left boundary nodes
      if(x_periodic) numGlobalNodes -= (2*ny - 1);
      // if y-periodic, subtract number of top boundary nodes
      if(y_periodic) numGlobalNodes -= (2*nx + 1);
      // if x- and y-periodic, subtract an additional node for bottom right corner
      if(x_periodic && y_periodic) numGlobalNodes -= 1;
    }
  }
  else if (ndim==3) {
    if (polyOrder == 1) {
      // there are two identical z planes: z=-1 and z=1
      // each z plane has nglobal_2d nodes
      numGlobalNodes = 2*getNumPerpGlobalNodes(nx,ny,2,polyOrder,periodicFlgs);
    }
    else if (polyOrder ==2)
    {
      // there are three z planes: z=-1, z=0, and z=1

      // the z=-1 and z=1 planes are identical
      // each has nglobal_2d nodes with p=2 in the x-y direction
      numGlobalNodes = 2*getNumPerpGlobalNodes(nx,ny,2,polyOrder,periodicFlgs);

      // also need to add nglobal_2d for z=0 plane, 
      // which effectively has p=1 in x-y direction
      numGlobalNodes += getNumPerpGlobalNodes(nx,ny,2,1,periodicFlgs);
    }
  }
  return numGlobalNodes;
}

int FemPerpPoisson::getNumLocalNodes(int ndim, int p) 
{
  int numNodes = -1;
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

void FemPerpPoisson::getPerpLocalToGlobalInteriorBLRT(std::vector<int>& lgMap, int ix, int iy, int nx, int ny, int ndim, int polyOrder, bool periodicFlgs[2])
{
  int ninterior = (2*polyOrder-1)*nx*ny-polyOrder*(nx+ny)+1;

  if(ndim==2) {
    if (polyOrder == 1)
    {
      lgMap[0] = F1_func(nx,ny,ninterior,ix,iy,periodicFlgs); // 0
      lgMap[1] = F1_func(nx,ny,ninterior,ix+1,iy,periodicFlgs); // 1
      lgMap[2] = F1_func(nx,ny,ninterior,ix,iy+1,periodicFlgs); // 2
      lgMap[3] = F1_func(nx,ny,ninterior,ix+1,iy+1,periodicFlgs); // 3
    }
    else if (polyOrder == 2)
    {
      lgMap[0] = F2_func(nx,ny,ninterior,ix,iy,periodicFlgs); // 0
      lgMap[1] = ix==0 ? 
         F2_func(nx,ny,ninterior,ix+1,iy,periodicFlgs)-1 : 
         F2_func(nx,ny,ninterior,ix,iy,periodicFlgs)+1; // 1
      lgMap[2] = F2_func(nx,ny,ninterior,ix+1, iy,periodicFlgs); // 2
      lgMap[3] = G2_func(nx,ny,ninterior,ix, iy,periodicFlgs); // 3
      lgMap[4] = G2_func(nx,ny,ninterior,ix+1, iy,periodicFlgs); // 4
      lgMap[5] = F2_func(nx,ny,ninterior,ix, iy+1,periodicFlgs); // 5
      lgMap[6] = ix==0 ? 
         F2_func(nx,ny,ninterior,ix+1,iy+1,periodicFlgs)-1 : 
         F2_func(nx,ny,ninterior,ix,iy+1,periodicFlgs)+1; // 6
      lgMap[7] = F2_func(nx,ny,ninterior,ix+1,iy+1,periodicFlgs); // 7
    }
  }
  else if (ndim==3) {
    // note: this creates a global mapping in x and y, but not in z.
    // this should be used only for perpendicular poisson solvers (which do not require continuity in z).
    // the mapping has a block structure, with each block corresponding to one of the local z nodes.
    if (polyOrder == 1)
    {
      // z = -1 block
      lgMap[0] = F1_func(nx,ny,ninterior,ix,iy,periodicFlgs); // 0
      lgMap[1] = F1_func(nx,ny,ninterior,ix+1,iy,periodicFlgs); // 1
      lgMap[2] = F1_func(nx,ny,ninterior,ix,iy+1,periodicFlgs); // 2
      lgMap[3] = F1_func(nx,ny,ninterior,ix+1,iy+1,periodicFlgs); // 3
      // z = 1 block
      // offset = nglobal in 2d for p=1
      int offset = getNumPerpGlobalNodes(nx,ny,2,1,periodicFlgs);
      lgMap[4] = offset + F1_func(nx,ny,ninterior,ix,iy,periodicFlgs); // 0
      lgMap[5] = offset + F1_func(nx,ny,ninterior,ix+1,iy,periodicFlgs); // 1
      lgMap[6] = offset + F1_func(nx,ny,ninterior,ix,iy+1,periodicFlgs); // 2
      lgMap[7] = offset + F1_func(nx,ny,ninterior,ix+1,iy+1,periodicFlgs); // 3
    }
    else if (polyOrder == 2)
    {
      // z = -1 block
      lgMap[0] = F2_func(nx,ny,ninterior,ix,iy,periodicFlgs); // 0
      lgMap[1] = ix==0 ? 
         F2_func(nx,ny,ninterior,ix+1,iy,periodicFlgs)-1 : 
         F2_func(nx,ny,ninterior,ix,iy,periodicFlgs)+1; // 1
      lgMap[2] = F2_func(nx,ny,ninterior,ix+1, iy,periodicFlgs); // 2
      lgMap[3] = G2_func(nx,ny,ninterior,ix, iy,periodicFlgs); // 3
      lgMap[4] = G2_func(nx,ny,ninterior,ix+1, iy,periodicFlgs); // 4
      lgMap[5] = F2_func(nx,ny,ninterior,ix, iy+1,periodicFlgs); // 5
      lgMap[6] = ix==0 ? 
         F2_func(nx,ny,ninterior,ix+1,iy+1,periodicFlgs)-1 : 
         F2_func(nx,ny,ninterior,ix,iy+1,periodicFlgs)+1; // 6
      lgMap[7] = F2_func(nx,ny,ninterior,ix+1,iy+1,periodicFlgs); // 7
      // z = 0 block (effectively p=1)
      int ninterior_p1 = (nx-1)*(ny-1);
      // offset = nglobal in 2d for p=2
      int offset = getNumPerpGlobalNodes(nx,ny,2,2,periodicFlgs);
      lgMap[8] = offset + F1_func(nx,ny,ninterior_p1,ix,iy,periodicFlgs); // 0
      lgMap[9] = offset + F1_func(nx,ny,ninterior_p1,ix+1,iy,periodicFlgs); // 1
      lgMap[10] = offset + F1_func(nx,ny,ninterior_p1,ix,iy+1,periodicFlgs); // 2
      lgMap[11] = offset + F1_func(nx,ny,ninterior_p1,ix+1,iy+1,periodicFlgs); // 3
      // z = 1 block
      // offset2 = nglobal in 2d for p=1
      int offset2 = getNumPerpGlobalNodes(nx,ny,2,1,periodicFlgs);
      lgMap[12] = F2_func(nx,ny,ninterior,ix,iy,periodicFlgs); // 0
      lgMap[13] = ix==0 ? 
         F2_func(nx,ny,ninterior,ix+1,iy,periodicFlgs)-1 : 
         F2_func(nx,ny,ninterior,ix,iy,periodicFlgs)+1; // 1
      lgMap[14] = F2_func(nx,ny,ninterior,ix+1, iy,periodicFlgs); // 2
      lgMap[15] = G2_func(nx,ny,ninterior,ix, iy,periodicFlgs); // 3
      lgMap[16] = G2_func(nx,ny,ninterior,ix+1, iy,periodicFlgs); // 4
      lgMap[17] = F2_func(nx,ny,ninterior,ix, iy+1,periodicFlgs); // 5
      lgMap[18] = ix==0 ? 
         F2_func(nx,ny,ninterior,ix+1,iy+1,periodicFlgs)-1 : 
         F2_func(nx,ny,ninterior,ix,iy+1,periodicFlgs)+1; // 6
      lgMap[19] = F2_func(nx,ny,ninterior,ix+1,iy+1,periodicFlgs); // 7
      for(int i=12; i<20; i++) {
        lgMap[i] += (offset + offset2);
      }
    }
  }
}

int FemPerpPoisson::F1_func(int nx, int ny, int ninterior, int ix, int iy, bool periodicFlgs[2]) 
{
  // when periodic in the x (resp. y) direction, we will map the right boundary
  // to the left boundary (resp. top boundary to the bottom boundary)
  bool x_periodic = periodicFlgs[0];
  bool y_periodic = periodicFlgs[1];

  // first handle corners, which are put at the end of the 2d global mapping
  // index of top right corner is always last
  int topright = getNumPerpGlobalNodes(nx,ny,2,1,periodicFlgs)-1;
  int bottomleftoffset = 3;
  if(x_periodic && y_periodic) bottomleftoffset = 0;
  else if (x_periodic || y_periodic) bottomleftoffset = 1;
  if(ix==0 && iy==0) return topright - bottomleftoffset;
  else if(ix==nx && iy==0) return topright - bottomleftoffset + !x_periodic;
  else if(ix==0 && iy==ny) return topright - !x_periodic;
  else if(ix==nx && iy==ny) return topright;
  // now handle edges
  else if(iy==0 || (iy==ny && y_periodic)) return ix+ninterior-1;
  else if(iy==ny) {
    //if mapping right side to left side because x-periodic, offset top so that there isn't a gap in the mapping where the right side usually is
    int res = ix+(nx+1)*ny-3;
    if(x_periodic) res -= ny-1; // offset by -(ny-1)
    return res;
  }
  else if(ix==0 || (ix==nx && x_periodic)) return iy+nx+ninterior-2;
  else if(ix==nx) return iy+nx+ny+ninterior-3;
  // interior
  else return ix-1+(nx-1)*(iy-1);
}

int FemPerpPoisson::F2_func(int nx, int ny, int ninterior, int ix, int iy, bool periodicFlgs[2]) 
{
  // when periodic in the x (resp. y) direction, we will map the right boundary
  // to the left boundary (resp. top boundary to the bottom boundary)
  bool x_periodic = periodicFlgs[0];
  bool y_periodic = periodicFlgs[1];

  // first handle corners, which are put at the end of the 2d global mapping
  // index of top right corner is always last
  int topright = getNumPerpGlobalNodes(nx,ny,2,2,periodicFlgs)-1;
  int bottomleftoffset = 3;
  if(x_periodic && y_periodic) bottomleftoffset = 0;
  else if (x_periodic || y_periodic) bottomleftoffset = 1;
  if(ix==0 && iy==0) return topright - bottomleftoffset;
  else if(ix==nx && iy==0) return topright - bottomleftoffset + !x_periodic;
  else if(ix==0 && iy==ny) return topright - !x_periodic;
  else if(ix==nx && iy==ny) return topright;
  // now handle edges
  if(iy==0 || (iy==ny && y_periodic)) return 2*ix+ninterior-1;
  else if(iy==ny) {
    //if mapping right side to left side because x-periodic, offset top so that there isn't a gap in the mapping where the left side usually is
    int res = 2*ix+(3*nx+2)*ny-3;
    if(x_periodic) res -= 2*ny-1; // offset by -(2*ny-1)
    return res;
  }
  else if(ix==0 || (ix==nx && x_periodic)) return 2*(iy+nx)+ninterior-2;
  else if(ix==nx) return 2*(iy+nx+ny)+ninterior-3;
  else return 2*(ix-1)+(3*nx-2)*(iy-1)+nx;
}

int FemPerpPoisson::G2_func(int nx, int ny, int ninterior, int ix, int iy, bool periodicFlgs[2]) 
{
  bool x_periodic = periodicFlgs[0];
  if(ix==0 || (ix==nx && x_periodic)) return 2*(iy+nx)+ninterior-1;
  else if(ix==nx) return 2*(iy+nx+ny)+ninterior-2;
  else return ix+(3*nx-2)*iy-1;
}

void FemPerpPoisson::getPerpStiffnessMatrix(Eigen::MatrixXd& stiff, int ndim, int polyOrder, double dx, double dy)
{
  double dfacx2 = 4.0/(dx*dx);
  double dfacy2 = 4.0/(dy*dy);

  if(ndim==2) {
    if(polyOrder == 1) 
    {
      stiff << 0.3333333333333333*dfacx2+0.3333333333333333*dfacy2,0.1666666666666667*dfacy2-0.3333333333333333*dfacx2,0.1666666666666667*dfacx2-0.3333333333333333*dfacy2,-0.1666666666666667*dfacx2-0.1666666666666667*dfacy2,0.1666666666666667*dfacy2-0.3333333333333333*dfacx2,0.3333333333333333*dfacx2+0.3333333333333333*dfacy2,-0.1666666666666667*dfacx2-0.1666666666666667*dfacy2,0.1666666666666667*dfacx2-0.3333333333333333*dfacy2,0.1666666666666667*dfacx2-0.3333333333333333*dfacy2,-0.1666666666666667*dfacx2-0.1666666666666667*dfacy2,0.3333333333333333*dfacx2+0.3333333333333333*dfacy2,0.1666666666666667*dfacy2-0.3333333333333333*dfacx2,-0.1666666666666667*dfacx2-0.1666666666666667*dfacy2,0.1666666666666667*dfacx2-0.3333333333333333*dfacy2,0.1666666666666667*dfacy2-0.3333333333333333*dfacx2,0.3333333333333333*dfacx2+0.3333333333333333*dfacy2;
    }
    else if (polyOrder == 2)
    {
      stiff << (0.5777777777777778*dfacx2)+(0.5777777777777778*dfacy2),(0.06666666666666667*dfacy2)-(0.888888888888889*dfacx2),(0.3111111111111111*dfacx2)+(0.1888888888888889*dfacy2),(0.06666666666666667*dfacx2)-(0.888888888888889*dfacy2),-(0.06666666666666667*dfacx2)-(0.4444444444444445*dfacy2),(0.1888888888888889*dfacx2)+(0.3111111111111111*dfacy2),-(0.4444444444444445*dfacx2)-(0.06666666666666667*dfacy2),(0.2555555555555556*dfacx2)+(0.2555555555555556*dfacy2),(0.06666666666666667*dfacy2)-(0.888888888888889*dfacx2),(1.777777777777778*dfacx2)+(0.5333333333333333*dfacy2),(0.06666666666666667*dfacy2)-(0.888888888888889*dfacx2),0.0,0.0,-(0.4444444444444445*dfacx2)-(0.06666666666666667*dfacy2),(0.888888888888889*dfacx2)-(0.5333333333333333*dfacy2),-(0.4444444444444445*dfacx2)-(0.06666666666666667*dfacy2),(0.3111111111111111*dfacx2)+(0.1888888888888889*dfacy2),(0.06666666666666667*dfacy2)-(0.888888888888889*dfacx2),(0.5777777777777778*dfacx2)+(0.5777777777777778*dfacy2),-(0.06666666666666667*dfacx2)-(0.4444444444444445*dfacy2),(0.06666666666666667*dfacx2)-(0.888888888888889*dfacy2),(0.2555555555555556*dfacx2)+(0.2555555555555556*dfacy2),-(0.4444444444444445*dfacx2)-(0.06666666666666667*dfacy2),(0.1888888888888889*dfacx2)+(0.3111111111111111*dfacy2),(0.06666666666666667*dfacx2)-(0.888888888888889*dfacy2),0.0,-(0.06666666666666667*dfacx2)-(0.4444444444444445*dfacy2),(0.5333333333333333*dfacx2)+(1.777777777777778*dfacy2),(0.888888888888889*dfacy2)-(0.5333333333333333*dfacx2),(0.06666666666666667*dfacx2)-(0.888888888888889*dfacy2),0.0,-(0.06666666666666667*dfacx2)-(0.4444444444444445*dfacy2),-(0.06666666666666667*dfacx2)-(0.4444444444444445*dfacy2),0.0,(0.06666666666666667*dfacx2)-(0.888888888888889*dfacy2),(0.888888888888889*dfacy2)-(0.5333333333333333*dfacx2),(0.5333333333333333*dfacx2)+(1.777777777777778*dfacy2),-(0.06666666666666667*dfacx2)-(0.4444444444444445*dfacy2),0.0,(0.06666666666666667*dfacx2)-(0.888888888888889*dfacy2),(0.1888888888888889*dfacx2)+(0.3111111111111111*dfacy2),-(0.4444444444444445*dfacx2)-(0.06666666666666667*dfacy2),(0.2555555555555556*dfacx2)+(0.2555555555555556*dfacy2),(0.06666666666666667*dfacx2)-(0.888888888888889*dfacy2),-(0.06666666666666667*dfacx2)-(0.4444444444444445*dfacy2),(0.5777777777777778*dfacx2)+(0.5777777777777778*dfacy2),(0.06666666666666667*dfacy2)-(0.888888888888889*dfacx2),(0.3111111111111111*dfacx2)+(0.1888888888888889*dfacy2),-(0.4444444444444445*dfacx2)-(0.06666666666666667*dfacy2),(0.888888888888889*dfacx2)-(0.5333333333333333*dfacy2),-(0.4444444444444445*dfacx2)-(0.06666666666666667*dfacy2),0.0,0.0,(0.06666666666666667*dfacy2)-(0.888888888888889*dfacx2),(1.777777777777778*dfacx2)+(0.5333333333333333*dfacy2),(0.06666666666666667*dfacy2)-(0.888888888888889*dfacx2),(0.2555555555555556*dfacx2)+(0.2555555555555556*dfacy2),-(0.4444444444444445*dfacx2)-(0.06666666666666667*dfacy2),(0.1888888888888889*dfacx2)+(0.3111111111111111*dfacy2),-(0.06666666666666667*dfacx2)-(0.4444444444444445*dfacy2),(0.06666666666666667*dfacx2)-(0.888888888888889*dfacy2),(0.3111111111111111*dfacx2)+(0.1888888888888889*dfacy2),(0.06666666666666667*dfacy2)-(0.888888888888889*dfacx2),(0.5777777777777778*dfacx2)+(0.5777777777777778*dfacy2);
    }
  }
  else if(ndim==3) {
    if(polyOrder==1) {
      // for 3d p=1 assume the perp stiffness matrix has diagonal block structure.
      // each diagonal block corresponds to the 2d stiffness matrix 
      // in the plane of one of the z nodes
      int nlocal_2d = getNumLocalNodes(2, polyOrder);
      MatrixXd stiff2d_p1(nlocal_2d, nlocal_2d);
      MatrixXd zeros = MatrixXd::Zero(nlocal_2d, nlocal_2d);
      getPerpStiffnessMatrix(stiff2d_p1, 2, polyOrder, dx, dy);
      // initialize stiffness matrix with block structure
      stiff << stiff2d_p1, zeros,
               zeros,      stiff2d_p1;
    }
    else if(polyOrder==2) {
      stiff << 0.362962962962963*dfacx2+0.362962962962963*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.2296296296296296*dfacx2+0.1703703703703704*dfacy2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.1703703703703704*dfacx2+0.2296296296296296*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.1259259259259259*dfacx2+0.1259259259259259*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,0.1703703703703704*dfacx2+0.1703703703703704*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.1259259259259259*dfacx2+0.1074074074074074*dfacy2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.1074074074074074*dfacx2+0.1259259259259259*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.04074074074074074*dfacx2+0.04074074074074074*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,1.185185185185185*dfacx2+0.3555555555555556*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.0,0.0,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.5925925925925926*dfacx2-0.3555555555555556*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.2222222222222222*dfacy2,0.2222222222222222*dfacy2,-0.2222222222222222*dfacy2,-0.2222222222222222*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.5925925925925926*dfacx2+0.1777777777777778*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.0,0.0,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.2962962962962963*dfacx2-0.1777777777777778*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.2296296296296296*dfacx2+0.1703703703703704*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.362962962962963*dfacx2+0.362962962962963*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.1259259259259259*dfacx2+0.1259259259259259*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.1703703703703704*dfacx2+0.2296296296296296*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,0.1259259259259259*dfacx2+0.1074074074074074*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.1703703703703704*dfacx2+0.1703703703703704*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.04074074074074074*dfacx2+0.04074074074074074*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.1074074074074074*dfacx2+0.1259259259259259*dfacy2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.0,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.3555555555555556*dfacx2+1.185185185185185*dfacy2,0.5925925925925926*dfacy2-0.3555555555555556*dfacx2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.0,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.2222222222222222*dfacx2,-0.2222222222222222*dfacx2,0.2222222222222222*dfacx2,-0.2222222222222222*dfacx2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.0,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.1777777777777778*dfacx2+0.5925925925925926*dfacy2,0.2962962962962963*dfacy2-0.1777777777777778*dfacx2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.0,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.0,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.5925925925925926*dfacy2-0.3555555555555556*dfacx2,0.3555555555555556*dfacx2+1.185185185185185*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.0,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,-0.2222222222222222*dfacx2,0.2222222222222222*dfacx2,-0.2222222222222222*dfacx2,0.2222222222222222*dfacx2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.0,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.2962962962962963*dfacy2-0.1777777777777778*dfacx2,0.1777777777777778*dfacx2+0.5925925925925926*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.0,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.1703703703703704*dfacx2+0.2296296296296296*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.1259259259259259*dfacx2+0.1259259259259259*dfacy2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.362962962962963*dfacx2+0.362962962962963*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.2296296296296296*dfacx2+0.1703703703703704*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.1074074074074074*dfacx2+0.1259259259259259*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.04074074074074074*dfacx2+0.04074074074074074*dfacy2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.1703703703703704*dfacx2+0.1703703703703704*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.1259259259259259*dfacx2+0.1074074074074074*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.5925925925925926*dfacx2-0.3555555555555556*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.0,0.0,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,1.185185185185185*dfacx2+0.3555555555555556*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,-0.2222222222222222*dfacy2,-0.2222222222222222*dfacy2,0.2222222222222222*dfacy2,0.2222222222222222*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.2962962962962963*dfacx2-0.1777777777777778*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.0,0.0,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.5925925925925926*dfacx2+0.1777777777777778*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.1259259259259259*dfacx2+0.1259259259259259*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.1703703703703704*dfacx2+0.2296296296296296*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.2296296296296296*dfacx2+0.1703703703703704*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.362962962962963*dfacx2+0.362962962962963*dfacy2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.04074074074074074*dfacx2+0.04074074074074074*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.1074074074074074*dfacx2+0.1259259259259259*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.1259259259259259*dfacx2+0.1074074074074074*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.1703703703703704*dfacx2+0.1703703703703704*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.2222222222222222*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.2222222222222222*dfacx2,-0.2222222222222222*dfacx2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,-0.2222222222222222*dfacy2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,0.3555555555555556*dfacx2+0.3555555555555556*dfacy2,0.1777777777777778*dfacy2-0.3555555555555556*dfacx2,0.1777777777777778*dfacx2-0.3555555555555556*dfacy2,-0.1777777777777778*dfacx2-0.1777777777777778*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.2222222222222222*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.2222222222222222*dfacx2,-0.2222222222222222*dfacx2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,-0.2222222222222222*dfacy2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.2222222222222222*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,-0.2222222222222222*dfacx2,0.2222222222222222*dfacx2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,-0.2222222222222222*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,0.1777777777777778*dfacy2-0.3555555555555556*dfacx2,0.3555555555555556*dfacx2+0.3555555555555556*dfacy2,-0.1777777777777778*dfacx2-0.1777777777777778*dfacy2,0.1777777777777778*dfacx2-0.3555555555555556*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.2222222222222222*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,-0.2222222222222222*dfacx2,0.2222222222222222*dfacx2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,-0.2222222222222222*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,-0.2222222222222222*dfacy2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,0.2222222222222222*dfacx2,-0.2222222222222222*dfacx2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.2222222222222222*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.1777777777777778*dfacx2-0.3555555555555556*dfacy2,-0.1777777777777778*dfacx2-0.1777777777777778*dfacy2,0.3555555555555556*dfacx2+0.3555555555555556*dfacy2,0.1777777777777778*dfacy2-0.3555555555555556*dfacx2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,-0.2222222222222222*dfacy2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,0.2222222222222222*dfacx2,-0.2222222222222222*dfacx2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.2222222222222222*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,-0.2222222222222222*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,-0.2222222222222222*dfacx2,0.2222222222222222*dfacx2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.2222222222222222*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,-0.1777777777777778*dfacx2-0.1777777777777778*dfacy2,0.1777777777777778*dfacx2-0.3555555555555556*dfacy2,0.1777777777777778*dfacy2-0.3555555555555556*dfacx2,0.3555555555555556*dfacx2+0.3555555555555556*dfacy2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,-0.2222222222222222*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,-0.2222222222222222*dfacx2,0.2222222222222222*dfacx2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.2222222222222222*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.1703703703703704*dfacx2+0.1703703703703704*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.1259259259259259*dfacx2+0.1074074074074074*dfacy2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.1074074074074074*dfacx2+0.1259259259259259*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.04074074074074074*dfacx2+0.04074074074074074*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,0.362962962962963*dfacx2+0.362962962962963*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.2296296296296296*dfacx2+0.1703703703703704*dfacy2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.1703703703703704*dfacx2+0.2296296296296296*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.1259259259259259*dfacx2+0.1259259259259259*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.5925925925925926*dfacx2+0.1777777777777778*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.0,0.0,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.2962962962962963*dfacx2-0.1777777777777778*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.2222222222222222*dfacy2,0.2222222222222222*dfacy2,-0.2222222222222222*dfacy2,-0.2222222222222222*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,1.185185185185185*dfacx2+0.3555555555555556*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.0,0.0,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.5925925925925926*dfacx2-0.3555555555555556*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.1259259259259259*dfacx2+0.1074074074074074*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.1703703703703704*dfacx2+0.1703703703703704*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.04074074074074074*dfacx2+0.04074074074074074*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.1074074074074074*dfacx2+0.1259259259259259*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,0.2296296296296296*dfacx2+0.1703703703703704*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.362962962962963*dfacx2+0.362962962962963*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.1259259259259259*dfacx2+0.1259259259259259*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.1703703703703704*dfacx2+0.2296296296296296*dfacy2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.0,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.1777777777777778*dfacx2+0.5925925925925926*dfacy2,0.2962962962962963*dfacy2-0.1777777777777778*dfacx2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.0,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.2222222222222222*dfacx2,-0.2222222222222222*dfacx2,0.2222222222222222*dfacx2,-0.2222222222222222*dfacx2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.0,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.3555555555555556*dfacx2+1.185185185185185*dfacy2,0.5925925925925926*dfacy2-0.3555555555555556*dfacx2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.0,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.0,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.2962962962962963*dfacy2-0.1777777777777778*dfacx2,0.1777777777777778*dfacx2+0.5925925925925926*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.0,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,-0.2222222222222222*dfacx2,0.2222222222222222*dfacx2,-0.2222222222222222*dfacx2,0.2222222222222222*dfacx2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.0,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.5925925925925926*dfacy2-0.3555555555555556*dfacx2,0.3555555555555556*dfacx2+1.185185185185185*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.0,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.1074074074074074*dfacx2+0.1259259259259259*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.04074074074074074*dfacx2+0.04074074074074074*dfacy2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,0.1703703703703704*dfacx2+0.1703703703703704*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.1259259259259259*dfacx2+0.1074074074074074*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,0.1703703703703704*dfacx2+0.2296296296296296*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.1259259259259259*dfacx2+0.1259259259259259*dfacy2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,0.362962962962963*dfacx2+0.362962962962963*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.2296296296296296*dfacx2+0.1703703703703704*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.2962962962962963*dfacx2-0.1777777777777778*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.0,0.0,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.5925925925925926*dfacx2+0.1777777777777778*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,-0.2222222222222222*dfacy2,-0.2222222222222222*dfacy2,0.2222222222222222*dfacy2,0.2222222222222222*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.5925925925925926*dfacx2-0.3555555555555556*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.0,0.0,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,1.185185185185185*dfacx2+0.3555555555555556*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.04074074074074074*dfacx2+0.04074074074074074*dfacy2,0.08888888888888889*dfacy2-0.1481481481481481*dfacx2,0.1074074074074074*dfacx2+0.1259259259259259*dfacy2,0.08888888888888889*dfacx2-0.1481481481481481*dfacy2,-0.08888888888888889*dfacx2-0.2962962962962963*dfacy2,0.1259259259259259*dfacx2+0.1074074074074074*dfacy2,-0.2962962962962963*dfacx2-0.08888888888888889*dfacy2,0.1703703703703704*dfacx2+0.1703703703703704*dfacy2,0.08888888888888889*dfacx2+0.08888888888888889*dfacy2,0.06666666666666667*dfacy2-0.08888888888888889*dfacx2,0.06666666666666667*dfacx2-0.08888888888888889*dfacy2,-0.06666666666666667*dfacx2-0.06666666666666667*dfacy2,0.1259259259259259*dfacx2+0.1259259259259259*dfacy2,0.06666666666666667*dfacy2-0.2962962962962963*dfacx2,0.1703703703703704*dfacx2+0.2296296296296296*dfacy2,0.06666666666666667*dfacx2-0.2962962962962963*dfacy2,-0.06666666666666667*dfacx2-0.5925925925925926*dfacy2,0.2296296296296296*dfacx2+0.1703703703703704*dfacy2,-0.5925925925925926*dfacx2-0.06666666666666667*dfacy2,0.362962962962963*dfacx2+0.362962962962963*dfacy2;
    }
  }
}

void FemPerpPoisson::getMassMatrix(Eigen::MatrixXd& mass, int ndim, int polyOrder)
{
  if(ndim==2) {
    if(polyOrder==1) {
      mass << 0.4444444444444444,0.2222222222222222,0.2222222222222222,0.1111111111111111,0.2222222222222222,0.4444444444444444,0.1111111111111111,0.2222222222222222,0.2222222222222222,0.1111111111111111,0.4444444444444444,0.2222222222222222,0.1111111111111111,0.2222222222222222,0.2222222222222222,0.4444444444444444;
    }
    else if(polyOrder==2) {
      mass << 0.1333333333333333,-0.1333333333333333,0.04444444444444445,-0.1333333333333333,-0.1777777777777778,0.04444444444444445,-0.1777777777777778,0.06666666666666667,-0.1333333333333333,0.7111111111111111,-0.1333333333333333,0.4444444444444444,0.4444444444444444,-0.1777777777777778,0.3555555555555556,-0.1777777777777778,0.04444444444444445,-0.1333333333333333,0.1333333333333333,-0.1777777777777778,-0.1333333333333333,0.06666666666666667,-0.1777777777777778,0.04444444444444445,-0.1333333333333333,0.4444444444444444,-0.1777777777777778,0.7111111111111111,0.3555555555555556,-0.1333333333333333,0.4444444444444444,-0.1777777777777778,-0.1777777777777778,0.4444444444444444,-0.1333333333333333,0.3555555555555556,0.7111111111111111,-0.1777777777777778,0.4444444444444444,-0.1333333333333333,0.04444444444444445,-0.1777777777777778,0.06666666666666667,-0.1333333333333333,-0.1777777777777778,0.1333333333333333,-0.1333333333333333,0.04444444444444445,-0.1777777777777778,0.3555555555555556,-0.1777777777777778,0.4444444444444444,0.4444444444444444,-0.1333333333333333,0.7111111111111111,-0.1333333333333333,0.06666666666666667,-0.1777777777777778,0.04444444444444445,-0.1777777777777778,-0.1333333333333333,0.04444444444444445,-0.1333333333333333,0.1333333333333333;
    }
  }
  else if(ndim==3) {
    if(polyOrder==1) {
      // for 3d p=1 assume the mass matrix has diagonal block structure.
      // each diagonal block corresponds to the 2d mass matrix 
      // in the plane of one of the z nodes
      int nlocal_2d = getNumLocalNodes(2, polyOrder);
      MatrixXd mass2d_p1(nlocal_2d, nlocal_2d);
      MatrixXd zeros = MatrixXd::Zero(nlocal_2d, nlocal_2d);
      getMassMatrix(mass2d_p1, 2, polyOrder);
      // initialize stiffness matrix with block structure
      mass << mass2d_p1, zeros,
               zeros,      mass2d_p1;
    }
    else if(polyOrder==2) {
      mass << 0.2074074074074074,-0.237037037037037,0.162962962962963,-0.237037037037037,-0.1925925925925926,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.237037037037037,-0.1925925925925926,-0.1925925925925926,-0.1333333333333333,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,-0.1333333333333333,0.1481481481481481,-0.1333333333333333,0.1259259259259259,-0.237037037037037,0.4740740740740741,-0.237037037037037,0.2962962962962963,0.2962962962962963,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.2962962962962963,0.2962962962962963,0.1481481481481481,0.1481481481481481,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.1481481481481481,0.1481481481481481,-0.1333333333333333,0.1185185185185185,-0.1333333333333333,0.162962962962963,-0.237037037037037,0.2074074074074074,-0.1925925925925926,-0.237037037037037,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1925925925925926,-0.237037037037037,-0.1333333333333333,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1333333333333333,-0.1925925925925926,0.1259259259259259,-0.1333333333333333,0.1481481481481481,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.4740740740740741,0.237037037037037,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.2962962962962963,0.1481481481481481,0.2962962962962963,0.1481481481481481,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.237037037037037,0.1185185185185185,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.237037037037037,0.4740740740740741,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.1481481481481481,0.2962962962962963,0.1481481481481481,0.2962962962962963,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.1185185185185185,0.237037037037037,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.237037037037037,-0.1925925925925926,0.2074074074074074,-0.237037037037037,0.162962962962963,-0.1925925925925926,-0.1333333333333333,-0.237037037037037,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.1259259259259259,-0.1925925925925926,-0.1333333333333333,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.2962962962962963,0.2962962962962963,-0.237037037037037,0.4740740740740741,-0.237037037037037,0.1481481481481481,0.1481481481481481,0.2962962962962963,0.2962962962962963,-0.1333333333333333,0.1185185185185185,-0.1333333333333333,0.1481481481481481,0.1481481481481481,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1925925925925926,-0.237037037037037,0.162962962962963,-0.237037037037037,0.2074074074074074,-0.1333333333333333,-0.1925925925925926,-0.1925925925925926,-0.237037037037037,0.1259259259259259,-0.1333333333333333,0.1481481481481481,-0.1333333333333333,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.2962962962962963,0.1481481481481481,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.4740740740740741,0.237037037037037,0.237037037037037,0.1185185185185185,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.2962962962962963,0.1481481481481481,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.1481481481481481,0.2962962962962963,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.237037037037037,0.4740740740740741,0.1185185185185185,0.237037037037037,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.1481481481481481,0.2962962962962963,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.2962962962962963,0.1481481481481481,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.237037037037037,0.1185185185185185,0.4740740740740741,0.237037037037037,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.2962962962962963,0.1481481481481481,-0.237037037037037,0.2962962962962963,-0.1925925925925926,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.1481481481481481,0.2962962962962963,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.1185185185185185,0.237037037037037,0.237037037037037,0.4740740740740741,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.1481481481481481,0.2962962962962963,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,-0.1333333333333333,0.1481481481481481,-0.1333333333333333,0.1259259259259259,-0.237037037037037,-0.1925925925925926,-0.1925925925925926,-0.1333333333333333,0.2074074074074074,-0.237037037037037,0.162962962962963,-0.237037037037037,-0.1925925925925926,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.1481481481481481,0.1481481481481481,-0.1333333333333333,0.1185185185185185,-0.1333333333333333,0.2962962962962963,0.2962962962962963,0.1481481481481481,0.1481481481481481,-0.237037037037037,0.4740740740740741,-0.237037037037037,0.2962962962962963,0.2962962962962963,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1333333333333333,-0.1925925925925926,0.1259259259259259,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,-0.237037037037037,-0.1333333333333333,-0.1925925925925926,0.162962962962963,-0.237037037037037,0.2074074074074074,-0.1925925925925926,-0.237037037037037,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.237037037037037,0.1185185185185185,-0.1925925925925926,0.1481481481481481,-0.1333333333333333,0.2962962962962963,0.1481481481481481,0.2962962962962963,0.1481481481481481,-0.237037037037037,0.2962962962962963,-0.1925925925925926,0.4740740740740741,0.237037037037037,-0.237037037037037,0.2962962962962963,-0.1925925925925926,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.1185185185185185,0.237037037037037,-0.1333333333333333,0.1481481481481481,-0.1925925925925926,0.1481481481481481,0.2962962962962963,0.1481481481481481,0.2962962962962963,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.237037037037037,0.4740740740740741,-0.1925925925925926,0.2962962962962963,-0.237037037037037,0.1481481481481481,-0.1333333333333333,0.1259259259259259,-0.1925925925925926,-0.1333333333333333,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,-0.1333333333333333,-0.237037037037037,-0.1925925925925926,0.162962962962963,-0.1925925925925926,0.1481481481481481,-0.237037037037037,-0.1925925925925926,0.2074074074074074,-0.237037037037037,0.162962962962963,-0.1333333333333333,0.1185185185185185,-0.1333333333333333,0.1481481481481481,0.1481481481481481,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.1481481481481481,0.1481481481481481,0.2962962962962963,0.2962962962962963,-0.1925925925925926,0.237037037037037,-0.1925925925925926,0.2962962962962963,0.2962962962962963,-0.237037037037037,0.4740740740740741,-0.237037037037037,0.1259259259259259,-0.1333333333333333,0.1481481481481481,-0.1333333333333333,-0.1925925925925926,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1333333333333333,-0.1925925925925926,-0.1925925925925926,-0.237037037037037,0.1481481481481481,-0.1925925925925926,0.162962962962963,-0.1925925925925926,-0.237037037037037,0.162962962962963,-0.237037037037037,0.2074074074074074;
    }
  }
}

void FemPerpPoisson::getNodToModMatrix(Eigen::MatrixXd& mat, int ndim, int polyOrder)
{
  if(ndim==2) {
    if(polyOrder==1) {
      mat << 0.5,0.5,0.5,0.5,-0.2886751345948129,0.2886751345948129,-0.2886751345948129,0.2886751345948129,-0.2886751345948129,-0.2886751345948129,0.2886751345948129,0.2886751345948129,0.1666666666666667,-0.1666666666666667,-0.1666666666666667,0.1666666666666667;
    }
    else if (polyOrder==2) {
      mat << -0.166667, 0.666667, -0.166667, 0.666667, 0.666667, -0.166667, 0.666667, -0.166667,-0.096225, 0.0, 0.096225, -0.3849, 0.3849, -0.096225, 0.0, 0.096225,-0.096225, -0.3849, -0.096225, 0.0, 0.0, 0.096225, 0.3849, 0.096225,0.166667, 0.0, -0.166667, 0.0, 0.0, -0.166667, 0.0, 0.166667,0.149071, -0.298142, 0.149071, 0.0, 0.0, 0.149071, -0.298142, 0.149071,0.149071, 0.0, 0.149071, -0.298142, -0.298142, 0.149071, 0.0, 0.149071,-0.0860663, 0.172133, -0.0860663, 0.0, 0.0, 0.0860663, -0.172133, 0.0860663,-0.0860663, 0.0, 0.0860663, 0.172133, -0.172133, -0.0860663, 0.0, 0.0860663 ;
    }
  }
  else if (ndim==3) {
    if(polyOrder==1) {
      mat << 0.3535533905932737,0.3535533905932737,0.3535533905932737,0.3535533905932737,0.3535533905932739,0.3535533905932739,0.3535533905932739,0.3535533905932739,-0.2041241452319315,0.2041241452319315,-0.2041241452319315,0.2041241452319315,-0.2041241452319317,0.2041241452319316,-0.2041241452319316,0.2041241452319316,-0.2041241452319315,-0.2041241452319315,0.2041241452319315,0.2041241452319315,-0.2041241452319316,-0.2041241452319317,0.2041241452319315,0.2041241452319317,-0.2041241452319315,-0.2041241452319315,-0.2041241452319315,-0.2041241452319315,0.2041241452319316,0.2041241452319316,0.2041241452319316,0.2041241452319315,0.1178511301977579,-0.1178511301977579,-0.1178511301977579,0.1178511301977579,0.117851130197758,-0.117851130197758,-0.117851130197758,0.117851130197758,0.1178511301977579,-0.1178511301977579,0.1178511301977579,-0.1178511301977579,-0.117851130197758,0.117851130197758,-0.117851130197758,0.117851130197758,0.1178511301977579,0.1178511301977579,-0.1178511301977579,-0.1178511301977579,-0.117851130197758,-0.117851130197758,0.117851130197758,0.117851130197758,-0.0680413817439772,0.0680413817439772,0.0680413817439772,-0.06804138174397718,0.0680413817439772,-0.06804138174397721,-0.0680413817439772,0.0680413817439772;
    }
    else if (polyOrder==2) {
      mat << -0.3535533905932737,0.4714045207910317,-0.3535533905932737,0.4714045207910317,0.4714045207910317,-0.3535533905932737,0.4714045207910317,-0.3535533905932737,0.4714045207910317,0.4714045207910317,0.4714045207910317,0.4714045207910317,-0.3535533905932737,0.4714045207910317,-0.3535533905932737,0.4714045207910317,0.4714045207910317,-0.3535533905932737,0.4714045207910317,-0.3535533905932737,0.06804138174397717,0.0,-0.06804138174397717,-0.2721655269759087,0.2721655269759087,0.06804138174397717,0.0,-0.06804138174397717,-0.2721655269759087,0.2721655269759087,-0.2721655269759087,0.2721655269759087,0.06804138174397717,0.0,-0.06804138174397717,-0.2721655269759087,0.2721655269759087,0.06804138174397717,0.0,-0.06804138174397717,0.06804138174397717,-0.2721655269759087,0.06804138174397717,0.0,0.0,-0.06804138174397717,0.2721655269759087,-0.06804138174397717,-0.2721655269759087,-0.2721655269759087,0.2721655269759087,0.2721655269759087,0.06804138174397717,-0.2721655269759087,0.06804138174397717,0.0,0.0,-0.06804138174397717,0.2721655269759087,-0.06804138174397717,0.06804138174397717,-0.2721655269759087,0.06804138174397717,-0.2721655269759087,-0.2721655269759087,0.06804138174397717,-0.2721655269759087,0.06804138174397717,0.0,0.0,0.0,0.0,-0.06804138174397717,0.2721655269759087,-0.06804138174397717,0.2721655269759087,0.2721655269759087,-0.06804138174397717,0.2721655269759087,-0.06804138174397717,0.0392837100659193,0.0,-0.0392837100659193,0.0,0.0,-0.0392837100659193,0.0,0.0392837100659193,0.1571348402636772,-0.1571348402636772,-0.1571348402636772,0.1571348402636772,0.0392837100659193,0.0,-0.0392837100659193,0.0,0.0,-0.0392837100659193,0.0,0.0392837100659193,0.0392837100659193,0.0,-0.0392837100659193,0.1571348402636772,-0.1571348402636772,0.0392837100659193,0.0,-0.0392837100659193,0.0,0.0,0.0,0.0,-0.0392837100659193,0.0,0.0392837100659193,-0.1571348402636772,0.1571348402636772,-0.0392837100659193,0.0,0.0392837100659193,0.0392837100659193,0.1571348402636772,0.0392837100659193,0.0,0.0,-0.0392837100659193,-0.1571348402636772,-0.0392837100659193,0.0,0.0,0.0,0.0,-0.0392837100659193,-0.1571348402636772,-0.0392837100659193,0.0,0.0,0.0392837100659193,0.1571348402636772,0.0392837100659193,0.105409255338946,-0.210818510677892,0.105409255338946,0.0,0.0,0.105409255338946,-0.210818510677892,0.105409255338946,0.0,0.0,0.0,0.0,0.105409255338946,-0.210818510677892,0.105409255338946,0.0,0.0,0.105409255338946,-0.210818510677892,0.105409255338946,0.105409255338946,0.0,0.105409255338946,-0.210818510677892,-0.210818510677892,0.105409255338946,0.0,0.105409255338946,0.0,0.0,0.0,0.0,0.105409255338946,0.0,0.105409255338946,-0.210818510677892,-0.210818510677892,0.105409255338946,0.0,0.105409255338946,0.105409255338946,0.0,0.105409255338946,0.0,0.0,0.105409255338946,0.0,0.105409255338946,-0.210818510677892,-0.210818510677892,-0.210818510677892,-0.210818510677892,0.105409255338946,0.0,0.105409255338946,0.0,0.0,0.105409255338946,0.0,0.105409255338946,-0.06804138174397717,0.0,0.06804138174397717,0.0,0.0,0.06804138174397717,0.0,-0.06804138174397717,0.0,0.0,0.0,0.0,0.06804138174397717,0.0,-0.06804138174397717,0.0,0.0,-0.06804138174397717,0.0,0.06804138174397717,-0.06085806194501844,0.1217161238900369,-0.06085806194501844,0.0,0.0,0.06085806194501844,-0.1217161238900369,0.06085806194501844,0.0,0.0,0.0,0.0,-0.06085806194501844,0.1217161238900369,-0.06085806194501844,0.0,0.0,0.06085806194501844,-0.1217161238900369,0.06085806194501844,-0.06085806194501844,0.0,0.06085806194501844,0.1217161238900369,-0.1217161238900369,-0.06085806194501844,0.0,0.06085806194501844,0.0,0.0,0.0,0.0,-0.06085806194501844,0.0,0.06085806194501844,0.1217161238900369,-0.1217161238900369,-0.06085806194501844,0.0,0.06085806194501844,-0.06085806194501844,0.1217161238900369,-0.06085806194501844,0.0,0.0,-0.06085806194501844,0.1217161238900369,-0.06085806194501844,0.0,0.0,0.0,0.0,0.06085806194501844,-0.1217161238900369,0.06085806194501844,0.0,0.0,0.06085806194501844,-0.1217161238900369,0.06085806194501844,-0.06085806194501844,0.0,-0.06085806194501844,0.1217161238900369,0.1217161238900369,-0.06085806194501844,0.0,-0.06085806194501844,0.0,0.0,0.0,0.0,0.06085806194501844,0.0,0.06085806194501844,-0.1217161238900369,-0.1217161238900369,0.06085806194501844,0.0,0.06085806194501844,-0.06085806194501844,0.0,0.06085806194501844,0.0,0.0,-0.06085806194501844,0.0,0.06085806194501844,0.1217161238900369,-0.1217161238900369,0.1217161238900369,-0.1217161238900369,-0.06085806194501844,0.0,0.06085806194501844,0.0,0.0,-0.06085806194501844,0.0,0.06085806194501844,-0.06085806194501844,0.0,-0.06085806194501844,0.0,0.0,0.06085806194501844,0.0,0.06085806194501844,0.1217161238900369,0.1217161238900369,-0.1217161238900369,-0.1217161238900369,-0.06085806194501844,0.0,-0.06085806194501844,0.0,0.0,0.06085806194501844,0.0,0.06085806194501844,0.03513641844631532,-0.07027283689263064,0.03513641844631532,0.0,0.0,-0.03513641844631532,0.07027283689263064,-0.03513641844631532,0.0,0.0,0.0,0.0,-0.03513641844631532,0.07027283689263064,-0.03513641844631532,0.0,0.0,0.03513641844631532,-0.07027283689263064,0.03513641844631532,0.03513641844631532,0.0,-0.03513641844631532,-0.07027283689263064,0.07027283689263064,0.03513641844631532,0.0,-0.03513641844631532,0.0,0.0,0.0,0.0,-0.03513641844631532,0.0,0.03513641844631532,0.07027283689263064,-0.07027283689263064,-0.03513641844631532,0.0,0.03513641844631532,0.03513641844631532,0.0,-0.03513641844631532,0.0,0.0,-0.03513641844631532,0.0,0.03513641844631532,-0.07027283689263064,0.07027283689263064,0.07027283689263064,-0.07027283689263064,0.03513641844631532,0.0,-0.03513641844631532,0.0,0.0,-0.03513641844631532,0.0,0.03513641844631532;
    }
  }
}

// C wrappers for interfacing with FemPerpPoisson class
extern "C" void* new_FemPerpPoisson(int nx, int ny, int ndim, int polyOrder, double dx, double dy, bool periodicFlgs[2], bcdata_t bc[2][2], bool writeMatrix, double laplacianWeight, double modifierConstant)
{
  FemPerpPoisson *f = new FemPerpPoisson(nx, ny, ndim, polyOrder, dx, dy, periodicFlgs, bc, writeMatrix, laplacianWeight, modifierConstant);
  return reinterpret_cast<void*>(f);
}

extern "C" void delete_FemPerpPoisson(FemPerpPoisson* f)
{
  delete f;
}

extern "C" void createGlobalSrc(FemPerpPoisson* f, double* localSrcPtr, int idx, int idy, double intSrcVol)
{
  f->createGlobalSrc(localSrcPtr, idx, idy, intSrcVol);
}

extern "C" void zeroGlobalSrc(FemPerpPoisson* f)
{
  f->zeroGlobalSrc();
} 

extern "C" void allreduceGlobalSrc(FemPerpPoisson* f, MPI_Comm comm)
{
  f->allreduceGlobalSrc(comm);
}

extern "C" void solve(FemPerpPoisson* f)
{
  f->solve();
}

extern "C" void getSolution(FemPerpPoisson* f, double* ptr, int idx, int idy)
{
  f->getSolution(ptr, idx, idy);
}

extern "C" void getNodalSolution(FemPerpPoisson* f, double* ptr, int idx, int idy)
{
  f->getNodalSolution(ptr, idx, idy);
}

