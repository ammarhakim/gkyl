// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for 2D FEM Poisson solver
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <FemGyroaverageImpl.h>
#include <FemMatrices.h>

// std includes
#include <string>
#include <vector>

using namespace Eigen;
static const int DIRICHLET_BC = 0;
static const int DIRICHLET_VARIABLE_BC = 2;
static const int DX = 0;
static const int DY = 1;
static const int LO = 0;
static const int HI = 1;

double take_lastGy(const double &in,const double &b) { return b; }
void vectorSumGy(double *in, double *inout, int *len, MPI_Datatype *dptr)
{
  int i;
  for(i=0; i< *len; ++i) {
    inout[i] += in[i];
  }
}

FemGyroaverage::FemGyroaverage(int nx_, int ny_, int ndim_, int polyOrder_, 
                       bool periodicFlgs_[2], bcdata_t bc_[2][2], bool writeMatrix_)
  : nx(nx_), ny(ny_), ndim(ndim_), 
    polyOrder(polyOrder_), 
    writeMatrix(writeMatrix_)
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
    }
  }

  allPeriodic = false;
  if(periodicFlgs[0] && periodicFlgs[1]) allPeriodic = true;
  cornerval = 0.0;

  int nglobal = getNumPerpGlobalNodes(nx, ny, ndim, polyOrder, periodicFlgs);
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  int nlocal_z = 1;
  if(ndim==3) nlocal_z += polyOrder;

  analyzed_ = false; // flag so that stiffness matrix only analyzed once

// initialize boundary condition indexing
  setupBoundaryIndices(bc, ndim, polyOrder);

// initialize gyaverage matrix, which has dependence on idx
  localGyavgModToNod = new MatrixXd[nx]; 
  for(int i=0; i<nx; i++) {
    localGyavgModToNod[i] = MatrixXd::Zero(nlocal, nlocal);
  }

// initialize global source vector
  globalSrc = VectorXd(nglobal);

// initialize modal-nodal transformation matrices
  localNodToMod = MatrixXd::Zero(nlocal, nlocal);
  localModToNod = MatrixXd::Zero(nlocal, nlocal);

  // nodToMod = V
  getPerpNodToModMatrix(localNodToMod, ndim, polyOrder);
  
  // modToNod = U = inv(V)
  localModToNod = localNodToMod.transpose().inverse().transpose();

  if (writeMatrix)
  {
    std::string outName = "poisson-nodtomod"; 
    outName += std::to_string(ndim) + "d";
    saveMarket(localNodToMod, outName);
    outName = "poisson-modtonod"; 
    outName += std::to_string(ndim) + "d";
    saveMarket(localModToNod, outName);
  }

  // create MPI type for eigen triplet vector
  MPI_Type_contiguous(sizeof(Triplet<double>), MPI_BYTE, &MPI_triplet_t);
  MPI_Type_commit(&MPI_triplet_t);

  // create MPI operation for vector sum
  MPI_Op_create((MPI_User_function *) vectorSumGy, true, &MPI_vectorSumGy_op);
}

FemGyroaverage::~FemGyroaverage() 
{
  stiffMat.resize(0,0);
  globalSrc.resize(0);
  sourceModVec.resize(0);
  x.resize(0);
}

// get start and end indices for each boundary
// note: corners not included!
void FemGyroaverage::setupBoundaryIndices(bcdata_t bc[2][2], int ndim, int polyOrder)
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

void FemGyroaverage::makeGlobalPerpStiffnessMatrix(double *modifierWeight, int idx, int idy)
{
  int nglobal = getNumPerpGlobalNodes(nx, ny, ndim, polyOrder, periodicFlgs);
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  int nlocal_z = 1;
  if(ndim==3) nlocal_z += polyOrder;

  int nonzeros = nlocal*(std::pow(2.0, 1.0*2)+1); // estimate number of nonzeros per row

  MatrixXd localMass = MatrixXd::Zero(nlocal, nlocal);

  std::vector<int> lgMap(nlocal);
  stiffTripletList.reserve(nonzeros*nglobal); // estimate number of nonzero elements

  getMassMatrix(localMass, modifierWeight, ndim, polyOrder);

  getPerpLocalToGlobalInteriorBLRT(lgMap,idx,idy,nx,ny,ndim,polyOrder,periodicFlgs);
      
// make triplets for constructing Eigen SparseMatrix
  for (unsigned k=0; k<nlocal; ++k)
  {
    for (unsigned m=0; m<nlocal; ++m) {
      double val = localMass(k,m);
      unsigned globalIdx_k = lgMap[k];
      unsigned globalIdx_m = lgMap[m];

      stiffTripletList.push_back(Triplet<double>(globalIdx_k, globalIdx_m, val));
    }
  }
}

void FemGyroaverage::finishGlobalPerpStiffnessMatrix()
{
  int nglobal = getNumPerpGlobalNodes(nx, ny, ndim, polyOrder, periodicFlgs);
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
  for (int d=0; d<2; ++d)
  {
    for (int side=0; side<2; ++side)
    {
      if ((bc[d][side].isSet && bc[d][side].type == DIRICHLET_BC) ||
          (bc[d][side].isSet && bc[d][side].type == DIRICHLET_VARIABLE_BC))
      {
        for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
          int start = bc[d][side].istart[zblock];
          int end = bc[d][side].iend[zblock];
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
  if((bc[DX][LO].type==DIRICHLET_BC || bc[DY][LO].type==DIRICHLET_BC) ) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][LO].cornerstart[zblock];
      stiffMatRowMajor.middleRows(corner, 1) = SparseMatrix<double,RowMajor>(1, stiffMatRowMajor.cols());
      nDirichlet += 1;
    }
  }
  // bottom right
  if((bc[DX][HI].type==DIRICHLET_BC || bc[DY][LO].type==DIRICHLET_BC)) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][LO].cornerend[zblock];
      stiffMatRowMajor.middleRows(corner, 1) = SparseMatrix<double,RowMajor>(1, stiffMatRowMajor.cols());
      nDirichlet += 1;
    }
  }
  // top left
  if(bc[DX][LO].type==DIRICHLET_BC || bc[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][HI].cornerstart[zblock];
      stiffMatRowMajor.middleRows(corner, 1) = SparseMatrix<double,RowMajor>(1, stiffMatRowMajor.cols());
      nDirichlet += 1;
    }
  }
  // top right
  if(bc[DX][HI].type==DIRICHLET_BC || bc[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][HI].cornerend[zblock];
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
  std::vector<Triplet<double> > identityTripletList;
  identityTripletList.reserve(nDirichlet); // estimate number of nonzero elements

  for (int d=0; d<2; ++d)
  {
    for (int side=0; side<2; ++side)
    {
      if ((bc[d][side].isSet && bc[d][side].type == DIRICHLET_BC) ||
          (bc[d][side].isSet && bc[d][side].type == DIRICHLET_VARIABLE_BC))
      {
        for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
          int start = bc[d][side].istart[zblock];
          int end = bc[d][side].iend[zblock];
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
  if((bc[DX][LO].type==DIRICHLET_BC || bc[DY][LO].type==DIRICHLET_BC) ) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][LO].cornerstart[zblock];
// copy cols corresponding to Dirichlet BCs to sourceModMat
      sourceModMat.middleCols(corner, 1) = stiffMat.middleCols(corner, 1);
// in stiffMat, set cols corresponding to Dirichlet BCs to zero
      stiffMat.middleCols(corner, 1) = SparseMatrix<double,ColMajor>(stiffMat.rows(), 1);
// set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
      identityTripletList.push_back(Triplet<double>(corner, corner, 1.));
    }
  }
  // bottom right
  if((bc[DX][HI].type==DIRICHLET_BC || bc[DY][LO].type==DIRICHLET_BC)) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][LO].cornerend[zblock];
// copy cols corresponding to Dirichlet BCs to sourceModMat
      sourceModMat.middleCols(corner, 1) = stiffMat.middleCols(corner, 1);
// in stiffMat, set cols corresponding to Dirichlet BCs to zero
      stiffMat.middleCols(corner, 1) = SparseMatrix<double,ColMajor>(stiffMat.rows(), 1);
// set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
      identityTripletList.push_back(Triplet<double>(corner, corner, 1.));
    }
  }
  // top left
  if(bc[DX][LO].type==DIRICHLET_BC || bc[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][HI].cornerstart[zblock];
// copy cols corresponding to Dirichlet BCs to sourceModMat
      sourceModMat.middleCols(corner, 1) = stiffMat.middleCols(corner, 1);
// in stiffMat, set cols corresponding to Dirichlet BCs to zero
      stiffMat.middleCols(corner, 1) = SparseMatrix<double,ColMajor>(stiffMat.rows(), 1);
// set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
      identityTripletList.push_back(Triplet<double>(corner, corner, 1.));
    }
  }
  // top right
  if(bc[DX][HI].type==DIRICHLET_BC || bc[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][HI].cornerend[zblock];
// copy cols corresponding to Dirichlet BCs to sourceModMat
      sourceModMat.middleCols(corner, 1) = stiffMat.middleCols(corner, 1);
// in stiffMat, set cols corresponding to Dirichlet BCs to zero
      stiffMat.middleCols(corner, 1) = SparseMatrix<double,ColMajor>(stiffMat.rows(), 1);
// set diagonal of rows/cols corresponding to Dirichlet BCs to 1 (identity block)
      identityTripletList.push_back(Triplet<double>(corner, corner, 1.));
    }
  }

  // this ensures no duplicates, so that entries for these rows are at most equal to 1
  dirichletIdentity.setFromTriplets(identityTripletList.begin(), identityTripletList.end(), take_lastGy);

  //if (modifierConstant!=0.0 && laplacianConstant==0.0) {
  //  stiffMat+=modifierConstant*dirichletIdentity;
  //} else {
    stiffMat+=dirichletIdentity;
  //}

// create vector of Dirichlet values
  SparseVector<double> dirichletVec(nglobal);
  dirichletVec.reserve(nDirichlet);
  for (int d=0; d<2; ++d)
  {
    for (int side=0; side<2; ++side)
    {
      if (bc[d][side].isSet && bc[d][side].type == DIRICHLET_BC)
      {
        for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
          for(unsigned i=bc[d][side].istart[zblock]; i<=bc[d][side].iend[zblock]; i++) {
            dirichletVec.coeffRef(i) = bc[d][side].value;
          }
        }
      }
    }
  }

// handle Dirichlet corners
  // bottom left
  // also make this corner effectively Dirichlet if all periodic BCs
  if((bc[DX][LO].type==DIRICHLET_BC || bc[DY][LO].type==DIRICHLET_BC) ) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][LO].cornerstart[zblock];
      dirichletVec.coeffRef(corner) = bc[DX][LO].type==DIRICHLET_BC ? bc[DX][LO].value : bc[DY][LO].value;
    }
  }
  // bottom right
  if((bc[DX][HI].type==DIRICHLET_BC || bc[DY][LO].type==DIRICHLET_BC)) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][LO].cornerend[zblock];
      dirichletVec.coeffRef(corner) = bc[DX][HI].type==DIRICHLET_BC ? bc[DX][HI].value : bc[DY][LO].value;
    }
  }
  // top left
  if(bc[DX][LO].type==DIRICHLET_BC || bc[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][HI].cornerstart[zblock];
      dirichletVec.coeffRef(corner) = bc[DX][LO].type==DIRICHLET_BC ? bc[DX][LO].value : bc[DY][HI].value;
    }
  }
  // top right
  if(bc[DX][HI].type==DIRICHLET_BC || bc[DY][HI].type==DIRICHLET_BC) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][HI].cornerend[zblock];
      dirichletVec.coeffRef(corner) = bc[DX][HI].type==DIRICHLET_BC ? bc[DX][HI].value : bc[DY][HI].value;
    }
  }

// calculate vector to subtract from source
  sourceModVec = sourceModMat*dirichletVec;

  DMSG("Finished initializing stiffness matrices");

  stiffMat.makeCompressed();

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
    outName += std::to_string(ndim) + "d";
    saveMarket(stiffMat, outName);
  }
}

// clear out existing stuff in source vector: this is required
// otherwise successive calls to advance() will accumulate into source
// from prevous calls, which is of course not what we want.
void FemGyroaverage::zeroGlobalSrcGy() {
  int nglobal = getNumPerpGlobalNodes(nx, ny, ndim, polyOrder, periodicFlgs);
  globalSrc.setZero(nglobal);
}

void FemGyroaverage::makeGyavgMatrix(double *rho1, double *rho2, double *rho3, int idx)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  MatrixXd localGyavg = MatrixXd::Zero(nlocal, nlocal);

  getGyavgMatrix(localGyavg, rho1, rho2, rho3, ndim, polyOrder);
  
  localGyavgModToNod[idx] = localGyavg*localModToNod;
}

// called within an indexer loop over idx and idy
void FemGyroaverage::createGlobalSrcGy(double* localSrcPtr, int idx, int idy)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);

  std::vector<int> lgMap(nlocal);
  std::vector<double> localSrc(nlocal), localGyavgSrc(nlocal);

  // copy data from input
  for (unsigned k=0; k<nlocal; ++k)
    localSrc[k] = localSrcPtr[k];
  
  getPerpLocalToGlobalInteriorBLRT(lgMap,idx,idy,nx,ny,ndim,polyOrder,periodicFlgs);

// evaluate local mass matrix times local source
  for (unsigned k=0; k<nlocal; ++k)
  {
    localGyavgSrc[k] = 0.0;
    for (unsigned m=0; m<nlocal; ++m)
    {
      localGyavgSrc[k] += localGyavgModToNod[idx](k,m)*localSrc[m];
    }
  }
// map into global source vector
  for(unsigned i=0; i<nlocal; ++i)
  {
    globalSrc.coeffRef(lgMap[i]) += localGyavgSrc[i];
  }
  if (writeMatrix)
  {
    std::string outName = "poisson-src-beforeBCs-";
    outName += std::to_string(ndim) + "d";
    saveMarket(globalSrc, outName);
  }
}

void FemGyroaverage::allreduceGlobalSrcGy(MPI_Comm comm)
{
  int nglobal = getNumPerpGlobalNodes(nx, ny, ndim, polyOrder, periodicFlgs);
  // all reduce (sum) globalSrc
  MPI_Allreduce(MPI_IN_PLACE, globalSrc.data(), nglobal, MPI_DOUBLE, MPI_vectorSumGy_op, comm);
}

void FemGyroaverage::allgatherGlobalStiffGy(MPI_Comm comm)
{
  // all gather (concatenate) stiffTripletList
  std::vector<Triplet<double> > stiffTripletListGathered;
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  stiffTripletListGathered.resize(stiffTripletList.size()*nprocs); // this resize is required.. just need to make sure to reserve enough space
  MPI_Allgather(stiffTripletList.data(), stiffTripletList.size(), MPI_triplet_t, stiffTripletListGathered.data(), stiffTripletList.size(), MPI_triplet_t, comm);
  stiffTripletList = stiffTripletListGathered;
}

void FemGyroaverage::solveGy()
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
  if((bc[DX][LO].type==DIRICHLET_BC || bc[DY][LO].type==DIRICHLET_BC) ) {
    for(unsigned zblock=0; zblock<nlocal_z; ++zblock) {
      unsigned corner = bc[DY][LO].cornerstart[zblock];
      globalSrc.coeffRef(corner) = bc[DX][LO].type==DIRICHLET_BC ? bc[DX][LO].value : bc[DY][LO].value;
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
  
// solve linear system(s)
// modify non-dirichlet rows of source by subtracting sourceModVec
  x = solver.solve(globalSrc-sourceModVec);

  if (writeMatrix)
  {
    std::string outName = "poisson-sol";
    outName += std::to_string(ndim) + "d";
    saveMarket(x, outName);
  }
}

void FemGyroaverage::getSolutionGy(double* solPtr, int idx, int idy)
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

void FemGyroaverage::getNodalSolutionGy(double* solPtr, int idx, int idy)
{
  int nlocal = getNumLocalNodes(ndim, polyOrder);
  std::vector<int> lgMap(nlocal);
  getPerpLocalToGlobalInteriorBLRT(lgMap, idx,idy,nx,ny,ndim,polyOrder, periodicFlgs);

  // transform global nodal x to local nodal sol
  for (unsigned k=0; k<nlocal; ++k) {
    solPtr[k] = x.coeffRef(lgMap[k]);
  }
}

int FemGyroaverage::getNumPerpGlobalNodes(int nx, int ny, int ndim, int polyOrder, bool periodicFlgs[2])
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

int FemGyroaverage::getNumLocalNodes(int ndim, int p) 
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

void FemGyroaverage::getPerpLocalToGlobalInteriorBLRT(std::vector<int>& lgMap, int ix, int iy, int nx, int ny, int ndim, int polyOrder, bool periodicFlgs[2])
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

int FemGyroaverage::F1_func(int nx, int ny, int ninterior, int ix, int iy, bool periodicFlgs[2]) 
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

int FemGyroaverage::F2_func(int nx, int ny, int ninterior, int ix, int iy, bool periodicFlgs[2]) 
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

int FemGyroaverage::G2_func(int nx, int ny, int ninterior, int ix, int iy, bool periodicFlgs[2]) 
{
  bool x_periodic = periodicFlgs[0];
  if(ix==0 || (ix==nx && x_periodic)) return 2*(iy+nx)+ninterior-1;
  else if(ix==nx) return 2*(iy+nx+ny)+ninterior-2;
  else return ix+(3*nx-2)*iy-1;
}

// C wrappers for interfacing with FemGyroaverage class
extern "C" void* new_FemGyroaverage(int nx, int ny, int ndim, int polyOrder, bool periodicFlgs[2], bcdata_t bc[2][2], bool writeMatrix)
{
  FemGyroaverage *f = new FemGyroaverage(nx, ny, ndim, polyOrder, periodicFlgs, bc, writeMatrix);
  return reinterpret_cast<void*>(f);
}

extern "C" void delete_FemGyroaverage(FemGyroaverage* f)
{
  delete f;
}

extern "C" void makeGlobalStiffGy(FemGyroaverage* f, double *modifierWeight, int idx, int idy)
{
  f->makeGlobalPerpStiffnessMatrix(modifierWeight, idx, idy);
}

extern"C" void makeGyavgMatrix(FemGyroaverage* f, double *rho1, double *rho2, double *rho3, int idx)
{
  f->makeGyavgMatrix(rho1, rho2, rho3, idx);
}

extern "C" void finishGlobalStiffGy(FemGyroaverage* f)
{
  f->finishGlobalPerpStiffnessMatrix();
}

extern "C" void createGlobalSrcGy(FemGyroaverage* f, double* localSrcPtr, int idx, int idy)
{
  f->createGlobalSrcGy(localSrcPtr, idx, idy);
}

extern "C" void zeroGlobalSrcGy(FemGyroaverage* f)
{
  f->zeroGlobalSrcGy();
} 

extern "C" void allreduceGlobalSrcGy(FemGyroaverage* f, MPI_Comm comm)
{
  f->allreduceGlobalSrcGy(comm);
}

extern "C" void allgatherGlobalStiffGy(FemGyroaverage* f, MPI_Comm comm)
{
  f->allgatherGlobalStiffGy(comm);
}

extern "C" void solveGy(FemGyroaverage* f)
{
  f->solveGy();
}

extern "C" void getSolutionGy(FemGyroaverage* f, double* ptr, int idx, int idy)
{
  f->getSolutionGy(ptr, idx, idy);
}

extern "C" void getNodalSolutionGy(FemGyroaverage* f, double* ptr, int idx, int idy)
{
  f->getNodalSolutionGy(ptr, idx, idy);
}

