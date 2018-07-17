#include "MappedPoissonImpl.h"
using namespace Eigen;

/////////////////////////////INITIALIZATION/////////////////////////////////////

//initialize pointer func
void* new_MapPoisson(int nx, int ny, double xl, double yl, double xu, double yu,
  int xbctypel, int ybctypel, int xbctypeu, int ybctypeu, double xbcu,
  double xbcl, double ybcu, double ybcl)
{
  MapPoisson *d = new MapPoisson(nx, ny, xl, yl, xu, yu, xbctypel, ybctypel, xbctypeu,
    ybctypeu, xbcu, xbcl, ybcu, ybcl);
  return reinterpret_cast<void*>(d);
}

//constructor
MapPoisson::MapPoisson(int nx_, int ny_, double xl_, double yl_, double xu_, double yu_,
int xbctypel_, int ybctypel_, int xbctypeu_, int ybctypeu_, double xbcu_,
double xbcl_, double ybcu_, double ybcl_)
    : nx(nx_), ny(ny_), xl(xl_), yl(yl_), xu(xu_), yu(yu_), xbctypel(xbctypel_),
      ybctypel(ybctypel_), xbctypeu(xbctypeu_), ybctypeu(ybctypeu_), xbcu(xbcu_),
      xbcl(xbcl_), ybcu(ybcu_), ybcl(ybcl_)
{

  dxc = (xu-xl)/nx;
  dyc = (yu-yl)/ny;
  xca = MatrixXd::Zero(nx, ny);
  yca = MatrixXd::Zero(nx, ny);

  //create xcomp mesh - should be monotonically increasing
  for (int j = 0; j < ny; j++){
    for (int i = 0; i < nx; i++){
      xca(i,j) = xl+0.5*dxc+i*dxc;
    }
  }
  //creat ycomp mesh - should be monotonically increasing
  for (int i = 0; i < nx; i++){
    for (int j = 0; j < ny; j++){
      yca(i,j) = yl+0.5*dyc+j*dyc;
    }
  }

  jsource = VectorXd::Zero(nx*ny);
}


//////////////////////////MATRIX CONSTRUCTORS//////////////////////////////////


//BLOCK CREATOR
blkstor MapPoisson::blkrow(int j)
{
  //initialize blocks
  MatrixXd blkjp = MatrixXd::Zero(nx,nx+2);
  MatrixXd blkj = MatrixXd::Zero(nx,nx+2);
  MatrixXd blkjm = MatrixXd::Zero(nx,nx+2);

  for (int i = 0; i < nx; i++){
   //center block
    blkj(i,i+1) = MapPoisson::cij(xca(i,j),yca(i,j),i,j);
    blkj(i,i) = MapPoisson::cimj(xca(i,j),yca(i,j),i,j);
    blkj(i,i+2) = MapPoisson::cipj(xca(i,j),yca(i,j),i,j);
   //upper block (upper in logical space)
    blkjp(i,i+1) = MapPoisson::cijp(xca(i,j),yca(i,j),i,j);
    blkjp(i,i) = MapPoisson::cimjp(xca(i,j),yca(i,j),i,j);
    blkjp(i,i+2) = MapPoisson::cipjp(xca(i,j),yca(i,j),i,j);
   //lower block (lower in logical space)
    blkjm(i,i+1) = MapPoisson::cijm(xca(i,j),yca(i,j),i,j);
    blkjm(i,i) = MapPoisson::cimjm(xca(i,j),yca(i,j),i,j);
    blkjm(i,i+2) = MapPoisson::cipjm(xca(i,j),yca(i,j),i,j);
  }

 // edit blocks according to xbcs
  MatrixXd bjp = blkjp.block(0,1,nx,nx);
  MatrixXd bj = blkj.block(0,1,nx,nx);
  MatrixXd bjm = blkjm.block(0,1,nx,nx);
  MatrixXd m = MatrixXd::Zero(nx,nx);

  //lower x values
  if (xbctypel == 0){ //Dirichlet
    m.col(0) = blkjp.col(0);
    bjp = bjp - m;
    m.col(0) = blkj.col(0);
    bj = bj - m;
    m.col(0) = blkjm.col(0);
    bjm = bjm - m;
  }
  else if (xbctypel == 1){ //neumann
    m.col(0) = blkjp.col(0);
    bjp = bjp + m;
    m.col(0) = blkj.col(0);
    bj = bj + m;
    m.col(0) = blkjm.col(0);
    bjm = bjm + m;
  }
  else if (xbctypel == 2){ //periodic
    m.col(nx-1) = blkjp.col(0);
    bjp = bjp + m;
    m.col(nx-1) = blkj.col(0);
    bj = bj + m;
    m.col(nx-1) = blkjm.col(0);
    bjm = bjm + m;
  }

  m = MatrixXd::Zero(nx,nx);
  //upper x values
  if (xbctypeu == 0){ //Dirichlet
    m.col(nx-1) = blkjp.col(nx+1);
    bjp = bjp - m;
    m.col(nx-1) = blkj.col(nx+1);
    bj = bj - m;
    m.col(nx-1) = blkjm.col(nx+1);
    bjm = bjm - m;
  }
  else if (xbctypeu == 1){ //neumann (extra from slope added on source side)
    m.col(nx-1) = blkjp.col(nx+1);
    bjp = bjp + m;
    m.col(nx-1) = blkj.col(nx+1);
    bj = bj + m;
    m.col(nx-1) = blkjm.col(nx+1);
    bjm = bjm + m;
  }
  else if (xbctypeu == 2){ //periodic
    m.col(0) = blkjp.col(nx+1);
    bjp = bjp + m;
    m.col(0) = blkj.col(nx+1);
    bj = bj + m;
    m.col(0) = blkjm.col(nx+1);
    bjm = bjm + m;
  }

  blkstor triblock = {bjp, bj, bjm};

  return triblock;
}



//LAPLACIAN CREATOR (BLOCK ORGANIZER)
void MapPoisson::laplace()
{
  //triplet dtype               //modify for y bcs (ybctyl/u = 3)
  typedef Triplet<double> T;
  std::vector<T> trips;

  for (int j = 1; j< ny-1; j++){
    blkstor threes = MapPoisson::blkrow(j);
    MatrixXd up = threes.upper;
    MatrixXd mid = threes.middle;
    MatrixXd down = threes.lower;

    //placement  of indices not at y bounds
    for (int k = 0; k < nx; k++){
      for (int i = 0; i < nx; i++){
        if (up(k,i)*up(k,i) != 0.0){ //sometimes it gives -0.0
          trips.push_back(T(k+j*nx, i+(j+1)*nx, up(k,i)));
        }
        if (mid(k,i)*mid(k,i) != 0.0){
          trips.push_back(T(k+j*nx, i+j*nx, mid(k,i)));
        }
        if (down(k,i)*down(k,i) != 0.0){
          trips.push_back(T(k+j*nx, i+(j-1)*nx, down(k,i)));
        }
      }
    }
  }

  blkstor threes = MapPoisson::blkrow(0);
  MatrixXd up = threes.upper;
  MatrixXd mid = threes.middle;
  MatrixXd down = threes.lower;
  ///make correct for triplet construction
  if (ybctypel == 0){ //dirichlet
    for (int k = 0; k < nx; k++){
      for (int i = 0; i < nx; i++){
        if (up(k,i)*up(k,i) != 0.0){ //upper block normal
          trips.push_back(T(k, i+nx, up(k,i)));
        }
        MatrixXd comb = mid-down; //mid block diffed with lower
        if (comb(k,i)*comb(k,i) != 0.0){
          trips.push_back(T(k, i, comb(k,i)));
        }
      }
    }
  }
  else if (ybctypel == 1){ //neumannn (extra from slope added on source side)
    for (int k = 0; k < nx; k++){
      for (int i = 0; i < nx; i++){
        if (up(k,i)*up(k,i) != 0.0){ //upper block normal
          trips.push_back(T(k, i+nx, up(k,i)));
        }
        MatrixXd comb = mid+down; //mid block added with lower
        if (comb(k,i)*comb(k,i) != 0.0){
          trips.push_back(T(k, i, comb(k,i)));
        }
      }
    }
  }
  else if (ybctypel == 2){ //periodic
    for (int k = 0; k < nx; k++){
      for (int i = 0; i < nx; i++){
        if (up(k,i)*up(k,i) != 0.0){
          trips.push_back(T(k, i+nx, up(k,i)));
        }
        if (mid(k,i)*mid(k,i) != 0.0){
          trips.push_back(T(k, i, mid(k,i)));
        }
        if (down(k,i)*down(k,i) != 0.0){ //lower in upper right corner
          trips.push_back(T(k, i+(ny-1)*nx, down(k,i)));
        }
      }
    }
  }
  else if (ybctypel == 3){ //nonsingular sphere map
    down = down.rowwise().reverse().eval();
    mid = mid+down;
    for (int k = 0; k < nx; k++){
      for (int i = 0; i < nx; i++){
        if (up(k,i)*up(k,i) != 0.0){
          trips.push_back(T(k, i+nx, up(k,i)));
        }
        if (mid(k,i)*mid(k,i) != 0.0){
          trips.push_back(T(k, i, mid(k,i)));
        }
      }
    }
  }


  threes = MapPoisson::blkrow(ny-1);
  up = threes.upper;
  mid = threes.middle;
  down = threes.lower;
  //y upper values
  if (ybctypeu == 0){ //dirichlet
    for (int k = 0; k < nx; k++){
      for (int i = 0; i < nx; i++){
        if (down(k,i)*down(k,i) != 0.0){ //lower block normal
          trips.push_back(T(k+(ny-1)*nx, i+(ny-2)*nx, down(k,i)));
        }
        MatrixXd comb = mid-up; //mid block diffed with upper
        if (comb(k,i)*comb(k,i) != 0.0){
          trips.push_back(T(k+(ny-1)*nx, i+(ny-1)*nx, comb(k,i)));
        }
      }
    }
  }
  else if (ybctypeu == 1){ //neumannn (extra from slope added on source side)
    for (int k = 0; k < nx; k++){
      for (int i = 0; i < nx; i++){
        if (down(k,i)*down(k,i) != 0.0){ //lower block normal
          trips.push_back(T(k+(ny-1)*nx, i+(ny-2)*nx, down(k,i)));
        }
        MatrixXd comb = mid+up; //mid block added with upper
        if (comb(k,i)*comb(k,i) != 0.0){
          trips.push_back(T(k+(ny-1)*nx, i+(ny-1)*nx, comb(k,i)));
        }
      }
    }
  }
  else if (ybctypeu == 2){ //periodic
    for (int k = 0; k < nx; k++){
      for (int i = 0; i < nx; i++){
        if (up(k,i)*up(k,i) != 0.0){ //upper in lower left corner
          trips.push_back(T(k+(ny-1)*nx, i, up(k,i)));
        }
        if (mid(k,i)*mid(k,i) != 0.0){
          trips.push_back(T(k+(ny-1)*nx, i+(ny-1)*nx, mid(k,i)));
        }
        if (down(k,i)*down(k,i) != 0.0){
          trips.push_back(T(k+(ny-1)*nx, i+(ny-2)*nx, down(k,i)));
        }
      }
    }
  }
  else if (ybctypeu == 3){ //nonsingular sphere map
    up = up.rowwise().reverse().eval();
    mid = mid+up;
    for (int k = 0; k < nx; k++){
      for (int i = 0; i < nx; i++){
        if (mid(k,i)*mid(k,i) != 0.0){
          trips.push_back(T(k+(ny-1)*nx, i+(ny-1)*nx, mid(k,i)));
        }
        if (down(k,i)*down(k,i) != 0.0){
          trips.push_back(T(k+(ny-1)*nx, i+(ny-2)*nx, down(k,i)));
        }
      }
    }
  }

  SparseMatrix<double, ColMajor> lhstrip(nx*ny,nx*ny);
  lhstrip.setFromTriplets(trips.begin(), trips.end());

  //deal with underdefinition (neumann/periodic) by elim last row, rep last w/ 1
  if (xbctypel != 0 and xbctypeu != 0 and ybctypel != 0 and ybctypeu != 0){
    for (int l = 0; l < nx*ny; l++){
      lhstrip.coeffRef(nx*ny-1,l) = 0.0;
    }
    lhstrip.prune(0.0);
    lhstrip.coeffRef(nx*ny-1,nx*ny-1) = 1.0;
  }


  solver.compute(lhstrip);  //compute factorization

  if(solver.info()!=Success) { //check if decomp worked
    std::cout << "ERROR: DECOMPOSITION FAILURE" << std::endl;// decomp failed
  }


    //std::cout << MatrixXd(lhstrip) << std:: endl;

}


///////////////////////////SOURCE CONVERTER/////////////////////////////////////

VectorXd MapPoisson::srcconv()
{
  //check if both dirs periodic - integrate out the mean
  double integral = 0.0;
  if (xbctypel == 2 and xbctypeu == 2 and ybctypel == 2 and ybctypeu == 2){
    int s = 0;
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        double gdetk = MapPoisson::gdet(MapPoisson::gij(xca(i,j),yca(i,j)));
        integral += sqrt(gdetk)*jsource(s)*dxc*dyc;
        s = s+1;
      }
    }
  }

  VectorXd rhst = VectorXd::Zero(nx*ny);

  //BC Part (currently set up for x dir bcs):
  int k = 0;
  if (ybctypel == 0){  //y dirichlet lower
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        if (j == 0){//lower bound
          rhst(k)+= -(MapPoisson::cijm(xca(i,j),yca(i,j),i,j)+MapPoisson::cimjm(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cipjm(xca(i,j),yca(i,j),i,j))*2*ybcl;
        }
        k = k+1;
      }
    }
  }

  k = 0;
  if (ybctypeu == 0){  //y dirichlet upper
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        if (j == ny-1){//upper bound
          rhst(k)+= -(MapPoisson::cijp(xca(i,j),yca(i,j),i,j)+MapPoisson::cipjp(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cimjp(xca(i,j),yca(i,j),i,j))*2*ybcu;
        }
        k = k+1;
      }
    }
  }

  k = 0;
  if (xbctypel == 0){  //x dirichlet lower
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        if (i == 0 and j != 0 and j != ny-1){//lower bound
          rhst(k)+= -(MapPoisson::cimj(xca(i,j),yca(i,j),i,j)+MapPoisson::cimjm(xca(i,j),yca(i,j),i,j)
                      +MapPoisson::cimjp(xca(i,j),yca(i,j),i,j))*2*xbcl;
        } else if (i == 0 and j ==0 and ybctypel!= 2){
          rhst(k)+= -(MapPoisson::cimj(xca(i,j),yca(i,j),i,j)
                      +MapPoisson::cimjp(xca(i,j),yca(i,j),i,j))*2*xbcl;
        } else if (i == 0 and j == ny-1 and ybctypel != 2){
          rhst(k)+= -(MapPoisson::cimj(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cimjm(xca(i,j),yca(i,j),i,j))*2*xbcl;
        } else if (i == 0 and ybctypel == 2){
          rhst(k)+= -(MapPoisson::cimj(xca(i,j),yca(i,j),i,j)+MapPoisson::cimjm(xca(i,j),yca(i,j),i,j)
                      +MapPoisson::cimjp(xca(i,j),yca(i,j),i,j))*2*xbcl;
        }
        k = k+1;
      }
    }
  }

  k = 0;
  if (xbctypeu == 0){  //x dirichlet upper
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        if (i == nx-1 and j != 0 and j!= ny-1){//upper bound
          rhst(k)+= -(MapPoisson::cipj(xca(i,j),yca(i,j),i,j)+MapPoisson::cipjm(xca(i,j),yca(i,j),i,j)
                      +MapPoisson::cipjp(xca(i,j),yca(i,j),i,j))*2*xbcu;
        } else if (i == nx-1 and j ==0 and ybctypel!= 2){
          rhst(k)+= -(MapPoisson::cipj(xca(i,j),yca(i,j),i,j)
                      +MapPoisson::cipjp(xca(i,j),yca(i,j),i,j))*2*xbcu;
        } else if (i == nx-1 and j == ny-1 and ybctypel!= 2){
          rhst(k)+= -(MapPoisson::cipj(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cipjm(xca(i,j),yca(i,j),i,j))*2*xbcu;
        } else if (i == nx-1 and ybctypel == 2){
          rhst(k)+= -(MapPoisson::cipj(xca(i,j),yca(i,j),i,j)+MapPoisson::cipjm(xca(i,j),yca(i,j),i,j)
                      +MapPoisson::cipjp(xca(i,j),yca(i,j),i,j))*2*xbcu;
        }
        k = k+1;
      }
    }
  }

  k = 0; //y neumann lower
  if (ybctypel == 1){  //y neumann
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        if (j == 0){//lower bound
          rhst(k)+= (MapPoisson::cijm(xca(i,j),yca(i,j),i,j)+MapPoisson::cipjm(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cimjm(xca(i,j),yca(i,j),i,j))*ybcl*dyc;
        }
        k = k+1;
      }
    }
  }

  k = 0; //y neumann upper
  if (ybctypeu == 1){  //y neumann
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        if (j == ny-1){//upper bound
          rhst(k)+= -(MapPoisson::cijp(xca(i,j),yca(i,j),i,j)+MapPoisson::cimjp(xca(i,j),yca(i,j),i,j)
                     +MapPoisson::cipjp(xca(i,j),yca(i,j),i,j))*ybcu*dyc;
        }
        k = k+1;
      }
    }
  }

  k = 0;
  if (xbctypel == 1){  //x neumann lower
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        if (i == 0 and j != 0 and j != ny-1){//lower bound
          rhst(k)+= (MapPoisson::cimj(xca(i,j),yca(i,j),i,j)+MapPoisson::cimjm(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cimjp(xca(i,j), yca(i,j),i,j))*xbcl*dxc;
        } else if (i == 0 and j ==0 and ybctypel!= 2){
          rhst(k)+= (MapPoisson::cimj(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cimjp(xca(i,j), yca(i,j),i,j))*xbcl*dxc;
        } else if (i == 0 and j == ny-1 and ybctypel!= 2){
          rhst(k)+= (MapPoisson::cimj(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cimjm(xca(i,j),yca(i,j),i,j))*xbcl*dxc;
        } else if (i == 0 and ybctypel == 2){
          rhst(k)+= (MapPoisson::cimj(xca(i,j),yca(i,j),i,j)+MapPoisson::cimjm(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cimjp(xca(i,j), yca(i,j),i,j))*xbcl*dxc;
        }
        k = k+1;
      }
    }
  }

  k = 0;
  if (xbctypeu == 1){  //x neumann upper
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        if (i == nx-1 and j != 0 and j!= ny-1){//upper bound
          rhst(k)+= -(MapPoisson::cipj(xca(i,j),yca(i,j),i,j)+MapPoisson::cipjp(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cipjm(xca(i,j),yca(i,j),i,j))*xbcu*dxc;
        } else if (i == nx-1 and j ==0 and ybctypel!= 2){
          rhst(k)+= -(MapPoisson::cipj(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cipjp(xca(i,j),yca(i,j),i,j))*xbcu*dxc;
        } else if (i == nx-1 and j == ny-1 and ybctypel!= 2){
          rhst(k)+= -(MapPoisson::cipj(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cipjm(xca(i,j),yca(i,j),i,j))*xbcu*dxc;
        } else if (i == nx-1 and ybctypel == 2){
          rhst(k)+= -(MapPoisson::cipj(xca(i,j),yca(i,j),i,j)+MapPoisson::cipjp(xca(i,j),yca(i,j),i,j)
                    +MapPoisson::cipjm(xca(i,j),yca(i,j),i,j))*xbcu*dxc;
        }
        k = k+1;
      }
    }
  }

  //add in contribution from source
  k = 0;
  for (int j = 0; j < ny; j++){
    for (int i = 0; i < nx; i++){
      double gdetij = MapPoisson::gdet(MapPoisson::gij(xca(i,j),yca(i,j)));
      rhst(k) += jsource(k)*sqrt(gdetij)*dxc*dyc - integral/(nx*ny);
      k = k+1;
    }
  }

  //deal with underdefinition (neumann/periodic) by elim last rhs term
  if (xbctypel != 0 and xbctypeu != 0 and ybctypel != 0 and ybctypeu != 0){
    rhst(nx*ny-1) = 0; //set last corner phi value to 0
    std::cout << "corner fixed to 0" << std::endl;
  }

  return rhst;
}


//////////////////////////METRIC CALCULATOR////////////////////////////////////


//metric lua func wrapper
void setMetricFuncPointer_MapPoisson(MapPoisson *d, void (*gfunc)(double *xcl, double *gl))
{
  d->setMetricFuncPointer(gfunc);
}

//lua metric converter to eigen dtype
Vector3d MapPoisson::gij(double xc, double yc){

  double xcp[3] = {0.0, xc, yc};  //oversized to acct for index difference
  double g[4];
  //c wrapped lua function for metric
  MapPoisson::gfunc(xcp, g);
  Vector3d gconv;
  gconv << g[1], g[2], g[3];

  if ( g[1] != g[1] or g[2] != g[2] or g[3] != g[3]){
    std::cout << gconv << std::endl;
    std::cout << "WARNING: CALCULATION WILL FAIL DUE TO NANS IN METRIC" << std::endl;
  }
  if (gconv(0)*gconv(2)-gconv(1)*gconv(1) < 0){
    std::cout << "determinant negative at:" << std::endl;
    std::cout << xc << "  " << yc << std::endl;
  }
  return gconv;
}

//simple 2x2 determinant
double MapPoisson::gdet(Vector3d gvec){
  double det = gvec(0)*gvec(2)-gvec(1)*gvec(1);
  if ( det != det ){
    //std::cout << det << std::endl;
    std::cout << "WARNING: CALCULATION WILL FAIL DUE TO NANS IN METRIC DETERMINANT" << std::endl;
  }
  if ( det < 0 ){
    //std::cout << det << std::endl;
    std::cout << "WARNING: CALCULATION WILL FAIL DUE NEGATIVE METRIC DETERMINANT" << std::endl;
  }
  return det;
}

//simple 2x2 inverse
Vector3d MapPoisson::ginv(Vector3d gvec){
  double det = gvec(0)*gvec(2)-gvec(1)*gvec(1);
  Vector3d inv;
  inv << gvec(2)/det, -gvec(1)/det, gvec(0)/det;
  if (inv(0) != inv(0) or inv(1) != inv(1) or inv(2) != inv(2)){
    //std::cout << inv << std::endl;
    std::cout << "WARNING: CALCULATION WILL FAIL DUE TO NANS IN INVERSE METRIC" << std::endl;
  }
  return inv;
}

//////////////////////////COEFFICIENT CALCULATORS//////////////////////////////

double MapPoisson::cij(double xc, double yc, int i, int j){//cij
  //calculate coefficients
  Vector3d g = MapPoisson::gij(xc-0.5*dxc,yc);
  Vector3d gi = MapPoisson::ginv(g);
  double h1 = -dyc*sqrt(MapPoisson::gdet(g))*gi(0)/dxc;
  g = MapPoisson::gij(xc,yc+0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h2 = -dxc*sqrt(MapPoisson::gdet(g))*gi(2)/dyc;
  g = MapPoisson::gij(xc+0.5*dxc,yc);
  gi = MapPoisson::ginv(g);
  double h3 = -dyc*sqrt(MapPoisson::gdet(g))*gi(0)/dxc;
  g = MapPoisson::gij(xc,yc-0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h4 = -dxc*sqrt(MapPoisson::gdet(g))*gi(2)/dyc;
  return h1+h2+h3+h4;
}

double MapPoisson::cimj(double xc, double yc, int i, int j){//ci-1,j
  //calculate coefficients
  Vector3d g = MapPoisson::gij(xc-0.5*dxc,yc);
  Vector3d gi = MapPoisson::ginv(g);
  double h1 = dyc*sqrt(MapPoisson::gdet(g))*gi(0)/dxc;
  g = MapPoisson::gij(xc,yc+0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h2 = -sqrt(MapPoisson::gdet(g))*gi(1)/4;
  g = MapPoisson::gij(xc,yc-0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h3 = sqrt(MapPoisson::gdet(g))*gi(1)/4;
  //dirichlet correctiom
  double coeff = 0;
  if (xbctypel==0 and xbctypeu==0 and ybctypel==0 and ybctypeu==0){ //currently only imp for all dirich
    if (i!=0 and i!=nx-1 and j!=0 and j!=ny-1){ //normal
      coeff = h1+h2+h3;
    } else if (j!=0 and j!=ny-1 and i==0){ //left side
      coeff = h1+h2+h3;
    } else if (j!=0 and j!=nx-1 and i==nx-1){ //right side
      coeff = h1+h2+h3;
    } else if (i!=0 and i!=nx-1 and j==0){ //bottom side
      coeff = h1+h2;
    } else if (i!=0 and i!=nx-1 and j==ny-1){ //top side
      coeff = h1+h3;
    } else if (j==0 and i==0){ //bottom left corner
      coeff = h1+h2;
    } else if (j==0 and i==nx-1){ //bottom right corner
      coeff = h1+h2;
    } else if (j==ny-1 and i==0){ //top left corner
      coeff = h1+h3;
    } else if (j==ny-1 and i==nx-1){ //top right corner
      coeff = h1+h3;
    }
  } else {
    coeff = h1+h2+h3;
  }
  return coeff;
}

double MapPoisson::cimjp(double xc, double yc, int i, int j){//ci-1,j+1
  //calculate coefficients
  Vector3d g = MapPoisson::gij(xc-0.5*dxc,yc);
  Vector3d gi = MapPoisson::ginv(g);
  double h1 = -sqrt(MapPoisson::gdet(g))*gi(1)/4;
  g = MapPoisson::gij(xc,yc+0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h2 = -sqrt(MapPoisson::gdet(g))*gi(1)/4;
  //dirichlet correction
  double coeff = 0;
  if (xbctypel==0 and xbctypeu==0 and ybctypel==0 and ybctypeu==0){ //currently only imp for all dirich
    if (i!=0 and i!=nx-1 and j!=0 and j!=ny-1){ //normal
      coeff = h1+h2;
    } else if (j!=0 and j!=ny-1 and i==0){ //left side
      coeff = h2;
    } else if (j!=0 and j!=nx-1 and i==nx-1){ //right side
      coeff = h1+h2;
    } else if (i!=0 and i!=nx-1 and j==0){ //bottom side
      coeff = h1+h2;
    } else if (i!=0 and i!=nx-1 and j==ny-1){ //top side
      coeff = h1;
    } else if (j==0 and i==0){ //bottom left corner
      coeff = h2;
    } else if (j==0 and i==nx-1){ //bottom right corner
      coeff = h1+h2;
    } else if (j==ny-1 and i==0){ //top left corner
      coeff = 0;
    } else if (j==ny-1 and i==nx-1){ //top right corner
      coeff = h1;
    }
  } else {
    coeff = h1+h2;
  }
  return coeff;
}

double MapPoisson::cijp(double xc, double yc, int i, int j){//ci,j+1
  //calculate coefficients
  Vector3d g = MapPoisson::gij(xc-0.5*dxc,yc);
  Vector3d gi = MapPoisson::ginv(g);
  double h1 = -sqrt(MapPoisson::gdet(g))*gi(1)/4;
  g = MapPoisson::gij(xc,yc+0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h2 = dxc*sqrt(MapPoisson::gdet(g))*gi(2)/dyc;
  g = MapPoisson::gij(xc+0.5*dxc,yc);
  gi = MapPoisson::ginv(g);
  double h3 = sqrt(MapPoisson::gdet(g))*gi(1)/4;
  //dirichlet correction
  double coeff = 0;
  if (xbctypel==0 and xbctypeu==0 and ybctypel==0 and ybctypeu==0){ //currently only imp for all dirich
    if (i!=0 and i!=nx-1 and j!=0 and j!=ny-1){ //normal
      coeff = h1+h2+h3;
    } else if (j!=0 and j!=ny-1 and i==0){ //left side
      coeff = h2+h3;
    } else if (j!=0 and j!=nx-1 and i==nx-1){ //right side
      coeff = h1+h2;
    } else if (i!=0 and i!=nx-1 and j==0){ //bottom side
      coeff = h1+h2+h3;
    } else if (i!=0 and i!=nx-1 and j==ny-1){ //top side
      coeff = h1+h2+h3;
    } else if (j==0 and i==0){ //bottom left corner
      coeff = h2+h3;
    } else if (j==0 and i==nx-1){ //bottom right corner
      coeff = h1+h2;
    } else if (j==ny-1 and i==0){ //top left corner
      coeff = h2+h3;
    } else if (j==ny-1 and i==nx-1){ //top right corner
      coeff = h1+h2;
    }
  } else {
    coeff = h1+h2+h3;
  }
  return coeff;
}

double MapPoisson::cipjp(double xc, double yc, int i, int j){//ci+1,j+1
  //calculate coefficients
  Vector3d g = MapPoisson::gij(xc+0.5*dxc,yc);
  Vector3d gi = MapPoisson::ginv(g);
  double h1 = sqrt(MapPoisson::gdet(g))*gi(1)/4;
  g = MapPoisson::gij(xc,yc+0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h2 = sqrt(MapPoisson::gdet(g))*gi(1)/4;
  //dirichlet correction
  double coeff = 0;
  if (xbctypel==0 and xbctypeu==0 and ybctypel==0 and ybctypeu==0){ //currently only imp for all dirich
    if (i!=0 and i!=nx-1 and j!=0 and j!=ny-1){ //normal
      coeff = h1+h2;
    } else if (j!=0 and j!=ny-1 and i==0){ //left side
      coeff = h1+h2;
    } else if (j!=0 and j!=nx-1 and i==nx-1){ //right side
      coeff = h2;
    } else if (i!=0 and i!=nx-1 and j==0){ //bottom side
      coeff = h1+h2;
    } else if (i!=0 and i!=nx-1 and j==ny-1){ //top side
      coeff = h1;
    } else if (j==0 and i==0){ //bottom left corner
      coeff = h1+h2;
    } else if (j==0 and i==nx-1){ //bottom right corner
      coeff = h2;
    } else if (j==ny-1 and i==0){ //top left corner
      coeff = h1;
    } else if (j==ny-1 and i==nx-1){ //top right corner
      coeff = 0;
    }
  } else {
    coeff = h1+h2;
  }
  return coeff;
}

double MapPoisson::cipj(double xc, double yc, int i, int j){//ci+1,j
  //calculate coefficients
  Vector3d g = MapPoisson::gij(xc,yc+0.5*dyc);
  Vector3d gi = MapPoisson::ginv(g);
  double h1 = sqrt(MapPoisson::gdet(g))*gi(1)/4;
  g = MapPoisson::gij(xc+0.5*dxc,yc);
  gi = MapPoisson::ginv(g);
  double h2 = dyc*sqrt(MapPoisson::gdet(g))*gi(0)/dxc;
  g = MapPoisson::gij(xc,yc-0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h3 = -sqrt(MapPoisson::gdet(g))*gi(1)/4;
  //dirichlet correction
  double coeff = 0;
  if (xbctypel==0 and xbctypeu==0 and ybctypel==0 and ybctypeu==0){ //currently only imp for all dirich
    if (i!=0 and i!=nx-1 and j!=0 and j!=ny-1){ //normal
      coeff = h1+h2+h3;
    } else if (j!=0 and j!=ny-1 and i==0){ //left side
      coeff = h1+h2+h3;
    } else if (j!=0 and j!=nx-1 and i==nx-1){ //right side
      coeff = h1+h2+h3;
    } else if (i!=0 and i!=nx-1 and j==0){ //bottom side
      coeff = h1+h2;
    } else if (i!=0 and i!=nx-1 and j==ny-1){ //top side
      coeff = h2+h3;
    } else if (j==0 and i==0){ //bottom left corner
      coeff = h1+h2;
    } else if (j==0 and i==nx-1){ //bottom right corner
      coeff = h1+h2;
    } else if (j==ny-1 and i==0){ //top left corner
      coeff = h2+h3;
    } else if (j==ny-1 and i==nx-1){ //top right corner
      coeff = h2+h3;
    }
  } else {
    coeff = h1+h2+h3;
  }
  return coeff;
}

double MapPoisson::cipjm(double xc, double yc, int i, int j){//ci+1,j-1
  //calculate coefficients
  Vector3d g = MapPoisson::gij(xc+0.5*dxc,yc);
  Vector3d gi = MapPoisson::ginv(g);
  double h1 = -sqrt(MapPoisson::gdet(g))*gi(1)/4;
  g = MapPoisson::gij(xc,yc-0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h2 = -sqrt(MapPoisson::gdet(g))*gi(1)/4;
  //dirichlet correction
  double coeff = 0;
  if (xbctypel==0 and xbctypeu==0 and ybctypel==0 and ybctypeu==0){ //currently only imp for all dirich
    if (i!=0 and i!=nx-1 and j!=0 and j!=ny-1){ //normal
      coeff = h1+h2;
    } else if (j!=0 and j!=ny-1 and i==0){ //left side
      coeff = h1+h2;
    } else if (j!=0 and j!=nx-1 and i==nx-1){ //right side
      coeff = h2;
    } else if (i!=0 and i!=nx-1 and j==0){ //bottom side
      coeff = h1;
    } else if (i!=0 and i!=nx-1 and j==ny-1){ //top side
      coeff = h1+h2;
    } else if (j==0 and i==0){ //bottom left corner
      coeff = h1;
    } else if (j==0 and i==nx-1){ //bottom right corner
      coeff = 0;
    } else if (j==ny-1 and i==0){ //top left corner
      coeff = h1+h2;
    } else if (j==ny-1 and i==nx-1){ //top right corner
      coeff = h2;
    }
  } else {
    coeff = h1+h2;
  }
  return coeff;
}

double MapPoisson::cijm(double xc, double yc, int i, int j){//ci,j-1
  //calculate coefficients
  Vector3d g = MapPoisson::gij(xc-0.5*dxc,yc);
  Vector3d gi = MapPoisson::ginv(g);
  double h1 = sqrt(MapPoisson::gdet(g))*gi(1)/4;
  g = MapPoisson::gij(xc+0.5*dxc,yc);
  gi = MapPoisson::ginv(g);
  double h2 = -sqrt(MapPoisson::gdet(g))*gi(1)/4;
  g = MapPoisson::gij(xc,yc-0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h3 = dxc*sqrt(MapPoisson::gdet(g))*gi(2)/dyc;
  //dirichlet correction
  double coeff = 0;
  if (xbctypel==0 and xbctypeu==0 and ybctypel==0 and ybctypeu==0){ //currently only imp for all dirich
    if (i!=0 and i!=nx-1 and j!=0 and j!=ny-1){ //normal
      coeff = h1+h2+h3;
    } else if (j!=0 and j!=ny-1 and i==0){ //left side
      coeff = h2+h3;
    } else if (j!=0 and j!=nx-1 and i==nx-1){ //right side
      coeff = h1+h3;
    } else if (i!=0 and i!=nx-1 and j==0){ //bottom side
      coeff = h1+h2+h3;
    } else if (i!=0 and i!=nx-1 and j==ny-1){ //top side
      coeff = h1+h2+h3;
    } else if (j==0 and i==0){ //bottom left corner
      coeff = h2+h3;
    } else if (j==0 and i==nx-1){ //bottom right corner
      coeff = h1+h3;
    } else if (j==ny-1 and i==0){ //top left corner
      coeff = h2+h3;
    } else if (j==ny-1 and i==nx-1){ //top right corner
      coeff = h1+h3;
    }
  } else {
    coeff = h1+h2+h3;
  }
  return coeff;
}

double MapPoisson::cimjm(double xc, double yc, int i, int j){//ci-1,j-1
  //calculate coefficients
  Vector3d g = MapPoisson::gij(xc-0.5*dxc,yc);
  Vector3d gi = MapPoisson::ginv(g);
  double h1 = sqrt(MapPoisson::gdet(g))*gi(1)/4;
  g = MapPoisson::gij(xc,yc-0.5*dyc);
  gi = MapPoisson::ginv(g);
  double h2 = sqrt(MapPoisson::gdet(g))*gi(1)/4;
  //dirichlet correction
  double coeff = 0;
  if (xbctypel==0 and xbctypeu==0 and ybctypel==0 and ybctypeu==0){ //currently only imp for all dirich
    if (i!=0 and i!=nx-1 and j!=0 and j!=ny-1){ //normal
      coeff = h1+h2;
    } else if (j!=0 and j!=ny-1 and i==0){ //left side
      coeff = h2;
    } else if (j!=0 and j!=nx-1 and i==nx-1){ //right side
      coeff = h1+h2;
    } else if (i!=0 and i!=nx-1 and j==0){ //bottom side
      coeff = h1;
    } else if (i!=0 and i!=nx-1 and j==ny-1){ //top side
      coeff = h1+h2;
    } else if (j==0 and i==0){ //bottom left corner
      coeff = 0;
    } else if (j==0 and i==nx-1){ //bottom right corner
      coeff = h1;
    } else if (j==ny-1 and i==0){ //top left corner
      coeff = h2;
    } else if (j==ny-1 and i==nx-1){ //top right corner
      coeff = h1+h2;
    }
  } else {
    coeff = h1+h2;
  }
  return coeff;
}


///////////////////////////SOLVER & FACTORIZER//////////////////////////////////

//FACTORIZE LAPLACIAN
void MapPoisson::factorize()
{
  MapPoisson::laplace();
}


//GET SOLUTION
void MapPoisson::phisolve()
{
  //convert/calculate source
  VectorXd rhs = MapPoisson::srcconv();

  //compute solution
  phi = solver.solve(rhs);

  if(solver.info()!=Success) { //check if solving failed
    std::cout << "ERROR: SOLVING FAILURE" << std::endl;// solving failed
  }

}


//factorizer wrap
void wrap_factorize(MapPoisson *d){
  return d->factorize();
}

//solver wrap
void wrap_phisolve(MapPoisson *d){
  return d->phisolve();
}


///////////////////////SOURCE AND OUTPUT CONVERTERS/////////////////////////////

//source
void MapPoisson::getSrcvalatIJ(int i, int j, double sitrij){
  int k = (j)*nx+(i);
  jsource(k) = sitrij;
}

//solution
double MapPoisson::getSolvalatIJ(int i, int j){
  int k = (j)*nx+(i);
  return phi(k);
}

//source wrapper
void wrap_getSrcvalatIJ(MapPoisson *d, int i, int j, double sitrij){
  return d->getSrcvalatIJ(i, j, sitrij);
}

//solution wrapper
double wrap_getSolvalatIJ(MapPoisson *d, int i, int j){
  return d->getSolvalatIJ(i, j);
}
