local ffi = require "ffi"

ffi.cdef [[ 
/* Basis function identifiers */
enum gkyl_basis_type {
  GKYL_BASIS_MODAL_SERENDIPITY,
  GKYL_BASIS_MODAL_TENSOR,
  GKYL_BASIS_MODAL_HYBRID,
  GKYL_BASIS_MODAL_GKHYBRID,
};

/**
 * Basis function object
 */
struct gkyl_basis {
  unsigned ndim, poly_order, num_basis;
  char id[64]; // "serendipity", "tensor", "hybrid, "gkhybrid"
  enum gkyl_basis_type b_type; // identifier for basis function
    
/**
 * Evaluate basis in unit cell (i.e. a hypercube with each side
 * [-1,1])
 *
 * @param z Location to evaluate basis. z \in [-1,1]^n
 * @param b On output, value of basis at 'z'
 */
  void (*eval)(const double *z, double *b);

/**
 * Evaluate expansion at point in the logical cell (hypercube)
 *
 * @param z Location to evaluate exansion. z \in [-1,1]^n
 * @param f Expansion coefficients
 * @return Expansion evaluated at z
 */
  double (*eval_expand)(const double *z, const double *f);

/**
 * Evaluate gradient, given expansion at point in the logical cell
 * (hypercube)
 *
 * @param dir Direction to compute gradient
 * @param z Location to evaluate exansion. z \in [-1,1]^n
 * @param f Expansion coefficients
 * @return Expansion evaluated at z
 */
  double (*eval_grad_expand)(int dir, const double *z, const double *f);

/**
 * Flip-sign function: changes signs of input expansion cofficients by
 * changing sign of odd monomial powers in specified direction. So if
 * dir=0, all odd powers of x appearing in the expansion will have
 * sign flipped.
 *
 * @param dir Direction to flip sign
 * @param f Input expansion
 * @param fout On output, flipped version of @a f
 */
  void (*flip_odd_sign)(int dir, const double *f, double *fout);

  /**
 * Flip-sign function: changes signs of input expansion cofficients by
 * changing sign of even monomial powers in specified direction. So if dir=0, all
 * even powers of x appearing in the expansion will have sign flipped.
 *
 * @param dir Direction to flip sign
 * @param f Input expansion
 * @param fout On output, flipped version of @a f
 */
  void (*flip_even_sign)(int dir, const double *f, double *fout);

/**
 * Construct list of nodes corresponding to this basis set. The nodes
 * coordinates are in the unit cell [-1,1]^ndim and stored such that
 * the coodinates of a node are contiguous, starting at index ndim*n,
 * n = 0, ... num_basis-1. The node_coords array must be pre-allocated
 * by the caller.
 */
  void (*node_list)(double *node_coords);

/**
 * Given expansion coefficients on nodal basis (nodes specified by the
 * node_list method), compute modal expansion coefficients.
 *
 * @param fnodal Coefficients of nodal expansion
 * @param fmodal On output, coefficients of modal expansion
 */
  void (*nodal_to_modal)(const double *fnodal, double *fmodal);

/**
 * Given expansion coefficients on nodal basis defined by Gauss-Legendre
 * quadrature points, compute modal expansion coefficients.
 *
 * @param fquad Coefficients of nodal expansion in quadrature node basis
 * @param fmodal On output, coefficients of modal expansion
 * @param linc2 Modal component being updated (allows for parallelization over basis functions)
 */
  void (*quad_nodal_to_modal)(const double *fquad, double *fmodal, long linc2);  

/**
 * Given expansion coefficients of DG modal basis, evaluate basis at Gauss-Legendre
 * quadrature points of order p+1 (the Gauss-Legendre nodal basis).
 * 
 * @param fmodal Coefficients of modal expansion
 * @param fquad On output, evaluation of modal expansion at Gauss-Legendre quadrature nodal basis
 * @param linc2 Quadrature node being updated (allows for parallelization over quadrature points)
 */
  void (*modal_to_quad_nodal)(const double *fmodal, double *fquad, long linc2);  
};
]]
