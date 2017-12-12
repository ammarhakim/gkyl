The file modal-basis.mac has Maxima code to work with modal basis
functions. These basis functions can be constructed in arbitrary
dimensions and wirh arbitray order. Maximal order, Serendipity and
Tensor-product basis sets are supported.

The first step to creating a basis set is to define a list of
variable. For example, to create a 2D basis set one can do:

vars : [x,y]$

This will make all basis sets be functions of x and y. There are three
methods to create basis functions:

makeMaxOrderBasis(vars, n)
makeSerendipBasis(vars, n)
makeTensorBasis(vars, n)

Here n is the polyOrder for the desired basis set. These methods will
create the basis as a set of mononmials. Usually, for use in creating
algorithms one may want to orthonormalize the basis. For this, the
following method is provided:

gsOrthoNorm(vars, basis)

where basis is the list of functions returned by the above three
methods. The list returned by this method consists of orthonormal
functions.

Important methods:

# innerProd(vars, w, f1, f2)

Computes the inner-product <w f1 f2> over the variables vars.

# calcMassMatrix(vars, w, basis)

Compute the matrix M_ij = <w basis[i] basis[j]>. Note that if w=1 and
basis are orthonormal, the mass matrix will be an identity matrix.

# calcWeightedGradStiffMatrix(v, vars, w, basis)

Compute the matrix S_ij[v] = <w diff(basis[i],v) basis[j]>. The
variable v must be in the list vars.

# calcGradStiffMatrix(v, vars, basis)

This is calcWeightedGradStiffMatrix with w=1

# calcStiffMatrix(v, vars, basis)

Computes the matrix S_ij[v] = <diff(basis[i],v) diff(basis[j],v)>. To
use this in computing FEM stiffness matrix one needs a further sum
over all variables.

# calcCrossMassMatrix(varsC, basisC, varsP, basisP, w)

Computes mass matrix between two sets of basis functions:

M_ij = <w basisC[i] basisP[i]>

The inner product is defined over the variables varsP. For this to
work correctly, varsP must be a the same or a superset of varsC.