# tprod
Matlab library for simple tensor operations using Einstein summation convention semantics

This function computes a generalized multi-dimensional matrix product based upon the einstein summation convention (plus extras). This means given 2 n-d inputs:

X = [ A x B x C x E .... ]
Y = [ D x F x G x H .... ]

we define the result, Z, to be (in ESC)
Z_{c,e,d,f} = X_{a,b,c,e} Y_{d,f,a,b}

(N.B. This syntax can be used directly with the etprod wrapper script).

This translates into tprod syntax as:
Z = tprod(X,[-1 -2 1 2],Y,[3 4 -1 -2])

N.B. if Y==[], then it is assumed to be a copy of X.

This result is produced by forming an inner-product (multiply+sum) for the pairs of dimensions which have the same (negative) label (i.e. -1 => X dim 1 and Y dim 3, and -2 => X dim 2 and Y dim 4) and forming an outer product for all the remaining dimensions, with the positive label determining where this dimension is placed in the output.

# Download from Mathworks File-Exchange
https://www.mathworks.com/matlabcentral/fileexchange/16275-tprod-arbitary-tensor-products-between-n-d-arrays?s_tid=mwa_osa_a
