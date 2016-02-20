# `dft`

A simple Matlab implementation of Hartree-Fock and Kohn-Sham calculations.

Dependencies:

* Maple/Fortran framework [dfauto](http://dx.doi.org/10.1016/S0010-4655(01)00148-5) for easy implementation of exchange-correlation functionals
* Auto-generated C++ library [libint](https://github.com/evaleev/libint) for 2-electron integrals
* Fortran library [gen1int](https://gitlab.com/bingao/gen1int) for 1-electron integrals
* Matlab MEX library [tprod](http://www.mathworks.com/matlabcentral/fileexchange/16275-tprod-arbitary-tensor-products-between-n-d-arrays) for fast and simple tensor products

All these dependencies are automatically downloaded and compiled during build of `dft`. The tprod package is rather dated and does not seem to work on modern OS X installations out of the box.

For installation, run

``` bash
make
```

For running two simple test calculations,

``` bash
make test
```