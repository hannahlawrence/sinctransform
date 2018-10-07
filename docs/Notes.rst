Notes
=========================================

FINUFFT
---------
The current installation instructions do not use the parallel capability of FINUFFT, which would provide an additional speedup. To do so, follow the installation instructions of FINUFFT regarding the openmp library and the appropriate compilation flags.

Runtime
--------
The code presented here scales much better than direct calculation (one can experiment with the tests directory to verify this). However, note that the runtime is heavily dependent on the speed of the nonuniform Fourier transforms, which in turn depend on a variety of factors.

Precision
----------
The requested precision is not guaranteed, so it may be worth running a few test cases in the same regime (size and scale) as the intended application to pin down the optimal input precision. When necessary, it is also possible to change some of the default parameters, depending on choice of quadrature, to achieve better precision; see the `Quadrature`_ page for details.


