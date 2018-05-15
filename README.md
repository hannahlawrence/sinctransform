# Sinc Transform

A C++ package for computing the sinc and sinc-squared transforms (as described by [1] in Citations) in 1, 2, and 3 dimensions. **Add note here about speed! (See individual files for more detailed documentation.)

**Links to readthedocs here:
[to be inserted]

[interim documentation below, possibly keep?]

### Installation

1. Install the finufft library (https://github.com/ahbarnett/finufft)
2. Download contents of sinctransform package 
3. In the makefile, change CURRENT to the folder containing the library: CURRENT = /some/path/to/sinctransform
4. In the makefile, change FINUFFT to the path to finufft: FINUFFT = /some/path/to/finufft

### Example and Test Programs
To run one of the simple example programs:
```
make examples/example1d      .
examples/example1d           
```
To run one of the existing time and performance tests:  
```
make tests/test1d        
tests/test1d           
```

To use in other programs, following the example code: 
```
make libsinc.a
```
then compile with both sinctransform and finufft static libraries, e.g.:
```
g++ -std=c++11 -Wall -g -o myprogram myprogram.cpp libsinc.a /some/path/to/finufft/lib/libfinufft.a -lfftw3 -lm
```
It may be necessary to include a flag telling the compiler where to find the FFTW library, which is a prerequisite for the finufft library. To do so, add the flag "-L/some/dir" such that /some/dir contains the static FFTW library (.a file).

# Contents

The sinc and sinc-squared transform of specified dimension:

```
sinc1d.cpp
sinc2d.cpp
sinc3d.cpp
sinctransform.hpp
```
To compute the sinc transform directly:

```
directsinc.cpp
directsinc.hpp
```

Programs used in directsinc.cpp, as well as all example and testing code, to print arrays, generate random arrays, etc:
```
sincutil.cpp
sincutil.hpp
```

For computing Gauss-Legendre quadrature weights (see [2] in citations):
```
fastGL.cpp
fastGL.hpp
```

Simple usage examples for specified dimension:
```	
example1d.cpp
example2d.cpp
example3d.cpp
```

Accuracy testing at precision range for specified dimension:
```
test1d.cpp
test2d.cpp
test3d.cpp
```


Matlab code (not wrappers) for sinc transform:
```
sinc1d.m
sinc2d.m
sinc3d.m
sincsq1d.m
sincsq2d.m
sincsq3d.m
```

# Citations

Sinc transform:

[1] Greengard, L., Lee, J-Y., & Inati, S. (2006).
The fast sinc transform and image reconstruction from non-uniform samples in k-space,
Communications in Applied Mathematics and Computational Science, 1, 121-132.

Computation of Gauss-Legendre quadrature weights via fastGL.cpp and fastGL.hpp:

[2] Ignace Bogaert,
Iteration-free computation of Gauss-Legendre quadrature nodes and weights,
SIAM Journal on Scientific Computing, Volume 36, Number 3, 2014, pages A1008-1026.

