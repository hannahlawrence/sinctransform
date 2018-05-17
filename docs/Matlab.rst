Matlab
=========================================

The test codes are built into the matlab files themselves, in the form of self-contained test and direct calculation functions near the end of each file. To run, simply call the function without any input arguments. The test functions also demonstrate appropriate usage, which is more straightforward than in C++. The input vectors can be modified to determine what input precision is appropriate for a particular regime. 

Note that the programs are self-contained -- they are not Matlab interfaces to the C++, but pure Matlab code. As such, the only setup step is to ensure the FINUFFT matlab interfaces are on the matlab path.

To use the FINUFFT Matlab interfaces, follow their documentation_. Then, use Matlab's "addpath" to ensure the FINUFFT programs (e.g. finufft2d3) will be found. 

.. _documentation:
	
	http://finufft.readthedocs.io/en/latest/index.html

Finally, make sure lgwt.m is in the same directory, or that its directory is also included in Matlab's path.
