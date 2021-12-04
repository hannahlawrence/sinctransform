Package Contents
=========================================

The sinc and sinc-squared transform of specified dimension:

.. code::

	sinc1d.cpp
	sinc2d.cpp
	sinc3d.cpp
	sinctransform.hpp

To compute the sinc transform directly:

.. code::

	directsinc.cpp
	directsinc.hpp


Programs used in directsinc.cpp, as well as all example and testing code, to print arrays, generate random arrays, etc:

.. code::

	sincutil.cpp
	sincutil.hpp


For computing Gauss-Legendre quadrature nodes and weights (see [2] in citations):

.. code::

	fastGL.cpp
	fastGL.hpp

Simple usage examples for specified dimension:

.. code::
	
	example1d.cpp
	example2d.cpp
	example3d.cpp

Accuracy testing at precision range for specified dimension:

.. code::

	test1d.cpp
	test2d.cpp
	test3d.cpp

Matlab code (not C++ wrappers) for sinc transform:

.. code::

	sinc1d.m
	sinc2d.m
	sinc3d.m
	sincsq1d.m
	sincsq2d.m
	sincsq3d.m

Matlab file containing precomputed constants for trapezoidal constants:

.. code::

	newconstants.mat


Matlab to compute Gauss-Legendre nodes and weights, and to generate a 3d Shepp-Logan phantom (credit to Greg von Winckel and Matthias Schabel, respectively; see licenses):

.. code::
	
	lgwt.m
	phantom3d.m

As one application-based example, we also include annotated Matlab code to reconstruct an image from its nonuniformly sampled Fourier data, using the autoquad and sinc functions:

.. code::
	
	recon2d.m
	recon3d.m

These use the following to compute the "optimal" quadrature weights for reconstruction, as described in [1] of citations, via the sinc-squared programs:

.. code::

	autoquad1d.m
	autoquad2d.m
	autoquad3d.m




