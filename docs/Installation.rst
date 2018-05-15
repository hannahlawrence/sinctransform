Installation and Usage
=========================================

1. Install the finufft_ library, including all of its required libraries
2. Download contents of sinctransform package, i.e. "git clone https://github.com/hannahlawrence/sinctransform.git"

.. _finufft: https://github.com/ahbarnett/finufft

Configuration
---------------

1. In the makefile, change CURRENT to the folder containing the library: 

.. code::

	CURRENT = /some/path/to/sinctransform

2. In the makefile, change FINUFFT to the path to finufft: 

.. code::
	
	FINUFFT = /some/path/to/finufft

3. In the makefile, change FFTW to the directory containing the fftw library: 

.. code::

	FFTW = /some/path/to/fft 

(Note: on a Mac, this may be /usr/local/lib)

Usage
----------------

To use in other programs, begin by creating the static library:

.. code::

	make libsinc.a

Then compile with both the sinctransform and finufft static libraries, e.g.

.. code::

	g++ -std=c++11 -Wall -g -o myprog myprog.cpp /path/to/libsinc.a /path/to/finufft/lib/libfinufft.a -lfftw3 -lm

It may be necessary to include a flag telling the compiler where to find the FFTW library, which is a prerequisite for the finufft library. To do so, add the flag "-L/some/dir" such that /some/dir contains the static FFTW library (.a file).


