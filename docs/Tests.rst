Time and Performance Tests
=========================================

Test scripts in 1D, 2D, and 3D are included in the test directory. They compare running time of the sinc and sincsq programs to the brute force approach. In addition, accuracy is measured by the L2 error with respect to the result of the direct calculation. To enable larger-scale testing, on which the complete direct computation becomes intractable, the L2 error is approximated by the L2 error of a tractable subset of the requested outputs, as determined by the "numeval" variable.

To run:

.. code::
	
	make tests/test2d    
	tests/test2d 

It is also possible to make and run all tests simultaneously:

.. code::
	
	make tests

