Overview
=========================================

This is a C++ package to compute the sinc and sinc-squared transforms, defined as follows with respect to some input :math:`\mathbf{k_1},...,\mathbf{k_n} \in \mathbf{R}^d` and :math:`q_1,...,q_n \in \mathbf{C}`:

.. math::

	\sum_{j=1}^m q_j\text{sinc}(\mathbf{k_i}-\mathbf{k_j}) \text{  or  } \sum_{j=1}^m q_j\text{sinc}^{2}(\mathbf{k_i}-\mathbf{k_j})

where we have

.. math::
	
	\text{sinc}(\mathbf{x})=\prod_{i=1}^d \frac{\text{sin}(x_i)}{x_i} \: \: \: \mathbf{x} \in \mathbf{R}^d

Sometimes, :math:`\text{sinc}` is defined with an additional :math:`\pi` coefficient, i.e. :math:`\text{sinc}(x)=\frac{\sin(\pi x)}{\pi x}`. There is always an option to specify which convention you prefer as one of the input arguments.

This code relies on the FINUFFT library to efficiently compute the nonuniform Fourier transform, and contains separate C++ code for 1, 2, and 3 dimensions. For completeness, there is also some Matlab code to perform the same functions. It is slightly slower and not as well-documented, but may be more convenient or easy to understand.

.. figure:: SincGraphBasic.png
    :width: 70%
    :align: center

    The sinc function in 1D

.. figure:: basicsincsqplot.png
    :width: 70%
    :align: center

    The sinc-squared function in 1D

