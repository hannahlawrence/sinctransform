Overview
=========================================

This is a C++ package_ to compute the sinc and sinc-squared transforms (based on [1]), defined as follows. Consider :math:`2m` points in :math:`d` dimensional space, :math:`\mathbf{a_1},...,\mathbf{a_m},\mathbf{k_1},...,\mathbf{k_m} \in \mathbf{R}^d`, with weights :math:`q_1,...,q_m \in \mathbf{C}`. The task is to compute

.. math::

	\sum_{j=1}^m q_j\text{sinc}(\mathbf{a_i}-\mathbf{k_j}) \text{  or  } \sum_{j=1}^m q_j\text{sinc}^{2}(\mathbf{a_i}-\mathbf{k_j})   \qquad \mbox{ for all } i=1,\ldots,m

.. _package: https://github.com/hannahlawrence/sinctransform

where the :math:`d` dimensional sinc kernel is

.. math::
	
	\text{sinc}(\mathbf{x})=\prod_{i=1}^d \frac{\text{sin}(x_i)}{x_i} \: \: \: \mathbf{x} \in \mathbf{R}^d

Sometimes, :math:`\text{sinc}` is defined with an additional :math:`\pi` scaling, i.e. :math:`\text{sinc}(x)=\frac{\sin(\pi x)}{\pi x}`. Our code has an option to specify which convention you prefer.

The above task has naive complexity :math:`O(m^2)`. Our code uses a
fast algorithm to compute it in :math:`O(m \log 1/\epsilon + K^d \log K)` time,
where :math:`\epsilon>0` is the desired relative tolerance, and
:math:`K` is the maximal size of a cuboid containing the points (i.e. space-bandwidth product, or number of oscillations over such a cuboid).

      Our code is in C++, handles 1, 2, and 3 dimensions, as well as the :math:`\text{sinc}^2` kernel.
      It relies on the FINUFFT library to efficiently compute the type-3 nonuniform Fourier transform from the points to a set of quadrature nodes for the Fourier transform of the kernel.

For completeness, there is also some MATLAB code to perform the same functions. The latter is slightly slower and not as well-documented, but may be more convenient or easy to understand.

.. figure:: SincGraphBasic.png
    :width: 70%
    :align: center

    The sinc function in 1D

.. figure:: basicsincsqplot.png
    :width: 70%
    :align: center

    The sinc-squared function in 1D

