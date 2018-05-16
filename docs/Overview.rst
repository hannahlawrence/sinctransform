Overview
=========================================

This is a C++ package to compute the sinc and sinc-squared transforms, defined as follows with respect to some input :math:`k_1,...,k_n \in \mathbf{R}` and :math:`q_1,...,q_n \in \mathbf{C}`:

.. math:

	\sum_{j=1}^m}q_j\text{sinc}(\mathbf{k_i}-\mathbf{k_j})

	\sum_{j=1}^m}q_j\text{sinc}^2(\mathbf{k_i}-\mathbf{k_j})

where we define


.. math:
	
	\text{sinc}(\mathbf{x})=\prod_{i=1}^r \frac{\text{sin}(x_i)}{x_i} \: \: \: \mathbf{x} \in \mathbf{R}^r



.. figure:: SincGraphBasic.png
    :width: 70%
    :align: center

    The sinc function in 1D

.. figure:: basicsincsqplot.png
    :width: 70%
    :align: center

    The sinc function in 2d