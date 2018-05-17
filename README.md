# Sinc Transform

A C++ package for computing the fast sinc and sinc-squared transforms in 1, 2, and 3 dimensions according to the algorithm proposed by Greengard et. al. in 2006. The transforms are defined as follows: 

<img src="docs/eqns.png" width="350"/>

Please see the primary documentation [here](http://fast-sinc-transform.readthedocs.io/en/latest/) for installation instructions, example code, derivations, licenses, etc. Note that this code crucially requires the [FINUFFT](https://github.com/ahbarnett/finufft) library, as detailed in the installaiton instructions. Finally, matlab codes (following the same sinc algorithm) are included for convenience as well.
