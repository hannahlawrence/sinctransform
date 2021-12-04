Citations and Licenses
=========================================

Main author of sinctransform: Hannah Lawrence.

sinctranform was written while Hannah was a summer intern at the
Numerical Algorithms Group, Center for Computational Biology,
Flatiron Institute, NYC.

Collaborators: Leslie Greengard, Alex Barnett, Jeremy Magland
(Flatiron Institute)


Relevant journal articles
-------------------------

Sinc transform:

[1] Greengard, L., Lee, J-Y., & Inati, S. (2006).
The fast sinc transform and image reconstruction from non-uniform samples in k-space,
Communications in Applied Mathematics and Computational Science, 1, 121-132.

Computation of Gauss-Legendre quadrature weights via fastGL.cpp and fastGL.hpp:

[2] Ignace Bogaert,
Iteration-free computation of Gauss-Legendre quadrature nodes and weights,
SIAM Journal on Scientific Computing, Volume 36, Number 3, 2014, pages A1008-1026.

Corrected trapezoidal quadrature weights:

[3] Sharad Kapur & Vladimir Rokhlin (1997).
High-Order Corrected Trapezoidal Quadrature Rules for Singular Functions, 
SIAM Journal of Numerical Analysis, Vol. 34, No. 4, pp. 1331–1356


FINUFFT library:

[4] (coming soon)

3d k-space sampling (for 3d reconstruction example):

[5] Park J, Shin T, Yoon SH, Goo JM, Park J-Y. A radial sampling strategy for uniform k-space coverage with retrospective respiratory gating in 3D ultrashort-echo-time lung imaging. NMR in biomedicine. 2016;29(5):576-587. doi:10.1002/nbm.3494.

Licenses
---------
Sinc and FINUFFT
~~~~~~~~~~~~~~~~~~~
Copyright (C) 2017-2018 The Simons Foundation, Inc. - All Rights Reserved.

Main author: Hannah Lawrence.

In collaboration with: Leslie Greengard, Alex Barnett, Jeremy Magland
(Flatiron Institute)

sinctransform is licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

C++ functions for Gauss-Legendre quadrature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Code based on [2], by J. Burkardt. Code distributed under the BSD license, which can be found at

	https://people.sc.fsu.edu/~jburkardt/txt/bsd_license.txt

3d phantom
~~~~~~~~~~~~~~~~~~~~~

From MathWorks File Exchange, via BSD license (see above)

Matlab functions for Gauss-Legendre quadrature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Copyright (c) 2009, Greg von Winckel 
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

