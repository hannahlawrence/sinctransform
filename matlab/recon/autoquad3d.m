function wvec=autoquad3d(klocs_d1,klocs_d2,klocs_d3)

% The ith weight w_i = 1 / sum sinc^2(delta_d1 * (klocs_d1_i-klocs_d1_j)) *
%                           j
%                              sinc^2(delta_d2 * (klocs_d2_i-klocs_d2_j)) *
%
%                              sinc^2(delta_d3 * (klocs_d3_i-klocs_d3_j))
%
% where sinc(x)=sin(pi * x) / (pi * x)

ifl=1;
q=ones(size(klocs_d1));
wvec=1./sincsq3d(ifl,klocs_d1,klocs_d2,klocs_d3,klocs_d1,klocs_d2,klocs_d3,q,1e-6,'legendre'); %lower, for speed
