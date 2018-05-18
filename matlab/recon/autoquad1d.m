function wvec=autoquad1d(klocs)

% The ith weight w_i = 1 / sum sinc^2(delta * (klocs_i-klocs_j))
%                           j
% where sinc(x)=sin(pi * x) / (pi * x)

ifl=1;
q=ones(size(klocs));
wvec=1./sincsq1d(ifl,klocs,q,1e-16);