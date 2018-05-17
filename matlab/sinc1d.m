function wtrans= sinc1d(ifl,klocs,q,tol)

if(nargin<1), test_sinc1d; return; end

%  
% wtrans(j) = sum sinc(klocs(k)-klocs(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% klocs = (real) sample locations
% q = sample strengths
% tol = requested precision

% ensure vector inputs are column vectors
q=q(:);
klocs=klocs(:);

% more aggressive tolerance to compensate for repeated calls and quadrature
newtol=max(tol/1000,1e-16);

rkmax=max(bsxfun(@max,zeros(size(klocs)),abs(klocs)));

if ifl==1
    rkmax=pi*rkmax;
    klocs=pi*klocs;
end
rsamp=2;
nx=ceil(rsamp*round(rkmax+3)); % following fortran code

[xx,ww]=lgwt(nx,-1,1);
h_at_xx=finufft1d3(klocs,q,-1,newtol,xx);

% factor of 1/2 required for change of integration bounds
wtrans=0.5*finufft1d3(xx,h_at_xx.*ww,1,newtol,klocs);

% For real inputs, only return real values. Otherwise, return complex
% wtrans
if isreal(q)
    wtrans=real(wtrans);
end

function test_sinc1d
n=100;
precisions=[1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15];
for p=1:length(precisions)
    klocs=(-10)+(20*rand(n,1));
    q=complex(rand(1,n)*10,10*rand(1,n));
    pr=precisions(p);
    ifl=1;
    my_wtrans= sinc1d(ifl,klocs,q,pr);
    correct_wtrans=slowsinc1d(ifl,klocs,q);
    if(size(correct_wtrans,1)~=size(my_wtrans,1))
        correct_wtrans=correct_wtrans.';
    end
    err=abs(correct_wtrans-my_wtrans);
    err=sqrt(err.'*err);
    fprintf("Requested: %g Error: %g\n", pr, err);
end

function correct_wtrans=slowsinc1d(ifl,klocs,q) 
    [a,b]=ndgrid(klocs,klocs);
    if ifl==1
        sincmat=sin(pi*(a-b))./(pi*(a-b));
    else
        sincmat=sin(a-b)./(a-b);
    end
    sincmat(arrayfun(@isnan,sincmat))=1;
    correct_wtrans=sum(repmat(q,length(q),1).*sincmat,2);



