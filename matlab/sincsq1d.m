function wtrans= sincsq1d(ifl,klocs,q,tol)

if(nargin<1), test_sincsq1d; return; end

%  
% wtrans(j) = sum sinc^2(klocs(k)-klocs(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% klocs = (real) sample locations
% q = sample strengths
% tol = requested precision

rkmax=max(bsxfun(@max,zeros(size(klocs)),abs(klocs)));

if ifl==1
    rkmax=pi*rkmax;
    klocs=pi*klocs;
end

newtol=max(tol/1000,1e-16);
rsamp=2;
nx=ceil(rsamp*round(rkmax+3));
[xx_pre,ww_pre]=lgwt(nx,-1,1);

% quadrature points and weights to cover [-2,2]
xx=vertcat((xx_pre-1),(xx_pre+1));
ww=vertcat(ww_pre,ww_pre);

h_at_xx=finufft1d3(klocs,q,-1,newtol,xx);

trianglevec=zeros(size(ww));
for i=1:(2*nx)
    val=abs(xx(i));
    trianglevec(i)=2-val;
end

wtrans=.25*finufft1d3(xx,h_at_xx.*ww.*trianglevec,1,newtol,klocs);

% For real inputs, only return real values. Otherwise, return complex
% wtrans
if isreal(q)
    wtrans=real(wtrans);
end

function test_sincsq1d
n=10;
precisions=[1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15];
for p=1:length(precisions)
    pr=precisions(p);
    klocs=rand(n,1)*20;
    q=complex(rand(1,n)*30,rand(1,n)*30);
    ifl=1;
    correct_wtrans=slowsincsq1d(ifl,klocs,q);
    my_wtrans= sincsq1d(ifl,klocs,q,pr);
    err=abs(correct_wtrans-my_wtrans);
    err=sqrt(err.'*err);
    fprintf("Requested: %g Error: %g\n", pr, err);
end

function correct_wtrans=slowsincsq1d(ifl,klocs,q)
    [a,b]=ndgrid(klocs,klocs);
    repmat(q,length(q),1);
    if ifl==1
        sincmat=(sin(pi*(a-b))./(pi*(a-b))).^2;
    else
        sincmat=(sin(a-b)./(a-b)).^2;
    end
    sincmat(arrayfun(@isnan,sincmat))=1;
    correct_wtrans=sum(repmat(q,length(q),1).*sincmat,2);




