function wtrans= sinc1d(ifl,a1,klocs,q,tol,mode)

if(nargin<1), test_sinc1d; return; end

%  
% wtrans(j) = sum sinc(a1(k)-klocs(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% a1 = (real) evaluation locations
% klocs = (real) sample locations
% q = sample strengths
% tol = requested precision
% mode = quadrature method: `legendre` or `trapezoid`

% ensure vector inputs are column vectors
q=q(:);
klocs=klocs(:);

% more aggressive tolerance to compensate for repeated calls and quadrature
newtol=max(tol/1000,1e-16);

rkmax=max(bsxfun(@max,zeros(size(klocs)),abs(klocs)));
rkmax=max(rkmax,max(bsxfun(@max,zeros(size(a1)),abs(a1))));

if ifl==1
    rkmax=pi*rkmax;
    klocs=pi*klocs;
    a1=pi*a1;
end
rsamp=2; % increase to impose higher accuracy; will increase runtime
nx=ceil(rsamp*round(rkmax+3)); 

if isequal(mode,'legendre')
    [xx,ww]=lgwt(nx,-1,1);
    h_at_xx=finufft1d3(klocs,q,-1,newtol,xx);
    % factor of 1/2 required for change of integration bounds
    wtrans=0.5*finufft1d3(xx,h_at_xx.*ww,1,newtol,a1);
else

a=-1; b=1;
e=21; % increase (up to 60) to impose higher accuracy; will increase runtime
if mod(nx,2)~=0; nx=nx+1; end % ensure even so that 0 is a quadrature point
n=nx; h=(b-a)/n; aind=e+1; bind=aind+n;
load('newconstants.mat') % Precomputed corrections to trapezoidal rule
constants=constantcell{e};
xx=a-(e*h):h:a+(n+e)*h;
xx=xx(:);
ww=zeros(length(xx),1);
ww(aind)=0.5; ww(bind)=0.5; 
ww(aind+1:bind-1)=1;
for k=1:e
    ww(aind-k) = ww(aind-k) - constants(k);
    ww(aind+k) = ww(aind+k) + constants(k);
    ww(bind-k) = ww(bind-k) + constants(k);
    ww(bind+k) = ww(bind+k) - constants(k);
end
ww=ww*h;

% Transform sources and modes to proper form for finufft1d1 and finufft1d2
ms=length(xx); actual_unif_space=xx(2)-xx(1); L=xx(1);
if max(abs(klocs))*actual_unif_space>pi || max(abs(a1))*actual_unif_space>pi
    fprintf('Error: cannot use finufft type 1/2 here; outside [-pi,pi]\n');
end
if mod(ms,2)==0
    DU1=-ms/2;
else
    DU1=(-ms+1)/2;
end
translation=(DU1-(L/actual_unif_space));
h_at_xx=finufft1d1(klocs*actual_unif_space,q.*exp(1i*actual_unif_space*klocs*translation),-1,newtol,ms); % -ms/2 <= k1 <= (ms-1)/2 
translation=((-1*DU1)*actual_unif_space)+L;
temp=finufft1d2(a1*actual_unif_space,1,newtol,h_at_xx.*ww);
wtrans=0.5*exp(1i*translation*a1).*temp;
end

% For real inputs, only return real values. Otherwise, return complex
if isreal(q)
    wtrans=real(wtrans);
end
end

function test_sinc1d
n=1000; ifl=1;
klocs=(-10)+(20*rand(n,1));
q=complex(rand(1,n)*10,10*rand(1,n));
a1=100*rand(size(klocs));
tic;correct_wtrans=slowsinc1d(ifl,a1,klocs,q);t3=toc;
precisions=[1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15];
for p=1:length(precisions)
    pr=precisions(p);
    tic;my_wtrans= sinc1d(ifl,a1,klocs,q,pr,'legendre');t1=toc;
    tic;my_wtrans_new=sinc1d(ifl,a1,klocs,q,pr,'trap');t2=toc;
    if(size(correct_wtrans,1)~=size(my_wtrans,1))
        correct_wtrans=correct_wtrans.';
    end
    err=norm(correct_wtrans-my_wtrans,2);
    err2=norm(correct_wtrans-my_wtrans_new,2);
    fprintf("Requested: %g Error (Leg): %g (Trap): %g\n", pr, err,err2);
    fprintf("              Time  (Leg): %g s (Trap): %g s (Direct): %g s\n",t1,t2,t3);
end
end

function correct_wtrans=slowsinc1d(ifl,a1,klocs,q)  % wrong????
    [a,b]=ndgrid(a1,klocs);
    if ifl==1
        sincmat=sin(pi*(a-b))./(pi*(a-b));
    else
        sincmat=sin(a-b)./(a-b);
    end
    sincmat(arrayfun(@isnan,sincmat))=1;
    correct_wtrans=sum(repmat(q,length(q),1).*sincmat,2);
end

function correct=superslowsinc1d(ifl,a1,klocs_d1,q)
% Alternative brute force calculation; even slower
% May be substituted in timing tests
    results=zeros(size(a1)); 
    for ind=1:length(results)
        sm=0;
        for j=1:length(klocs_d1)
            if ifl==0
                val=q(j)*(sinc(a1(ind)-klocs_d1(j)));
            else
                val=q(j)*(sinc(pi*(a1(ind)-klocs_d1(j))));
            end
            sm=sm+val;
        end
        results(ind)=sm;
    end
    correct=results;
end

function val=sinc(x)
    if x==0
        val=1;
    else
        val=sin(x)/x;
    end
end