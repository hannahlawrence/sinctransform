function wtrans= sincsq1d(ifl,a1,klocs,q,tol,mode)

if(nargin<1), test_sincsq1d; return; end

%  
% wtrans(j) = sum sinc^2(a1(k)-klocs(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% a1 = (real) evaluation locations
% klocs = (real) sample locations
% q = sample strengths
% tol = requested precision

% ensure vector inputs are correct size
q=q(:).';
klocs=klocs(:);

rkmax=max(bsxfun(@max,zeros(size(klocs)),abs(klocs)));
rkmax=max(rkmax,max(bsxfun(@max,zeros(size(a1)),abs(a1))));

if ifl==1
    rkmax=pi*rkmax;
    klocs=pi*klocs;
    a1=pi*a1;
end

newtol=max(tol,1e-16);

if isequal(mode,'legendre')
    rsamp=2; % increase to impose higher accuracy; will increase runtime
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
    wtrans=.25*finufft1d3(xx,h_at_xx.*ww.*trianglevec,1,newtol,a1);
else

rsamp=3; % increase to impose higher accuracy; will increase runtime
nx=ceil(rsamp*round(rkmax+3));
a=-2; b=0;
e=25; % increase (up to 60) to impose higher accuracy
if mod(nx,2)~=0; nx=nx+1; end % ensure even so that 0 is a quadrature point
n=nx; h=(b-a)/n;
xx=a-(e*h):h:b+(n+e)*h; xx=xx(:);
aind=e+1; zind=aind+n; bind=zind+n; 

leftvec=zeros(size(xx));
rightvec=zeros(size(xx));
trianglevec=zeros(size(xx));
for i=1:length(leftvec)
    val=xx(i);
    leftvec(i)=2+val;
    rightvec(i)=2-val;
    trianglevec(i)=2-abs(val);
end

load('newconstants.mat')
constants=constantcell{e};

ww_trap=zeros(size(xx));
ww_trap(aind)=0.5; ww_trap(bind)=0.5;
ww_trap(aind+1:bind-1)=1; 
ww_trap(aind:bind)=ww_trap(aind:bind).*trianglevec(aind:bind); % incorporate real vals
% above: ready for regular trapezoid rule
ww_left=zeros(size(xx)); %corrections from left side
ww_right=zeros(size(xx)); %corrections from right side
for k=1:e
    ww_left(aind-k) = ww_left(aind-k) - leftvec(aind-k)*constants(k);
    ww_left(aind+k) = ww_left(aind+k) + leftvec(aind+k)*constants(k);
    ww_left(zind-k) = ww_left(zind-k) + leftvec(zind-k)*constants(k);
    ww_left(zind+k) = ww_left(zind+k) - leftvec(zind+k)*constants(k);
end
for k=1:e
    ww_right(zind-k) =  ww_right(zind-k)- rightvec(zind-k)*constants(k);
    ww_right(zind+k) =  ww_right(zind+k)+ rightvec(zind+k)*constants(k);
    ww_right(bind-k) =  ww_right(bind-k)+ rightvec(bind-k)*constants(k);
    ww_right(bind+k) =  ww_right(bind+k)- rightvec(bind+k)*constants(k);
end
ww=ww_trap+ww_left+ww_right; % Use corrections
ww=h*ww;

% Transform sources and modes to proper form for finufft1d1 and finufft1d2
ms=length(xx); actual_unif_space=xx(2)-xx(1);
if max(abs(klocs))*actual_unif_space>pi || max(abs(a1))*actual_unif_space>pi
    fprintf('Error: cannot use finufft type 1/2 here; outside [-pi,pi]\n');
end
L=xx(1);
if mod(ms,2)==0
    DU1=-ms/2;
else
    DU1=(-ms+1)/2;
end
translation=(DU1-(L/actual_unif_space));
h_at_xx=finufft1d1(klocs*actual_unif_space,q.'.*exp(1i*actual_unif_space*klocs*translation),-1,newtol,ms); % -ms/2 <= k1 <= (ms-1)/2 
translation=((-1*DU1)*actual_unif_space)+L;
temp=finufft1d2(a1*actual_unif_space,1,newtol,h_at_xx.*ww);
wtrans=0.25*exp(1i*translation*a1).*temp;
end

% For real inputs, only return real values. Otherwise, return complex
% wtrans
if isreal(q)
    wtrans=real(wtrans);
end
wtrans=wtrans(:);
end

function test_sincsq1d
n=10000; ifl=1;
numeval=1000;
klocs=rand(n,1)*200;
a1=rand(size(klocs));
q=complex(rand(1,n)*30,rand(1,n)*30);
tic;correct_wtrans=superslowsincsq1d(ifl,a1(1:numeval),klocs,q);t3=toc;
precisions=[1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15];
for p=1:length(precisions)
    pr=precisions(p);
    tic;my_wtrans=sincsq1d(ifl,a1,klocs,q,pr,'legendre');t1=toc;
    tic;my_wtrans_new=sincsq1d(ifl,a1,klocs,q,pr,'trap');t2=toc;
    err=norm(correct_wtrans-my_wtrans(1:numeval),2)/norm(correct_wtrans,2);
    err2=norm(correct_wtrans-my_wtrans_new(1:numeval),2)/norm(correct_wtrans,2);
    fprintf("Requested: %g Error (Leg): %g (Trap): %g\n", pr, err,err2);
    fprintf("              Time  (Leg): %g s (Trap): %g s (Direct): %g s\n",t1,t2,t3);
end
end

function correct_wtrans=slowsincsq1d(ifl,a1,klocs,q)
    [a,b]=ndgrid(a1,klocs);
    repmat(q,length(q),1);
    if ifl==1
        sincmat=(sin(pi*(a-b))./(pi*(a-b))).^2;
    else
        sincmat=(sin(a-b)./(a-b)).^2;
    end
    sincmat(arrayfun(@isnan,sincmat))=1;
    correct_wtrans=sum(repmat(q,length(q),1).*sincmat,2);
end

function correct=superslowsincsq1d(ifl,a1,klocs_d1,q) 
% Alternative brute force calculation; even slower
% May be substituted in timing tests
    results=zeros(size(a1)); 
    for ind=1:length(results)
        sm=0;
        for j=1:length(klocs_d1)
            if ifl==0
                val=q(j)*(sinc(a1(ind)-klocs_d1(j)))^2;
            else
                val=q(j)*(sinc(pi*(a1(ind)-klocs_d1(j))))^2;
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


