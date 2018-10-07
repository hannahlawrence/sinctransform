function wtrans=sinc2d(ifl,a1,a2,klocs_d1,klocs_d2,q,tol,mode)

if(nargin<1), test_sinc2d; return; end

%  
% wtrans(j) = sum sinc(a1(k)-klocs_d1(j)) * sinc(a2(k)-klocs_d2(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% a1 = (real) evaluation locations in dimension 1
% a2 = (real) evaluation locations in dimension 2
% klocs_d1 = (real) sample locations in dimension 1
% klocs_d2 = (real) sample locations in dimension 2
% q = sample strengths
% tol = requested precision

% ensure vector inputs are correct size
q=q(:).';
klocs_d1=klocs_d1(:);
klocs_d2=klocs_d2(:);

rkmaxx=max(bsxfun(@max,zeros(size(klocs_d1)),abs(klocs_d1)));
rkmaxx=max(rkmaxx,max(bsxfun(@max,zeros(size(a1)),abs(a1))));
rkmaxy=max(bsxfun(@max,zeros(size(klocs_d2)),abs(klocs_d2)));
rkmaxy=max(rkmaxy,max(bsxfun(@max,zeros(size(a2)),abs(a2))));

if ifl==1
    rkmaxx=pi*rkmaxx;
    rkmaxy=pi*rkmaxy;
    klocs_d1=pi*klocs_d1;
    klocs_d2=pi*klocs_d2;
    a1=pi*a1;
    a2=pi*a2;
end


newtol=max(tol,1e-16);

if isequal(mode,'legendre')
    rsamp=2; % increase to impose higher accuracy; will increase runtime
    nx=ceil(rsamp*round(rkmaxx+3));
    ny=ceil(rsamp*round(rkmaxy+3));
    [xx,wwx]=lgwt(nx,-1,1);
    [yy,wwy]=lgwt(ny,-1,1);
    % create 2D grid of points and corresponding weights
    [c,d]=ndgrid(xx,yy);
    allxx=c(:);
    allyy=d(:);
    [e,f]=ndgrid(wwx,wwy);
    allww=e(:).*f(:);
    h_at_xxyy=finufft2d3(klocs_d1,klocs_d2,q,-1,newtol,allxx,allyy);
    wtrans=0.25*finufft2d3(allxx,allyy,h_at_xxyy.*allww,1,newtol,a1,a2);
else
rsamp=3; % increase to impose higher accuracy; will increase runtime
nx=ceil(rsamp*round(rkmaxx+3));
ny=ceil(rsamp*round(rkmaxy+3));
a=-1; b=1; 
e=25; % increase (up to 60) to impose higher accuracy
load('newconstants.mat')
constants=constantcell{e};
if mod(nx,2)~=0; nx=nx+1; end % ensure even so that 0 is a quadrature point
n=nx; h=(b-a)/n;
aind=e+1; bind=aind+n;
xx=a-(e*h):h:a+(n+e)*h; xx=xx(:);
ww_x=zeros(length(xx),1);
ww_x(aind)=0.5; ww_x(bind)=0.5; 
ww_x(aind+1:bind-1)=1;
for k=1:e
    ww_x(aind-k) = ww_x(aind-k) - constants(k);
    ww_x(aind+k) = ww_x(aind+k) + constants(k);
    ww_x(bind-k) = ww_x(bind-k) + constants(k);
    ww_x(bind+k) = ww_x(bind+k) - constants(k);
end
ww_x=ww_x*h;

if mod(ny,2)~=0; ny=ny+1; end % ensure even so that 0 is a quadrature point
n=ny; h=(b-a)/n;
aind=e+1; bind=aind+n;
yy=a-(e*h):h:a+(n+e)*h; yy=yy(:);
ww_y=zeros(length(yy),1);
ww_y(aind)=0.5; ww_y(bind)=0.5; 
ww_y(aind+1:bind-1)=1;
for k=1:e
    ww_y(aind-k) = ww_y(aind-k) - constants(k);
    ww_y(aind+k) = ww_y(aind+k) + constants(k);
    ww_y(bind-k) = ww_y(bind-k) + constants(k);
    ww_y(bind+k) = ww_y(bind+k) - constants(k);
end
ww_y=ww_y*h;

% create 2D grid of points and corresponding weights
[c,d]=ndgrid(xx,yy);
allxx=c(:);
allyy=d(:);
[e,f]=ndgrid(ww_x,ww_y);
allww=e(:).*f(:);

% Transform sources and modes to proper form for finufft2d1 and finufft2d2
msx=length(xx); msy=length(yy);
actual_unif_spacex=xx(2)-xx(1); actual_unif_spacey=yy(2)-yy(1);
if max(abs(klocs_d1))*actual_unif_spacex>pi || max(abs(a1))*actual_unif_spacex>pi || max(abs(klocs_d2))*actual_unif_spacey>pi || max(abs(a2))*actual_unif_spacey>pi
    fprintf('Error: cannot use finufft type 1/2 here; outside [-pi,pi]\n');
end
Lx=xx(1); Ly=yy(1);
if mod(msx,2)==0
    DU1x=-msx/2;
else
    DU1x=(-msx+1)/2;
end
if mod(msy,2)==0
    DU1y=-msy/2;
else
    DU1y=(-msy+1)/2;
end
translationx=(DU1x-(Lx/actual_unif_spacex));
translationy=(DU1y-(Ly/actual_unif_spacey));

h_at_xxyy=finufft2d1(klocs_d1*actual_unif_spacex,klocs_d2*actual_unif_spacey,q.'.*exp(1i*((actual_unif_spacex*klocs_d1*translationx) + (actual_unif_spacey*klocs_d2*translationy))),-1,newtol,msx,msy);
h_at_xxyy=h_at_xxyy(:);
translationx=((-1*DU1x)*actual_unif_spacex)+Lx;
translationy=((-1*DU1y)*actual_unif_spacey)+Ly;
temp=finufft2d2(a1*actual_unif_spacex,a2*actual_unif_spacey,1,newtol,reshape(h_at_xxyy.*allww,msx,msy));
wtrans=0.25*exp(1i*(translationx*a1+translationy*a2)).*temp;
end

% For real inputs, only return real values. Otherwise, return complex
if isreal(q)
    wtrans=real(wtrans);
end
wtrans=wtrans(:);
end

function test_sinc2d
n=1000; ifl=1;
numeval=100;
klocs_d1=-10+(20*rand(n,1));
klocs_d2=-10+(20*rand(n,1));
a1=-10+(200*rand(n,1));
a2=-10+(200*rand(n,1));
q=complex(rand(1,n)*30,rand(1,n)*30);
tic;correct=superslowsinc2d(ifl,a1(1:numeval),a2(1:numeval),klocs_d1,klocs_d2,q); tt3=toc;
precisions=[1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15];
for p=1:length(precisions)
    pr=precisions(p);
    tic;myresult=sinc2d(ifl,a1,a2,klocs_d1,klocs_d2,q,pr,'legendre');tt=toc;
    tic;myresult2=sinc2d(ifl,a1,a2,klocs_d1,klocs_d2,q,pr,'trap');tt2=toc;
    err=norm(correct-myresult(1:numeval),2)/norm(correct,2);
    err2=norm(correct-myresult2(1:numeval),2)/norm(correct,2);
    fprintf("Requested: %g Error (Leg): %g Error (Trap): %g\n", pr,err,err2);
    fprintf("              Time  (Leg): %g s (Trap): %g s (Direct): %g s\n",tt,tt2,tt3);
end
end

function correct = slowsinc2d(ifl,a1,a2,klocs_d1,klocs_d2,q)
    q=q(:).';
    [a1,b1]=ndgrid(a1,klocs_d1);
    [a2,b2]=ndgrid(a2,klocs_d2);
    if ifl==1
        x=sin(pi*(a1-b1))./(pi*(a1-b1));
        y=sin(pi*(a2-b2))./(pi*(a2-b2));
    else
        x=sin(a1-b1)./(a1-b1);
        y=sin(a2-b2)./(a2-b2);
    end
    x(arrayfun(@isnan,x))=1;
    y(arrayfun(@isnan,y))=1;
    sincmat=x.*y;
    correct=sum(repmat(q,length(klocs_d1),1).*sincmat,2);
end

function correct=superslowsinc2d(ifl,a1,a2,klocs_d1,klocs_d2,q) % wrong here?
% Alternative brute force calculation; even slower
% May be substituted in timing tests    
    results=zeros(size(a1)); 
    for ind=1:length(results)
        sm=0;
        for j=1:length(klocs_d1)
            if ifl==0
                val=q(j)*(sinc(a1(ind)-klocs_d1(j)))*(sinc(a2(ind)-klocs_d2(j)));
            else
                val=q(j)*(sinc(pi*(a1(ind)-klocs_d1(j))))*(sinc(pi*(a2(ind)-klocs_d2(j))));
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




