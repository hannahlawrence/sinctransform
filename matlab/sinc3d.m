function wtrans=sinc3d(ifl,a1,a2,a3,klocs_d1,klocs_d2,klocs_d3,q,tol,mode)

if(nargin<1), test_sinc3d; return; end

%  
% wtrans(j) = sum sinc(a1(k)-klocs_d1(j)) * sinc(a2(k)-klocs_d2(j)) * sinc(a3(k)-klocs_d3(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% a1 = (real) evaluation locations in dimension 1
% a2 = (real) evaluation locations in dimension 2
% a3 = (real) evaluation locations in dimension 3
% klocs_d1 = (real) sample locations in dimension 1
% klocs_d2 = (real) sample locations in dimension 2
% klocs_d3 = (real) sample locations in dimension 3
% q = sample strengths
% tol = requested precision

% ensure vector inputs are correct size
q=q(:).';
klocs_d1=klocs_d1(:);
klocs_d2=klocs_d2(:);
klocs_d3=klocs_d3(:);

newtol=max(tol,1e-16);

rkmaxx=max(bsxfun(@max,zeros(size(klocs_d1)),abs(klocs_d1)));
rkmaxx=max(rkmaxx,max(bsxfun(@max,zeros(size(a1)),abs(a1))));
rkmaxy=max(bsxfun(@max,zeros(size(klocs_d2)),abs(klocs_d2)));
rkmaxy=max(rkmaxy,max(bsxfun(@max,zeros(size(a2)),abs(a2))));
rkmaxz=max(bsxfun(@max,zeros(size(klocs_d3)),abs(klocs_d3)));
rkmaxz=max(rkmaxz,max(bsxfun(@max,zeros(size(a3)),abs(a3))));

if ifl==1
    rkmaxx=pi*rkmaxx;
    rkmaxy=pi*rkmaxy;
    rkmaxz=pi*rkmaxz;
    klocs_d1=pi*klocs_d1;
    klocs_d2=pi*klocs_d2;
    klocs_d3=pi*klocs_d3;
    a1=pi*a1;
    a2=pi*a2;
    a3=pi*a3;
end
rsamp=3; % increase to impose higher accuracy; will increase runtime
nx=ceil(rsamp*round(rkmaxx+3)); 
ny=ceil(rsamp*round(rkmaxy+3));
nz=ceil(rsamp*round(rkmaxz+3));

if isequal(mode,'legendre')
    [xx,wwx]=lgwt(nx,-1,1); 
    [yy,wwy]=lgwt(ny,-1,1);
    [zz,wwz]=lgwt(nz,-1,1);
    [c,d,e]=ndgrid(xx,yy,zz);
    allxx=c(:);
    allyy=d(:);
    allzz=e(:);
    [f,g,h]=ndgrid(wwx,wwy,wwz);
    allww=f(:).*g(:).*h(:);
    h_at_xxyyzz=finufft3d3(klocs_d1,klocs_d2,klocs_d3,q,-1,newtol,allxx,allyy,allzz);
    wtrans=(1/8)*finufft3d3(allxx,allyy,allzz,h_at_xxyyzz.*allww,1,newtol,a1,a2,a3);
else

a=-1; b=1;
%e=min(round(21+(-1*log10(newtol)-2)*25/14),60); % increase (up to 60) to impose higher accuracy; will increase runtime
e=25; % increase (up to 60) to impose higher accuracy
load('newconstants.mat'); 
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

if mod(nz,2)~=0; nz=nz+1; end % ensure even so that 0 is a quadrature point
n=nz; h=(b-a)/n;
aind=e+1; bind=aind+n;
zz=a-(e*h):h:a+(n+e)*h; zz=zz(:);
ww_z=zeros(length(zz),1);
ww_z(aind)=0.5; ww_z(bind)=0.5; 
ww_z(aind+1:bind-1)=1;
for k=1:e
    ww_z(aind-k) = ww_z(aind-k) - constants(k);
    ww_z(aind+k) = ww_z(aind+k) + constants(k);
    ww_z(bind-k) = ww_z(bind-k) + constants(k);
    ww_z(bind+k) = ww_z(bind+k) - constants(k);
end
ww_z=ww_z*h;

[c,d,e]=ndgrid(xx,yy,zz);
allxx=c(:);
allyy=d(:);
allzz=e(:);
[f,g,h]=ndgrid(ww_x,ww_y,ww_z);
allww=f(:).*g(:).*h(:);

% Transform sources and modes to proper form for finufft3d1 and finufft3d2
msx=length(xx); msy=length(yy); msz=length(zz);
actual_unif_spacex=xx(2)-xx(1); 
actual_unif_spacey=yy(2)-yy(1);
actual_unif_spacez=zz(2)-zz(1);
if max(abs(klocs_d1))*actual_unif_spacex>pi || max(abs(a1))*actual_unif_spacex>pi || max(abs(klocs_d2))*actual_unif_spacey>pi || max(abs(a2))*actual_unif_spacey>pi || max(abs(klocs_d3))*actual_unif_spacey>pi || max(abs(a3))*actual_unif_spacey>pi
    fprintf('Error: cannot use finufft type 1/2 here; outside [-pi,pi]\n');
end
Lx=xx(1); Ly=yy(1); Lz=zz(1);
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
if mod(msz,2)==0
    DU1z=-msz/2;
else
    DU1z=(-msz+1)/2;
end
translationx=(DU1x-(Lx/actual_unif_spacex));
translationy=(DU1y-(Ly/actual_unif_spacey));
translationz=(DU1z-(Lz/actual_unif_spacez));

h_at_xxyyzz=finufft3d1(klocs_d1*actual_unif_spacex,klocs_d2*actual_unif_spacey,klocs_d3*actual_unif_spacez,q.'.*exp(1i*((actual_unif_spacex*klocs_d1*translationx) + (actual_unif_spacey*klocs_d2*translationy) + (actual_unif_spacez*klocs_d3*translationz))),-1,newtol,msx,msy,msz);
h_at_xxyyzz=h_at_xxyyzz(:);

translationx=((-1*DU1x)*actual_unif_spacex)+Lx;
translationy=((-1*DU1y)*actual_unif_spacey)+Ly;
translationz=((-1*DU1z)*actual_unif_spacez)+Lz;
temp=finufft3d2(a1*actual_unif_spacex,a2*actual_unif_spacey,a3*actual_unif_spacez,1,newtol,reshape(h_at_xxyyzz.*allww,msx,msy,msz));
wtrans=(1/8)*exp(1i*(translationx*a1+translationy*a2+translationz*a3)).*temp;

% For real inputs, only return real values. Otherwise, return complex
if isreal(q)
    wtrans=real(wtrans);
end
wtrans=wtrans(:);

end

% For real inputs, only return real values. Otherwise, return complex
% wtrans
if isreal(q)
    wtrans=real(wtrans);
end
end

function test_sinc3d
n=10000; ifl=1;
numeval=100;
precisions=[1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15];
precisions=precisions(1:3:length(precisions));
klocs_d1=-10+(50*rand(n,1));
klocs_d2=-10+(20*rand(n,1));
klocs_d3=-10+(20*rand(n,1));
q=complex(rand(1,n)*30,rand(1,n)*30); 
a1=-10+(20*rand(n,1));
a2=-10+(20*rand(n,1));
a3=-10+(20*rand(n,1));
tic;correct=superslowsinc3d(ifl,a1(1:numeval),a2(1:numeval),a3(1:numeval),klocs_d1,klocs_d2,klocs_d3,q);t3=toc;
for p=1:length(precisions)
    pr=precisions(p);
    tic;myresult=sinc3d(ifl,a1,a2,a3,klocs_d1,klocs_d2,klocs_d3,q,pr,'legendre');t1=toc;
    tic;myresult2=sinc3d(ifl,a1,a2,a3,klocs_d1,klocs_d2,klocs_d3,q,pr,'trap');t2=toc;
    err1=norm(correct-myresult(1:numeval),2)/norm(correct,2);
    err2=norm(correct-myresult2(1:numeval),2)/norm(correct,2);
    fprintf("Requested: %g Error (Leg): %g (Trap): %g\n", pr, err1,err2);
    fprintf("              Time  (Leg): %g s (Trap): %g s (Direct): %g s\n",t1,t2,t3);
end
end

function correct=slowsinc3d(ifl,a1,a2,a3,klocs_d1,klocs_d2,klocs_d3,q)
    q=q(:).';
    [a1,b1]=ndgrid(a1,klocs_d1);
    [a2,b2]=ndgrid(a2,klocs_d2);
    [a3,b3]=ndgrid(a3,klocs_d3);
    if ifl==1
        x=sin(pi*(a1-b1))./(pi*(a1-b1));
        y=sin(pi*(a2-b2))./(pi*(a2-b2));
        z=sin(pi*(a3-b3))./(pi*(a3-b3));
    else
        x=sin(a1-b1)./(a1-b1);
        y=sin(a2-b2)./(a2-b2);
        z=sin(a3-b3)./(a3-b3);
    end
    x(arrayfun(@isnan,x))=1;
    y(arrayfun(@isnan,y))=1;
    z(arrayfun(@isnan,z))=1;
    sincmat=x.*y.*z;
    correct=sum(repmat(q,length(klocs_d1),1).*sincmat,2);
end

function correct=superslowsinc3d(ifl,a1,a2,a3,klocs_d1,klocs_d2,klocs_d3,q) % wrong here??
% Alternative brute force calculation; even slower
% May be substituted in timing tests     
    results=zeros(size(a1)); 
    for ind=1:length(results)
        sm=0;
        for j=1:length(klocs_d1)
            if ifl==0
                val=q(j)*(sinc(a1(ind)-klocs_d1(j)))*(sinc(a2(ind)-klocs_d2(j)))*(sinc(a3(ind)-klocs_d3(j)));
            else
                val=q(j)*(sinc(pi*(a1(ind)-klocs_d1(j))))*(sinc(pi*(a2(ind)-klocs_d2(j))))*(sinc(pi*(a3(ind)-klocs_d3(j))));
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
