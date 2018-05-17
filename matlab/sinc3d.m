function wtrans=sinc3d(ifl,klocs_d1,klocs_d2,klocs_d3,q,tol)

if(nargin<1), test_sinc3d; return; end

%  
% wtrans(j) = sum sinc(klocs_d1(k)-klocs_d1(j)) * sinc(klocs_d2(k)-klocs_d2(j)) * sinc(klocs_d3(k)-klocs_d3(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% klocs_d1 = (real) sample locations in dimension 1
% klocs_d2 = (real) sample locations in dimension 2
% klocs_d3 = (real) sample locations in dimension 3
% q = sample strengths
% tol = requested precision

newtol=max(tol/1000,1e-16);

rkmaxx=max(bsxfun(@max,zeros(size(klocs_d1)),abs(klocs_d1)));
rkmaxy=max(bsxfun(@max,zeros(size(klocs_d2)),abs(klocs_d2)));
rkmaxz=max(bsxfun(@max,zeros(size(klocs_d3)),abs(klocs_d3)));

if ifl==1
    rkmaxx=pi*rkmaxx;
    rkmaxy=pi*rkmaxy;
    rkmaxz=pi*rkmaxz;
    klocs_d1=pi*klocs_d1;
    klocs_d2=pi*klocs_d2;
    klocs_d3=pi*klocs_d3;
end
rsamp=2;
nx=ceil(rsamp*round(rkmaxx+3)); 
ny=ceil(rsamp*round(rkmaxy+3));
nz=ceil(rsamp*round(rkmaxz+3));

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

wtrans=(1/8)*finufft3d3(allxx,allyy,allzz,h_at_xxyyzz.*allww,1,newtol,klocs_d1,klocs_d2,klocs_d3);

% For real inputs, only return real values. Otherwise, return complex
% wtrans
if isreal(q)
    wtrans=real(wtrans);
end

function test_sinc3d
n=10;
precisions=[1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15];
for p=1:length(precisions)
    pr=precisions(p);
    klocs_d1=-10+(20*rand(n,1));
    klocs_d2=-10+(20*rand(n,1));
    klocs_d3=-10+(20*rand(n,1));
    q=complex(rand(1,n)*30,rand(1,n)*30); 
    ifl=0;
    correct=slowsinc3d(ifl,klocs_d1,klocs_d2,klocs_d3,q);
    myresult=sinc3d(ifl,klocs_d1,klocs_d2,klocs_d3,q,pr);
    err=abs(correct-myresult);
    err=sqrt(err.'*err);
    fprintf("Requested: %g Error: %g\n", pr, err);
end

function correct=slowsinc3d(ifl,klocs_d1,klocs_d2,klocs_d3,q)
    [a1,b1]=ndgrid(klocs_d1,klocs_d1);
    [a2,b2]=ndgrid(klocs_d2,klocs_d2);
    [a3,b3]=ndgrid(klocs_d3,klocs_d3);
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