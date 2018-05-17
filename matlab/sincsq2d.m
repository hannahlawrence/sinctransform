function wtrans=sincsq2d(ifl,klocs_d1,klocs_d2,q,tol)

if(nargin<1), test_sincsq2d; return; end

%  
% wtrans(j) = sum sinc^2(klocs_d1(k)-klocs_d1(j)) * sinc^2(klocs_d2(k)-klocs_d2(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% klocs_d1 = (real) sample locations in dimension 1
% klocs_d2 = (real) sample locations in dimension 2
% q = sample strengths
% tol = requested precision

newtol=max(tol/1000,1e-16);
rkmaxx=max(bsxfun(@max,zeros(size(klocs_d1)),abs(klocs_d1)));
rkmaxy=max(bsxfun(@max,zeros(size(klocs_d2)),abs(klocs_d2)));

if ifl==1
    rkmaxx=pi*rkmaxx;
    rkmaxy=pi*rkmaxy;
    klocs_d1=pi*klocs_d1;
    klocs_d2=pi*klocs_d2;
end
rsamp=2;
nx=ceil(rsamp*round(rkmaxx+3)); 
ny=ceil(rsamp*round(rkmaxy+3));

[xx,wwx]=lgwt(nx,-1,1);
[yy,wwy]=lgwt(ny,-1,1);

% quadrature points and weights to cover [-2,2] in each dimension

xx=vertcat(xx-1,xx+1);
wwx=vertcat(wwx,wwx);
wwx=wwx.*(2-abs(xx));
yy=vertcat(yy-1,yy+1);
wwy=vertcat(wwy,wwy);
wwy=wwy.*(2-abs(yy));
[c,d]=ndgrid(xx,yy);
allxx=c(:);
allyy=d(:);
[e,f]=ndgrid(wwx,wwy);
allww=(e(:)).*(f(:));

h_at_xxyy=finufft2d3(klocs_d1,klocs_d2,q,-1,newtol,allxx,allyy);

wtrans=(1/16)*finufft2d3(allxx,allyy,h_at_xxyy.*allww,1,newtol,klocs_d1,klocs_d2);

% For real inputs, only return real values. Otherwise, return complex
% wtrans
if isreal(q)
    wtrans=real(wtrans);
end

function test_sincsq2d
n=50;
precisions=[1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14 1e-15];
for p=1:length(precisions)
    pr=precisions(p);
    klocs_d1=-pi+(2*pi*rand(n,1));
    klocs_d2=-pi+(2*pi*rand(n,1));
    q=complex(rand(1,n)*30,rand(1,n)*30);
    ifl=0;
    myresult=sincsq2d(ifl,klocs_d1,klocs_d2,q,pr);
    correct=slowsincsq2d(ifl,klocs_d1,klocs_d2,q);
    err=abs(correct-myresult);
    err=sqrt(err.'*err);
    fprintf("Requested: %g Error: %g\n", pr, err);
end

function correct=slowsincsq2d(ifl,klocs_d1,klocs_d2,q)
    [a1,b1]=ndgrid(klocs_d1,klocs_d1);
    [a2,b2]=ndgrid(klocs_d2,klocs_d2);
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
    sincmat=sincmat.^2;
    correct=sum(repmat(q,length(klocs_d1),1).*sincmat,2); 