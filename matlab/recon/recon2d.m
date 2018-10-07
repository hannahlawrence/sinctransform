function recon2d

% Test image reconstruction using autoquad2d with radial and spiral
% sampling in k-space
% 
% Given measurements m and quadrature weights (corresponding to sampling
% scheme) W:
%
% Reconstruction scheme 1: given measurements m, estimate data to be F*(Wm)
% Reconstruction scheme 2 (PCG-right): 
%   Solve W^(1/2)FF*W^(1/2) x = W^(1/2) m            via conjugate gradient
%   y=W^(1/2) x
%   F*y is estimate of data
% Reconstruction scheme 3 (PCG-left):
%   Solve F*WF x = F*W m                             via conjugate gradient
%   x is estimate of data

%%%% Set uniform points (desired function locations in real space)
N1=200; N2=N1;
unifx=(0:N1-1)/N1;
unify=(0:N2-1)/N2;
[a,b]=ndgrid(unifx,unify);
unif_d1=a(:);
unif_d2=b(:);

%%%% Generate phantom
phant=phantom('Modified Shepp-Logan',N1); %can be changed
f=phant(:);

%%%% Cartesian sampling in k-space
kx=0:N1-1;
ky=0:N2-1;
[c,d]=ndgrid(kx,ky);
k_d1=c(:);
k_d2=d(:); % One could use finufft types 2 and 1 here, but for clarity 
% (and since type 3 is reasonably fast on this many points), type 3 is used
measurements=finufft2d3(unif_d1,unif_d2,f,-1,1e-10,2*pi*k_d1,2*pi*k_d2); 
retrieved_cartesian=(1/(N1*N2))*finufft2d3(k_d1,k_d2,measurements,1,1e-15,2*pi*unif_d1,2*pi*unif_d2);

%%%% Radial sampling in k-space
function x=polarx(r,th)
    x=r*cos(th);
end
function y=polary(r,th)
    y=r*sin(th);
end
r=45; 
numsectors=ceil(2*(2*pi*r+1)+1);
numsectors=round(numsectors/2); % TO DELETE
theta=2*pi/numsectors;
thetavals=0:theta:(2*pi);
rvals=0:(r*length(thetavals))/(N1*N2):r;
xp=bsxfun(@polarx,rvals',thetavals);
yp=bsxfun(@polary,rvals',thetavals);
k_d1=xp(:);
k_d2=yp(:);
weights=autoquad2d(k_d1,k_d2); %quadrature weights

    function vec=rF(strengths)
        vec=finufft2d3(unif_d1,unif_d2,strengths,-1,1e-15,2*pi*k_d1,2*pi*k_d2);
    end
    function vec=rFstar(strengths)
        vec=finufft2d3(k_d1,k_d2,strengths,1,1e-15,2*pi*unif_d1,2*pi*unif_d2);
    end

measurements=rF(f);
retrieved_radial_weights=rFstar(measurements.*weights);
retrieved_radial_weights=retrieved_radial_weights/(N1*N2); % Method 1 [1]

direct_radial=rFstar(measurements)/(N1*N2);

% Preconditioning Method 1 (Right)
sqw=sqrt(weights);
    function vec=rmat(x)
        %FFstar=sinc2d(1,k_d1,k_d2,sqw.*x,1e-14);
        FFstar=rF(rFstar(x.*sqw));
        vec=sqw.*FFstar;
    end
x=pcg(@rmat,sqw.*measurements);
y=sqw.*x;
retrieved_radial_PCGright=rFstar(y);

% Preconditioning Method 2 (Left)
    function vec=rlmat(x)
        vec=rFstar(weights.*rF(x));
    end
retrieved_radial_PCGleft=pcg(@rlmat,rFstar(weights.*measurements));

%%%% Archimedean Spiral Sampling (see [1])
N=length(unif_d1);
N=round(N*2/3); % TO RETURN
x=zeros(N,1);
y=zeros(N,1);
kmax=50; %64; %TO RETURN
for n=1:N
    val=kmax*sqrt(n/N);
    x(n)=val*cos(3*pi*val);
    y(n)=val*sin(3*pi*val);
end
k_d1=x(:);
k_d2=y(:);
weights=autoquad2d(k_d1,k_d2); %quadrature weights

    function vec=F(strengths)
        vec=finufft2d3(unif_d1,unif_d2,strengths,-1,1e-15,2*pi*k_d1,2*pi*k_d2);
    end
    function vec=Fstar(strengths)
        vec=finufft2d3(k_d1,k_d2,strengths,1,1e-15,2*pi*unif_d1,2*pi*unif_d2);
    end
measurements=F(f);
retrieved_spiral_weights=Fstar(weights.*measurements);
retrieved_spiral_weights=retrieved_spiral_weights/(N1*N2); % Method 1 [1]

direct_spiral=Fstar(F(f))/(N1*N2);

% Preconditioning Method 1 (Right)
sqw=sqrt(weights);
    function vec=mat(x)
        %FFstar=sinc2d(1,k_d1,k_d2,sqw.*x,1e-14);
        FFstar=F(Fstar(x.*sqw));
        vec=sqw.*FFstar;
    end
x=pcg(@mat,sqw.*measurements);
y=sqw.*x;
retrieved_spiral_PCGright=Fstar(y);

% Preconditioning Method 2 (Left)
    function vec=lmat(x)
        vec=Fstar(weights.*F(x));
    end
retrieved_spiral_PCGleft=pcg(@lmat,Fstar(weights.*measurements));

%%%% Take a look at the results
close all
figure('pos',[50 330 1350 200]); suptitle('Radial'); % radial
subplot(1,4,1);
%imagesc(reshape(real(retrieved_cartesian),N1,N2));
imagesc(reshape(real(rFstar(measurements)),N1,N2));
colorbar();
title('Cartesian')
subplot(1,4,2);
imagesc(reshape(real(retrieved_radial_weights),N1,N2));
colorbar();
title('Weights')
subplot(1,4,3);
imagesc(reshape(real(retrieved_radial_PCGright),N1,N2));
colorbar();
title('PCG (Right)')
subplot(1,4,4);
imagesc(reshape(real(retrieved_radial_PCGleft),N1,N2));
colorbar();
title('PCG (Left)')

figure('pos',[50 1500 1350 200]); suptitle('Spiral'); % spiral
subplot(1,4,1);
%imagesc(reshape(real(retrieved_cartesian),N1,N2));
imagesc(reshape(real(Fstar(measurements)),N1,N2));
colorbar();
title('Cartesian')
subplot(1,4,2);
imagesc(reshape(real(retrieved_spiral_weights),N1,N2));
colorbar();
title('Weights')
subplot(1,4,3);
imagesc(reshape(real(retrieved_spiral_PCGright),N1,N2));
colorbar();
title('PCG (Right)')
subplot(1,4,4);
imagesc(reshape(real(retrieved_spiral_PCGleft),N1,N2));
colorbar();
title('PCG (Left)')

figure('pos',[455 74 494 183]); suptitle('Direct Application of F*')
subplot(1,2,1);
imagesc(reshape(real(direct_radial),N1,N2)); title('Radial')
subplot(1,2,2);
imagesc(reshape(real(direct_spiral),N1,N2)); title('Spiral')

end
