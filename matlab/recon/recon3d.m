function recon3d % Note: still being developed

% Test image reconstruction using autoquad3d with spiral or nested
% spherical sampling in k-space
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
%
% Note: autoquad3d may be slow on large numbers of points; you can manually
% reduce the number of quadrature points or requested precision in sincsq3d
% (whcih is called by autoquad3d), if necessary

%%%% Set uniform points (desired function locations in real space)
N1=20;
N2=N1;
N3=N1;
unifx=(0:N1-1)/N1;
unify=(0:N2-1)/N2;
unifz=(0:N3-1)/N3;
[a,b,c]=ndgrid(unifx,unify,unifz);
unif_d1=a(:);
unif_d2=b(:);
unif_d3=c(:);

%%%% Generate 3d phantom
phant=phantom3d('Modified Shepp-Logan',N1);
f=phant(:);

%%%% Cartesian
kx=0:N1-1;
ky=0:N2-1;
kz=0:N3-1;
[c,d,e]=ndgrid(kx,ky,kz);
k_d1=c(:);
k_d2=d(:); 
k_d3=e(:); % One could use finufft types 2 and 1 here, but for clarity 
% (and since type 3 is reasonably fast on this many points), type 3 is used
measurements=finufft3d3(unif_d1,unif_d2,unif_d3,f,-1,1e-10,2*pi*k_d1,2*pi*k_d2,2*pi*k_d3); 
retrieved_cartesian=(1/(N1*N2*N3))*finufft3d3(k_d1,k_d2,k_d3,measurements,1,1e-15,2*pi*unif_d1,2*pi*unif_d2,2*pi*unif_d3);

%%%% Generate interesting k_d1,k_d2,k_d3 (sample points in k-space, 3d)
N=N1*N2*N3;
x=zeros(N,1);
y=zeros(N,1);
z=zeros(N,1);
kmax=35;
for n=1:N % Something somewhat like the Archimedean spiral in 3d
    r=kmax*sqrt(n/N);
    theta=3*pi*r;
    phi=3*pi*r;
    x(n)=r*sin(theta)*cos(phi);
    y(n)=r*sin(theta)*sin(phi);
    z(n)=r*cos(theta);
end
k_d1=x(:); k_d2=y(:); k_d3=z(:);

%[k_d1,k_d2,k_d3]=interleaved3d(50,10); % Another option: nested shells

%weights=autoquad3d(k_d1,k_d2,k_d3); % may be slow for large # points
weights=ones(size(k_d1)); % benchmark
    function vec=F(strengths)
        vec=finufft3d3(unif_d1,unif_d2,unif_d3,strengths,-1,1e-15,2*pi*k_d1,2*pi*k_d2,2*pi*k_d3);
    end
    function vec=Fstar(strengths)
        vec=finufft3d3(k_d1,k_d2,k_d3,strengths,1,1e-15,2*pi*unif_d1,2*pi*unif_d2,2*pi*unif_d3);
    end

measurements=F(f);
retrieved_weights=Fstar(measurements.*weights);

% Preconditioning Method 1 (Right)
sqw=sqrt(weights);
  function vec=rmat(x)
        %FFstar=sinc2d(1,k_d1,k_d2,sqw.*x,1e-14);
        FFstar=F(Fstar(x.*sqw));
        vec=sqw.*FFstar;
    end
x=pcg(@rmat,sqw.*measurements);
y=sqw.*x;
retrieved_PCGright=Fstar(y);

% Preconditioning Method 2 (Left)
    function vec=rlmat(x)
        vec=Fstar(weights.*F(x));
    end
retrieved_PCGleft=pcg(@rlmat,Fstar(weights.*measurements));

% Look at the sampling scheme in 3d
figure;
pt=scatter3(k_d1(1),k_d2(1),k_d3(1));
title('3D K-space Sampling')
xlim([min(k_d1) max(k_d1)]);
ylim([min(k_d2) max(k_d2)]);
zlim([min(k_d3) max(k_d3)]);

for i=2:length(k_d1)
    pt.XData=k_d1(1:i);
    pt.YData=k_d2(1:i);
    pt.ZData=k_d3(1:i);
    pause(0.001);
end

% Look at the results via 2D slices
figure('pos',[50 1500 1350 200]);
for slice=1:N3
    subplot(1,4,1);
    
    mretrieved_cartesian=reshape(retrieved_cartesian,N1,N2,N3);
    imagesc(real(mretrieved_cartesian(:,:,slice)));
    title('Cartesian');
    
%     imagesc(phant(:,:,slice));
%     colorbar();
%     title('True')

    subplot(1,4,2);
    mretrieved_weights=reshape(retrieved_weights,N1,N2,N3);
    imagesc(real(mretrieved_weights(:,:,slice)));
    title('Weights');
    colorbar();
    
    subplot(1,4,3);
    mretrieved_PCGright=reshape(retrieved_PCGright,N1,N2,N3);
    imagesc(real(mretrieved_PCGright(:,:,slice)));
    colorbar();
    title('PCG (Right)')
    
    subplot(1,4,4);
    mretrieved_PCGleft=reshape(retrieved_PCGleft,N1,N2,N3);
    imagesc(real(mretrieved_PCGleft(:,:,slice)));
    colorbar();
    title('PCG (Left)')
    
    suptitle(strcat('Slice ',num2str(slice)))
    pause(0.01);
end
end