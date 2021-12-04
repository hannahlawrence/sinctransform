function sincinterp2d() % IN PROGRESS
% Interpolate from uniform samples to nonuniform targets
% With sufficient samples, reconstruction should be exact (Nyquist) 
unifx=-5:0.01:5;
unify=-5:0.01:5;
[X,Y]=ndgrid(unifx,unify);
T=unifx(2)-unifx(1);
nonunifx=-5+10*rand(100,1); nonunifx=sort(nonunifx);
nonunify=-5+10*rand(100,1); nonunify=sort(nonunify);
[NX,NY]=ndgrid(nonunifx,nonunify);
a=0.5;
b=1;
%fU=exp(unif.^2);
fU=2*a^2*b^2*sinc(2*pi*a*X).*sinc(2*pi*a*Y);%sin(2*pi*a*X)./(2*pi*a*X).*sin(2*pi*a*Y)./(2*pi*a*Y); % Example function
interp=sinc2d(1,NX(:)/T,NY(:)/T,X(:)/T,Y(:)/T,fU(:).',1e-10,'trap');
interp=reshape(interp,size(NX,1),size(NX,2));
close all; subplot(1,2,1);mesh(X,Y,fU); title('True');
subplot(1,2,2); mesh(NX,NY,interp); title('Interpolated');
title('Uniform --> Nonuniform')

% OTHER METHODS!
% Interpolate from nonuniform samples to nonuniform targets 
% Option 1 for Weights: space between samples (see Choi and Munson)
% Option 2 for Weights: sinc-squared according to bandwidth (see Choi and Munson)
nonunifx=-5+10*rand(800,1); nonunifx=sort(nonunifx);
nonunify=-5+10*rand(800,1); nonunify=sort(nonunify);
[X,Y]=ndgrid(nonunifx,nonunify);
a=0.5;
b=1;
fU=2*a^2*b^2*sinc(2*pi*a*X).*sinc(2*pi*a*Y); % Example function
%fU=exp(-(nonunif.^2)); % Gaussian
sigma=a*7;%*4;%1.1; % band limit, with some clearance...
targetx=-5+10*rand(200,1); targetx=sort(targetx);
targety=-5+10*rand(200,1); targety=sort(targety);
[NX,NY]=ndgrid(targetx,targety);
weightsx=diff(nonunifx,1); % Option 1
weightsx=vertcat(weightsx,0); % Boundary case
weightsy=diff(nonunify,1); % Option 1
weightsy=vertcat(weightsy,0); % Boundary case
weightsx=weightsx(:);
weightsy=weightsy(:); weightsy=weightsy';
weightmat=weightsx*weightsy;
interp1=(sigma/pi)*sinc2d(0,sigma*NX(:),sigma*NY(:),sigma*X(:),sigma*Y(:),(weightmat(:).*fU(:)).',1e-10,'trap'); 
interp=reshape(interp1,size(NX,1),size(NX,2));
figure; subplot(1,2,1);mesh(X,Y,fU); title('True');
subplot(1,2,2); mesh(NX,NY,interp); title('Interpolated');
title('Nonuniform --> Nonuniform, Differences')

weights2=(pi/sigma)./sincsq2d(0,sigma*X(:),sigma*Y(:),sigma*X(:),sigma*Y(:),ones(1,length(X(:))),1e-10,'trap'); %Option 2
interp2=(sigma/pi)*sinc2d(0,sigma*NX(:),sigma*NY(:),sigma*X(:),sigma*Y(:),(weights2(:).*fU(:)).',1e-10,'trap'); 
interp=reshape(interp2,size(NX,1),size(NX,2));
figure; subplot(1,2,1);mesh(X,Y,fU); title('True');
subplot(1,2,2); mesh(NX,NY,interp); title('Interpolated');
title('Nonuniform --> Nonuniform, Sincsq Weights')

end

function y=sinc(x)
y=zeros(size(x));
y(x==0)=1;
y(x~=0)=sin(x(x~=0))./x(x~=0);
end