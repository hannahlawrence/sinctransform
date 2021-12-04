function sincinterp1d() % IN PROGRESS
% Interpolate from uniform samples to nonuniform targets
% With sufficient samples, reconstruction should be exact (Nyquist) 
unif=-5:0.01:5;
T=unif(2)-unif(1);
nonunif=-5+10*rand(400,1); nonunif=sort(nonunif);
a=0.5;
b=1;
%fU=exp(unif.^2);
fU=2*a*b*sin(2*pi*a*unif)./(2*pi*a*unif); % Example function
fU(unif==0)=2*a*b;
fU=exp(-(unif.^2));
%interp=asinc1d(1,nonunif/T,-100:100,fU,1e-10,'trap');
interp=sinc1d(1,nonunif/T,-500:500,fU,1e-10,'trap');
close all; plot(unif,fU,'k-'); hold on
plot(nonunif,interp,'ro-')
legend('True','Interpolated')
title('Uniform --> Nonuniform')

% Interpolate from nonuniform samples to nonuniform targets 
% Option 1 for Weights: space between samples (see Choi and Munson)
% Option 2 for Weights: sinc-squared according to bandwidth (see Choi and Munson)
nonunif=-5+10*rand(800,1); nonunif=sort(nonunif);
a=0.5;
b=1;
fU=2*a*b*sin(2*pi*a*nonunif)./(2*pi*a*nonunif); % Example function
fU(nonunif==0)=2*a*b;
fU=exp(-(nonunif.^2)); % Gaussian
%fU=exp(sin(nonunif)/10);
sigma=3;%*4;%1.1; % band limit, with some clearance...
target=-5+10*rand(200,1); target=sort(target);
weights=diff(nonunif,1); % Option 1
weights=vertcat(weights,0); % Boundary case
weights2=(pi/sigma)./sincsq1d(0,sigma*nonunif,sigma*nonunif,ones(1,length(nonunif)),1e-10,'trap'); %Option 2
interp1=(sigma/pi)*sinc1d(0,sigma*target,sigma*nonunif,weights.*fU,1e-10,'trap'); 
interp2=(sigma/pi)*sinc1d(0,sigma*target,sigma*nonunif,weights2.*fU,1e-10,'trap'); 
figure; plot(nonunif,fU,'k-'); hold on
plot(target,interp1,'ro')
plot(target,interp2,'go')
legend('True','Interpolated (Diff)','Interpolated (Sinc-Sq)')
title('Nonuniform --> Nonuniform')
end

function mx=maxconsecdist(vec)
    mx=0;
    for i=1:length(vec)-1
        if abs(vec(i+1)-vec(i))>mx
            mx=abs(vec(i+1)-vec(i));
        end
    end
end


