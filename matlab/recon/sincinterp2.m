function sincinterp2() % IN PROGRESS
% uniform --> anywhere, e.g. nonuniform. should be exact? proof of nyquist theorem?
unif=-5:0.01:5;%linspace(-5,5,400);
T=unif(2)-unif(1);
nonunif=-5+10*rand(400,1); nonunif=sort(nonunif);
a=0.5;
b=1;
%fU=exp(unif.^2);
fU=2*a*b*sin(2*pi*a*unif)./(2*pi*a*unif);
fU(unif==0)=2*a*b;
%interp=asinc1d(1,nonunif/T,-100:100,fU,1e-10,'trap');
interp=asinc1d(1,nonunif/T,-500:500,fU,1e-10,'trap');
if 0
    close all
    plot(unif,fU,'k-'); hold on
    plot(nonunif,interp,'r-')
    legend('True','Interpolated')
end


% nonuniform --> anywhere, e.g. nonuniform. optimal weights from
% choi-munson?
nonunif=-5+10*rand(800,1); nonunif=sort(nonunif);
a=0.5;
b=1;
fU=2*a*b*sin(2*pi*a*nonunif)./(2*pi*a*nonunif);
fU(nonunif==0)=2*a*b;
%fU=exp(nonunif.^2);
sigma=a*7;%1.1; % band limit, with some clearance...WHY?
target=-5+10*rand(200,1); target=sort(target);
weights=diff(nonunif,1);
weights=vertcat(weights,0); % what to do at boundary??
weights2=(pi/sigma)./asincsq1d(0,nonunif,nonunif,ones(size(nonunif)),1e-10,'legendre'); %ERROR FOR TRAP??
% FIGURE OUT WHAT TRAP ERROR IS ABOVE!!!!

%weights=mean(weights);
%weights=ones(size(weights));
interp=(sigma/pi)*asinc1d(0,sigma*target,sigma*nonunif,weights.*fU,1e-10,'trap'); % asinc1d(ifl,a1,klocs,q,tol,mode)
close all
plot(nonunif,fU,'k-'); hold on
plot(target,interp,'ro')
legend('True','Interpolated')

nonunif=linspace(-5,10,401); nonunif=nonunif(:);
display(maxconsecdist(nonunif))
%fNU=2*sin(nonunif)+cos(nonunif/10)+cos(nonunif/10)+sin(4*nonunif)./(4*nonunif);
a=0.5;
b=1;
fNU=2*a*b*sin(2*pi*a*nonunif)./(2*pi*a*nonunif);
% need: bandlimited function? and to know what highest freq is?
% highest frequency is 1/(2*pi)
% need max. distance lower than Nyquist rate?
%sigma=1/(2*pi);
sigma=1.4*pi;
sigma=1.1*a;
%sigma=1/2;
%weights=autoquad1d(nonunif); % this is not quite what paper suggests
%weights=1./sincsq1d(0,sigma*nonunif,ones(size(nonunif)),1e-15); % one of these is NaN
weights=ones(size(fNU))/length(fNU); 
length(find(isnan(weights)))
q=weights.*fNU;
length(find(isnan(q)))
unif=linspace(-5,5,400);
% need dist < 1/2fmax = pi
interp=asinc1d(1,sigma*unif',sigma*nonunif,q,1e-10,'legendre');
close all
hold on
plot(nonunif,fNU,'bo');
plot(unif,interp,'ro')
legend('True','Interpolated')
end

function mx=maxconsecdist(vec)
    %mat=bsxfun(@(x,y) abs(x-y),reshape(vec,length(vec),1),reshape(vec,1,length(vec)));
    %val=max(max(mat));
    mx=0;
    for i=1:length(vec)-1
        if abs(vec(i+1)-vec(i))>mx
            mx=abs(vec(i+1)-vec(i));
        end
    end
end


