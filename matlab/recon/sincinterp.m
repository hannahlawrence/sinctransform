function sincinterp() % IN PROGRESS
%nonunif=-5+10*rand(400,1); %need: max consecutive distance?
%nonunif=sort(nonunif);
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


