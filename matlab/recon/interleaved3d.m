function [k_d1,k_d2,k_d3]=interleaved3d(rmax,numr)
rvals=0.001:rmax/numr:rmax;
imaxvals= ceil(rvals/2); % imax for each radius
pmaxvals= ceil(rvals/2); % pmax for each radius

k_d1=[];
k_d2=[];
k_d3=[];
for ind=1:length(rvals)
    radius=rvals(ind);
    imax=imaxvals(ind); % number interleaves
    pmax=pmaxvals(ind); % points per interleave
    for i=1:imax
        for p=1:pmax
            zval=(2*p-pmax-1)/pmax; % move this up? no?
            k_d3=[k_d3 radius*zval];
            inval=sqrt(pmax*pi/imax)*asin(zval) + (2*pi*i)/imax;
            coeff=sqrt(1-(zval)^2);
            k_d1=[k_d1 radius*cos(inval)*coeff];
            k_d2=[k_d2 radius*sin(inval)*coeff];
        end
    end
end

show=0;
if show
    figure;
    pt=scatter3(k_d1(1),k_d2(1),k_d3(1));
    xlim([min(k_d1) max(k_d1)]);
    ylim([min(k_d2) max(k_d2)]);
    zlim([min(k_d3) max(k_d3)]);

    for v=2:length(k_d1)
        pt.XData=k_d1(1:v);
        pt.YData=k_d2(1:v);
        pt.ZData=k_d3(1:v);
        pause(0.001);
    end
end
k_d1=k_d1(:);
k_d2=k_d2(:);
k_d3=k_d3(:);
