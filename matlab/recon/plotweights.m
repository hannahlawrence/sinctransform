function plotweights(k_d1,k_d2,weights) % IN PROGRESS
if nargin<1
    testplotweights;
    return
end
mn=min(min(k_d1),min(k_d2));
mx=max(max(k_d1),max(k_d2));
vec=mn:0.5:mx;
[v,w]=meshgrid(vec,vec);
temp=griddata(k_d1,k_d2,weights,v,w);
mesh(vec,vec,temp);
hold on
plot3(k_d1,k_d2,weights,'ro');
end

function testplotweights()
N=1000;
x=zeros(N,1);
y=zeros(N,1);
kmax=64;
for n=1:N
    val=kmax*sqrt(n/N);
    if mod(3*pi*val,2*pi)>=(pi/4)
        x(n)=val*cos(3*pi*val); %make a gap
        y(n)=val*sin(3*pi*val);
    else
        x(n)=0;
        y(n)=0;
    end
end
k_d1=x(:);
k_d2=y(:);
weights=autoquad2d(k_d1,k_d2); %quadrature weights

end