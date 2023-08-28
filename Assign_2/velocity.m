close all
clear
clc
% read xc
load xc.dat
nim1=length(xc);
% nim1 = ni-1 = number of grid lines. Number of cell nodes = ni
% For a 10x10 mesh, ni=no of nodes=12,nim1=no of grid lines=11
ni=nim1+1;
%
% read yc
load yc.dat
njm1=length(yc);
nj=njm1+1;
%
% read u
load u.dat
u2d=reshape(u,ni,nj);
% read v
load v.dat
v2d=reshape(v,ni,nj);
%
% xc and yc are the coordinates of the grid
%
%       o---------------o  yc(j)
%       |               |
%       |               |
%       |       P       |
%       |               |
%       |               |
%       o---------------o  yc(j-1)
%
%     xc(i-1)          xc(i)          
%
% The cell-node (P) has the coordinates xp(i), yp(j)
%
% compute the x-coordinates of the cell centres
for i=2:nim1
   xp(i)=0.5*(xc(i)+xc(i-1));
end
xp(1)=xc(1);
xp(ni)=xc(nim1);
%
% take the transpose of x
xp=xp';

% compute the y-coordinates of the cell centres
for j=2:njm1
   yp(j)=0.5*(yc(j)+yc(j-1));
end
yp(1)=yc(1);
yp(nj)=yc(njm1);
%
% take the transpose of y
yp=yp';

% plot the velocity field
% length of vectors = vec
vec= 5
figure
quiver(xp,yp,u2d,v2d,vec)
axis('equal');
xlabel('x'); ylabel('y')
%print vectxy.ps -deps

%To plot the mesh
figure
[x,y]= meshgrid(xc,yc);
plot(x,yc);
hold on;
plot(xc,y');
xlabel('x'); ylabel('y')

%To compute the delta_xe:dist b/w node 
for i=1:nim1
d_x(i)= xp(i+1)-xp(i);
end

for i=1:njm1
d_y(i)= yp(i+1)-yp(i);
end

%To compute Peclet number, Pe
rho=1;gamma=1/500;
for i=2:nim1
    Pe=rho*u2d(i)*d_x(i)/gamma;
end