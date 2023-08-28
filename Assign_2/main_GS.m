close all
clear
clc
max_e = 0.00001;
%% read xc
load xc.dat
niml=length(xc);
% nim1 = ni-1 = number of grid lines. Number of cell nodes = ni
% For a 10x10 mesh, ni=no of nodes=12,nim1=no of grid lines=11
% niml = 26, ni = 27, number of cells = 25.
ni=niml+1;

%% read yc
load yc.dat
njml=length(yc);
nj=njml+1;

%% compute the x-coordinates of the cell centres
for i=2:niml
   xp(i)=0.5*(xc(i)+xc(i-1));
end
xp(1)=xc(1);
xp(ni)=xc(niml);
%
% take the transpose of x
xp=xp';

%% compute the y-coordinates of the cell centres
for j=2:njml
   yp(j)=0.5*(yc(j)+yc(j-1));
end
yp(1)=yc(1);
yp(nj)=yc(njml);
%
% take the transpose of y
yp=yp';

%% read u
load u.dat
u2d=reshape(u,ni,nj);

%% read v
load v.dat
v2d=reshape(v,ni,nj);

%% u2df, v2df values on faces niml, njml:
u2d = u2d';
v2d = v2d';
for i=2:niml-1
    for j=1:njml-1
        u2df(i,j) = u2d(i,j+1)+(u2d(i+1,j+1)-u2d(i,j+1))*(xc(i)-xp(i))/(xp(i+1)-xp(i));
    end
end
for j=1:njml-1
    u2df(1,j) = u2d(1,j+1);
    u2df(niml,j) = u2d(niml+1,j+1);
end
for i=1:niml-1
    for j=2:njml-1
        v2df(i,j) = v2d(i+1,j)+(v2d(i+1,j+1)-v2d(i+1,j))*(yc(j)-yp(j))/(yp(j+1)-yp(j));
    end
end
for i=1:niml-1
    v2df(i,1) = v2d(i+1,1);
    v2df(i,njml) = v2d(i+1,njml+1);
end

%% To compute the delta_xe:dist b/w node 
for i=1:niml
d_x(i)= xp(i+1)-xp(i);
end

for i=1:njml
d_y(i)= yp(i+1)-yp(i);
end

%% plot the velocity field
u2d = u2d';
v2d = v2d';
% length of vectors = vec
vec = 5;
figure
quiver(xp,yp,u2d,v2d,vec)
axis('equal');
xlabel('x'); 
ylabel('y');
%print vectxy.ps -deps

%To plot the mesh
figure
[x,y]= meshgrid(xc,yc);
plot(x,yc);
hold on;
plot(xc,y');
xlabel('x'); 
ylabel('y');

%% Temperature data:
Nx = ni-2;
Ny = nj-2;
lambda = 1/500;
rho = 1;
% Initialization:
T = ones(Nx+2,Ny+2);

% Boundary conditions:
for j = 7:Ny+2;
    T(Nx+2,j) = 10;
end
for j=Ny-3:Ny+1
    T(1,j)= 20;
end

%%
% u2d, v2d values are at nodes E,W,P,S,N
% xc, yc values are at cell centers or on faces e,w,s,n
% xp, yp values are at nodes E,W,P,S,N
% d_x, d_y vlaues are nodal distances (not cell sizes).

%% coefficients:

for i=2:Nx+1
    for j=2:Ny+1
        Fw(i,j) = rho*u2df(i-1,j-1)*(yc(j)-yc(j-1));
        Fe(i,j) = rho*u2df(i,j-1)*(yc(j)-yc(j-1));
        Fs(i,j) = rho*v2df(i-1,j-1)*(xc(i)-xc(i-1));
        Fn(i,j) = rho*v2df(i-1,j)*(xc(i)-xc(i-1));
        Dw(i,j) = lambda*(yc(j)-yc(j-1))/(xp(i)-xp(i-1));
        De(i,j) = lambda*(yc(j)-yc(j-1))/(xp(i+1)-xp(i));
        Ds(i,j) = lambda*(xc(i)-xc(i-1))/(yp(i)-yp(i-1));
        Dn(i,j) = lambda*(xc(i)-xc(i-1))/(yp(i+1)-yp(i));
        W = [Fw(i,j),(Dw(i,j)+Fw(i,j)/2),0];
        E = [-Fe(i,j),(De(i,j)-Fe(i,j)/2),0];
        S = [Fs(i,j),(Ds(i,j)+Fs(i,j)/2),0];
        N = [-Fn(i,j),(Dn(i,j)-Fn(i,j)/2),0];
        aW(i,j) = max(W);
        aE(i,j) = max(E);
        aS(i,j) = max(S);
        aN(i,j) = max(N);
        aP(i,j) = aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j);
    end
end

%% Solving:
% Gauss-Siedel Method:
for k=1:100000
    sum(k) = 0;
    % loop starts:
    for i=2:Nx+1
        for j=2:Ny+1
            T(i,j) = (aE(i,j)*T(i+1,j)+aW(i,j)*T(i-1,j)+aN(i,j)*T(i,j+1)+aS(i,j)*T(i,j-1))/aP(i,j);
        end
    end
    % Boundaries update:
    for j=1:Ny-4
        T(1,j) = T(2,j);
    end
    for j=1:6
        T(Nx+2,j) = T(Nx+1,j);
    end
    for i=1:Nx+2
        T(i,1) = T(i,2);
        T(i,Ny+2) = T(i,Ny+1);
    end
    % For error calculation:
    for i=2:Nx+1
        for j=2:Ny+1
            sum(k) = sum(k) + abs(aE(i,j)*T(i+1,j)+aW(i,j)*T(i-1,j)+aN(i,j)*T(i,j+1)+aS(i,j)*T(i,j-1)-aP(i,j)*T(i,j));
        end
    end
    f = rho*1*2*0.068*10;
    sum(k) = sum(k)/f;
    iteration(k) = k;
    % residual check:
    if(sum(k)<max_e)
        break
    end
    
end

%% rotating the T matrix
T_rot = T';

%% Temperature contour plot:
[X,Y] = meshgrid(xp,yp);
figure();
contourf(X,Y,T_rot,30);
colormap(jet);
colorbar;
title('Temperature contour', 'FontSize', 18);

%% Residue plot:
figure();
plot(iteration,sum);
title('Residue plot');
xlabel('Iteration number');
ylabel('Residue error');



%% To compute Peclet number, Pe
rho=1;gamma=1/500;
for i=2:niml
    Pe=rho*u2d(i)*d_x(i)/gamma;
end