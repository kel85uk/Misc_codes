%
clear all; clc; close all;

L = 100;
UN = 1;
US = 0;
UE = 0;
UW = 0;
Nx = 20; Ny = 15;
tol = 1e-6;
maxit = 100;
restart = 100;
% Create coordinates
dx = 1;%L/Nx;
dy = 1;%L/Ny;
xi = 0.0; xf = Nx;
yi = 0.0; yf = Ny;
xp = xi+dx/2:dx:xf-dx/2;
yp = yi+dy/2:dy:yf-dy/2;
xu = xi:dx:xf;
yu = yi+dy/2:dy:yf-dy/2;
xv = xi+dx/2:dx:xf-dx/2;
yv = yi:dy:yf;
[XP,YP] = meshgrid(xp,yp);
XP = XP'; YP = YP';
[XU,YU] = meshgrid(xu,yu);
XU = XU'; YU = YU';
[XV,YV] = meshgrid(xv,yv);
XV = XV'; YV = YV';
% Set primitive variables
P = zeros(Nx,Ny);
U = zeros(Nx,Ny);
V = zeros(Nx,Ny);
Ub = zeros(Nx+1,Ny);
Vb = zeros(Nx,Ny+1);

% Create solutions vector
sol = [reshape(U,[],1);reshape(V,[],1);reshape(P,[],1)];

ex = ones(Nx,1)*1/dx^2;
DxxU = spdiags([ex -2*ex ex], [-1 0 1], Nx, Nx); %1D discrete Laplacian in the x-direction ;
ey = ones(Ny,1)*1/dy^2;
DyyU = spdiags([ey, -2*ey ey], [-1 0 1], Ny, Ny); %1D discrete Laplacian in the y-direction ;
DyyU(1,1) = -3/dy^2;
DyyU(end,end) = -3/dy^2;
Au = kron(speye(Ny),DxxU) + kron(DyyU,speye(Nx));
bu = zeros(Nx*Ny,1); %zeros(Nx*Ny,1);
for i=1:Nx
	for j=1:Ny
		index=i+(Nx)*(j-1);
		if (j == Ny)
			bu(index) = -2*UN*1/dy^2;
		end
		if (i == Nx)
			Au(index,:) = 0;
			Au(index,index) = 1;
		end
	end
end
DxxV = spdiags([ex -2*ex ex], [-1 0 1], Nx, Nx); %1D discrete Laplacian in the x-direction ;
DyyV = spdiags([ey, -2*ey ey], [-1 0 1], Ny, Ny); %1D discrete Laplacian in the y-direction ;
DxxV(1,1) = -3/dx^2;
DxxV(end,end) = -3/dx^2;
Av = kron(speye(Ny),DxxV) + kron(DyyV,speye(Nx));
bv = zeros(Nx*Ny,1); %zeros(Nx*Ny,1);
for i=1:Nx
	for j=1:Ny
		index=i+(Nx)*(j-1);
		if (j == Ny)
			Av(index,:) = 0;
			Av(index,index) = 1;
		end
	end
end
%full(Av)
Delplusx = spdiags([ones(Nx*Ny,1)*1/dx^2, -ones(Nx*Ny,1)*1/dx^2],[0 1],Nx*Ny,Nx*Ny);
Delplusy = spdiags([ones(Nx*Ny,1)*1/dy^2, -ones(Nx*Ny,1)*1/dy^2],[0 Nx],Nx*Ny,Nx*Ny);

for i=1:Nx
	for j=1:Ny
		index=i+(Nx)*(j-1);
		if (i == Nx)
			Delplusx(index,:) = 0;
			Delplusx(index,index) = 1;
		end
	end
end
Bp = [Delplusx;Delplusy];
BpT = Bp';
% Modify for boundary conditions
for i=1:Nx
	for j=1:Ny
		index=i+(Nx)*(j-1);
		if (i == Nx)
			Delplusx(index,:) = 0;
		end
		if (j == Ny)
			Delplusy(index,:) = 0;
		end
	end
end
Bp = [Delplusx;Delplusy];
bp = zeros(Nx*Ny,1);

Aall = [kron([1,0;0,0],Au) + kron([0,0;0,1],Av),Bp;BpT,zeros(size(BpT,1),size(Bp,2))];
Aall(end,end-Nx*Ny+1:end) = ones(1,Nx*Ny);
b_all = [bu;bv;bp];
%sol = Aall\b_all;
sol = gmres(Aall,b_all,restart,tol,maxit);

Ub(2:end,1:end) = reshape(sol(1:Nx*Ny),Nx,Ny);
Uplot = interp2(XU',YU',Ub',XP',YP');
figure(1)
surfc(XP',YP',Uplot);
view([0 90])

Vb(1:end,2:end) = reshape(sol(Nx*Ny+1:2*Nx*Ny),Nx,Ny);
Vplot = interp2(XV',YV',Vb',XP',YP');
figure(2)
surfc(XP',YP',Vplot);
view([0 90])

Pplot = reshape(sol(2*Nx*Ny+1:3*Nx*Ny),Nx,Ny);
figure(3)
surfc(XP,YP,Pplot);
view([0 90])

figure(4);
quiver(XP',YP',Uplot,Vplot);

figure(5)
spy(Aall)
%format compact
%full(Aall)
