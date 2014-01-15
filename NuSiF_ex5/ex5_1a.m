%
dx = 1;
dy = 1;
Nx = 3;
Ny = 3;

ex = ones(Nx,1)*1/dx^2;
DxxU = spdiags([ex -2*ex ex], [-1 0 1], Nx, Nx); %1D discrete Laplacian in the x-direction ;
ey = ones(Ny,1)*1/dy^2;
DyyU = spdiags([ey, -2*ey ey], [-1 0 1], Ny, Ny); %1D discrete Laplacian in the y-direction ;
%DyyU(1,1) = -3/dy^2;
%DyyU(end,end) = -3/dy^2;
Au = kron(speye(Ny),DxxU) + kron(DyyU,speye(Nx));
bu = zeros(Nx*Ny,1); %zeros(Nx*Ny,1);
for i=1:Nx
	for j=1:Ny
		index=i+(Nx)*(j-1);
		if (i == Nx)
			Au(index,:) = 0;
			Au(index,index) = 1;
		end
	end
end
DxxV = spdiags([ex -2*ex ex], [-1 0 1], Nx, Nx); %1D discrete Laplacian in the x-direction ;
DyyV = spdiags([ey, -2*ey ey], [-1 0 1], Ny, Ny); %1D discrete Laplacian in the y-direction ;
%DxxV(1,1) = -3/dx^2;
%DxxV(end,end) = -3/dx^2;
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

Aall = [kron([1,0;0,0],Au) + kron([0,0;0,1],Av),Bp;BpT,zeros(size(BpT,1),size(Bp,2))];
Aall(end,end-Nx*Ny+1:end) = ones(1,Nx*Ny); % For normalization of pressure

full(Aall)
dlmwrite('matrix_ohneBC',Aall,'\t')
