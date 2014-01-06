	function elintsubr(ndim,npelm,x,nunk_pel,solution,itype,icheli)
	implicit none
	integer ndim, npelm, nunk_pel, itype, icheli
	double precision x(npelm,ndim), solution(nunk_pel), elintsubr
	double precision e(3,2), be(2), delta, gradphi(3,2), h, Uinf, C, gamma, rho_0, rho, dphidx, dphidy, normphisq
	integer i,j
	C=340.d0
	gamma=1.4d0
	rho_0=1.2d0
			
	e(1,1) = x(2,2) - x(3,2)
	e(2,1) = x(3,2) - x(1,2)
	e(3,1) = x(1,2) - x(2,2)
	e(1,2) = x(3,1) - x(2,1)
	e(2,2) = x(1,1) - x(3,1)
	e(3,2) = x(2,1) - x(1,1)
	delta = e(3,1) * e(1,2) - e(3,2) * e(1,1)
	
	gradphi(1:3,1:2) = e(1:3,1:2)/delta
	
	! calulate derivatives of solution and new density
	dphidx = 0.d0
	dphidy = 0.d0
	do i = 1,3
	dphidx = dphidx + solution(i)*gradphi(i,1)
	dphidy = dphidy + solution(i)*gradphi(i,2)
	end do
	
	normphisq =  dphidx**2 + dphidy**2
	rho = rho_0*(1-(gamma-1)/(gamma+1)*1/C**2*normphisq)**(1/(gamma-1))
	
	elintsubr = 0d0
	do i = 1, 3
	elintsubr = elintsubr + rho*abs(delta)/6d0
	end do ! i = 1, 3
	end