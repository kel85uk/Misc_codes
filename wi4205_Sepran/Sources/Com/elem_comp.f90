      subroutine elemsubr ( ndim, npelm, x, nunk_pel, elem_mat,  &
                            elem_vec, elem_mass, uold, itype )
      implicit none
      integer ndim, npelm, nunk_pel, itype
      double precision x(npelm,ndim), elem_mat(nunk_pel,nunk_pel),  &
                       elem_vec(nunk_pel), elem_mass(nunk_pel),     &
                       uold(nunk_pel)

!     elem_mass   o   In this two-dimensional array the student must store the
!                     element mass matrix, provided the mass matrix must be
!                     computed, in the following way:
!                     elem_mass(i,j) = s_{ij} ; i,j = 1(1)nunk_pel
!                     This matrix should only be filled if a mass matrix is
!                     required, for example for time-dependent problems.
!     elem_mat    o   In this two-dimensional array the student must store the
!                     element matrix, in the following way:
!                     elem_mat(i,j) = s_{ij} ; i,j = 1(1)nunk_pel.
!                     The degrees of freedom in an element are stored
!                     sequentially, first all degrees of freedom corresponding
!                     to the first point, then to the second, etcetera.
!     elem_vec    o   In this array the student must store the element vector,
!                     in the following way:
!                     elem_vec(i) = f_i; i = 1(1)nunk_pel
!                     It concerns the derived quantity that must be computed
!     itype       i   Type number of element.
!                     This parameter is defined in the input block PROBLEM
!     ndim        i   Dimension of the space.
!     npelm       i   Number of points per element
!     nunk_pel    i   number of degrees of freedom in the element
!     uold        i   In this array the old solution, as indicated by V1,
!                     is stored. This solution may contain the boundary
!                     conditions only, if the array has been created by
!                     prescribe_boundary_conditions, but also a starting vector
!                     if V1 has been created by create or even the previous
!                     solution in an iteration process if nonlinear_equations
!                     is used.
!     x           i   Contains the coordinates of the nodes of the element using
!                     the local node numbering
!                     x(i,1) contains the x-coordinate of the i-th node in the
!                     element and x(i,2) the y-coordinate of this node
! **********************************************************************
!
!                       LOCAL PARAMETERS

      double precision e(3,2), be(2), delta, gradphi(3,2), h, Uinf, C, gamma, rho_0, rho, dphidx, dphidy, normphisq
      integer i,j
      C=340.d0
      Uinf = C*0.5d0
      gamma=1.4d0
      rho_0=1.2d0
!
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!    none
! **********************************************************************
!
!                       I/O
!
!    none
! **********************************************************************
!
!                       ERROR MESSAGES
!
!    none
! **********************************************************************
!
!                       PSEUDO CODE
!
!    The element matrix, element right-hand side and if the problem so
!    requires the element mass matrix are filled by the user, depending on
!    the parameter itype
!
!    The element matrix and element vector are defined in the Lecture Notes
!    "Numerical Methods in Scientific Computing"
!    See also the description in the manual for the formulas
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      if (itype == 1) then
! Type 1 = internal element      
	e(1,1) = x(2,2) - x(3,2)
	e(2,1) = x(3,2) - x(1,2)
	e(3,1) = x(1,2) - x(2,2)
	
	e(1,2) = x(3,1) - x(2,1)
	e(2,2) = x(1,1) - x(3,1)
	e(3,2) = x(2,1) - x(1,1)
	
	delta = e(3,1)*e(1,2) - e(3,2)*e(1,1)
	
	gradphi(1:3,1:2) = e(1:3,1:2)/delta
	
	! calulate derivatives of solution and new density
	    dphidx = 0.d0
	    dphidy = 0.d0
	    do i = 1,3
	      dphidx = dphidx + uold(i)*gradphi(i,1)
	      dphidy = dphidy + uold(i)*gradphi(i,2)
	    end do
	    
	normphisq =  dphidx**2 + dphidy**2
	rho = rho_0*(1-(gamma-1)/(gamma+1)*1/C**2*normphisq)**(1/(gamma-1))
	
	
	do j = 1,3
	  do i = 1,3
	    elem_mat(i,j) = 0.5d0 *abs(delta) *rho* &
	    (gradphi(i,1)*gradphi(j,1) + gradphi(i,2)*gradphi(j,2))
	  end do
	end do
	
	elem_vec = 0.d0
      
      end if
    if (itype == 2) then
! Type 2 = Boundary element Gamma 2
	be(1)=x(2,1) - x(1,1)
	be(2)=x(2,2) - x(1,2)
	h = sqrt( (be(1))**2 + (be(2))**2)
	
	gradphi(1,1) = be(1)/h
	gradphi(1,2) = be(2)/h
	gradphi(2,1) = be(1)/h
	gradphi(2,2) = be(2)/h
	
	
	! calulate derivatives of solution and new density
	    dphidx = 0.d0
	    dphidy = 0.d0
	    do i = 1,2
	      dphidx = dphidx + uold(i)*gradphi(i,1)
	      dphidy = dphidy + uold(i)*gradphi(i,2)
	    end do
	    
	normphisq =  dphidx**2 + dphidy**2
	rho = rho_0 !*(1-(gamma-1)/(gamma+1)*1/C**2*normphisq)**(1/(gamma-1))
	
	elem_mat = 0.d0
	elem_vec(1:2) = h*0.5d0*Uinf*rho
      end if
     if (itype == 3) then
 ! Type 3 = Boundary element Gamma 4
 	be(1)=x(2,1) - x(1,1)
	be(2)=x(2,2) - x(1,2)
	h = sqrt( (be(1))**2 + (be(2))**2)
	
	gradphi(1,1) = be(1)/h
	gradphi(1,2) = be(2)/h
	gradphi(2,1) = be(1)/h
	gradphi(2,2) = be(2)/h
	
	
	! calulate derivatives of solution and new density
	    dphidx = 0.d0
	    dphidy = 0.d0
	    do i = 1,2
	      dphidx = dphidx + uold(i)*gradphi(i,1)
	      dphidy = dphidy + uold(i)*gradphi(i,2)
	    end do
	    
	normphisq =  dphidx**2 + dphidy**2
	rho = rho_0 !*(1-(gamma-1)/(gamma+1)*1/C**2*normphisq)**(1/(gamma-1))
 	    
	elem_mat = 0.d0
	elem_vec(1:2) = -h*0.5d0*Uinf*rho
      end if    

      end
