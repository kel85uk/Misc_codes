      subroutine eldervsubr ( ndim, npelm, x, nunk_pel, elem_vec,  &
                              solution, itype, icheld, len_outvec )
! ======================================================================
!
!        programmer     Guus Segal
!        version 1.0    date 02-10-2009
!
!
!   copyright (c) 2009-2009  "Ingenieursbureau SEPRA"
!   permission to copy or distribute this software or documentation
!   in hard copy or soft copy granted only by written license
!   obtained from "Ingenieursbureau SEPRA".
!   all rights reserved. no part of this publication may be reproduced,
!   stored in a retrieval system ( e.g., in memory, disk, or core)
!   or be transmitted by any means, electronic, mechanical, photocopy,
!   recording, or otherwise, without written permission from the
!   publisher.
! **********************************************************************
!
!                       DESCRIPTION
!
!      User element subroutine to compute derivatives to be used by SEPRAN
!      This is the special version for
!      the Numerical Analysis Lab of Delft University of Technology
! **********************************************************************
!
!                       KEYWORDS
!
!    derivatives
!    vector
!    element
! **********************************************************************
!
!                       MODULES USED
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      implicit none
      integer ndim, npelm, nunk_pel, itype, icheld, len_outvec
      double precision x(npelm,ndim), solution(nunk_pel), &
                       elem_vec(len_outvec)

!     elem_vec    o   In this array the student must store the element vector,
!                     in the following way:
!                     elem_vec(i) = f_i; i = 1(1)nunk_pel
!                     It concerns the derived quantity that must be computed
!     icheld      i   Choice parameter to distinguish between the various
!     itype       i   Type number of element.
!                     This parameter is defined in the input block PROBLEM
!     ndim        i   Dimension of the space.
!     npelm       i   Number of points per element
!     nunk_pel    i   number of degrees of freedom in the element
!     solution    i   In this array the solution from which the derived
!                     quantities must be computed, as indicated by Vi, is stored
!     x           i   Contains the coordinates of the nodes of the element using
!                     the local node numbering
!                     x(i,1) contains the x-coordinate of the i-th node in the
!                     element and x(i,2) the y-coordinate of this node
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
!      Dit gedeelte moet door de practicant worden ingevuld
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
! **********************************************************************
!
!                       PSEUDO CODE
!
!    The element vector is filled by the user, depending on the parameters
!    itype and icheld
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      double precision e(3,2), delta, gradphi(3,2), dudx
      integer i, j
	e(1,1) = x(2,2) - x(3,2)
	e(2,1) = x(3,2) - x(1,2)
	e(3,1) = x(1,2) - x(2,2)
	
	e(1,2) = x(3,1) - x(2,1)
	e(2,2) = x(1,1) - x(3,1)
	e(3,2) = x(2,1) - x(1,1)
	
	delta = e(3,1)*e(1,2) - e(3,2)*e(1,1)
	
	gradphi(1:3,1:2) = e(1:3,1:2)/delta
	
	if (icheld == 1) then
	dudx = 0.d0
	do i = 1,3
	  dudx = dudx + solution(i)*gradphi(i,1)
	end do
	
	elem_vec(1:3) = -dudx
	else
	dudx = 0.d0
	do i = 1,3
	  dudx = dudx + solution(i)*gradphi(i,2)
	end do
	
	elem_vec(1:3) = -dudx
	
	end if
      end
