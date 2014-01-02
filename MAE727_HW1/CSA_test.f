      program CSA_test
      integer i, ij, iT, Ni, Nn, NT, Nmax
      integer, parameter :: n = 500
      integer, parameter :: long = selected_real_kind(9,99)
      integer, dimension(n) :: x, s, xi = 0
      real(kind=long), dimension(n) :: J, E, QE, Ji, Ei, QEi = 0.
      real(kind=long) :: S_QE, S_J, boltz, g, D0, DI, Mo, PI, h, k = 0.
      real(kind=long) :: S_QEI, S_JI, Ne, Ntot, deb, dm, elec, a0 = 0.
      real(kind=long) :: C, Na, Me, Ke, Kei, Ked, QTe, QP, QA = 0.
      real(kind=long) :: eps_i, eps_a, eps_m, eps_t, S_QE2, S_QE1, 
     &                   SI, SII = 0.
      real temp, mol_parti, mol_partf, dtemp 
      open(51,file='CSAtest.out')
c      open(53,file='SAHA_stuff.out')
      open(52,file='H_atom2.in')
      PI = 4.d0*atan(1.d0)
      Me = 9.1093897d-31
      Na = 6.02214129d23
      elec = 1.60217733d-19
      k = 1.3806488d-23
      C = 299792458.d0
      h = 6.62606957d-34
c     boltzmann's constant 1/cm.K
      boltz = k/(100.d0*h*C)
      read (52,*) Ni, Nn, Mo, D0, Ntot !Ni = number of energy levels for atom, Nn = number of levels for ion, Mo = mass of 1 mol of atom, D0 = dissociation energy of molecule (cm-1), Ntot = Total number density (cm-3)
      read (52,*) temp, dtemp, NT !temp = initial temperature, dtemp = temperature increment, NT = number of temperature steps
      write(51,21) 'T (K)', 'qA (cm-3)', 'qP (cm-3)',
     & 'qTe (cm-3)', 'qTA (cm-3)', 'QM_int (cm-3)',
     & 'QM (cm-3)', 'eps_A (cm-1)', 'eps_i (cm-1)', 'eps_t (cm-1)',
     & 'Kei (cm-3)', 'Ked (cm-3)',
     & 'Ke (cm-6)', 'Log(Ke)'
c     temperature input at eV to K      
c      temp = temp*11604.505
      do i=1,Ni+1
        read(52,*) x(i), s(i), J(i), E(i)
      end do
      do i=1,Nn
        read(52,*) xi(i), Ji(i), Ei(i)
      end do
c     Calculate first Bohr radius
      a0 = 0.529177249d-8 !in cm
      Ne = 100. !Placeholder number for electron number density   
      do iT = 1,NT
c     Calculate density cutoff radius      
      dm = (1.d0/Ntot)**(1.d0/3.d0) !in cm
c     Calculate Debye shielding length
      deb = 0.69*sqrt(temp/Ne) !in cm
      if(dm .gt. deb) then
      dm = deb
      end if 
      Nmax = floor(sqrt(dm/a0))
      ij = 0
c     Search for the corresponding energy level of index Nmax = s
      do i=1,Ni
        if(Nmax .lt. s(i)) then
        Nmax = x(i) - 1
        ij = 1
        end if
      end do
      if (ij .eq. 0) then
        Nmax = x(Ni)
      end if
      do i=1,Nmax
        QE(i) = (J(i)+1.d0)*exp(-E(i)/(boltz*temp))
      end do
      do i=1,Nn
        QEi(i) = (Ji(i)+1.d0)*exp(-Ei(i)/(boltz*temp))
      end do      
c     Calculate the atom internal partition function
      do i=1,Nmax
      S_QE = S_QE + QE(i)
      end do
c      do i=2,4
c      S_QE2 = S_QE2 + QE(i)
c      end do
c      do i=1,2
c      S_QE1 = S_QE1 + QE(i)
c      end do
c     Calculate ion internal partition function
      do i=1,Nn
      S_QEI = S_QEI +QEi(i)
      end do
c     To account for hydrogen ion internal partition function (proton)
      if (Nn == 0) then
      S_QEI = 1.0d0
      end if
c     Obtain the ionization energy of the atom measured from the ground state (subject to cutoff)         
      DI = E(Nmax+1)
c     Atom and ion translational partition function (Assuming negligible mass difference)
      QT = (2.d0*PI*Mo*k*temp/(1000.d0*Na*h*h))**(1.5)*1.d-6 !multiplication factor to convert to 1/cm^3
c      QT = 1.878e20*sqrt(Mo*temp)*Mo*temp
      QTe = (2.d0*PI*Me*k*temp/(h*h))**(1.5)*1.d-6 !electron translation partition function
c     Electron total partition function (QEe*QTe) QEe = 2.0
      QTe = 2.d0*QTe
c     Atom total partition function
      QA = S_QE*QT      
c     Ion total partition function (QTH*QEH)
      QP = QT*S_QEI
c     Calculate molecular partition function
      mol_parti = 0.0
      mol_partf = 0.0
      call vibrot(temp,mol_parti,mol_partf)
c     Calculate internal energy of atoms
      do i=1,Nmax
        eps_a = eps_a + E(i)*(J(i)+1.d0)*exp(-E(i)/(boltz*temp))
      end do
      do i=1,Nn
        eps_i = eps_i + Ei(i)*(Ji(i)+1.d0)*exp(-Ei(i)/(boltz*temp))
      end do      
      eps_t = (3.d0)*boltz*temp/(2.d0)
      eps_a = eps_a/S_QE! + (3.d0)*boltz*temp/(2.d0)
      eps_i = eps_i/S_QEI
c     Calculate equilibrium constant Ke
      Ke = exp(-(DI+D0)/(boltz*temp))*QP*QTe*QA/(mol_partf*Na*Na)
      Kei = exp(-(DI)/(boltz*temp))*QP*QTe/(QA*Na)
      Ked = exp(-(D0)/(boltz*temp))*QA*QA/(mol_partf*Na)
      write(51,22) temp, S_QE, S_QEI, QTe, QT, mol_parti, 
     & mol_partf, eps_a, eps_i, eps_t, Kei, Ked, Ke, log(Ke)
c      write(*,'(1p4E15.6)') eps_t
c      write(53,20) temp, QP*QP*(QTe*QTe)*exp(-(DI+D0)/(boltz*temp))/
c     &             (mol_partf*Na*Na*Na)
c      write(*,*) Nmax
   20 format (2(1p4E10.3,3X))
   21 format (14(A,5X))
   22 format (14(1p4E15.6,3X))
      temp = temp + dtemp
      S_QE = 0.d0
      S_QEI = 0.d0
      S_QE1 = 0.d0
      S_QE2 = 0.d0
      end do
c     Obtain the coordinates for the potential function based on Jarmain's method   
      call turnpt
      end
