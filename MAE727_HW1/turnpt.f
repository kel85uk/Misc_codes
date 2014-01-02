c***********************************************************************
      subroutine turnpt
c calculates turning points using jarmain's method
      real k,m,n
      dimension r(2,80),a(5),b(5),rbyre(2),rbohr(2,80),ev(80),
     1 title(20)
      open (5,file='turnptcoa.inp',status='old',form='formatted')
      open(6,file='turnpt_kklloh.out')
      gamma=0.
      delta=0.
      beta=0.
      nf=20
   10 continue
      nf=nf+1
      write(6,40)
   40 format(/)
      read(5,11) jvmax,jvmax1, (title(i),i=1,8)
      write(6,11) jvmax,jvmax1, (title(i),i=1,8)
   11 format(2i5,8a4)
      read(5,30) (title(i),i=1,19)
   30 format(20a4)
      write(6,30) (title(i),i=1,19)
      read(5,1)d,re,we,wexe,weye,weze,be,alpha
    1 format(8e10.3)
      write(6,2)d,re,we,wexe,weye,weze,be,alpha
    2 format(3h d=,f10.2,2x,3hre=,f10.5,2x,3hwe=,f10.3,2x,5hwexe=,
     1 f10.2,2x,5hweye=,f10.4,2x,5hweze=,f10.6,2x,3hbe=,f10.5,2x,
     2 6halpha=,f10.6)
      if(d.lt.1.) stop
      k=wexe/we
      m=weye/we
      n=weze/we
c
c calculation starts
      jj=0
      do 4 jv=1,jvmax
      v=jv-0.5
      ev(jv)=we*v-wexe*v**2+weye*v**3+weze*v**4
      a(1)=(5./6.)*k
      a(2)=(43./40.)*k**2-(11./10.)*m
      a(3)=(177./122.)*k**2-(81./28.)*k*m+(93./70.)*n
      a(4)=(2867./1152.)*k**4-(311./48.)*k**2*m+(1487./420.)*k*n
     >+(331./168.)*m**2
      a(5)=(11531./2816.)*k**5-(101527./2464.)*k**3*m+(29627./3696.)
     1 *k**2*n+(19205./528.)*k*m**2-(47821./4620.)*m*n
      s1=1.
      do 5 jr=1,5
      s1=s1+a(jr)*v**jr
    5 continue
c
      alphab=alpha/be
      gammab=gamma/be
      b(1)=-(2./3.)*alphab
      b(2)=-(3./5.)*k*alphab+(8./15.)*gammab
      b(3)=(-(23./28.)*k**2+(29./35.)*m)*alphab+(52./105.)*k*gammab
      b(4)=(-(4721./2520.)*k**3+(481./210.)*k*m-(65./63.)*n)*alphab
     1 +((73./105.)*k**2-(44./63.)*m)*alphab
      b(5)=(-(1451./704.)*k**4+(9847./1848.)*k**2*m-(1327./462.)*
     1 k*n-(1061./660.)*m**2)*alphab+((1513./1386.)*k**3-(2278./
     21155.)*k*m+(68./77.)*n)*gammab
      s2=0.
      do 6 jr=1,5
      s2=s2+b(jr)*v**jr
    6 continue
      fbyre=2.0*sqrt(be/we)*sqrt(v)*s1
      rbyre(1)=sqrt(s1/(s1+s2)+fbyre**2)+fbyre
      rbyre(2)=sqrt(s1/(s1+s2)+fbyre**2)-fbyre
      do 7 i=1,2
    7 rbohr(i,jv)=rbyre(i)*re/0.529177
    4 continue
c
c write out results
      do 8 jv=1,jvmax
      j=jv-1
      write(6,9) j,ev(jv),rbohr(2,jv),rbohr(1,jv)
    9 format(i5,f15.2,10x,2f10.5)
    8 continue
c
      call coordf(jvmax,d,ev,re,rbohr)
c
c determine the value and slope of exponential curves at maxima
c
c lower limit
      rl=rbohr(2,jvmax)
      ul=ev(jvmax)/219474.6
      dudrl=((ev(jvmax)-ev(jvmax-1))/219474.6)/(rbohr(2,
     1 jvmax)-rbohr(2,jvmax-1))
      bl=-dudrl/ul
      al=ul/exp(-bl*rl)
c
c upper limit
      ru=rbohr(1,jvmax)
      uu=ev(jvmax)/219474.6
      dudru=((ev(jvmax)-ev(jvmax-1))/219474.6)/(rbohr(1,
     1 jvmax)-rbohr(1,jvmax-1))
      bu=dudru/uu
      au=uu/exp(bu*ru)
      write(6,22) rl,al,bl,ru,au,bu
   22 format(4h rl=1pe15.7,3x,2ha=,e15.7,3x,2hb=,e15.7/
     1       4h ru=e15.7,3x,2ha=,e15.6,3x,2hb=,e15.7)
c
c generate data for spline program of langhoff
      write(nf,25) (title(i),i=1,8)
   25 format(8a4)
      do 12 jv1=1,jvmax
      jv=jvmax-jv1+1
      evbyu=ev(jv)/219474.6
      write(nf,13) rbohr(2,jv),evbyu
   13 format(2f10.6)
   12 continue
      rebohr=(rbohr(1,1)+rbohr(2,1))/2.
      zero=0.
      write(nf,13) rebohr,zero
      do 14 jv=1,jvmax
      evbyu=ev(jv)/219474.6
      write(nf,13) rbohr(1,jv),evbyu
   14 continue
c
c generate data for intensity program
      nf1=nf+20
c      do 20 jv=1,jvmax1
c      j2=jv-1
c      write(nf1,21) j2,ev(jv)
c   21 format(i5,f15.3)
c   20 continue
      write(nf1,23) rl,al,bl,ru,au,bu
   23 format(1p6e15.7)
      go to 10
      close(5)
      close(6)
      return
      end
c***********************************************************************
      subroutine coordf(jvmax,d,ev,re,rbohr)
c prints coordinates required by Charlie's Franck-Condon program
c input parameters
c   jvmax=number of vibrational quantum numbers (vmax+1)
c   ev(40)=vibrational energy level, cm-1
c   rbohr(1,80)=outer turning point, rbohr(2,80)=inner turning point, Bohr
c output parameters
c   x(301)=r-coordinates in Bohr
c   y(301)=energy level in atomic units
      dimension ev(80),rbohr(2,80),x(161),y(161)
c
c make x and y
      i=0
      do jv=jvmax,1,-1
        i=i+1
        x(i)=rbohr(2,jv)
        y(i)=ev(jv)/219474.6
      enddo
      i=i+1
      x(i)=re/0.52917715
      y(i)=0.  
      do jv=1,jvmax
        i=i+1
        x(i)=rbohr(1,jv)
        y(i)=ev(jv)/219474.6
      enddo         
      imax=i
c
c write out, Charlie's format, Bohr
      write(6,40)
  40  format(' Charlies format, Bohr')
      write(6,10) (x(i),i=1,5)
  10  format(' x=',5(f12.7,','))
      write(6,20) (x(i),i=6,imax)
  20  format(3x,5(f12.7,','))
      write(6,30) (y(i),i=1,5)
  30  format(' y=',5(f12.7,','))
      write(6,20) (y(i),i=6,imax)
c
c write out, Winifred's format, Angstrom
c      write(6,50)
c  50  format(/' Winiffred format, Angstrom, Hartree'/
c     > '      KNOTS           PUNCH')
c      do i=1,imax
c        write(6,60) x(i),y(i)-d/219474.6
c  60    format(2f15.8)
c      enddo
c      write(6,70)
c  70  format('END')
      return
      end
