c***********************************************************************
      subroutine vibrot(temp,parti,partf)
c calculates vibrational-rotational state energies
c double-precision accuracy needed for 32-bit machines
c uses the method of Liu and Vinokur in calculating partition functions: "Equilibrium
c   Gas Flow Computations. I. Accurate and Efficient Calculation of Equilibrium Gas
c   Properties," AIAA paper 89-1736, June 1989
c written by Chul Park, KAIST, 2005
      dimension dum(20),T(5)
c
      open(5,file='kklloh_h2_vibrot.inp')
      open(6,file='vibrot.out')
c--------------------------------------------------------------------
c read in thermochemical data inputs
      write(6,50)
   50 format(/'INPUT'/)
      call thrinp(3,temp,parti,partf)
c--------------------------------------------------------------------
      close(5)
      close(6)
      return
      end
c***********************************************************************
      subroutine barier (n,iprt,de,we,wexe,re,alphae,redums,nprt)
      real ke,kej
      common /barr1/jrepul
      common /barr/ xrepul,urepul,x(0:600,2),u(0:600,2),
     &              du1(0:600,2),du2(0:600,2)
      common /iter1/  nit,iupopt
      common /iter2/  tol
      common /const/ avogad,boltzm,unigas,plank,spdlgt,pi
      common /spetro/du3(0:600),du4(0:600),bej(0:600),
     &               wej(0:600),xej(0:600),cej(0:600)
c
      h=plank*1.0e7
      c=spdlgt*100.
      cm2erg=h*c
      recm=re*1.0e-8
      ifun =1
c      ifun = 2
      if (wexe .eq. 0.  .or.  alphae .eq. 0.) ifun=2
      ke=4.0*(pi*c*we)**2*redums
      bere2=h/(8.0*pi*pi*c*redums)*1.0e16
      be=bere2/(re*re)
      rcmurt=1./(2.*pi*c*sqrt(redums))
      if (ifun .eq. 1) then
c       hulburt-hirschfelder potential
          f=alphae*we/(6.0*be*be)
          g=8.0*wexe/be
          a=sqrt(ke/(2.0*de*cm2erg))*1.0e-8
          c=1.0-(1.0+f)/(a*re)
          b=2.-(1./c)*(7./12.-1./(a*re)**2*(1.25*(1.+f*(2.0+f))-g/12.0))
      else
c       lippincott potential
          delta=ke*recm*recm/(2.0*de*cm2erg)
          b=1.065
          bsqrtd=b*sqrt(delta)
          a=0.8*(1.-1./bsqrtd)
          c=bsqrtd/re
          alphae=6.*be*be/we*a*bsqrtd
          wexe=0.125*be*(3.+12.*a*bsqrtd+delta*(6.+a*b*b*(15.*b-12.)))
      endif
c compute first repulsive rotational quantum no.
      xrepul=0.1*re
10    continue
      do 20 i=1,nit
      id=3
      if(nprt.gt.3) write(6,40) i,ifun,a,b,c
40    format(' in barier, entering dp2. i,ifun,a,b,c=',2i3,4e10.3)
         call dp2 (id,ifun,a,b,c,de,re,delta,bere2,
     &                        xrepul,db1,db2,db3,dv1,dv2,dv3)
         dxr=-(db2*dv1-db1*dv2)/(db3*dv1-db1*dv3)
         wf=1.0
         if (i .lt. iupopt) wf=float(i)/float(iupopt)
         xrepul=xrepul+wf*dxr
         if (i .gt. 3  .and.  abs(dxr)/xrepul .le. tol) goto 30
20    continue
      if(nprt.gt.2) write (65,610) n,nit,tol,xrepul
610   format (///' ***** fail to find first repulsive quantum no. for'
     &           ' electronic state ',i3,' in ',i5,' iterations *****'/
     &           ' ***** tol = ',e14.7,5x,'xrepul = ',e14.7,' *****'///)
      stop
30    continue
      id=3
      if(nprt.gt.2) write(6,40) i,ifun,a,b,c
      call dp2 (id,ifun,a,b,c,de,re,delta,bere2,
     &                       xrepul,db1,db2,db3,dv1,dv2,dv3)
      repul=0.5*(sqrt(1.+2.*dv1*(xrepul+re)**3/bere2)-1.)
      jrepul=repul
      if (iprt .ge. 1) then
         urp =dv1+repul*(repul+1)*db1
         urpp=dv2+repul*(repul+1)*db2
	 izero=0
	 if(nprt.gt.2) write(6,50) ifun,a,b,c
   50 format(' in barier, entering potent. ifun,a,b,c=',i3,1p4e10.3)
         call potent(izero,izero,izero,ifun,a,b,c,de,re,delta,bere2,
     &     xrepul,urepul)
         urepul=urepul-0.5*repul*(repul+1)*(xrepul+re)*db1
         if(nprt.gt.2) write (65,601)n,jrepul,xrepul,urepul,urp,urpp,a,
     &     b,c,wexe,alphae
601      format (//5x,'electronic state ',i5/
     &    10x,'jrepul=',i14,10x,'xrepul=',1pe14.7,10x,'urepul=',1pe14.7/
     &    10x,'   urp=',1pe14.7,10x,'  urpp=',1pe14.7/
     &    10x,'     a=',1pe14.7,10x,'     b=',1pe14.7,10x,'     c=',
     &    1pe14.7/10x,'  wexe=',1pe14.7,10x,'alphae=',1pe14.7/)
      endif
      if (repul .gt. 600) then
         if(nprt.gt.2) write (65,600) jrepul,600-1
600      format (///' ***** jrepul = ',i5,', it is reset to ',
     &                     '600-1 = ',i5,' *****'///)
         jrepul=600-1
      endif
c compute minimum and maximum of intermolecular potential
      m=1
      x(0,m)=0.0
      u(0,m)=0.0
      do 200 j=1,jrepul
         x(j,m)=x(j-1,m)
         jjp1=j*(j+1)
         do 160 i=1,nit
	 id=2
            call dp2 (id,ifun,a,b,c,de,re,delta,bere2,
     &                       x(j,m),db1,db2,db3,dv1,dv2,dv3)
            dx=-(jjp1*db1+dv1)/(jjp1*db2+dv2)
            wf=1.0
            if (i .lt. iupopt) wf=float(i)/float(iupopt)
            x(j,m)=x(j,m)+wf*dx
            if (i .lt. 3) goto 160
            if (abs(dx/x(j,m)) .lt. tol) goto 170
160      continue
         if(nprt.gt.2) write (65,611) nit,n,j,x(j,m),dx
611   format (' *** fail to find minimum of potential in ',i5,
     &    ' iterations *** n=',i3,' j=',i3,' x=',1pe14.7,' dx=',1pe14.7)
170      continue
200   continue
      m=2
      x(jrepul+1,m)=xrepul+0.1
      x(0,m)=10.0
      u(0,m)=de
      do 300 j=jrepul,1,-1
         x(j,m)=x(j+1,m)
         jjp1=j*(j+1)
	 itwo=2
         do 260 i=1,nit
            call dp2 (itwo,ifun,a,b,c,de,re,delta,bere2,
     &                       x(j,m),db1,db2,db3,dv1,dv2,dv3)
            dx=-(jjp1*db1+dv1)/(jjp1*db2+dv2)
            wf=1.0
            if (i .lt. iupopt) wf=float(i)/float(iupopt)
            x(j,m)=x(j,m)+wf*dx
            if (x(j,m) .lt. xrepul  .or.  x(j,m) .gt. 20.)
     &                                              x(j,m)=x(j+1,m)+1.0
            if (i .lt. 3) goto 260
            if (abs(dx/x(j,m)) .lt. tol) goto 270
260      continue
         if(nprt.gt.2) write (65,612) nit,n,j,x(j,m),dx
612   format (' *** fail to find maximum of potential in ',i5,
     &    ' iterations *** n=',i3,' j=',i3,' x=',1pe14.7,' dx=',1pe14.7)
270      continue
300   continue
      sign=1.
      do 320 m=1,2
         call potent (600,1,jrepul,ifun,a,b,c,de,re,delta,bere2,
     &                                                x(0,m),u(0,m))
         call dp1    (600,0,jrepul,ifun,a,b,c,de,re,delta,bere2,
     &                                     x(0,m),du1(0,m),du2(0,m))
         if (m .eq. 2) sign=-1.
         do 310 j=1,jrepul
            if (abs(du1(j,m)) .lt. 0.1 .and. sign*du2(j,m) .ge. 0.)
     &          goto 310
            if (j .lt. 30) then
                u(j,m)=u(0,m)
                goto 310
            endif
            if(nprt.gt.2) write (65,613) n,m,j,x(j,m),u(j,m),du1(j,m),
     &        du2(j,m)
613   format (' *** incorrect value of min. or max. ***  n ='
     &         i3,' m=',i1,' j=',i3,' x=',1pe12.5,', u=',1pe12.5,
     &         ', du1=',1pe12.5,', du2=',1pe12.5)
310      continue
320   continue
      call dp4 (0,jrepul,ifun,a,b,c,de,re,delta,bere2,x(0,1))
      do 330 j=0,jrepul
         rej=x(j,1)+re
         kej=du2(j,1)*(cm2erg*1.0e16)
         a1j=du3(j)*rej/(3.*du2(j,1))
         a1j2=a1j*a1j
         a2j=du4(j)*rej*rej/(12.*du2(j,1))
         bej(j)=bere2/(rej*rej)
         wej(j)=rcmurt*sqrt(kej)
         xej(j)=1.5*bej(j)/wej(j)*(1.25*a1j2-a2j)
         cej(j)=0.125*bej(j)*(3.*a2j-1.75*a1j2)
c         write (1,700) j,rej,kej,wej(j),wej(j)*xej(j),bej(j),cej(j)
330   continue
700   format (i5,7(e14.6))
      if (iprt .ge. 2) then
         if(nprt.gt.2) write (65,602)
         if(nprt.gt.2) write (65,603) (j,(x(j,m),u(j,m),du1(j,m),
     &     du2(j,m),m=1,2),
     &                    j=0,jrepul)
602   format (///7x,'j',11x,'xmin',10x,'umin',9x,'du1min',8x,'du2min',
     &                  15x,'xmax',10x,'umax',9x,'du1max',8x,'du2max'//)
603   format (1x,i7,6x,1pe14.6,1pe14.6,1pe14.6,1pe14.6,
     &              6x,1pe14.6,1pe14.6,1pe14.6,1pe14.6)
      endif
      return
      end
c*******************************************************************************
      subroutine cspln (np,idrv,f,x,df,dx,b,c,d)
      dimension f(np),x(np),df(2),dx(np),b(np),c(np),d(np)
c
      npm1=np-1
      if (idrv .eq. 1) then
         b( 1)=df(1)/3.
         b(np)=df(2)/3.
         do 110 k=1,npm1
            d(k)=1.0/(x(k+1)-x(k))
            dx(k+1)=d(k)
110         c(k)=(f(k+1)-f(k))*d(k)*d(k)
         do 120 k=2,npm1
120         b(k)=c(k)+c(k-1)
            b(2)=b(2)-b(1)*d(1)
            b(npm1)=b(npm1)-b(np)*d(npm1)
         do 130 k=2,npm1
130         c(k)=2.0*(d(k-1)+d(k))
         call trid (np,dx,c,d,b,2,npm1)
         do 140 k=1,npm1
           ck=(f(k+1)-f(k))*dx(k+1)*dx(k+1)
           c(k)=3.*(ck-dx(k+1)*(b(k+1)+2.*b(k)))
140        d(k)=dx(k+1)*( 3.*dx(k+1)*(b(k+1)+b(k))-2.*ck )
         do 150 k=1,npm1
150         b(k)=3.0*b(k)
      elseif (idrv .eq. 2) then
         c( 1)=df(1)/6.
         c(np)=df(2)/6.
         do 210 k=1,npm1
            d(k)=x(k+1)-x(k)
            dx(k+1)=d(k)
210         b(k)=(f(k+1)-f(k))/d(k)
         do 220 k=2,npm1
220         c(k)=b(k)-b(k-1)
            c(2)=c(2)-c(1)*d(1)
            c(npm1)=c(npm1)-c(np)*d(npm1)
         do 230 k=2,npm1
230         b(k)=2.*(d(k-1)+d(k))
         call trid (np,dx,b,d,c,2,npm1)
         do 240 k=1,npm1
            dxinv=1.0/dx(k+1)
            b(k)=(f(k+1)-f(k))*dxinv-(c(k+1)+2.*c(k))*dx(k+1)
            d(k)=(c(k+1)-c(k))*dxinv
240      continue
         do 250 k=1,npm1
250         c(k)=3.*c(k)
      else
         write (*,600)
600      format (///' ***** you must provide either 1st or 2nd '
     &              'derivatives at the boundaries *****'///)
      endif
      return
      end

c***********************************************************************
      subroutine dpf (iprt,rmass,symfac,g,eel,dzero,we,wexe,re,alphae,
     > t,qi,nprt,  jmax,vmax,gjx,ejx,evx,evjx)
      integer v,vmax
      common /barr1/jrepul
      common /barr/ xrepul,urepul,x(0:600,2),u(0:600,2),
     &              du1(0:600,2),du2(0:600,2)
      common /const/ avogad,boltzm,unigas,plank,spdlgt,pi
      common /iter1/  nit,iupopt
      common /iter2/  tol
      common /spetro/du3(0:600),du4(0:600),bej(0:600),
     &               wej(0:600),xej(0:600),cej(0:600)
      dimension vmax(0:400),gjx(0:400),ejx(0:400),evx(0:400,
     & 0:70),evjx(0:400,0:70)
      data avogad/6.022169e26/,boltzm/1.380622e-23/,plank/6.626196e-34/,
     & spdlgt/2.997925e8/,pi/3.141592/
      unigas=avogad*boltzm*1.0e-3
      np=1
      nel=1
      n=1
c      iprt=1
      nit=500
      iupopt=5
      tol=1.0e-5
      redums=rmass/avogad*1.0e3
c
c  calculate partition function and derivatives for diatomic molecules
c  all units are in si unit, except spectroscopic data in cm-1
c
      e0=we*0.5-wexe*0.25
      de=dzero+e0
      d1=dzero/8067.5
c
c initialize partition functions
         qi    =0.0
c   internal partition function
      cm2k=plank*spdlgt/boltzm*100.
      level=0
c
         gn=symfac*g
	 if(nprt.gt.2) write(6,10) de,we,wexe
   10 format(' in dpf, entering barier. de,we,wexe=',1p3e10.3)
         call barier (n,iprt,de,we,wexe,re,alphae,redums,nprt)
c  zero point energy
cl       e0=cej(0)+wej(0)*(0.5-0.25*xej(0))
clv         e0=cej(0)+wej(0)*(0.5-0.25*xej(0))
         jmax=jrepul
         do 50  j=0,jrepul
          er=bej(0)*j*(j+1)
         cqi=we-alphae*j*(j+1)
            gj=gn*(2.0*j+1.0)
            do 30 v=0,500
               vphalf=v+0.5
c  vibrational energy
            ev=(cqi*vphalf)-wexe*(vphalf*vphalf)
            devdv=cqi-2.0*wexe*vphalf
cl             ev=cej(j)+wej(j)*vphalf*(1.0-xej(j)*vphalf)
cl             devdv=wej(j)*(1.0-2.*xej(j)*vphalf)
clv               ev=cej(j)+wej(j)*vphalf*(1.0-xej(j)*vphalf)
clv               devdv=wej(j)*(1.0-2.*xej(j)*vphalf)
c  rotational-vibrational energy
cd             evj=ev+er
clj            evj=ev+u(j,1)
               evj=ev+u(j,1)
               ejx(j)=u(j,1)
               evx(j,v)=ev
               evjx(j,v)=evj
               if (evj .gt. u(j,2) .or. devdv .lt. 0. .or. evj .lt. 0.)
     &           goto 40
c              converts energy unit from cm-1 to k
               tj=cm2k*(eel+evj-e0)
               level=level+1
c
c partition function 
                  qij=gj*exp(-tj/t)
                  gjx(j)=gj
                  qi=qi+qij
30          continue
40          continue
            vmax(j)=v
50    continue
80    continue
c
      return
      end
c***********************************************************************
      subroutine dp1 (nj,js,je,ifun,a,b,c,de,re,delta,bere2,x,du1,du2)
      dimension x(0:nj),du1(0:nj),du2(0:nj)
c
      if (ifun .eq. 1) then
c       derivatives of hulburt-hirschfelder potential
         du1cst=a*de
         du2cst=-2.*a*du1cst
         b2=2.*b
         b4m2=4.*b-2.
         b6m6=6.*b-6.
         b8m2=8.*b-2.
         do 10 j=js,je
            rinv=1.0/(x(j)+re)
            db1=-2.*j*(j+1)*bere2*rinv*rinv*rinv
            db2=-3.*rinv*db1
            ax=a*x(j)
            ax2=ax*ax
            eax=exp(-ax)
            eax2=eax*eax
            du1(j)=db1+du1cst*(2.*eax
     &                +eax2*(c*ax2*(3.+ax*(b4m2-b2*ax))-2.))
            du2(j)=db2+du2cst*(eax
     &                -eax2*(c*ax*(3.+ax*(b6m6-ax*(b8m2-b2*ax)))+2.))
10       continue
      else
c       derivatives of lippincott potential
         yconst=sqrt(delta/re)
         ab=a*b
         c2=c*c
         do 20 j=js,je
            rinv=1./(x(j)+re)
            db1=-2.*j*(j+1)*bere2*rinv*rinv*rinv
            db2=-3.*rinv*db1
            ybyx=yconst*sqrt(rinv)
            y=ybyx*x(j)
            yrinv=y*rinv
            yp=ybyx-0.5*yrinv
            ypp=rinv*(0.75*yrinv-ybyx)
            y2=y*y
            ey2=exp(-y2)
            aby=ab*y
            z=aby*(1.-ey2)
            zp=ab*(1.-(1.-2.*y2)*ey2)
            zpp=aby*(6.-4.*y2)*ey2
            ecx=exp(-c*x(j))
            zecx=z*ecx
            zpecx=zp*ecx
            dl1=2.*y*ey2-zpecx
            dl2=((2.-4.*y2)*ey2-zpp*ecx)*yp
            du1(j)=db1+de*( dl1*yp+c*zecx )
            du2(j)=db2+de*( dl1*ypp+(dl2+2.*c*zpecx)*yp-c2*zecx )
20       continue
      endif
      return
      end
c***********************************************************************
      subroutine dp2 (id,ifun,a,b,c,de,re,delta,bere2,
     &                      x,db1,db2,db3,dv1,dv2,dv3)
c
      if (ifun .eq. 1) then
c       derivatives of hulburt-hirschfelder potential
         du1cst=a*de
         du2cst=-2.*a*du1cst
         b2=2.*b
         b4m2=4.*b-2.
         b6m6=6.*b-6.
         b8m2=8.*b-2.
          rinv=1./(x+re)
          db1=-2.*bere2*rinv*rinv*rinv
          db2=-3.*rinv*db1
          ax=a*x
          ax2=ax*ax
          eax=exp(-ax)
          eax2=eax*eax
          dv1=du1cst*(2.*eax+eax2*(c*ax2*(3.+ax*(b4m2-b2*ax))-2.))
          dv2=du2cst*(eax-eax2*(c*ax*(3.+ax*(b6m6-ax*(b8m2-b2*ax)))+2.))
          if (id .eq. 3) then
             db3=-4.*rinv*db2
             dv3=-a*du2cst*(eax+eax2*(c*(3.+ax*((12.*b-18.)
     &                 -ax*((36.*b-18.)-ax*((24.*b-4.)-4.*b*ax))))-4.))
          endif
      else
c       derivatives of lippincott potential
         yconst=sqrt(delta/re)
         ab=a*b
         c2=c*c
          rinv=1./(x+re)
          db1=-2.*bere2*rinv*rinv*rinv
          db2=-3.*rinv*db1
          ybyx=yconst*sqrt(rinv)
          y=ybyx*x
          yrinv=y*rinv
          yp=ybyx-0.5*yrinv
          ypp=rinv*(0.75*yrinv-ybyx)
          y2=y*y
          ey2=exp(-y2)
          aby=ab*y
          z=aby*(1.-ey2)
          zp=ab*(1.-(1.-2.*y2)*ey2)
          zpp=aby*(6.-4.*y2)*ey2
          ecx=exp(-c*x)
          zecx=z*ecx
          zpecx=zp*ecx
          dl1=2.*y*ey2-zpecx
          dl2=((2.-4.*y2)*ey2-zpp*ecx)*yp
          dv1=de*( dl1*yp+c*zecx )
          dv2=de*( dl1*ypp+(dl2+2.*c*zpecx)*yp-c2*zecx )
          if (id .eq. 3) then
             db3=-4.*rinv*db2
             yppp=rinv*rinv*(2.25*ybyx-1.875*yrinv)
             dl3=( (8.*y2-12.)*y-ab*(6.-y2*(24.-8.*y2))*ecx )*ey2*yp*yp
             dv3=de*( dl1*yppp +3.*(dl2+c*zpecx)*ypp
     &               +yp*(dl3+3.*c*(zpp*ecx*yp-c*zpecx))+c*c2*zecx )
          endif
      endif
      return
      end
c***********************************************************************
      subroutine dp4(js,je,ifun,a,b,c,de,re,delta,bere2,x)
      common /spetro/du3(0:600),du4(0:600),bej(0:600),
     &               wej(0:600),xej(0:600),cej(0:600)
      dimension x(0:600)
c
      if (ifun .eq. 1) then
c       derivatives of hulburt-hirschfelder potential
         du1cst=a*de
         du2cst=-2.*a*du1cst
         du3cst=-a*du2cst
         du4cst=-a*du3cst
         b2=2.*b
         b4=4.*b
         b8=8.*b
         b4m2=4.*b-2.
         b6m6=6.*b-6.
         b8m2=8.*b-2.
         b12m18=12.*b-18.
         b12m24=12.*b-24.
         b24m4=24.*b-4.
         b36m18=36.*b-18.
         b64m8=64.*b-8.
         b96m72=96.*b-72.
         b14448=144.*b-48.
         do 10 j=js,je
            rinv=1.0/(x(j)+re)
            db1=-2.*j*(j+1)*bere2*rinv*rinv*rinv
            db2=-3.*rinv*db1
            db3=-4.*rinv*db2
            db4=-5.*rinv*db3
            ax=a*x(j)
            ax2=ax*ax
            eax=exp(-ax)
            eax2=eax*eax
c           du1(j)=db1+du1cst*(2.*eax
c    &                +eax2*(c*ax2*(3.+ax*(b4m2-b2*ax))-2.))
c           du2(j)=db2+du2cst*(eax
c    &                -eax2*(c*ax*(3.+ax*(b6m6-ax*(b8m2-b2*ax)))+2.))
            du3(j)=db3+du3cst*(eax+eax2*(c*(3.+ax*(b12m18
     &                 -ax*(b36m18-ax*(b24m4-b4*ax))))-4.))
            du4(j)=db4+du4cst*(eax-eax2*(c*(b12m24-ax*(b96m72
     &                 -ax*(b14448-ax*(b64m8-b8*ax))))+8.))
10       continue
      else
c       derivatives of lippincott potential
         yconst=sqrt(delta/re)
         ab=a*b
         c2=c*c
         c3=c2*c
         c4=c3*c
         do 20 j=js,je
            rinv=1./(x(j)+re)
            db1=-2.*j*(j+1)*bere2*rinv*rinv*rinv
            db2=-3.*rinv*db1
            db3=-4.*rinv*db2
            db4=-5.*rinv*db3
            ybyx=yconst*sqrt(rinv)
            y=ybyx*x(j)
            yrinv=y*rinv
            yp1=       ybyx-0.5000*yrinv
            yp2=(-     ybyx+0.7500*yrinv)*rinv
            yp3=( 2.25*ybyx-1.8750*yrinv)*rinv**2
            yp4=(-7.50*ybyx+6.5625*yrinv)*rinv**3
            y2=y*y
            ey2=exp(-y2)
            aby=ab*y
            abey2=ab*ey2
            abyey2=aby*ey2
            zp0=aby*(1.-ey2)
            zp1=ab-abey2*(1.-y2*  2.        )
            zp2=abyey2*(  6.-y2*  4.        )
            zp3=abey2 *(  6.-y2*(24.- 8.*y2))
            zp4=abyey2*(-60.+y2*(80.-16.*y2))
            ecx=exp(-c*x(j))
            zp0ecx=zp0*ecx
            zp1ecx=zp1*ecx
            zp2ecx=zp2*ecx
            zp3ecx=zp3*ecx
            zp4ecx=zp4*ecx
            yey2=y*ey2
            dl1=   2.       *yey2-zp1ecx
            dl2=(  2.-4.*y2)* ey2-zp2ecx
            dl3=(-12.+8.*y2)*yey2-zp3ecx
            dl4=(-12.+y2*(48.-16.*y2))*ey2-zp4ecx
c           du1(j)=db1+de*(dl1*yp1+c*zp0ecx)
c           du2(j)=db2+de*(dl1*yp2+dl2*yp1**2+2.*c*zp1ecx*yp1-c2*zp0ecx)
            du3(j)=db3+de*(dl1*yp3+3.*(dl2*yp1*yp2+c*zp1ecx*yp2)
     &                    +dl3*yp1**3+3.*c*zp2ecx*yp1**2
     &                    -3.*c2*zp1ecx*yp1+c3*zp0ecx)
            du4(j)=db4+de*(dl1*yp4+dl2*(4.*yp3*yp1+3.*yp2**2)
     &                    +6.*dl3*yp1**2*yp2+dl4*yp1**4
     &                    +4.*c*zp1ecx*yp3+10.*c*zp2ecx*yp1*yp2
     &                    +4.*c*zp3ecx*yp1**3-6.*c2*zp1ecx*yp2
     &                    -6.*c2*zp2ecx*yp1**2+4.*c3*zp1ecx*yp1
     &                    -c4*zp0ecx)
20       continue
      endif
      return
      end
c***********************************************************************
      subroutine partfx(m,specie,amass,t,tv,te,nprt,partf,parti) 
c evaluates partition function
c input parameters: 
c   m=species index: 
c   specie=specie
c   amass=atomic or molecular weight, gram/mol
c   spect(i,j,m)=spectral constants
c     i: 1=degen,2=term,3=we,4=wexe,5=weye,6=re,7=dzero
c        8=be,9=alphae,10=de,11=betae
c     j: max number of electronic levels is 12
c   t=heavy particle translational and rotational temperature, k
c   tv=vibrational temperature, k
c   te=electron temperature, k
c   nprt=print index. print if nprt.ge.2
c calculated parameter: 
c   factr=sysmmetry factor; 0.5 for homo diatomics, 1.0 otherwise
c output parameter: 
c   partf=partition function, including translational component
c   parti=internal partition function, excluding translational component
c
      integer v,vmax
      character*4 spnm,specie,stnm
      common/coma1/spnm(5),nelec(5),spect(20,12,5),spwt(5),rmass(5),
     1 factr(5),stnm(3,20)
      dimension vmax(0:400),gjx(0:400),ejx(0:400),evx(0:400,
     & 0:70),evjx(0:400,0:70)
c
c translational component
      qt=1.878e20*sqrt(amass*t)*amass*t 
c
c diatomic molecule
c      if((im(m).eq.1).and.(ih(m).eq.2)) then
      iprt=1
      q=0.
      do 20 iel=1,nelec(m)
      g=spect(1,iel,m)
      eel=spect(2,iel,m)
      we=spect(3,iel,m)
      wexe=spect(4,iel,m)
      re=spect(6,iel,m)
      dzero=spect(7,iel,m)
      alphae=spect(9,iel,m)
      call dpf(iprt,rmass(m),factr(m),g,eel,dzero,we,wexe,re,alphae,t,
     & qi,nprt,  jmax,vmax,gjx,ejx,evx,evjx)
      if(nprt.gt.2) write(6,40) m,iel,t,qi
   40 format('in partfx, diatomic, m,iel=',2i3,'. T,qi=',1p2e10.3)
c
c write
      write(6,10) (stnm(i,iel),i=1,3),jmax
   10 format(3a4,i10,18x,       '   ! state,jmax')
      do j=0,jmax
        write(6,30) j,vmax(j),gjx(j),ejx(j)
   30   format(2i10,f10.0,f10.2,'   ! j,vmax,gj,ej')
        do v=0,vmax(j)
          write(6,50) v,evx(j,v),evjx(j,v)
   50     format(i5,2f10.2,15x,'   ! v,ev,evj')
        enddo
      enddo
      q=q+qi
   20 continue
c      endif
c
      parti=q
      partf=q*qt
      return
      end 
c***********************************************************************
      subroutine potent (nj,js,je,ifun,a,b,c,de,re,delta,bere2,x,u)
      dimension x(0:nj),u(0:nj)
c
      if (ifun .eq. 1) then
c       hulburt-hirschfelder potential
         do 10 j=js,je
            ax=a*x(j)
            eax=exp(-ax)
            u(j)=de*((1.-eax)**2+c*ax**3*(1.+b*ax)*eax*eax)
     &          +j*(j+1)*bere2/(x(j)+re)**2
10    continue
      else
c       lippincott potential
         yconst=sqrt(delta/re)
         do 20 j=js,je
            y=yconst*x(j)/sqrt(x(j)+re)
            u(j)=de*(1.-exp(-y*y))*(1.-a*b*y*exp(-c*x(j)))
     &          +j*(j+1)*bere2/(x(j)+re)**2
20    continue
      endif
      return
      end
c***********************************************************************
      subroutine thrinp(nprt,temp,parti,partf)
c thermodynamic input subroutine
c input: nprt=print index.
c    nprt=0 no print
c        =1 chemical rate formulas only printed
c        =2 inputs and chemical rate formulas printed
c        =3 further details printed
      character*4 spnm,stnm
      common/coma1/spnm(5),nelec(5),spect(20,12,5),spwt(5),rmass(5),
     1 factr(5),stnm(3,20)
      dimension dum(20)
      character*4 dum
      jsp=1
      read(5,60) (dum(i),i=1,19),spnm(jsp),spwt(jsp),nelec(jsp)
   60 format(19a4/a4,6x,f10.6,9i5)
      if(nprt.gt.0) write(6,120) (dum(i),i=1,19),spnm(jsp),spwt(jsp),
     1 nelec(jsp)
  120 format(19a4/a4,6x,f10.7,1pe10.3,9i5)
   80 format(19a4)
   90 format(8e10.3)
      read(5,80) (dum(i),i=1,19)
      if(nprt.gt.1) write(6,80) (dum(i),i=1,19)
      read(5,90) rmass(jsp),factr(jsp)
      if(nprt.gt.1) write(6,90) rmass(jsp),factr(jsp)
      read(5,80) (dum(i),i=1,19)
      if(nprt.gt.1) write(6,80) (dum(i),i=1,19)
      read(5,80) (dum(i),i=1,19)
      if(nprt.gt.1) write(6,80) (dum(i),i=1,19)
      read(5,80) (dum(i),i=1,19)
      if(nprt.gt.1) write(6,80) (dum(i),i=1,19)
      nlev=nelec(jsp)
      do j=1,nlev
        read(5,80) (stnm(i,j),i=1,3)
        if(nprt.gt.1) write(6,80) (dum(i),i=1,19)
        read(5,130) (spect(i,j,jsp),i=1,7)
  130   format(f10.0,f10.2,f10.4,f10.5,f10.7,f10.5,f10.0)
        if(nprt.gt.1) write(6,130) (spect(i,j,jsp),i=1,7)
        read(5,140) (spect(i,j,jsp),i=8,11)
  140   format(10x,f10.3,f10.6,e10.3,e10.3)
        if(nprt.gt.1) write(6,150) (spect(i,j,jsp),i=8,11)
  150   format(10x,f10.3,f10.6,1pe10.3,e10.3)
      enddo
c
c      write(*,*) 'Temperature=' 
c      read(*,*) temp
      temp=temp
      jsp=1
      spwt1=spwt(jsp)*1.0e3
      call partfx(jsp,spnm(jsp),spwt1,temp,temp,temp,nprt,partf,
     1 parti)
      write(6,320) parti,partf
  320 format(' parti,partf=',1p4e10.3)
      return
      end
c***********************************************************************
      subroutine trid (np,a,b,c,f,nl,nu)
      dimension a(np),b(np),c(np),f(np)
c
      b(nl)=1./b(nl)
      c(nl)=c(nl)*b(nl)
      f(nl)=f(nl)*b(nl)
      do 10 j=nl+1,nu
         b(j)=1./(b(j)-a(j)*c(j-1))
         c(j)=c(j)*b(j)
         f(j)=(f(j)-a(j)*f(j-1))*b(j)
10    continue
c
      do 20 j=nu-1,nl,-1
         f(j)=f(j)-c(j)*f(j+1)
20    continue
      return
      end
