!
!   Программа считает поле в заданном направлении (вперед или назад) на заданное
!   расстояние от источника. Источник может располагаться на любом расстоянии
!   от начала координат. Введен файл прямого доступа и используются функции
!   MS библиотеки для двоичного обмкна
!
!========================================================================================================
MODULE RAMMod

   SAVE
!
!  Mr=bathymetry points, Mz=depth grid, mp=pade terms, MNL=Max N Layers.
!
   Integer,   Parameter :: Mr=2000, Mz=1000000, Mp=10, Mzr=1000, MNL=10

   Integer    Nz, Np, Mdr, Ndr, Ndz, iZ, Nzplt, Lz, iBt, ir( Mzr ),  irot, Nzr
   Integer    NLayers, Nprofile, flHandle, Direction, nOut, nBath;

   Real (Kind=4) ::    Dir( Mzr ), rs, zs, dr, dz, zmplt, eta, eps, omega,rmax, &
                       c0, k0, rDist, rRange, rdr, rp, theta

   Real (Kind=4) ::    RhoB (Mz), rB(Mr), zB(Mr), cW(Mz), cp(Mz), cs(Mz),  &
                       attnp(Mz), attns(Mz), tlg(Mz), RealP(Mz), ImagP(Mz)

   Complex (Kind=4) :: lamW(Mz), lamB(Mz),  muB(Mz),   u(Mz),     v(Mz)

END MODULE RAMMod

!========================================================================================================

MODULE WrkMod

  USE RAMMod, ONLY : Mr,Mz,Mp
  PRIVATE Mr, Mz, Mp

  SAVE
!
!  Mr=bathymetry points, Mz=depth grid, mp=pade terms.
!
   Complex (Kind=4) :: pd1 (Mp),  pd2(Mp),                                                          &
                       t1(6, Mp), t2(6, Mp), t3(6, Mp), t4(6, Mp), t5(6, Mp), t6(6, Mp), t7(6, Mp), &
                       s1(Mz,Mp), s2(Mz,Mp), s3(Mz,Mp), s4(Mz,Mp), s5(Mz,Mp), s6(Mz,Mp), s7(Mz,Mp), &
                       r1(Mz,Mp), r2(Mz,Mp), r3(Mz,Mp), r4(Mz,Mp), r5(Mz,Mp), r6(Mz,Mp), r7(Mz,Mp)

END MODULE WrkMod

!========================================================================================================

MODULE ParMem

   SAVE
!
!  Mr=bathymetry points, Mz=depth grid, mp=pade terms, MNL=Max N Layers.
!
   Integer,   Parameter :: MaxPrfs=100, MaxPrfPnts=300

   Real (Kind=4)  ::   zRhoBSt(MaxPrfPnts,MaxPrfs),   zcWSt(MaxPrfPnts,MaxPrfs),     &
                       zcpSt(MaxPrfPnts,MaxPrfs),     zcsSt(MaxPrfPnts,MaxPrfs),     &
                       zattnpSt(MaxPrfPnts,MaxPrfs),  zattnsSt(MaxPrfPnts,MaxPrfs)
   Real (Kind=4)  ::   RhoBSt(MaxPrfPnts,MaxPrfs),    cWSt(MaxPrfPnts,MaxPrfs),      &
                       cpSt(MaxPrfPnts,MaxPrfs),      csSt(MaxPrfPnts,MaxPrfs),      &
                       attnpSt(MaxPrfPnts,MaxPrfs),   attnsSt(MaxPrfPnts,MaxPrfs),   &
                       RSt(MaxPrfs)

   Integer    NPrfs,   NCwPnts(MaxPrfs),NRhoPnts(MaxPrfs),NCpPnts(MaxPrfs),          &
                       NCsPnts(MaxPrfs),NAtpPnts(MaxPrfs),NAtsPnts(MaxPrfs)

END MODULE ParMem


!========================================================================================================

MODULE MathConst

   SAVE

   REAL    (KIND=4), PARAMETER :: PI  = 3.1415926535898
   COMPLEX (KIND=4), PARAMETER :: CI  = ( 0.0, 1.0 )
   COMPLEX (KIND=8), PARAMETER :: CID = ( 0.0D0, 1.0D0 )

END MODULE MathConst

!========================================================================================================
!
!      Dramatic personae
!
!      All above in modulae
!
!      SUBROUTINE setup  ( g0 )
!      SUBROUTINE profl  ( nprof )
!      SUBROUTINE zread  ( mz,nz,dz,prof,npnts,zPatt,Patt )
!      SUBROUTINE pread  ( zprof,prof,npnt )
!      SUBROUTINE rdProfls
!      SUBROUTINE outpt
!      SUBROUTINE updat
!      SUBROUTINE selfs  ( nu )
!      SUBROUTINE matrc  (mz,nz,mp,np,iz,jz,dz)
!      SUBROUTINE solve  (mz,nz,mp,np,iz,g0)
!      SUBROUTINE rpade  (mp,np,k0,dr,theta,g0)
!      SUBROUTINE epade  (mp,np,ns,ip,k0,c0,dr,irot,theta,g0,nu8)
!      SUBROUTINE deriv  (m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
!      SUBROUTINE gauss  (m,n,a,b)
!      SUBROUTINE pivot  (m,n,i,a,b)
!      SUBROUTINE fndrt  (a,n,z,m)
!      SUBROUTINE guerre (a,n,m,z,err,nter)
!

      PROGRAM rams

	  USE RamMod
	  USE MathConst
	  Use WrkMod

      complex g0
      integer closeFlag

    open(unit=1,status='old',file='rams.3.in')
    open(unit=2,status='unknown',file='tl.nLine.Txt', recl=1024)
	open(unit=3,status='unknown',file='TLrz',form='unformatted',access='stream')
	open(unit=9,status='unknown',file='RePrz',form='unformatted',access='stream')
	open(unit=10,status='unknown',file='ImPrz',form='unformatted',access='stream')


      nOut = 0
      call setup( g0 )

!     March the seismo-acoustic field out in range.

      DO While( rDist.lt.rmax .and. rRange.gt.0 )
         call updat
         call solve(mz,nz,mp,np,iz,g0)
         rDist=rDist+dr
         if( Direction .gt. 0 ) rRange=rRange+dr
         if( Direction .lt. 0 ) rRange=rRange-dr
         call outpt
!        nOut=nOut+1
      END DO

      close(1)
      close(2)
      close(3)
	close(9)
	close(10)

      open(unit=4,status='unknown',file='dbnds.txt')
      write( 4,101 ) nOut, rmax, ndr*dr, lz, zmplt, dz*ndz
      close(4)

      stop
101   Format(1X,'Num Range Records    : ', I8,/,    &
             1X,'Max Range            : ', F9.0,/,  &
             1X,'Step Along Range     : ', F12.3,/, &
             1X,'Num Horizontal Lines : ', I8,/,    &
             1X,'Max Depth            : ', F9.0,/,  &
             1X,'Step Along Depth     : ', F12.3    &
            )
      END

!========================================================================================================
!
!     Initialize the parameters, seismo-acoustic field, and matrices.
!
      SUBROUTINE setup( g0 )

	  USE MathConst
	  USE RamMod
	  USE ParMem

	  Real zr( Mzr )
      complex g0,nu

      zr(1:Mzr) = -1
      read(1,*)
      read(1,*)  freq, rs, zs, Direction
      read(1,*) ( zr(i),i=1,Mzr )
      i = 0
      DO While (.true.)
         i=i+1
         IF(zr(i).lt.0.0) EXIT
      END DO
      Nzr = i-1

      read(1,*)rmax,dr,ndr
      read(1,*)zmax,dz,ndz,zmplt
      read(1,*)c0,np,irot,theta
!
      i=0
      DO While (.true.)
         i=i+1
         read(1,*)rb(i),zb(i)
         IF(rb(i).lt.0.0) EXIT
      END DO

      zb(i)=zb(i-1)
      rb(i)=rb(i-1)+rs+rmax+2.0*dr
      nBath=i;
!
      nu=ci
      eta=1.0/(40.0*pi*alog10(exp(1.0)))
      eps=1.0e-20
      mdr=0
      rDist=dr
      if( Direction .gt. 0 ) rRange=rs+dr
      if( Direction .lt. 0 ) rRange=rs-dr
      omega=2.0*pi*freq

      DO i=1, Nzr
         ri=1.0+zr(i)/dz
         ir(i)=ifix(ri)
         dir(i)=ri-float( ir(i) )
      END DO

      k0=omega/c0
      nz=zmax/dz-0.5
      nzplt=zmplt/dz-0.5

!     Define bottom index
      iBt=1
      Do while( rb(iBt) .lt. rs .and. rb(iBt+1) .le. rs .and. iBt+1 .lt. nBath )
         iBt=iBt+1
      End Do
!     z=zb(1)
!     iZ=z/dz
!     write (*,*) 'primary z(1), iBz and iZ: ',z, iBt, iZ

      z=zb(iBt)+(rs-rb(iBt))*(zb(iBt+1)-zb(iBt))/(rb(iBt+1)-rb(iBt))
!     z=zb(iBt)
      iZ=z/dz
!     write (*,*) 'then z(iBt), iBt and iZ: ',z, iBt, iZ
!
      DO i=1,2*nz+4
         u(i)=0.0
         v(i)=0.0
      END DO

      lz = (nzplt-ndz-1)/ndz + 1
!     write(3)lz
!
!     The initial profiles and starting field.
!
      call rdProfls

!     Define profile index
      Nprofile=1
      Do while( rSt(Nprofile) .lt. rs .and. rSt(Nprofile+1) .le. rs .and. Nprofile .lt. NPrfs )
         Nprofile=Nprofile+1
      End Do

      call profl(Nprofile)
      call selfs( nu )
      call outpt
!
!     The propagation matrices.
!
      call epade(mp,np,1,1,k0,c0,dr,irot,theta,g0,nu)
      call matrc(mz,nz,mp,np,iz,iz,dz)
!
      RETURN
      END

!========================================================================================================
!
!     Set up the profiles.
!
      SUBROUTINE profl(nprof)
	  USE MathConst
	  USE RamMod
	  USE ParMem
	  integer nprof

      call zread(mz,nz,dz,cw, NcwPnts(nprof), zcwSt(:,nprof),cwSt(:,nprof))
      call zread(mz,nz,dz,cp, NcpPnts(nprof), zcpSt(:,nprof),cpSt(:,nprof))
      call zread(mz,nz,dz,cs, NcsPnts(nprof), zcsSt(:,nprof),csSt(:,nprof))
      call zread(mz,nz,dz,rhob,  NRhoPnts(nprof), zRhoBSt(:,nprof),RhoBSt(:,nprof))
      call zread(mz,nz,dz,attnp, NAtpPnts(nprof), zattnpSt(:,nprof),attnpSt(:,nprof))
      call zread(mz,nz,dz,attns, NAtsPnts(nprof), zattnsSt(:,nprof),attnsSt(:,nprof))

      if ( nprof.eq.NPrfs ) then
         if( Direction .gt. 0 )then
             rp=2.0*rb(nBath)+rs
         else
             rp=Rst(nprof)
         endif
      else
         if( Direction .gt. 0 ) then
             rp=Rst(nprof+1)+rs
         else
             rp=Rst(nprof)
         endif
      End If
!
    1 Continue

      DO i=1,nz+2
         lamw(i)=cw(i)**2
         lamb(i)=rhob(i)*((cp(i)/(1.0+ci*eta*attnp(i)))**2-  &
                 2.0*(cs(i)/(1.0+ci*eta*attns(i)))**2)
         mub(i)=rhob(i)*(cs(i)/(1.0+ci*eta*attns(i)))**2
      END DO
!
      RETURN
      END

!========================================================================================================
!
!     Read up all the profiles.
!

      SUBROUTINE rdProfls
      USE ParMem
      integer i

      Do i=1, MaxPrfs
         Rst(i)=0
      End Do

      i=0
      DO While(.true.)
         i=i+1
         call pread(zcwSt(:,i),cwSt(:,i),NcwPnts(i))
         call pread(zcpSt(:,i),cpSt(:,i),NcpPnts(i))
         call pread(zcsSt(:,i),csSt(:,i),NcsPnts(i))
         call pread(zRhoBSt(:,i),RhoBSt(:,i),NRhoPnts(i))
         call pread(zattnpSt(:,i),attnpSt(:,i),NAtpPnts(i))
         call pread(zattnsSt(:,i),attnsSt(:,i),NAtsPnts(i))
         read(1,*,END=1) Rst(i+1)
      END DO
!
    1 Continue
      NPrfs = i
!

      RETURN
      END

!========================================================================================================
!
!     Profile reader.
!
      SUBROUTINE pread(zprof,prof,npnt)
      USE ParMem, ONLY: MaxPrfs, MaxPrfPnts
      real prof(MaxPrfPnts),zprof(MaxPrfPnts)

      DO i=1,MaxPrfPnts
         prof(i)  =-1.0
         zprof(i) =-1.0
      END DO
      npnt = 0

      DO While(.true.)
         read(1,*)zi,profi
         IF(zi.lt.0.0) EXIT
         npnt = npnt+1
         zprof(npnt)=zi
         prof (npnt)=profi
      END DO

      RETURN
      END



!========================================================================================================
!
!     Profile reader and interpolator.
!
      SUBROUTINE zread(mz,nz,dz,prof,npnts,zPatt,Patt)
      USE ParMem, ONLY: MaxPrfs, MaxPrfPnts
      real Patt(MaxPrfPnts),zPatt(MaxPrfPnts)
      real prof(mz)

      DO i=1,nz+2
         prof(i)=-1.0
      END DO

      zi=zPatt(1)
      profi=Patt(1)
      prof(1)=profi
      i=1.5+zi/dz
      prof(i)=profi
      iold=i

      DO j=2,npnts
         zi=zPatt(j)
         profi=Patt(j)
         i=1.5+zi/dz
         IF(i.eq.iold)i=i+1
         prof(i)=profi
         iold=i
      END DO

      prof(nz+2)=prof(i)
      i=1
      j=1

      DO While ( j.lt.nz+2 )
         i=i+1
         IF (prof(i).lt.0.0) Cycle
         IF(i-j.ne.1) THEN
            DO k=j+1,i-1
               prof(k)=prof(j)+float(k-j)*(prof(i)-prof(j))/float(i-j)
            END DO
         END IF
         j=i
      END DO

      RETURN
      END


!========================================================================================================
!
!     Output transmission loss.
!
      SUBROUTINE outpt

	  USE RamMod

	  complex ur
	  real*8  uNorm
	  real    tl(Nzr), vTL
	  integer*4 i,j, nb
      character(1) bt

	  bt=char(13)

      mdr=mdr+1

      IF(mdr.eq.ndr)THEN

	     uNorm = 0

		 DO i=1,Nz
		    uNorm = uNorm + abs(u(i)*conjg(u(i)))
		    IF( i .EQ. iZ ) vTL = uNorm*dz
		 END DO
		 uNorm = DSQRT( uNorm )
		 IF( uNorm  > 10000 .OR. uNorm  .EQ. Inf .OR.  uNorm  .EQ. NaN ) then
		    STOP
		 ENDIF
		 DO i=1,Nzr
            ur=(1.0-dir(i))*u(2*ir(i)-1)+dir(i)*u(2*ir(i)+1)
            tl(i)=20.0*alog10(cabs(ur)+eps)-10.0*alog10(rDist+eps)
         END DO
         vTL=10.0*alog10(vTL+eps)-10.0*alog10(rDist+eps)

         write(2,100) rDist,rRange, ( tl(i),i=1,Nzr ),vTL,uNorm, dz*iz
         write(*,101) rDist,rRange, tl(1), vTL, uNorm, bt
         mdr=0
!
         j=0
         iflag=1
         DO i=1+ndz,nzplt,ndz
         ur=u(2*i-1)
         j=j+1
         tlg(j)=20.0*alog10(cabs(ur)+eps)-10.0*alog10(rDist+eps)
	 realP(j)=real(ur)/sqrt(rDist)
	 imagP(j)=aimag(ur)/sqrt(rDist)
!
!        Mark the ocean bottom.
!
!         IF((i.gt.iz).and.(iflag.eq.1))THEN
!         tlg(j)=0.0
!         iflag=0
!         END IF
!
         END DO

         nOut=nOut+1
         write(3)(tlg(j),j=1,lz)
	 write(9)(realP(j),j=1,lz)
	 write(10)(imagP(j),j=1,lz)
      END IF
!
      RETURN
100   Format(2x,2F12.2,100F8.2,4x,F8.2,F16.4)
101   Format(2x,2F12.2,F8.2,4x,F8.2,F16.4,A1$)
!100   Format(2x,2F12.2,F8.2,F16.4,A1$)
!101   Format(2x,2F12.2,<Nzr>F8.2,F16.4)
      END

!========================================================================================================
!
!     Matrix updates and energy conservation.
!
      SUBROUTINE updat
      USE RamMod
!
!     Varying bathymetry.
!
      jz=iz
      if( Direction .gt. 0 ) z=zb(iBt)+(rRange+0.5*dr-rb(iBt))*(zb(iBt+1)-zb(iBt))/(rb(iBt+1)-rb(iBt))
      if( Direction .lt. 0 ) z=zb(iBt)+(rRange-0.5*dr-rb(iBt))*(zb(iBt+1)-zb(iBt))/(rb(iBt+1)-rb(iBt))
!     z=zb(iBt)+(rDist+0.5*dr-rb(iBt))*(zb(iBt+1)-zb(iBt))/(rb(iBt+1)-rb(iBt))
      iz=z/dz

      IF(iz.ne.jz)THEN

		 Call matrc(mz,nz,mp,np,iz,jz,dz)
!
!        An approximate energy-flux correction.
!
         IF(iz.lt.jz)THEN
           fact=sqrt(cw(iz)**3)/sqrt(rhob(iz)*cp(iz)**3)
           DO i=iz+1,jz
              u(2*i+1)=u(2*i+1)*fact
           END DO
         Else
           fact=sqrt(rhob(iz)*cp(iz)**3)/sqrt(cw(iz)**3)
           DO i=jz+1,iz
              u(2*i+1)=u(2*i+1)*fact
           END DO
         END IF
!
      END IF
!
      if( Direction .gt. 0 .and. rRange.ge.rb(iBt+1) ) iBt=iBt+1
      if( Direction .lt. 0 .and. rRange.le.rb(iBt) .and. iBt.gt.1 ) iBt=iBt-1
!     IF(rDist.ge.rb(iBt+1))iBt=iBt+1
!
!     Varying profiles.
!
!      IF(rDist.ge.rp)THEN
!         Nprofile=Nprofile+1
!         call profl(Nprofile)
!         call matrc(mz,nz,mp,np,iz,iz,dz)
!      END IF
      IF(Direction .gt. 0 .and. rRange.ge.rp )THEN
         Nprofile=Nprofile+1
         call profl(Nprofile)
         call matrc(mz,nz,mp,np,iz,iz,dz)
      END IF
      IF(Direction .lt. 0 .and. rRange.le.rp .and. Nprofile.gt.1 )THEN
         Nprofile=Nprofile-1
         call profl(Nprofile)
         call matrc(mz,nz,mp,np,iz,iz,dz)
      END IF
!
      RETURN
      END

!========================================================================================================
!
!     The self starter.
!
      SUBROUTINE selfs( nu )
      USE MathConst
	  USE RamMod
	  USE WrkMod
      complex g0,nu
!
!     Conditions for the delta function.
!
      si=1.0+zs/dz
      is=ifix(si)
      dis=si-float(is)
      u(2*is-1)=(1.0-dis)*sqrt(2.0*pi/k0)/dz
      u(2*is+1)=dis*sqrt(2.0*pi/k0)/dz
!
!     Divide the delta function by (1-nu*X)**2 to get a smooth rhs.
!
      pd1(1)=0.0
      pd2(1)=-nu
      g0=(1.0,0.0)
      call matrc(mz,nz,mp,1,iz,iz,dz)
      call solve(mz,nz,mp,1,iz,g0)
      call solve(mz,nz,mp,1,iz,g0)
!
!     Apply the operator (1-nu*X)**2*(1+X)**(-1/4)*exp(ci*k0*r*sqrt(1+X)).
!
      call epade(mp,np,1,0,k0,c0,dr,0,0.0,g0,nu)
      call matrc(mz,nz,mp,np,iz,iz,dz)
      call solve(mz,nz,mp,np,iz,g0)
!
      RETURN
      END

!========================================================================================================
!
!     The heptadiagonal matrices.
!
      SUBROUTINE matrc(mz,nz,mp,np,iz,jz,dz)
	  USE WrkMod
	  USE RamMod, ONLY : k0, omega, rhob, lamw, lamb, mub

      complex l1,l2,l3,l4,l5,l6,l7,m1,m2,m3,m4,m5,m6,m7,a11,a12,a21,  &
              a22,a23,a31,a32,a33,a34,a35   !  ,lamb(mz),mub(mz)

!      real k0,rhob(mz),lamw(mz)
!
      If(iz.eq.jz)Then  !   New matrices when iz.eq.jz.
         ia=2
         ib=nz-2
      Else              !   Updated matrices when iz.ne.jz.
         ia=min(iz,jz)
         ib=max(iz,jz)
      End If
!
      DO i=ia-1,iz      !   Eq. (4) of jasa 86, 1459-1464 in fluid layer.

         l1=0.0
         l2=(lamw(i)+lamw(i+1))/12.0
         l3=0.0
         l4=(lamw(i)+6.0*lamw(i+1)+lamw(i+2))/12.0
         l5=0.0
         l6=(lamw(i+1)+lamw(i+2))/12.0
         l7=0.0
         m1=0.0
         m2=lamw(i+1)/dz**2+omega**2/6.0+  &
            0.5*(lamw(i)-lamw(i+1))/dz**2+  &
            0.5*(lamw(i)-lamw(i+1))/dz**2
         m3=0.0
         m4=-2.0*lamw(i+1)/dz**2+2.0*omega**2/3.0+  &
             0.5*(2.0*lamw(i+1)-lamw(i)-lamw(i+2))/dz**2+  &
             0.5*(lamw(i)-2.0*lamw(i+1)+lamw(i+2))/dz**2
         m5=0.0
         m6=lamw(i+1)/dz**2+omega**2/6.0+  &
            0.5*(lamw(i+2)-lamw(i+1))/dz**2+  &
            0.5*(lamw(i+2)-lamw(i+1))/dz**2
         m7=0.0
!
         DO j=1,np

            r1(2*i-1,j)=k0**2*l1+pd2(j)*(m1-k0**2*l1)
            r2(2*i-1,j)=k0**2*l2+pd2(j)*(m2-k0**2*l2)
            r3(2*i-1,j)=k0**2*l3+pd2(j)*(m3-k0**2*l3)
            r4(2*i-1,j)=k0**2*l4+pd2(j)*(m4-k0**2*l4)
            r5(2*i-1,j)=k0**2*l5+pd2(j)*(m5-k0**2*l5)
            r6(2*i-1,j)=k0**2*l6+pd2(j)*(m6-k0**2*l6)
            r7(2*i-1,j)=k0**2*l7+pd2(j)*(m7-k0**2*l7)
!
            s1(2*i-1,j)=k0**2*l1+pd1(j)*(m1-k0**2*l1)
            s2(2*i-1,j)=k0**2*l2+pd1(j)*(m2-k0**2*l2)
            s3(2*i-1,j)=k0**2*l3+pd1(j)*(m3-k0**2*l3)
            s4(2*i-1,j)=k0**2*l4+pd1(j)*(m4-k0**2*l4)
            s5(2*i-1,j)=k0**2*l5+pd1(j)*(m5-k0**2*l5)
            s6(2*i-1,j)=k0**2*l6+pd1(j)*(m6-k0**2*l6)
            s7(2*i-1,j)=k0**2*l7+pd1(j)*(m7-k0**2*l7)

         END DO  !  j=1,np

      END DO  !  i=ia-1,iz
!
      DO i=ia-1,iz          !   Eq. (2) of jasa 86, 1459-1464 in fluid layer.

         m1=(-1.0/6.0)*(lamw(i)+2.0*lamw(i+1))/dz+  &
             (1.0/6.0)*(-lamw(i)+lamw(i+1))/dz
         m2=0.0
         m3=(1.0/6.0)*(lamw(i)-lamw(i+2))/dz+  &
            (1.0/3.0)*(-lamw(i)+lamw(i+2))/dz
         m4=omega**2
         m5=(1.0/6.0)*(lamw(i+2)+2.0*lamw(i+1))/dz+  &
            (1.0/6.0)*(lamw(i+2)-lamw(i+1))/dz
         m6=0.0
!
         DO j=1,np

            r1(2*i,j)=m1
            r2(2*i,j)=m2
            r3(2*i,j)=m3
            r4(2*i,j)=m4
            r5(2*i,j)=m5
            r6(2*i,j)=m6
            r7(2*i,j)=m7
            s1(2*i,j)=0.0
            s2(2*i,j)=0.0
            s3(2*i,j)=0.0
            s4(2*i,j)=0.0
            s5(2*i,j)=0.0
            s6(2*i,j)=0.0
            s7(2*i,j)=0.0

         END DO  !  j=1,np

      END DO  !  i=ia-1,iz
!
!     Eq. (4) of jasa 86, 1459-1464 in solid layer.
!
      DO i=iz+1,ib+2

         l1=0.0
         l2=(lamb(i)+2.0*mub(i)+lamb(i+1)+2.0*mub(i+1))/12.0
         l3=(2.0/6.0)*(-mub(i)+mub(i+1))/dz
         l4=(lamb(i)+2.0*mub(i)+6.0*lamb(i+1)+12.0*mub(i+1)+  &
             lamb(i+2)+2.0*mub(i+2))/12.0
         l5=(2.0/3.0)*(-mub(i)+mub(i+2))/dz
         l6=(lamb(i+1)+2.0*mub(i+1)+lamb(i+2)+2.0*mub(i+2))/12.0
         l7=(2.0/6.0)*(mub(i+2)-mub(i+1))/dz
         m1=0.0
         m2=(lamb(i+1)+2.0*mub(i+1))/dz**2+  &
             omega**2*(rhob(i)+rhob(i+1))/12.0+  &
             0.5*(lamb(i)+2.0*mub(i)-lamb(i+1)-2.0*mub(i+1))/dz**2+  &
             0.5*(lamb(i)-lamb(i+1))/dz**2
         m3=(omega**2/6.0)*(-rhob(i)+rhob(i+1))/dz+  &
             2.0*(-mub(i)+mub(i+1))/dz**3
         m4=-2.0*(lamb(i+1)+2.0*mub(i+1))/dz**2+  &
             omega**2*(rhob(i)+6.0*rhob(i+1)+rhob(i+2))/12.0+  &
             0.5*(2.0*lamb(i+1)+4.0*mub(i+1)-lamb(i)-2.0*mub(i)-  &
             lamb(i+2)-2.0*mub(i+2))/dz**2+  &
             0.5*(lamb(i)-2.0*lamb(i+1)+lamb(i+2))/dz**2
         m5=(omega**2/3.0)*(-rhob(i)+rhob(i+2))/dz+  &
             2.0*(mub(i)-mub(i+2))/dz**3
         m6=(lamb(i+1)+2.0*mub(i+1))/dz**2+  &
             omega**2*(rhob(i+1)+rhob(i+2))/12.0+  &
             0.5*(lamb(i+2)+2.0*mub(i+2)-lamb(i+1)-2.0*mub(i+1))/dz**2+  &
             0.5*(lamb(i+2)-lamb(i+1))/dz**2
         m7=(omega**2/6.0)*(rhob(i+2)-rhob(i+1))/dz+  &
             2.0*(mub(i+2)-mub(i+1))/dz**3
!
         DO j=1,np

            r1(2*i-1,j)=k0**2*l1+pd2(j)*(m1-k0**2*l1)
            r2(2*i-1,j)=k0**2*l2+pd2(j)*(m2-k0**2*l2)
            r3(2*i-1,j)=k0**2*l3+pd2(j)*(m3-k0**2*l3)
            r4(2*i-1,j)=k0**2*l4+pd2(j)*(m4-k0**2*l4)
            r5(2*i-1,j)=k0**2*l5+pd2(j)*(m5-k0**2*l5)
            r6(2*i-1,j)=k0**2*l6+pd2(j)*(m6-k0**2*l6)
            r7(2*i-1,j)=k0**2*l7+pd2(j)*(m7-k0**2*l7)
!
            s1(2*i-1,j)=k0**2*l1+pd1(j)*(m1-k0**2*l1)
            s2(2*i-1,j)=k0**2*l2+pd1(j)*(m2-k0**2*l2)
            s3(2*i-1,j)=k0**2*l3+pd1(j)*(m3-k0**2*l3)
            s4(2*i-1,j)=k0**2*l4+pd1(j)*(m4-k0**2*l4)
            s5(2*i-1,j)=k0**2*l5+pd1(j)*(m5-k0**2*l5)
            s6(2*i-1,j)=k0**2*l6+pd1(j)*(m6-k0**2*l6)
            s7(2*i-1,j)=k0**2*l7+pd1(j)*(m7-k0**2*l7)
         END DO  !  j=1,np

      END DO  !  i=iz+1,ib+2
!
!     Eq. (2) of jasa 86, 1459-1464 in solid layer.
!
      Do i=iz+1,ib+2

                 l1=0.0
         l2=(mub(i)+mub(i+1))/12.0
         l3=0.0
         l4=(mub(i)+6.0*mub(i+1)+mub(i+2))/12.0
         l5=0.0
         l6=(mub(i+1)+mub(i+2))/12.0
         l7=0.0
         m1=(-1.0/6.0)*(lamb(i)+mub(i)+2.0*lamb(i+1) + 2.0*mub(i+1))/dz+  &
             (1.0/6.0)*(-lamb(i)+lamb(i+1))/dz
         m2=mub(i+1)/dz**2 + omega**2*(rhob(i)+rhob(i+1))/12.0+  &
            (mub(i)-mub(i+1))/dz**2
         m3=(1.0/6.0)*(lamb(i)+mub(i)-lamb(i+2)-mub(i+2))/dz+  &
            (1.0/3.0)*(-lamb(i)+lamb(i+2))/dz
         m4=-2.0*mub(i+1)/dz**2+ omega**2*(rhob(i)+6.0*rhob(i+1)+rhob(i+2))/12.0+  &
            (2.0*mub(i+1)-mub(i)-mub(i+2))/dz**2
         m5=(1.0/6.0)*(lamb(i+2)+mub(i+2)+2.0*lamb(i+1) + 2.0*mub(i+1))/dz+  &
            (1.0/6.0)*(lamb(i+2)-lamb(i+1))/dz
         m6=mub(i+1)/dz**2+ omega**2*(rhob(i+1)+rhob(i+2))/12.0+  &
           (mub(i+2)-mub(i+1))/dz**2
         m7=0.0
!
         Do j=1,np
            r1(2*i,j)=k0**2*l1+pd2(j)*(m1-k0**2*l1)
            r2(2*i,j)=k0**2*l2+pd2(j)*(m2-k0**2*l2)
            r3(2*i,j)=k0**2*l3+pd2(j)*(m3-k0**2*l3)
            r4(2*i,j)=k0**2*l4+pd2(j)*(m4-k0**2*l4)
            r5(2*i,j)=k0**2*l5+pd2(j)*(m5-k0**2*l5)
            r6(2*i,j)=k0**2*l6+pd2(j)*(m6-k0**2*l6)
            r7(2*i,j)=k0**2*l7+pd2(j)*(m7-k0**2*l7)
!
            s1(2*i,j)=k0**2*l1+pd1(j)*(m1-k0**2*l1)
            s2(2*i,j)=k0**2*l2+pd1(j)*(m2-k0**2*l2)
            s3(2*i,j)=k0**2*l3+pd1(j)*(m3-k0**2*l3)
            s4(2*i,j)=k0**2*l4+pd1(j)*(m4-k0**2*l4)
            s5(2*i,j)=k0**2*l5+pd1(j)*(m5-k0**2*l5)
            s6(2*i,j)=k0**2*l6+pd1(j)*(m6-k0**2*l6)
            s7(2*i,j)=k0**2*l7+pd1(j)*(m7-k0**2*l7)
         End Do  ! j=1,np

      End Do  ! i=iz+1,ib+2
!
!     Continuity conditions for a fluid-solid interface.
!
      a11=lamw(iz)/lamw(iz+2)
      a12=-2.0*dz*omega**2/lamw(iz+2)
!
      a21=-dz*lamw(iz+1)/mub(iz+1)
      a22= dz*lamb(iz+1)/mub(iz+1)
      a23= 1.0
!
      a31=-2.0*(mub(iz+1)+mub(iz+2)) * lamw(iz+1)/(mub(iz+1)*lamb(iz))
      a32= 2.0*(mub(iz+1)+mub(iz+2)) * lamb(iz+1)/(mub(iz+1)*lamb(iz))
      a33= 2.0*rhob(iz+1)*dz*omega**2/lamb(iz) -  &
           2.0*(mub(iz)+2.0*mub(iz+1)+mub(iz+2))/(dz*lamb(iz))
      a34=lamb(iz+2)/lamb(iz)
      a35=2.0*(mub(iz)+2.0*mub(iz+1)+mub(iz+2))/(dz*lamb(iz))
!
      Do j=1,np
!
         i=2*iz-1
         r2(i,j)=r2(i,j)+a11*r6(i,j)
         s2(i,j)=s2(i,j)+a11*s6(i,j)
         r7(i,j)=a12*r6(i,j)
         s7(i,j)=a12*s6(i,j)
         r6(i,j)=0.0
         s6(i,j)=0.0
!
         i=2*iz
         r1(i,j)=r1(i,j)+a11*r5(i,j)
         s1(i,j)=s1(i,j)+a11*s5(i,j)
         r6(i,j)=a12*r5(i,j)
         s6(i,j)=a12*s5(i,j)
         r5(i,j)=0.0
         s5(i,j)=0.0
!
         i=2*iz+1
         r5(i,j)=r5(i,j)+a33*r2(i,j)
         s5(i,j)=s5(i,j)+a33*s2(i,j)
         r6(i,j)=r6(i,j)+a34*r2(i,j)
         s6(i,j)=s6(i,j)+a34*s2(i,j)
         r7(i,j)=r7(i,j)+a35*r2(i,j)
         s7(i,j)=s7(i,j)+a35*s2(i,j)
         r4(i,j)=r4(i,j)+a32*r2(i,j)
         s4(i,j)=s4(i,j)+a32*s2(i,j)
         r2(i,j)=a31*r2(i,j)
         s2(i,j)=a31*s2(i,j)
!
         r2(i,j)=r2(i,j)+a21*r3(i,j)
         s2(i,j)=s2(i,j)+a21*s3(i,j)
         r4(i,j)=r4(i,j)+a22*r3(i,j)
         s4(i,j)=s4(i,j)+a22*s3(i,j)
         r7(i,j)=r7(i,j)+a23*r3(i,j)
         s7(i,j)=s7(i,j)+a23*s3(i,j)
         r3(i,j)=0.0
         s3(i,j)=0.0
!
         i=2*iz+2
         r3(i,j)=r3(i,j)+a32*r1(i,j)
         s3(i,j)=s3(i,j)+a32*s1(i,j)
         r4(i,j)=r4(i,j)+a33*r1(i,j)
         s4(i,j)=s4(i,j)+a33*s1(i,j)
         r5(i,j)=r5(i,j)+a34*r1(i,j)
         s5(i,j)=s5(i,j)+a34*s1(i,j)
         r6(i,j)=r6(i,j)+a35*r1(i,j)
         s6(i,j)=s6(i,j)+a35*s1(i,j)
         r1(i,j)=a31*r1(i,j)
         s1(i,j)=a31*s1(i,j)
!
         r1(i,j)=r1(i,j)+a21*r2(i,j)
         s1(i,j)=s1(i,j)+a21*s2(i,j)
         r3(i,j)=r3(i,j)+a22*r2(i,j)
         s3(i,j)=s3(i,j)+a22*s2(i,j)
         r6(i,j)=r6(i,j)+a23*r2(i,j)
         s6(i,j)=s6(i,j)+a23*s2(i,j)
         r2(i,j)=0.0
         s2(i,j)=0.0

      End Do  !  j=1,np
!
!     The matrix decomposition.
!
      Do j=1,np
!
         Do i=2*ia-3,2*iz

            If(i.gt.3) Then
               r2(i,j)=r2(i,j)-r1(i,j)*r5(i-3,j)
               r3(i,j)=r3(i,j)-r1(i,j)*r6(i-3,j)
               r4(i,j)=r4(i,j)-r1(i,j)*r7(i-3,j)
            End If
            If(i.gt.2) Then
               r3(i,j)=r3(i,j)-r2(i,j)*r5(i-2,j)
               r4(i,j)=r4(i,j)-r2(i,j)*r6(i-2,j)
               r5(i,j)=r5(i,j)-r2(i,j)*r7(i-2,j)
            End If
            If(i.gt.1) Then
               r4(i,j)=r4(i,j)-r3(i,j)*r5(i-1,j)
               r5(i,j)=r5(i,j)-r3(i,j)*r6(i-1,j)
               r6(i,j)=r6(i,j)-r3(i,j)*r7(i-1,j)
            End If
            r4(i,j)=1.0/r4(i,j)
            r5(i,j)=r5(i,j)*r4(i,j)
            r6(i,j)=r6(i,j)*r4(i,j)
            r7(i,j)=r7(i,j)*r4(i,j)

         End Do ! i=2*ia-3,2*iz
!
         i1=2*nz-2
         i2=2*nz-1
         i3=2*nz

         Do i=2*ib+4,2*iz+1,-1

            If(i.LT.i1) Then
               r6(i,j)=r6(i,j)-r7(i,j)*r3(i+3,j)
               r5(i,j)=r5(i,j)-r7(i,j)*r2(i+3,j)
               r4(i,j)=r4(i,j)-r7(i,j)*r1(i+3,j)
            End If
            If(i.LT.i2) Then
               r5(i,j)=r5(i,j)-r6(i,j)*r3(i+2,j)
               r4(i,j)=r4(i,j)-r6(i,j)*r2(i+2,j)
               r3(i,j)=r3(i,j)-r6(i,j)*r1(i+2,j)
            End If
            If(i.LT.i3) Then
               r4(i,j)=r4(i,j)-r5(i,j)*r3(i+1,j)
               r3(i,j)=r3(i,j)-r5(i,j)*r2(i+1,j)
               r2(i,j)=r2(i,j)-r5(i,j)*r1(i+1,j)
            End If
            r4(i,j)=1.0/r4(i,j)
            r3(i,j)=r3(i,j)*r4(i,j)
            r2(i,j)=r2(i,j)*r4(i,j)
            r1(i,j)=r1(i,j)*r4(i,j)

         End Do  !  i=2*ib+4,2*iz+1,-1
!
!        The six-by-six matrix near interface.
!
         i0=2*iz-3
         Do i=1,3
            t4(i,j)=1.0
            t5(i,j)=r5(i0+i,j)
            t6(i,j)=r6(i0+i,j)
            t7(i,j)=r7(i0+i,j)
         End Do
!
         Do i=4,6
            t1(i,j)=r1(i0+i,j)
            t2(i,j)=r2(i0+i,j)
            t3(i,j)=r3(i0+i,j)
            t4(i,j)=1.0
            t5(i,j)=0.0
            t6(i,j)=0.0
            t7(i,j)=0.0
         End Do
!
         t2(4,j)=t2(4,j)-t1(4,j)*t5(1,j)
         t3(4,j)=t3(4,j)-t1(4,j)*t6(1,j)
         t4(4,j)=t4(4,j)-t1(4,j)*t7(1,j)
!
         t3(4,j)=t3(4,j)-t2(4,j)*t5(2,j)
         t4(4,j)=t4(4,j)-t2(4,j)*t6(2,j)
         t5(4,j)=t5(4,j)-t2(4,j)*t7(2,j)
!
         t4(4,j)=t4(4,j)-t3(4,j)*t5(3,j)
         t5(4,j)=t5(4,j)-t3(4,j)*t6(3,j)
         t6(4,j)=t6(4,j)-t3(4,j)*t7(3,j)
!
         t5(4,j)=t5(4,j)/t4(4,j)
         t6(4,j)=t6(4,j)/t4(4,j)
         t4(4,j)=1.0/t4(4,j)
!
         t2(5,j)=t2(5,j)-t1(5,j)*t5(2,j)
         t3(5,j)=t3(5,j)-t1(5,j)*t6(2,j)
         t4(5,j)=t4(5,j)-t1(5,j)*t7(2,j)
!
         t3(5,j)=t3(5,j)-t2(5,j)*t5(3,j)
         t4(5,j)=t4(5,j)-t2(5,j)*t6(3,j)
         t5(5,j)=t5(5,j)-t2(5,j)*t7(3,j)
!
         t4(5,j)=t4(5,j)-t3(5,j)*t5(4,j)
         t5(5,j)=t5(5,j)-t3(5,j)*t6(4,j)
!
         t5(5,j)=t5(5,j)/t4(5,j)
         t4(5,j)=1.0/t4(5,j)
!
         t2(6,j)=t2(6,j)-t1(6,j)*t5(3,j)
         t3(6,j)=t3(6,j)-t1(6,j)*t6(3,j)
         t4(6,j)=t4(6,j)-t1(6,j)*t7(3,j)
!
         t3(6,j)=t3(6,j)-t2(6,j)*t5(4,j)
         t4(6,j)=t4(6,j)-t2(6,j)*t6(4,j)
!
         t4(6,j)=t4(6,j)-t3(6,j)*t5(5,j)
!
         t4(6,j)=1.0/t4(6,j)
!
      End Do ! j=1,np
!
      RETURN
      END

!========================================================================================================
!
!     The heptadiagonal solver.
!
      SUBROUTINE solve(mz,nz,mp,np,iz,g0)

	  USE RamMod, ONLY : u,v
	  USE WrkMod
      complex g0

      eps=1.0e-30
!
      DO j=1,np
!
!        The right side.
!
         i=1
         v(i+2)=s2(i,j)*u(i)+s3(i,j)*u(i+1)+s4(i,j)*u(i+2)+s5(i,j)*u(i+3)+s6(i,j)*u(i+4)+s7(i,j)*u(i+5)+eps
         Do i=2,2*nz
            v(i+2)=s1(i,j)*u(i-1)+s2(i,j)*u(i)+s3(i,j)*u(i+1)+               &
                   s4(i,j)*u(i+2)+s5(i,j)*u(i+3)+s6(i,j)*u(i+4)+             &
            s7(i,j)*u(i+5)+eps
         End Do
!
         i1=2*iZ-1
!
!        The elimination steps.
!
         i=1
         v(i+2)=(v(i+2)-r2(i,j)*v(i)-r3(i,j)*v(i+1))*r4(i,j)+eps
         Do i=2,2*iZ
            v(i+2)=(v(i+2)-r1(i,j)*v(i-1)-r2(i,j)*v(i) - r3(i,j)*v(i+1))*r4(i,j)+eps
         End Do
         i=2*nz
         v(i+2)=(v(i+2)-r6(i,j)*v(i+4)-r5(i,j)*v(i+3))*r4(i,j)+eps
         Do i=2*nz-1,2*iZ+1,-1
            v(i+2)=(v(i+2)-r7(i,j)*v(i+5)-r6(i,j)*v(i+4) - r5(i,j)*v(i+3))*r4(i,j)+eps
         End Do
!
!        The six-by-six system near the ocean bottom is solved.
!
         v(4+i1)=(v(4+i1)-t1(4,j)*v(1+i1)-t2(4,j)*v(2+i1) - t3(4,j)*v(3+i1))*t4(4,j)+eps
         v(5+i1)=(v(5+i1)-t1(5,j)*v(2+i1)-t2(5,j)*v(3+i1) - t3(5,j)*v(4+i1))*t4(5,j)+eps
         u(6+i1)=(v(6+i1)-t1(6,j)*v(3+i1)-t2(6,j)*v(4+i1) - t3(6,j)*v(5+i1))*t4(6,j)+eps
!
         u(5+i1)=v(5+i1)-t5(5,j)*u(6+i1)+eps
         u(4+i1)=v(4+i1)-t5(4,j)*u(5+i1)-t6(4,j)*u(6+i1)+eps
         u(3+i1)=v(3+i1)-t5(3,j)*u(4+i1)-t6(3,j)*u(5+i1) - t7(3,j)*u(6+i1)+eps
         u(2+i1)=v(2+i1)-t5(2,j)*u(3+i1)-t6(2,j)*u(4+i1) - t7(2,j)*u(5+i1)+eps
         u(1+i1)=v(1+i1)-t5(1,j)*u(2+i1)-t6(1,j)*u(3+i1) - t7(2,j)*u(4+i1)+eps
!
!        The back substitution steps.
!
         Do i=2*iZ-3,1,-1
            u(i+2)=v(i+2)-r5(i,j)*u(i+3)-r6(i,j)*u(i+4) - r7(i,j)*u(i+5)+eps
         End Do
         Do i=2*iZ+4,2*nz
            u(i+2)=v(i+2)-r3(i,j)*u(i+1)-r2(i,j)*u(i) - r1(i,j)*u(i-1)+eps
         End Do

      End DO
!
      DO i=1,nz
         u(2*i+1)=g0*u(2*i+1)
         u(2*i+2)=g0*u(2*i+2)
      END DO
!
      RETURN
      END

!========================================================================================================
!
!     The rotated Pade coefficients [J. Acoust. Soc. Am. 101, 760-766
!     (1997)].
!
      SUBROUTINE rpade(mp,np,k0,dr,theta,g0)
      USE MathConst
	  USE WrkMod
!
      complex g0,tfact,rot0,rot1,rot2
      real k0
      tfact=cexp(-ci*theta*pi/360.0)
      den=float(2*np+1)
      rot0=1.0
!
      DO j=1,np
!
!        The Pade coefficients.
!
         pade1=(2.0/den)*sin(float(j)*pi/den)**2
         pade2=cos(float(j)*pi/den)**2
!
!        The rotated Pade coefficients.
!
         rot1=tfact*pade1/(1.0+pade2*(tfact**2-1.0))**2
         rot2=tfact**2*pade2/(1.0+pade2*(tfact**2-1.0))
         rot0=rot0+pade1*(tfact**2-1.0)/(1.0+pade2*(tfact**2-1.0))
!
!        The Crank-Nicolson coefficients.
!
         pd1(j)=rot2+0.5*ci*k0*dr*rot1
         pd2(j)=rot2-0.5*ci*k0*dr*rot1
!
      END DO

      rot0=rot0/tfact
      g0=cexp(ci*k0*dr*rot0)
!
      RETURN
      END

!========================================================================================================
!
!     The coefficients of the rational approximation.
!
      SUBROUTINE epade(mp,np,ns,ip,k0,c0,dr,irot,theta,g0,nu8)
	  USE MathConst
	  USE WrkMod
!
      implicit real*8 (a-h,o-z)
      complex*16 z1,z2,g,dg,dh1,dh2,dh3,a,b,nu
      complex*8 g0,nu8
      real*4 k0,c0,dr,theta
      parameter (m=40)
      dimension bin(m,m),a(m,m),b(m),dg(m),dh1(m),dh2(m),dh3(m),fact(m)

      sig=k0*dr
      n=2*np
      g0=cexp(ci*k0*dr)
!
      IF((ip.eq.1).and.(irot.eq.1))THEN
         call rpade(mp,np,k0,dr,theta,g0)
         RETURN
      END IF
!
      IF(ip.eq.1)THEN
         nu=0.0d0
         alp=0.0d0
      ELSE
         nu=nu8
         alp=-0.25d0
      END IF
!
!     The factorials.
!
      fact(1)=1.0d0
      DO i=2,n
         fact(i)=dfloat(i)*fact(i-1)
      END DO
!
!     The binomial coefficients.
!
      DO i=1,n+1
         bin(i,1)=1.0d0
         bin(i,i)=1.0d0
      END DO

      DO i=3,n+1
         DO j=2,i-1
            bin(i,j)=bin(i-1,j-1)+bin(i-1,j)
         END DO
      END DO
!
      DO i=1,n
         DO j=1,n
            a(i,j)=0.0d0
         END DO
      END DO
!
!     The accuracy constraints.
!
      call deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
!
      DO i=1,n
         b(i)=dg(i+1)
      END DO

      DO i=1,n
         IF(2*i-1.le.n)a(i,2*i-1)=fact(i)
         DO j=1,i
            IF(2*j.le.n)a(i,2*j)=-bin(i+1,j+1)*fact(j)*dg(i-j+1)
         END DO
      END DO
!
!     The stability constraints.
!
      IF(ns.ge.1)THEN
         z1=-3.0d0
         b(n)=-1.0d0
         DO j=1,np
            a(n,2*j-1)=z1**j
            a(n,2*j)=0.0d0
         END DO
      END IF
!
      IF(ns.ge.2)THEN
         z1=-1.5d0
         b(n-1)=-1.0d0
         DO j=1,np
            a(n-1,2*j-1)=z1**j
            a(n-1,2*j)=0.0d0
         END DO
      END IF
!
      call gauss(m,n,a,b)
!
      dh1(1)=1.0d0
      DO j=1,np
         dh1(j+1)=b(2*j-1)
      END DO

	  call fndrt(dh1,np,dh2,m)

	  DO j=1,np
         pd1(j)=-1.0d0/dh2(j)
      END DO
!
      dh1(1)=1.0d0
      DO j=1,np
         dh1(j+1)=b(2*j)
      END DO

	  call fndrt(dh1,np,dh2,m)

	  DO j=1,np
         pd2(j)=-1.0d0/dh2(j)
      END DO
!
      RETURN
      END

!========================================================================================================
!
!     The operator function.
!
!     FUNCTION g(sig,x,alp,nu)
!     USE MathConst
!
!     complex*16 g,nu
!     real*8 alp,sig,x
!     g=(1.0d0-nu*x)**2*cdexp(alp*dlog(1.0d0+x) + ciD*sig*(-1.0d0+dsqrt(1.0d0+x)))
!     RETURN
!     END
!
!========================================================================================================
!
!     The derivatives of the operator function at x=0.
!
      SUBROUTINE deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
      USE MathConst

      implicit real*8 (a-h,o-z)
      complex*16 dg(m),dh1(m),dh2(m),dh3(m),nu
      real*8 bin(m,m)
!
      dh1(1)=0.5d0*ciD*sig
      exp1=-0.5d0
      dh2(1)=alp
      exp2=-1.0d0
      dh3(1)=-2.0d0*nu
      exp3=-1.0d0
      DO i=2,n
         dh1(i)=dh1(i-1)*exp1
         exp1=exp1-1.0d0
         dh2(i)=dh2(i-1)*exp2
         exp2=exp2-1.0d0
         dh3(i)=-nu*dh3(i-1)*exp3
         exp3=exp3-1.0d0
      END DO
!
      dg(1)=1.0d0
      dg(2)=dh1(1)+dh2(1)+dh3(1)
      DO i=2,n
         dg(i+1)=dh1(i)+dh2(i)+dh3(i)
         DO j=1,i-1
            dg(i+1)=dg(i+1)+bin(i,j)*(dh1(j)+dh2(j)+dh3(j))*dg(i-j+1)
         END DO
      END DO
!
      RETURN
      END

!========================================================================================================
!
!     Gaussian elimination.
!
      SUBROUTINE gauss(m,n,a,b)

      implicit real*8 (a-h,o-z)
      complex*16 a(m,m),b(m)
!
!     Downward elimination.
!
      DO i=1,n
         IF(i.lt.n) Call pivot(m,n,i,a,b)
         a(i,i)=1.0d0/a(i,i)
         b(i)=b(i)*a(i,i)
         IF(i.lt.n)THEN
            DO j=i+1,n
               a(i,j)=a(i,j)*a(i,i)
            END DO
            DO k=i+1,n
               b(k)=b(k)-a(k,i)*b(i)
               DO j=i+1,n
                  a(k,j)=a(k,j)-a(k,i)*a(i,j)
               END DO
            END DO
         END IF
      END DO
!
!     Back substitution.
!
      DO i=n-1,1,-1
         DO j=i+1,n
            b(i)=b(i)-a(i,j)*b(j)
         END DO
      END DO
!
      RETURN
      END


!========================================================================================================
!
!     Rows are interchanged for stability.
!
      SUBROUTINE pivot(m,n,i,a,b)
      implicit real*8 (a-h,o-z)
      complex*16 temp,a(m,m),b(m)
!
      i0=i
      amp0=cdabs(a(i,i))
      DO j=i+1,n
         amp=cdabs(a(j,i))
         IF(amp.gt.amp0)THEN
            i0=j
            amp0=amp
         END IF
      END DO
      IF(i0.eq.i)RETURN
!
      temp=b(i)
      b(i)=b(i0)
      b(i0)=temp
      DO j=i,n
         temp=a(i,j)
         a(i,j)=a(i0,j)
         a(i0,j)=temp
      END DO
!
      RETURN
      END

!========================================================================================================
!
!     The root-finding subroutine.
!
      SUBROUTINE fndrt(a,n,z,m)
      USE MathConst
      complex*16 a(m),z(m),root
      real*8 err
!
      IF(n.eq.1)THEN
         z(1)=-a(1)/a(2)
         RETURN
      END IF
      IF(n.eq.2)Go To 1
!
      DO k=n,3,-1
!
!        Obtain an approximate root.
!
         root=0.0d0
         err=1.0d-12
         call guerre(a,k,m,root,err,1000)
!
!        Refine the root by iterating five more times.
!
         err=0.0d0
         call guerre(a,k,m,root,err,5)
         z(k)=root
!
!        Divide out the factor (z-root).
!
         DO i=k,1,-1
            a(i)=a(i)+root*a(i+1)
         END DO
         DO i=1,k
            a(i)=a(i+1)
         END DO
!
      END DO
!
!     Solve the quadratic equation.
!
    1 Continue
          z(2)=0.5*(-a(2)+sqrt(a(2)**2-4.0*a(1)*a(3)))/a(3)
      z(1)=0.5*(-a(2)-sqrt(a(2)**2-4.0*a(1)*a(3)))/a(3)
!
      RETURN
      END


!========================================================================================================
!
!     This subroutine finds a root of a polynomial of degree n > 2
!     by Laguerre's method.
!
      SUBROUTINE guerre(a,n,m,z,err,nter)
      USE MathConst

      complex*16 a(m),az(50),azz(50),z,dz,p,pz,pzz,f,g,h
      real*8 amp1,amp2,rn,eps,err
      eps=1.0d-20
      rn=real(n)
!
!     The coefficients of p'(z) and p''(z).
!
      DO i=1,n
         az(i)=float(i)*a(i+1)
      END DO
      DO i=1,n-1
         azz(i)=float(i)*az(i+1)
      END DO
!
      iter=0
      jter=0

      DO While(.true.)

         p=a(n)+a(n+1)*z

         DO i=n-1,1,-1
            p=a(i)+z*p
         END DO

         IF(abs(p).lt.eps)RETURN
!
         pz=az(n-1)+az(n)*z
         DO i=n-2,1,-1
            pz=az(i)+z*pz
         END DO

!
         pzz=azz(n-2)+azz(n-1)*z
         DO i=n-3,1,-1
            pzz=azz(i)+z*pzz
         END DO

!
!        The Laguerre perturbation.
!
         f=pz/p
         g=f**2-pzz/p
         h=sqrt((rn-1.0d0)*(rn*g-f**2))
         amp1=abs(f+h)
         amp2=abs(f-h)
         IF(amp1.gt.amp2)THEN
            dz=-rn/(f+h)
         else
            dz=-rn/(f-h)
         END IF
!
         iter=iter+1
!
!        Rotate by 90 degrees to avoid limit cycles.
!
         jter=jter+1
         IF(jter.eq.10)THEN
            jter=1
            dz=dz*ciD
         END IF
         z=z+dz
!
         IF(iter.eq.100)THEN
            write(*,*)' '
            write(*,*)'   Laguerre method not converging.'
            write(*,*)'   Try a different combination of DR and NP.'
            write(*,*)' '
            stop
         END IF

        IF((abs(dz).le.err).or.(iter.ge.nter)) EXIT

!
      END DO
!
      RETURN
      END
