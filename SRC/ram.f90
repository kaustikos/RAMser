!==============================================================================
MODULE RAMMod

   SAVE
!
!  Mr=bathymetry points, Mz=depth grid, mp=pade terms, Mzr=max receivers horizons for txt output
!
   INTEGER,   PARAMETER :: Mr=2000, Mz=1000000, Mp=10, Mzr=1000

   INTEGER    flHandle, nOut, Mdr, Nz, Nzr, Ndr
   INTEGER    ir(Mzr)

   REAL (Kind=4) ::    k0, c0, omega
   REAL (Kind=4) ::    alpw(Mz), alpb(Mz), cw(Mz), cb(Mz), rhob(Mz), attn(Mz)
   REAL (Kind=4) ::    rb(Mr),zb(Mr)
   REAL (Kind=4) ::    ksqw(Mz), f1(Mz), f2(Mz), f3(Mz), tlg(Mz), realP(Mz), imagP(Mz)

   REAL (Kind=4) ::    Dir( Mzr ), r, rDist, rp, rSrc
!                       cs(Mz),  attns(Mz)              ! Dummy arrays for consistancy with RAMS

   COMPLEX    u(Mz), v(Mz), ksqb(Mz), ksq(Mz), pd1(Mp), pd2(Mp)
   COMPLEX    r1(Mz,Mp), r2(Mz,Mp), r3(Mz,Mp), s1(Mz,Mp), s2(Mz,Mp), s3(Mz,Mp)

!   Real (Kind=4) ::    rs, zs, dr, dz, zmplt, eta, eps, omega,rmax, &
!                       rRange, rdr, theta

END MODULE RAMMod



!========================================================================================================

MODULE MathConst

   SAVE

   REAL    (KIND=4), PARAMETER :: PI   = 3.1415926535898
   COMPLEX (KIND=4), PARAMETER :: CI  = ( 0.0, 1.0 )
   COMPLEX (KIND=8), PARAMETER :: CID = ( 0.0D0, 1.0D0 )

END MODULE MathConst

!========================================================================================================

      PROGRAM RAM
  	  USE RAMMod

      INTEGER  closeFlag
!
      OPEN(unit=1,status='old',file='rams.3.in')
      OPEN(unit=2,status='unknown',file='tl.nLine.Txt', recl=1024)
      OPEN(unit=3,status='unknown',file='TLrz',form='unformatted',access='stream')
	OPEN(unit=9,status='unknown',file='RePrz',form='unformatted',access='stream')
	OPEN(unit=10,status='unknown',file='ImPrz',form='unformatted',access='stream')



      nOut = 0
      r    = 0
!
      CALL setup(np,ns,ndz,iz,nzplt,lz,ib,dr,dz,eta,eps,rmax,rs)
!
!     March the acoustic field out in range.
!
      DO WHILE( r.lt.rmax )
         CALL updat(np,iz,ib,dr,dz,eta,rmax, rs )
         CALL solve(np,iz)
         r=r+dr
         rDist=rDist+dr
         CALL outpt(ndz,iz,dz,nzplt,lz,eps)
!        nOut=nOut+1
      END DO

      CLOSE(1)
      CLOSE(2)
      close(3)
      close(9)
      close(10)

      OPEN(unit=4,status='unknown',file='DomainBounds.Info')
      WRITE( 4,101 ) nOut, rmax, ndr*dr, lz, dz*ndz, dz*ndz
      CLOSE(4)

      STOP
101   FORMAT(1X,'Num Range Records    : ', I8,/,    &
             1X,'Max Range            : ', F9.0,/,  &
             1X,'Step Along Range     : ', F12.3,/, &
             1X,'Num Horizontal Lines : ', I8,/,    &
             1X,'Max Depth            : ', F12.3,/,  &
             1X,'Step Along Depth     : ', F12.3    &
            )
      END
!
!========================================================================================================
!
!     Initialize the parameters, acoustic field, and matrices.
!
      SUBROUTINE setup(np,ns,ndz,iz,nzplt,lz,ib, dr,dz,eta,eps,rmax,rs)

  	  USE RAMMod
  	  USE MathConst

	  REAL zr( Mzr ) 
!
      zr(1:Mzr) = -1      
      
      READ(1,*)
      READ(1,*)freq, rSrc, zs, Direct   ! Direct - Dummy for consistency with RAMS
      READ(1,*) ( zr(i),i=1,Mzr )
      
      i = 0
      DO WHILE (.true.)
         i=i+1
         IF(zr(i).lt.0.0) EXIT
      END DO
      Nzr = i-1 

      READ(1,*)rmax,dr,ndr
      READ(1,*)zmax,dz,ndz,zmplt
      READ(1,*)c0, np,ns,rs
!
      i=0
      DO WHILE (.TRUE.)
         i=i+1
         READ(1,*)rb(i),zb(i)
         IF (rb(i).lt.0.0) EXIT
      END DO
      rb(i)=rb(i-1) + rSrc + 2.0*rmax
      zb(i)=zb(i-1)
      nBott = i
!
      eta=1.0/(40.0*PI*alog10(exp(1.0)))
      eps=1.0e-20
      ib=1
      mdr=0
      r=dr
      omega=2.0*PI*freq
      DO i=1, Nzr
         ri=1.0+zr(i)/dz
         ir(i)=ifix(ri)
         dir(i)=ri-float( ir(i) )
      END DO
!      
      k0=omega/c0
      nz=zmax/dz-0.5
      nzplt=zmplt/dz-0.5
      z=zb(1)
      iz=1.0+z/dz
      iz=max(2,iz)
      iz=min(nz,iz)
      IF( rs .lt. dr ) rs=2.0*rmax
!
      IF( nz+2 .gt. mz )THEN
         write(*,*)'   Need to increase parameter mz to ',nz+2
         stop
      END IF
      IF( np .gt. mp )THEN
         write(*,*)'   Need to increase parameter mp to ',np
         STOP
      END IF
      IF(i.gt.mr)THEN
         write(*,*)'   Need to increase parameter mr to ',i
         STOP
      END IF
!
      DO j=1, mp
         r3(1,j)=0.0
         r1(nz+2,j)=0.0
      END DO
      DO i=1, nz+2
         u(i)=0.0
         v(i)=0.0
      END DO
      lz=0
      DO i=ndz, nzplt, ndz
      lz=lz+1
      END DO
!     write(3)lz
!
!     The initial profiles and starting field.
!
      rp = 0
      DO WHILE( rp .le. rSrc )
         rDist = rp
         CALL profl(dz,eta,rmax) 
      END DO
      rDist = rSrc + dr

      DO WHILE( ib+1 .lt. nBott .AND. rDist.ge.rb(ib+1) )
        ib=ib+1
      END DO
      z=zb(ib)+(rDist+0.5*dr-rb(ib))*(zb(ib+1)-zb(ib))/(rb(ib+1)-rb(ib))
      iz=1.0+z/dz
      iz=max(2,iz)
      iz=min(nz,iz)
      
      CALL selfs(np,ns,iz,zs,dr,dz)
      CALL outpt(ndz,iz,nzplt,lz,eps)
!
!     The propagation matrices.
!
      CALL epade(np,ns,1,dr)
      CALL matrc(np,iz,iz,dz)
!
      RETURN
      END
!
!========================================================================================================
!
!     Set up the profiles.
!
      SUBROUTINE profl( dz,eta,rmax)
   	  USE RAMMod
   	  USE MathConst
!
      CALL zread(mz,nz,dz,cw)
      CALL zread(mz,nz,dz,cb)
      CALL dummyrd                  ! Dummy read cs for consistency with RAMS
      CALL zread(mz,nz,dz,rhob)
      CALL zread(mz,nz,dz,attn)
      CALL dummyrd                  ! Dummy read attns for consistency with RAMS
      rp=rSrc+2.0*rmax
      READ(1,*,END=1)rp
!
    1 CONTINUE
      DO i=1, nz+2
         ksqw(i)=(omega/cw(i))**2-k0**2
         ksqb(i)=((omega/cb(i))*(1.0+CI*eta*attn(i)))**2-k0**2
         alpw(i)=sqrt(cw(i)/c0)
         alpb(i)=sqrt(rhob(i)*cb(i)/c0)
      END DO
      
!
      RETURN
      END
!
!========================================================================================================
!
!     Profile reader and interpolator.
!
      SUBROUTINE zread(mz,nz,dz,prof)
      REAL prof(mz)
!
      DO i=1,nz+2
         prof(i)=-1.0
      END DO
      READ(1,*)zi,profi
	  prof(1)=profi
      i=1.5+zi/dz
      prof(i)=profi
      iold=i

      DO WHILE( .TRUE. )
         READ(1,*) zi,profi
         IF ( zi.lt.0.0 ) EXIT
         i=1.5+zi/dz
         IF(i.eq.iold)i=i+1
         prof(i)=profi
         iold=i
      END DO

      prof(nz+2)=prof(i)

      i=1
      j=1
      DO WHILE ( j.lt.nz+2 )
         i=i+1
         IF (prof(i).lt.0.0) CYCLE
         IF(i-j.ne.1) THEN
            DO k=j+1,i-1
               prof(k)=prof(j)+float(k-j)*(prof(i)-prof(j))/float(i-j)
            END DO
         END IF
         j=i
      END DO

      RETURN
      END
!
!========================================================================================================
!
!     Dummy read profile
!
      SUBROUTINE dummyrd
!
      DO WHILE( .TRUE. )
         READ(1,*) zi,profi
         IF ( zi.lt.0.0 ) EXIT
      END DO

      RETURN
      END
!
!========================================================================================================
!
!     The tridiagonal matrices.
!
      SUBROUTINE matrc(np,iz,jz,dz)
         
      USE RAMMod  
         
      COMPLEX d1,d2,d3,rfact
!
      a1=k0**2/6.0
      a2=2.0*k0**2/3.0
      a3=k0**2/6.0
      cfact=0.5/dz**2
      dfact=1.0/12.0
!
!     New matrices when iz.eq.jz.
!
      IF(iz.eq.jz)THEN
         i1=2
         i2=nz+1
         DO i=1,iz
            f1(i)=1.0/alpw(i)
            f2(i)=1.0
            f3(i)=alpw(i)
            ksq(i)=ksqw(i)
         END DO
         DO i=iz+1,nz+2
            f1(i)=rhob(i)/alpb(i)
            f2(i)=1.0/rhob(i)
            f3(i)=alpb(i)
            ksq(i)=ksqb(i)
         END DO
      END IF
!
!     Updated matrices when iz.ne.jz.
!
      IF(iz.gt.jz)THEN
         i1=jz
         i2=iz+1
         DO i=jz+1,iz
            f1(i)=1.0/alpw(i)
            f2(i)=1.0
            f3(i)=alpw(i)
            ksq(i)=ksqw(i)
         END DO
      END IF
!
      IF(iz.lt.jz)THEN
         i1=iz
         i2=jz+1
         DO i=iz+1,jz
            f1(i)=rhob(i)/alpb(i)
            f2(i)=1.0/rhob(i)
            f3(i)=alpb(i)
            ksq(i)=ksqb(i)
         END DO
      END IF
!
      DO i=i1,i2
!
!     Discretization by Galerkin's method.
!
         c1=cfact*f1(i)*(f2(i-1)+f2(i))*f3(i-1)
         c2=-cfact*f1(i)*(f2(i-1)+2.0*f2(i)+f2(i+1))*f3(i)
         c3=cfact*f1(i)*(f2(i)+f2(i+1))*f3(i+1)
         d1=c1+dfact*(ksq(i-1)+ksq(i))
         d2=c2+dfact*(ksq(i-1)+6.0*ksq(i)+ksq(i+1))
         d3=c3+dfact*(ksq(i)+ksq(i+1))
!
         DO j=1,np
            r1(i,j)=a1+pd2(j)*d1
            r2(i,j)=a2+pd2(j)*d2
            r3(i,j)=a3+pd2(j)*d3
            s1(i,j)=a1+pd1(j)*d1
            s2(i,j)=a2+pd1(j)*d2
            s3(i,j)=a3+pd1(j)*d3
         END DO

      END DO
!
!     The matrix decomposition.
!
      DO j=1,np
         
         DO i=i1,iz
            rfact=1.0/(r2(i,j)-r1(i,j)*r3(i-1,j))
            r1(i,j)=r1(i,j)*rfact
            r3(i,j)=r3(i,j)*rfact
            s1(i,j)=s1(i,j)*rfact
            s2(i,j)=s2(i,j)*rfact
            s3(i,j)=s3(i,j)*rfact
         END DO
!
         DO i=i2,iz+2,-1
            rfact=1.0/(r2(i,j)-r3(i,j)*r1(i+1,j))
            r1(i,j)=r1(i,j)*rfact
            r3(i,j)=r3(i,j)*rfact
            s1(i,j)=s1(i,j)*rfact
            s2(i,j)=s2(i,j)*rfact
            s3(i,j)=s3(i,j)*rfact
         END DO
!
         r2(iz+1,j)=r2(iz+1,j)-r1(iz+1,j)*r3(iz,j)
         r2(iz+1,j)=r2(iz+1,j)-r3(iz+1,j)*r1(iz+2,j)
         r2(iz+1,j)=1.0/r2(iz+1,j)
!
      END DO
!
      RETURN
      END
!
!========================================================================================================
!
!     The tridiagonal solver.
!
      SUBROUTINE solve(np,iz)
      USE RAMMod     
      
      eps=1.0e-30
!
      DO j=1,np
!
!     The right side.
!
         DO i=2,nz+1
            v(i)=s1(i,j)*u(i-1)+s2(i,j)*u(i)+s3(i,j)*u(i+1)+eps
         END DO
!
!     The elimination steps.
!
         DO i=3,iz
            v(i)=v(i)-r1(i,j)*v(i-1)+eps
         END DO
         DO i=nz,iz+2,-1
            v(i)=v(i)-r3(i,j)*v(i+1)+eps
         END DO
!
         u(iz+1)=(v(iz+1)-r1(iz+1,j)*v(iz)-r3(iz+1,j)*v(iz+2))*         &
                  r2(iz+1,j)+eps
!
!     The back substitution steps.
!
         DO i=iz,2,-1
            u(i)=v(i)-r3(i,j)*u(i+1)+eps
         END DO
         DO i=iz+2,nz+1
           u(i)=v(i)-r1(i,j)*u(i-1)+eps
         END DO
      
      END DO
!
      RETURN
      END
!
!=============================================================================
!
!     Matrix updates.
!
      SUBROUTINE updat(np,iz,ib,dr,dz,eta,rmax, rs)
      USE RAMMod   
      USE MathConst
!
!     Varying bathymetry.
!
      IF(rDist.ge.rb(ib+1))ib=ib+1
      jz=iz
      z=zb(ib)+(rDist+0.5*dr-rb(ib))*(zb(ib+1)-zb(ib))/(rb(ib+1)-rb(ib))
      iz=1.0+z/dz
      iz=max(2,iz)
      iz=min(nz,iz)
      IF(iz.ne.jz)CALL matrc(np,iz,jz,dz)
!
!     Varying profiles.
!
      IF(r.ge.rp)THEN
         CALL profl(dz,eta,rmax)
         CALL matrc(np,iz,iz,dz)
      END IF
!
!     Turn off the stability constraints.
!
      IF(r.ge.rs)THEN
         ns=0
         rs=2.0*rmax
         CALL epade(np,ns,1,dr)
         CALL matrc(np,iz,iz,dz)
      END IF
!
      RETURN
      END
!
!========================================================================================================
!
!     The self-starter.
!
      SUBROUTINE selfs(np,ns,iz,zs,dr,dz)
      USE RAMMod
      USE MathConst
!
!     Conditions for the delta function.
!
      si=1.0+zs/dz
      is=ifix(si)
      dis=si-float(is)
      u(is)=(1.0-dis)*sqrt(2.0*PI/k0)/(dz*alpw(is))
      u(is+1)=dis*sqrt(2.0*PI/k0)/(dz*alpw(is))
!
!     Divide the delta function by (1-X)**2 to get a smooth rhs.
!
      pd1(1)=0.0
      pd2(1)=-1.0
      CALL matrc(1,iz,iz,dz,alpw,alpb)
      CALL solve(1,iz)
      CALL solve(1,iz)
!
!     Apply the operator (1-X)**2*(1+X)**(-1/4)*exp(ci*k0*r*sqrt(1+X)).
!
      CALL epade(np,ns,2,dr)
      CALL matrc(np,iz,iz,dz)
      CALL solve(np,iz)
!
      RETURN
      END
!
!========================================================================================================
!
!     Output transmission loss.
!
      SUBROUTINE outpt(ndz,iz,dz,nzplt,lz,eps)

	  USE RAMMod
      
      COMPLEX ur
      REAL TL(Mzr)
      CHARACTER(1) bt
	  REAL*8  uNorm
	  INTEGER nb

      bt=CHAR(13)      
      mdr=mdr+1

      IF(mdr.eq.ndr)THEN

	     uNorm = 0
		 DO i=1,Nz
		    uNorm = uNorm + abs(u(i)*conjg(u(i)))
		 END DO
		 uNorm = DSQRT( uNorm )
		 
		 DO i=1,Nzr
            ur=(1.0-dir(i))*f3(ir(i))*u(ir(i))+dir(i)*f3(ir(i)+1)*u(ir(i)+1)             
            tl(i)=20.0*alog10(cabs(ur)+eps)-10.0*alog10(r+eps)
         END DO
         
	 write(2,100) r,rDist, ( tl(i),i=1,Nzr ),uNorm, uNorm, dz*iz
!         WRITE(2,'(2x,F10.1,F10.1, <Nzr>F8.2,4x,F16.4)') r, rDist, ( tl(i),i=1,Nzr ),uNorm
         WRITE(*,'(2x,F10.1,F10.1, F8.2,4x,F16.4,A1$)')  r, rDist,  tl(1), uNorm, bt
         mdr=0
!
         j=0
!        iflag=1
         DO i=1+ndz,nzplt,ndz
            ur=u(i)*f3(i)
            j=j+1
            tlg(j)=20.0*alog10(cabs(ur)+eps)-10.0*alog10(r+eps)
		 realP(j)=real(ur)/sqrt(r)
		 imagP(j)=aimag(ur)/sqrt(r)  
!
!        Mark the ocean bottom.
!
!           IF((i.gt.iz).and.(iflag.eq.1))THEN
  !         tlg(j)=0.0
!           iflag=0
!           END IF
!
         END DO

         nOut=nOut+1

         write(3)(tlg(j),j=1,lz)
	 write(9)(realP(j),j=1,lz)
	 write(10)(imagP(j),j=1,lz) 
      END IF
!
      RETURN

100   Format(2x,2F12.2,100F8.2,4x,F8.2,F2.4) 
      END
!
!========================================================================================================
!
!     The coefficients of the rational approximation.
!
      SUBROUTINE epade(np,ns,ip,dr)
      USE RAMMod
!
      IMPLICIT REAL*8 (a-h,o-z)
      COMPLEX*16 z1,z2,g,dg,dh1,dh2,dh3,a,b
!     COMPLEX*16 ci
      REAL*8 nu
      REAL*4 dr
      PARAMETER (m=40)
      DIMENSION bin(m,m),a(m,m),b(m),dg(m),dh1(m),dh2(m),dh3(m),fact(m)

!     PI=4.0d0*datan(1.0d0)    ! it isn't usable in this procedure
!     ci=dcmplx(0.0d0,1.0d0)
      sig=k0*dr
      n=2*np
!
      IF(ip.eq.1)THEN
         nu=0.0d0
         alp=0.0d0
      ELSE
         nu=1.0d0
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
      CALL deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
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
      CALL gauss(m,n,a,b)
!
      dh1(1)=1.0d0
      DO j=1,np
         dh1(j+1)=b(2*j-1)
      END DO
      CALL fndrt(dh1,np,dh2,m)
      DO j=1,np
         pd1(j)=-1.0d0/dh2(j)
      END DO
!
      dh1(1)=1.0d0
      DO j=1,np
         dh1(j+1)=b(2*j)
      END DO
      CALL fndrt(dh1,np,dh2,m)
      DO j=1,np
         pd2(j)=-1.0d0/dh2(j)
      END DO
!
      RETURN
      END
!
!========================================================================================================
!
!     The operator function.
!
!     FUNCTION g(ci,sig,x,alp,nu)
!     COMPLEX*16 ci,g
!     REAL*8 alp,sig,x,nu
!     g=(1.0d0-nu*x)**2*cdexp(alp*dlog(1.0d0+x)+      &
!        ci*sig*(-1.0d0+dsqrt(1.0d0+x)))
!     RETURN
!     END
!
!========================================================================================================
!
!     The derivatives of the operator function at x=0.
!
      SUBROUTINE deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
      USE MathConst
      IMPLICIT REAL*8 (a-h,o-z)
      COMPLEX*16 dg(m),dh1(m),dh2(m),dh3(m)
      REAL*8 bin(m,m),nu
!
      dh1(1)=0.5d0*CID*sig
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
!
!========================================================================================================
!
!     Gaussian elimination.
!
      SUBROUTINE gauss(m,n,a,b)
      IMPLICIT REAL*8 (a-h,o-z)
      COMPLEX*16 a(m,m),b(m)
!
!     Downward elimination.
!
      DO i=1,n

         IF(i.lt.n)CALL pivot(m,n,i,a,b)
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
!
!========================================================================================================
!
!     Rows are interchanged for stability.
!
      SUBROUTINE pivot(m,n,i,a,b)
      IMPLICIT REAL*8 (a-h,o-z)
      COMPLEX*16 temp,a(m,m),b(m)
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
!
!========================================================================================================
!
!     The root-finding subroutine. 
!
      SUBROUTINE fndrt(a,n,z,m)
      COMPLEX*16 a(m),z(m),root
      REAL*8 err
!
      IF(n.eq.1)THEN
         z(1)=-a(1)/a(2)
         RETURN
      END IF

      IF(n.ne.2)THEN
!
         DO k=n,3,-1
!
!     Obtain an approximate root.
!
            root=0.0d0
            err=1.0d-12
            CALL guerre(a,k,m,root,err,1000)
!
!     Refine the root by iterating five more times.
!
            err=0.0d0
            CALL guerre(a,k,m,root,err,5)
            z(k)=root
!
!     Divide out the factor (z-root).
!
            DO i=k,1,-1
               a(i)=a(i)+root*a(i+1)
            END DO
            DO i=1,k
               a(i)=a(i+1)
            END DO
!
         END DO
      END IF
      
!
!     Solve the quadratic equation.
!
      z(2)=0.5*(-a(2)+sqrt(a(2)**2-4.0*a(1)*a(3)))/a(3)
      z(1)=0.5*(-a(2)-sqrt(a(2)**2-4.0*a(1)*a(3)))/a(3)
!
      RETURN
      END
!
!========================================================================================================
!
!     This subroutine finds a root of a polynomial of degree n > 2
!     by Laguerre's method.
!
      SUBROUTINE guerre(a,n,m,z,err,nter)
      USE MathConst
      COMPLEX*16 a(m),az(50),azz(50),z,dz,p,pz,pzz,f,g,h
      REAL*8 amp1,amp2,rn,eps,err

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

      DO WHILE(.TRUE.)

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
!     The Laguerre perturbation.
!
         f=pz/p
         g=f**2-pzz/p
         h=sqrt((rn-1.0d0)*(rn*g-f**2))
         amp1=abs(f+h)
         amp2=abs(f-h)
         IF(amp1.gt.amp2)THEN
           dz=-rn/(f+h)
         ELSE
           dz=-rn/(f-h)
         END IF
!
         iter=iter+1
!
!     Rotate by 90 degrees to avoid limit cycles. 
!
         jter=jter+1
         IF(jter.eq.10)THEN
           jter=1
           dz=dz*CID
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
!
        IF((abs(dz).le.err).or.(iter.ge.nter)) EXIT

      END DO
!
      RETURN
      END
