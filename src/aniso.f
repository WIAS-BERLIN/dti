CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute all eigenvalues (lambda) of a 3x3 matrix and the corresponding EV (theta)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine eigen3(y,lambda,theta,ierr)
      implicit logical (a-z)
      integer ierr
      real*8 y(6),a(3,3),lambda(3),theta(3,3)
      integer i,j,l,ISUPPZ(6),lwork,iwork(50),liwork,n,m
      real*8 work(104),vl,vu,eps
      n=3
      m=3
      eps=1.d-50
      l=1
      DO i=1,3
         DO j=i,3
            a(i,j)=y(l)
            l=l+1
         END DO
      END DO
      lwork=104
      liwork=50
      call dsyevr('V','A','U',n,a,n,vl,vu,1,n,eps,m,lambda,
     1            theta,n,ISUPPZ,work,lwork,iwork,liwork,ierr)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Compute all eigenvalues (lambda) of a 3x3 matrix 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine eigen30(y,lambda,ierr)
      implicit logical (a-z)
      integer ierr
      real*8 y(6),a(3,3),lambda(3),theta(3,3)
      integer i,j,l,ISUPPZ(6),lwork,iwork(50),liwork,n,m
      real*8 work(104),vl,vu,eps
      n=3
      m=3
      eps=1.d-50
      l=1
      DO i=1,3
         DO j=i,3
            a(i,j)=y(l)
            l=l+1
         END DO
      END DO
      lwork=104
      liwork=50
      call dsyevr('N','A','U',n,a,n,vl,vu,1,n,eps,m,lambda,
     1            theta,n,ISUPPZ,work,lwork,iwork,liwork,ierr)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute anisotropic distance
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function adist(a,x,y,z,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer x,y,z
      real*8 a(6),xx,yy,zz,vext(3)
      xx=x*vext(1)
      yy=y*vext(2)
      zz=z*vext(3)      
      adist=a(1)*xx*xx+a(4)*yy*yy+a(6)*zz*zz+
     1               2.d0*(a(2)*xx*yy+a(3)*xx*zz+a(5)*yy*zz)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute statistical penalty for dti
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function dtidist2(th1,th2,Bcov)
      implicit logical (a-z)
      real*8 th1(6),th2(6),Bcov(6,6),z,zd,zd2
      integer i,j
      z=0.d0
      DO i=1,6
         zd=th1(i)-th2(i)
         z=z+zd*zd*Bcov(i,i)
         if(i.lt.6) THEN
            DO j=i+1,6
               zd2=th1(j)-th2(j)
               z=z+2.d0*zd*zd2*Bcov(i,j)
            END DO
         END IF
      END DO
      dtidist2=z
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for anisotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangex(a,h,ia,ie,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ia,ie
      real*8 a(6),h,vext(3)
      real*8 p1,p2,p3,p4,p5,p6,z,s,t,p55,p66,p44
      p1=a(1)
      p2=a(4)
      p3=a(6)
      p4=a(2)
      p5=a(3)
      p6=a(5)
      p44=p4*p4
      p55=p5*p5
      p66=p6*p6
      s=p2*p3-p6*p6
      t=p5*p4/p2/p3
      z=(p1-p44/p2-p55/p3+2.d0*p6*t)*s-p66*p55/p3-p66*p44/p2+
     1                                             2.d0*p66*p6*t
      if(abs(z).le.1e-40) THEN
C         call dblepr("denominator in rangex",21,z,1)
         z=1
      END IF
      z=s/z
      if(z.le.0) THEN
         z=0.d0
      ELSE 
         z=sqrt(z)
      END IF
      ia=-z*h/vext(1)
      ie=z*h/vext(1)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for anisotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangey(a,ix,h,ja,je,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ja,je,ix
      real*8 a(6),h,vext(3)
      real*8 p1,p2,p3,p4,p5,p6,z,s,t,p55,p66,p44,x
      x=ix/h*vext(1)
      p1=a(1)
      p2=a(4)
      p3=a(6)
      p4=a(2)
      p5=a(3)
      p6=a(5)
      p44=p4*p4
      p55=p5*p5
      p66=p6*p6
      s=p2*p3-p6*p6
      if(s.le.1e-10) THEN
C         call dblepr("denominator in rangey",21,s,1)
         z=1
      END IF
      t=-(p4*p3 - p6*p5)* x
      z=(p44*p3*p3-2.d0*p4*p3*p6*p5+p66*p55-p1*p3*s+p55*s)*x*x+p3*s
      if(z.le.0) THEN
         z=0.d0
      ELSE 
         z=sqrt(z)
      END IF
      ja=(t-z)*h/s/vext(2)
      je=(t+z)*h/s/vext(2)
      z=z*h/s/vext(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute range of x for anisotropic neighborhood
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rangez(a,ix,jy,h,ka,ke,vext)
C
C    a - diffusion tensor  ( a_11, a_12, a_13, a_22, a_23, a_33
C    ix - x-coordinate 
C    h - bandwidth
C    rh0 - regularization parameter    ( use  a+ rho/ni I  to define ellipsoid )
C    ni -  sum of weights (measures the variability of a
C    ia,ie -  rane of x values (restricted to the grid)
      implicit logical (a-z)
      integer ka,ke,ix,jy
      real*8 a(6),h,vext(3)
      real*8 p1,p2,p3,p4,p5,p6,z,t,x,y
      x=ix/h*vext(1)
      y=jy/h*vext(2)
      p1=a(1)
      p2=a(4)
      p3=a(6)
      p4=a(2)
      p5=a(3)
      p6=a(5)
      z=(p6*y + p5*x)
      t= - z/p3
      z=(-z*t - p2*y*y - p1*x*x - 2.d0*p4*x*y + 1.d0)/p3
      if(z.le.0) THEN
         z=0.d0
      ELSE 
         z=sqrt(z)
      END IF
      ka=(t-z)*h/vext(3)
      ke=(t+z)*h/vext(3)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute DTI-Indices for a slice
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti2Dfa(D,n1,n2,mask,fa,md,adir)
      implicit logical (a-z)
      integer n1,n2
      logical mask(n1,n2)
      real*8 D(6,n1,n2),fa(n1,n2),md(n1,n2),adir(3,n1,n2)
      integer i1,i2,ierr
      real*8 lambda(3),evec(3,3),trc,d1,d2,d3,a1,a2,a3,dd
      DO i1=1,n1
         DO i2=1,n2
            if(mask(i1,i2)) THEN
               call eigen3(D(1,i1,i2),lambda,evec,ierr)
               a1=lambda(1)
               a2=lambda(2)
               a3=lambda(3)
               trc=(a1+a2+a3)/3.d0
               adir(1,i1,i2)=evec(1,3)
               adir(2,i1,i2)=evec(2,3)
               adir(3,i1,i2)=evec(3,3)
               md(i1,i2)=trc
               d1=a1-trc
               d2=a2-trc
               d3=a3-trc
               dd=a1*a1+a2*a2+a3*a3
               IF(dd.gt.1.d-12) THEN
               fa(i1,i2)=sqrt(1.5d0*(d1*d1+d2*d2+d3*d3)/dd)
               ELSE
               fa(i1,i2)=0.d0
               ENDIF
            ELSE
               md(i1,i2)=0.d0
               fa(i1,i2)=0.d0
               adir(1,i1,i2)=1.d0
               adir(2,i1,i2)=0.d0
               adir(3,i1,i2)=0.d0
            END IF
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti2Dga(D,n1,n2,mask,ga,md,adir)
      implicit logical (a-z)
      integer n1,n2
      logical mask(n1,n2)
      real*8 D(6,n1,n2),ga(n1,n2),md(n1,n2),adir(3,n1,n2)
      integer i1,i2,ierr
      real*8 lambda(3),evec(3,3),trc,d1,d2,d3,a1,a2,a3,dd
      DO i1=1,n1
         DO i2=1,n2
            if(mask(i1,i2)) THEN
               call eigen3(D(1,i1,i2),lambda,evec,ierr)
               a1=log(lambda(1))
               a2=log(lambda(2))
               a3=log(lambda(3))
               trc=(a1+a2+a3)/3.d0
               adir(1,i1,i2)=evec(1,3)
               adir(2,i1,i2)=evec(2,3)
               adir(3,i1,i2)=evec(3,3)
               md(i1,i2)=trc
               d1=a1-trc
               d2=a2-trc
               d3=a3-trc
               dd=a1*a1+a2*a2+a3*a3
               IF(dd.gt.1.d-12) THEN
               ga(i1,i2)=sqrt(1.5d0*(d1*d1+d2*d2+d3*d3)/dd)
               ELSE
               ga(i1,i2)=0.d0
               ENDIF
            ELSE
               md(i1,i2)=0.d0
               ga(i1,i2)=0.d0
               adir(1,i1,i2)=1.d0
               adir(2,i1,i2)=0.d0
               adir(3,i1,i2)=0.d0
            END IF
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute DTI-Indices for a volume
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dtiind3D(D,n1,n2,n3,mask,fa,ga,md,adir,bary)
      implicit logical (a-z)
      integer n1,n2,n3
      logical mask(n1,n2,n3)
      real*8 D(6,n1,n2,n3),fa(n1,n2,n3),md(n1,n2,n3),adir(3,n1,n2,n3),
     1       bary(3,n1,n2,n3),ga(n1,n2,n3)
      integer i1,i2,i3,ierr
      real*8 lambda(3),evec(3,3),trc,d1,d2,d3,a1,a2,a3,dd
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
            if(mask(i1,i2,i3)) THEN
               call eigen3(D(1,i1,i2,i3),lambda,evec,ierr)
               a1=lambda(1)
               a2=lambda(2)
               a3=lambda(3)
               trc=(a1+a2+a3)/3.d0
               adir(1,i1,i2,i3)=evec(1,3)
               adir(2,i1,i2,i3)=evec(2,3)
               adir(3,i1,i2,i3)=evec(3,3)
               md(i1,i2,i3)=trc
               d1=a1-trc
               d2=a2-trc
               d3=a3-trc
               dd=a1*a1+a2*a2+a3*a3
               IF(dd.gt.1.d-12) THEN
               fa(i1,i2,i3)=sqrt(1.5d0*(d1*d1+d2*d2+d3*d3)/dd)
               bary(1,i1,i2,i3)=(a3-a2)/trc/3.d0
               bary(2,i1,i2,i3)=2.d0*(a2-a1)/trc/3.d0
               bary(3,i1,i2,i3)=a1/trc
               ELSE
               fa(i1,i2,i3)=0.d0
               bary(1,i1,i2,i3)=0.d0
               bary(2,i1,i2,i3)=0.d0
               bary(3,i1,i2,i3)=1.d0
               ENDIF
               d1=log(a1)
               d2=log(a2)
               d3=log(a3)
               dd=(d1+d2+d3)/3.d0
               d1=d1-dd
               d2=d2-dd
               d3=d3-dd
               ga(i1,i2,i3)=sqrt(d1*d1+d2*d2+d3*d3)
            ELSE
               md(i1,i2,i3)=0.d0
               fa(i1,i2,i3)=0.d0
               ga(i1,i2,i3)=0.d0
               adir(1,i1,i2,i3)=1.d0
               adir(2,i1,i2,i3)=0.d0
               adir(3,i1,i2,i3)=0.d0
               bary(1,i1,i2,i3)=0.d0
               bary(2,i1,i2,i3)=0.d0
               bary(3,i1,i2,i3)=1.d0
            END IF
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute DTI-Indices for a volume
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti3Dall(D,n1,n2,n3,mask,fa,ga,md,adir,ev)
      implicit logical (a-z)
      integer n1,n2,n3
      logical mask(n1,n2,n3)
      real*8 D(6,n1,n2,n3),fa(n1,n2,n3),md(n1,n2,n3),adir(3,n1,n2,n3),
     1       ev(3,n1,n2,n3),ga(n1,n2,n3)
      integer i1,i2,i3,ierr
      real*8 evec(3,3),trc,d1,d2,d3,a1,a2,a3,dd
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
            if(mask(i1,i2,i3)) THEN
               call eigen3(D(1,i1,i2,i3),ev(1,i1,i2,i3),evec,ierr)
               a1=ev(1,i1,i2,i3)
               a2=ev(2,i1,i2,i3)
               a3=ev(3,i1,i2,i3)
               trc=(a1+a2+a3)/3.d0
               adir(1,i1,i2,i3)=evec(1,3)
               adir(2,i1,i2,i3)=evec(2,3)
               adir(3,i1,i2,i3)=evec(3,3)
               md(i1,i2,i3)=trc
               d1=a1-trc
               d2=a2-trc
               d3=a3-trc
               dd=a1*a1+a2*a2+a3*a3
               IF(dd.gt.1.d-12) THEN
               fa(i1,i2,i3)=sqrt(1.5d0*(d1*d1+d2*d2+d3*d3)/dd)
               ELSE
               fa(i1,i2,i3)=0.d0
               ev(1,i1,i2,i3)=0.d0
               ev(2,i1,i2,i3)=0.d0
               ev(3,i1,i2,i3)=0.d0
               ENDIF
               d1=log(a1)
               d2=log(a2)
               d3=log(a3)
               dd=(d1+d2+d3)/3.d0
               d1=d1-dd
               d2=d2-dd
               d3=d3-dd
               ga(i1,i2,i3)=sqrt(d1*d1+d2*d2+d3*d3)
            ELSE
               md(i1,i2,i3)=0.d0
               fa(i1,i2,i3)=0.d0
               ga(i1,i2,i3)=0.d0
               adir(1,i1,i2,i3)=1.d0
               adir(2,i1,i2,i3)=0.d0
               adir(3,i1,i2,i3)=0.d0
               ev(1,i1,i2,i3)=0.d0
               ev(2,i1,i2,i3)=0.d0
               ev(3,i1,i2,i3)=0.d0
            END IF
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute DTI-eigenvalues for a volume
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti3Dev(D,n1,n2,n3,mask,ev)
      implicit logical (a-z)
      integer n1,n2,n3
      logical mask(n1,n2,n3)
      real*8 D(6,n1,n2,n3),ev(3,n1,n2,n3)
      integer i1,i2,i3,ierr
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(mask(i1,i2,i3)) THEN
                  call eigen30(D(1,i1,i2,i3),ev(1,i1,i2,i3),ierr)
               ELSE
                  ev(1,i1,i2,i3)=0.d0
                  ev(2,i1,i2,i3)=0.d0
                  ev(3,i1,i2,i3)=0.d0
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Compute DTI-eigenvectors for a volume
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine dti3Dand(D,n1,n2,n3,mask,andir)
      implicit logical (a-z)
      integer n1,n2,n3
      logical mask(n1,n2,n3)
      real*8 D(6,n1,n2,n3),andir(3,n1,n2,n3)
      integer i1,i2,i3,ierr
      real*8 lambda(3),evec(3,3)
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               if(mask(i1,i2,i3)) THEN
                  call eigen3(D(1,i1,i2,i3),lambda,evec,ierr)
                  andir(1,i1,i2,i3)=evec(1,3)
                  andir(2,i1,i2,i3)=evec(2,3)
                  andir(3,i1,i2,i3)=evec(3,3)
               ELSE
                  andir(1,i1,i2,i3)=0.d0
                  andir(2,i1,i2,i3)=0.d0
                  andir(3,i1,i2,i3)=0.d0
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   determine sum of location weights for a given geometry a(3) and given 
C   bandwidth
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Algorithmus zur Nullstellenbestimmung einer monotonen Funktion auf(0,\infty)
      subroutine gethani(x,y,value,a,vext,eps,bw)
      implicit logical(a-z)
      real*8 x,y,value,a(6),vext(3),eps,bw
      real*8 fw1,fw2,fw3,z
      real*8 sofw3D
      external sofw3D
      if(x.ge.y) RETURN
      fw1=sofw3D(a,x,vext)
      fw2=sofw3D(a,y,vext)
      DO WHILE(fw1.gt.value)
         x=x*x/y
         fw1=sofw3D(a,x,vext)
      END DO
      DO WHILE(fw2.le.value)
         y=y*y/x
         fw2=sofw3D(a,y,vext)
      END DO
      DO WHILE(min(fw2/value,value/fw1).gt.1.d0+eps)
         z=x+(value-fw1)/(fw2-fw1)*(y-x)
         fw3=sofw3D(a,z,vext)
         if(fw3.le.value) THEN
            x=z
            fw1=fw3
         ENDIF
         if(fw3.ge.value) THEN
            y=z
            fw2=fw3
         ENDIF
               call rchkusr()
      END DO
      if(fw2/value.gt.value/fw1) THEN
          bw=x+(value-fw1)/(fw2-fw1)*(y-x)
      ELSE
          bw=y-(fw2-value)/(fw2-fw1)*(y-x)
      ENDIF
      RETURN
      END  
      subroutine getvofh(a,bw,vext,vol)
      implicit logical(a-z)
      real*8 a(6),bw,vext(3),vol,sofw3D
      vol=sofw3D(a,bw,vext)
      RETURN
      END
C  compute sum of weights for anisotropic smoothing      
      real*8 function sofw3D(a,bw,vext)
      implicit logical(a-z)
      real*8 a(6),bw,vext(3)
      integer ia1,ie1,ia2,ie2,ia3,ie3,i1,i2,i3
      real*8 wij,sw,h2,adist
      external adist
      h2=bw*bw
      sw=1.d0
      call rangex(a,bw,ia1,ie1,vext)
      DO i1=1,ie1
         call rangey(a,i1,bw,ia2,ie2,vext)
         DO i2=ia2,ie2
            call rangez(a,i1,i2,bw,ia3,ie3,vext)
            DO i3=ia3,ie3
               wij=max(0.d0,1.d0-adist(a,i1,i2,i3,vext)/h2)
               sw=sw+2.d0*wij
            END DO
         END DO
      END DO
C  now case i1=0
      call rangey(a,0,bw,ia2,ie2,vext)
      DO i2=1,ie2
         call rangez(a,0,i2,bw,ia3,ie3,vext)
         DO i3=ia3,ie3
            wij=max(0.d0,1.d0-adist(a,0,i2,i3,vext)/h2)
            sw=sw+2.d0*wij
         END DO
      END DO
C  now case i1=i2=0
      call rangez(a,0,0,bw,ia3,ie3,vext)
         DO i3=1,ie3
            wij=max(0.d0,1.d0-adist(a,0,0,i3,vext)/h2)
            sw=sw+2.d0*wij
         END DO
      sofw3D=sw
      RETURN
      END
      subroutine sofw3Ds(a,bw,vext,sw)
      implicit logical(a-z)
      real*8 a(6),bw,vext(3)
      integer ia1,ie1,ia2,ie2,ia3,ie3,i1,i2,i3
      real*8 wij,sw,h2,adist
      external adist
      h2=bw*bw
      sw=1.d0
      call rangex(a,bw,ia1,ie1,vext)
      DO i1=1,ie1
         call rangey(a,i1,bw,ia2,ie2,vext)
         DO i2=ia2,ie2
            call rangez(a,i1,i2,bw,ia3,ie3,vext)
            DO i3=ia3,ie3
               wij=max(0.d0,1.d0-adist(a,i1,i2,i3,vext)/h2)
               sw=sw+2.d0*wij
            END DO
         END DO
      END DO
C  now case i1=0
      call rangey(a,0,bw,ia2,ie2,vext)
      DO i2=1,ie2
         call rangez(a,0,i2,bw,ia3,ie3,vext)
         DO i3=ia3,ie3
            wij=max(0.d0,1.d0-adist(a,0,i2,i3,vext)/h2)
            sw=sw+2.d0*wij
         END DO
      END DO
C  now case i1=i2=0
      call rangez(a,0,0,bw,ia3,ie3,vext)
         DO i3=1,ie3
            wij=max(0.d0,1.d0-adist(a,0,0,i3,vext)/h2)
            sw=sw+2.d0*wij
         END DO
      RETURN
      END
