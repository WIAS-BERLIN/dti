      subroutine ftensnl(par,nb,b,fv)
C
C  compute f(par) 
C
      implicit none
      integer nb
      double precision par(7),b(6,nb),fv(nb)
      double precision D(6)
      call rho2D(par(2),D)
      call sihat(par(1),D,b,fv,nb)
C
C   this gives vector  th0*exp(-b g_i^T D g_i) in fv
C
      RETURN
      END

      subroutine gtensnl(par,nb,b,fv,grad)
C
C  compute gradient of f(par)
C
      implicit none
      integer nb
      double precision par(7),b(6,nb),fv(nb),grad(nb,7)
      integer i
      double precision D(6),z,z1
      call rho2D(par(2),D)
      call sihat(par(1),D,b,fv,nb)
C
C   this gives vector  th0*exp(-b g_i^T D g_i) in fv
C
C
C     derivative with respect to theta0
C
      DO i=1,nb
         grad(i,1)=fv(i)/par(1)
      END DO
C
C     derivatives with respect to par
C
      DO i=1,nb
         z=fv(i)
         z1=2.d0*b(1,i)*par(2)+b(2,i)*par(3)+b(3,i)*par(4)
         grad(i,2)=-z*z1
         z1=2.d0*b(4,i)*par(3)+b(2,i)*par(2)+b(5,i)*par(4)
         grad(i,3)=-z*z1
         z1=2.d0*b(6,i)*par(4)+b(3,i)*par(2)+b(5,i)*par(3)
         grad(i,4)=-z*z1
         z1=2.d0*b(4,i)*par(5)+b(5,i)*par(6)
         grad(i,5)=-z*z1
         z1=2.d0*b(6,i)*par(6)+b(5,i)*par(5)
         grad(i,6)=-z*z1
         z1=2.d0*b(6,i)*par(7)
         grad(i,7)=-z*z1
      END DO
C
C     We now have the gradient in grad
C
      RETURN
      END
