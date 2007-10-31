      subroutine nlrdti(s,nb,n1,n2,n3,mask,b,th0,D,niter,eps,Varth,
     1                  res,rss)
      implicit logical (a-z)
      integer nb,n1,n2,n3,s(nb,n1,n2,n3),niter
      logical mask(n1,n2,n3)
      real*8 D(6,n1,n2,n3),b(6,nb),Varth(28,n1,n2,n3),res(nb,n1,n2,n3),
     1    th0(n1,n2,n3),F(nb),eps,rss(n1,n2,n3)
      integer i1,i2,i3,j
      real*8 theta(7)
      DO i3=1,n3
         call intpr("Nonlinear regression for slice No:",34,i3,1)
         DO i2=1,n2
            DO i1=1,n1
               if(mask(i1,i2,i3)) THEN
               call solvedtb(s(1,i1,i2,i3),nb,b,th0(i1,i2,i3),
     1                       D(1,i1,i2,i3),Varth(1,i1,i2,i3),
     2                       res(1,i1,i2,i3),niter,eps,
     3                       rss(n1,n2,n3))
               ELSE
                  DO j=1,6
                     D(j,i1,i2,i3)=0.d0
                  END DO
                  rss(i1,i2,i3)=0.d0
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine solvedti(s,nb,b,th0,D,Varth,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C
      implicit logical (a-z)
      integer nb,s(nb),niter
      real*8 D(6),b(6,nb),th0,Varth(28),F(nb),eps,
     1       theta(7)
      integer i,j,k,info,iter
      real*8 ntheta(7),Vtheta(7,7),z,thcorr,gamma,alpha,delta,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss,th1,th2,th3,th4,th5,th6,th7,
     3       nth1,nth2,nth3,nth4,nth5,nth6,nth7,res,X(7)
      alpha=0.5D0
      delta=0.25D0
      th1=th0
      th2=D(1)
      th3=D(2)
      th4=D(3)
      th5=D(4)
      th6=D(5)
      th7=D(6)
C      call dblepr("theta",5,theta,7)
      gamma=1.d0
      alpha=0.7d0
      oldrss=1.d50
      rss=0.d0
      DO i=1,nb
         z=b(1,i)*th2+b(2,i)*th3+b(3,i)*th4+b(4,i)*th5+
     1     b(5,i)*th6+b(6,i)*th7
         z=exp(-z)
         res=s(i)-th1*z
         rss=rss+res*res
         F(i)=res
      END DO
      call dblepr("rss",3,rss,1)
      DO iter=1,niter
         DO j=1,7
            dg(j)=0.d0
            DO k=j,7
               ak(j,k)=0.d0
            END DO
         END DO            
         DO i=1,nb
            z=b(1,i)*th2+b(2,i)*th3+b(3,i)*th4+
     1        b(4,i)*th5+b(5,i)*th6+b(6,i)*th7
            z=exp(-z)
            X(1)= -z
            z=z*th1
            X(2)=b(1,i)*z
            X(3)=b(2,i)*z
            X(4)=b(3,i)*z
            X(5)=b(4,i)*z
            X(6)=b(5,i)*z
            X(7)=b(6,i)*z
            DO j=1,7
               dg(j)=dg(j)+X(j)*F(i)
               DO k=j,7
                  ak(j,k)=ak(j,k)+X(j)*X(k)
               END DO
            END DO 
         END DO
         maxabsdg=abs(dg(1))
         DO j=2,7
            maxabsdg=max(maxabsdg,abs(dg(j)))
         END DO
         relrss = (oldrss-rss)/rss
         IF(maxabsdg.lt.eps.or.relrss.lt.1d-5) THEN
C  prepare things for return if gradient is close to 0
            th0=th1
            D(1)=th2
            D(2)=th3
            D(3)=th4
            D(4)=th5
            D(5)=th6
            D(6)=th7
            i=1
            DO j=1,7
               DO k=j,7
                  Varth(i)=ak(j,k)
                  i=i+1
               END DO
            END DO
            RETURN
         END IF
         gamma=min(gamma/alpha,1.d0)
C  End of step 3
         notacc=.TRUE.
         DO WHILE (notacc) 
            IF(gamma.lt.1.d0) THEN
               DO j=1,7
                  DO k=j,7
                     ck(j,k)=gamma*ak(j,k)
                  END DO
                  ck(j,j)=ck(j,j)+1.d0-gamma
               END DO
            ELSE
C   we may still need ak and dg so copy them to pk and ck
               DO j=1,7
                  DO k=j,7
                     ck(j,k)=ak(j,k)
                  END DO
               END DO
            END IF
            DO j=1,7
               pk(j)=dg(j)
            END DO
C   Now solve  ak%*%dtheta= dg
	    call dposv("U",7,1,ck,7,pk,7,info)
            call dblepr("pk",2,pk,7)
            call dblepr("dg",2,dg,7)
C  Step 4 we have pk 
            IF(info.ne.0) THEN
               gamma=alpha*gamma
               call dblepr("gamma1",6,gamma,1)
C  thats step 6
            ELSE
C  comute things needed for decision in step 5 
C  if successful F, nrss, and theta will be reused in the  
C  next iteration
               nth1=th1-gamma*pk(1)
               nth2=th2-gamma*pk(2)
               nth3=th3-gamma*pk(3)
               nth4=th4-gamma*pk(4)
               nth5=th5-gamma*pk(5)
               nth6=th6-gamma*pk(6)
               nth7=th7-gamma*pk(7)
               ntheta(1)=nth1
               ntheta(2)=nth2
               ntheta(3)=nth3
               ntheta(4)=nth4
               ntheta(5)=nth5
               ntheta(6)=nth6
               ntheta(7)=nth7
               call dblepr("ntheta",6,ntheta,7)
               call intpr("s",1,s,nb)
               nrss=0.d0
               DO i=1,nb
                 z=b(i,1)*nth2+b(i,2)*nth3+b(i,3)*nth4+
     1             b(i,4)*nth5+b(i,5)*nth6+b(i,6)*nth7
                  res=s(i)-nth1*exp(-z)
                  nrss=nrss+res*res
                  F(i)=res
               END DO
               crss=0.d0
               DO j=1,7
                  crss=crss+dg(j)*pk(j)
               END DO
               crss=rss-delta*gamma*crss
               call dblepr("crss",4,crss,1)
               call dblepr("nrss",4,nrss,1)
               call dblepr("F",1,F,nb)
               IF(nrss.le.crss) THEN
                  notacc=.FALSE.
C  accept new estimate, prepare for next iteration
               ELSE
                  gamma=alpha*gamma
                  IF(gamma.lt.0.1) return
               call dblepr("gamma2",6,gamma,1)
C  decrease gamma and try new regularization
               END IF
            END IF
         END DO
         th1=nth1
         th2=nth2
         th3=nth3
         th4=nth4
         th5=nth5
         th6=nth6
         th7=nth7
         oldrss=rss
         rss=nrss
      call dblepr("rss",3,rss,1)
         call rchkusr()
      call intpr("iter",4,iter,1)
      END DO
      th0=th1
      D(1)=th2
      D(2)=th3
      D(3)=th4
      D(4)=th5
      D(5)=th6
      D(6)=th7
      i=1
      DO j=1,7
         DO k=j,7
            Varth(i)=ak(j,k)
            i=i+1
         END DO
      END DO
      RETURN
      END
      subroutine solvedtb(s,nb,b,th0,D,Varth,F,niter,eps,rss)
C
C  Implements the regularized Gauss-Newton Algortithm (10.2.8)
C  from Schwetlick (1979)
C
      implicit logical (a-z)
      integer nb,s(nb),niter
      real*8 D(6),b(6,nb),th0,Varth(28),F(nb),eps
      integer i,j,k,info,iter
      real*8 Vtheta(7,7),z,thcorr,gamma,alpha,delta,
     1       dg(7),pk(7),ak(7,7),ck(7,7),rss,nrss,crss,maxabsdg,
     2       oldrss,relrss,theta(6),ntheta(6),res,X(7),nth0
      alpha=0.5D0
      delta=0.25D0
      DO j=1,6
         theta(j)=D(j)
      END DO
C      call dblepr("theta",5,theta,7)
      gamma=1.d0
      alpha=0.7d0
      oldrss=1.d50
      rss=0.d0
      DO i=1,nb
         z=0.d0
         DO j=1,6
            z=z+b(j,i)*theta(j)
         END DO
         z=exp(-z)
         res=s(i)-th0*z
         rss=rss+res*res
         F(i)=res
      END DO
C      call dblepr("rss",3,rss,1)
      DO iter=1,niter
         DO j=1,7
            dg(j)=0.d0
            DO k=j,7
               ak(j,k)=0.d0
            END DO
         END DO            
         DO i=1,nb
            z=0.d0
            DO j=1,6
               z=z+b(j,i)*theta(j)
            END DO
            z=exp(-z)
            X(7)= -z
            z=z*th0
            DO j=1,6
               X(j)=b(j,i)*z
            END DO
            DO j=1,7
               dg(j)=dg(j)+X(j)*F(i)
               DO k=j,7
                  ak(j,k)=ak(j,k)+X(j)*X(k)
               END DO
            END DO 
         END DO
         maxabsdg=0.d0
         DO j=1,7
            maxabsdg=max(maxabsdg,abs(dg(j)))
         END DO
         relrss = (oldrss-rss)/rss
         IF(maxabsdg.lt.eps.or.relrss.lt.1d-5) THEN
C  prepare things for return if gradient is close to 0
            DO j=1,6
               D(j)=theta(j)
            END DO
            i=1
            DO j=1,7
               DO k=j,7
                  Varth(i)=ak(j,k)
                  i=i+1
               END DO
            END DO
            RETURN
         END IF
         gamma=min(gamma/alpha,1.d0)
C  End of step 3
         notacc=.TRUE.
         DO WHILE (notacc) 
            IF(gamma.lt.1.d0) THEN
               DO j=1,7
                  DO k=j,7
                     ck(j,k)=gamma*ak(j,k)
                  END DO
                  ck(j,j)=ck(j,j)+1.d0-gamma
               END DO
            ELSE
C   we may still need ak and dg so copy them to pk and ck
               DO j=1,7
                  DO k=j,7
                     ck(j,k)=ak(j,k)
                  END DO
               END DO
            END IF
            DO j=1,7
               pk(j)=dg(j)
            END DO
C   Now solve  ak%*%dtheta= dg
	    call dposv("U",7,1,ck,7,pk,7,info)
C            call dblepr("pk",2,pk,7)
C            call dblepr("dg",2,dg,7)
C  Step 4 we have pk 
            IF(info.ne.0) THEN
               gamma=alpha*gamma
C               call dblepr("gamma1",6,gamma,1)
C  thats step 6
            ELSE
C  comute things needed for decision in step 5 
C  if successful F, nrss, and theta will be reused in the  
C  next iteration
               DO j=1,6
                  ntheta(j)=theta(j)-gamma*pk(j)
               END DO
               nth0=th0-gamma*pk(7)
C               call dblepr("ntheta",6,ntheta,6)
C               call dblepr("nth0",4,nth0,1)
C               call intpr("s",1,s,nb)
               nrss=0.d0
               DO i=1,nb
                  z=0.d0
                  DO j=1,6
                     z=z+b(j,i)*ntheta(j)
                  END DO
                  res=s(i)-nth0*exp(-z)
                  nrss=nrss+res*res
                  F(i)=res
               END DO
               crss=0.d0
               DO j=1,7
                  crss=crss+dg(j)*pk(j)
               END DO
               crss=rss-delta*gamma*crss
C               call dblepr("F",1,F,nb)
C               call dblepr("crss",4,crss,1)
C               call dblepr("nrss",4,nrss,1)
               IF(nrss.le.crss) THEN
                  notacc=.FALSE.
C  accept new estimate, prepare for next iteration
               ELSE
                  gamma=alpha*gamma
C               call dblepr("gamma2",6,gamma,1)
C  decrease gamma and try new regularization
               END IF
            END IF
         END DO
         th0=nth0
         DO j=1,6
            theta(j)=ntheta(j)
         END DO
         oldrss=rss
         rss=nrss
C      call dblepr("rss",3,rss,1)
         call rchkusr()
C      call intpr("iter",4,iter,1)
      END DO
      DO j=1,6
         D(j)=theta(j)
      END DO
      i=1
      DO j=1,7
         DO k=j,7
            Varth(i)=ak(j,k)
            i=i+1
         END DO
      END DO
      RETURN
      END
      subroutine replvar(x,ngrad,n1,n2,n3,rind,tind,lind,nind,sigma2,
     1                   mean)
C
C  estimate variances from replicated gradient images
C  returns  df*variance estimate which should be distributed 
C  chisq(df) * variance
C
      implicit logical (a-z)
      integer ngrad,n1,n2,n3,nind
      integer x(ngrad,n1,n2,n3),rind(ngrad),tind(nind),lind(nind)
      real*8 sigma2(n1,n2,n3),mean(n1,n2,n3)
      integer i,j,k,l,m,r
      real*8 z
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               sigma2(i,j,k)=0.d0
            END DO
         END DO
      END DO
      DO l=1,nind
         if(tind(l).le.1) CYCLE
         DO i=1,n1
            DO j=1,n2
               DO k=1,n3
                  mean(i,j,k)=0.d0
               END DO
            END DO
         END DO
         m=lind(l)
         DO r=1,ngrad
            if(rind(r).ne.m) CYCLE
            DO i=1,n1
               DO j=1,n2
                  DO k=1,n3
                     z=x(r,i,j,k)
                     mean(i,j,k)=mean(i,j,k)+z
                     sigma2(i,j,k)=sigma2(i,j,k)+z*z
                  END DO
               END DO
            END DO
         END DO
         DO i=1,n1
            DO j=1,n2
               DO k=1,n3
                  z=mean(i,j,k)
                  sigma2(i,j,k)=sigma2(i,j,k)-z*z/tind(l)
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
