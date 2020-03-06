      subroutine initdata(si,n1,n2,n3,nb,maxvalue)
C
C   project all values to (1,maxvalue) to avoid infinite estimates
C
      integer n1,n2,n3,nb
      double precision   si(n1,n2,n3,nb),sii,maxvalue
      integer i1,i2,i3,k
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               DO k=1,nb
                  sii=si(i1,i2,i3,k)
                  if(sii.le.0.0) si(i1,i2,i3,k)=1
                  if(sii.gt.maxvalue) si(i1,i2,i3,k)=maxvalue
               END DO
            END DO
         END DO
      END DO
      RETURN
      END
      subroutine outlier(si,n,nb,s0ind,siind,ls0,sinew,ind)
C
C   replace physically meaningless Si values by mean S0
C
      implicit none
      integer n,nb,ls0,s0ind(ls0),siind(*)
      double precision si(nb,n),sinew(nb,n)
      integer ind(n)
      integer i,j1,j,ls0m1
      double precision s0,sji
      integer changed
      ls0m1=ls0-1
      DO i=1,n
         s0=0
         DO j1=1,ls0
            j=s0ind(j1)
            sji=si(j,i)
            s0=s0+sji
            sinew(j,i)=sji
         END DO
         s0=(s0+ls0m1)/ls0
         changed=0
         DO j1=1,nb-ls0
            j=siind(j1)
            sji=si(j,i)
            if(sji.ge.s0) THEN
               sji=s0
               changed=1
            END IF
            sinew(j,i)=sji
         END DO
         ind(i)=changed
      END DO
      RETURN
      END
      subroutine outlierp(si,n,nb,s0ind,ls0,siind,lsi,sinew,nb1)
C
C   replace physically meaningless Si values by mean S0
C
      implicit none
      integer n,nb,nb1,ls0,lsi,s0ind(ls0),siind(lsi)
      double precision si(nb,n),sinew(nb1,n)
      integer i,j1,j,ls0m1,changed
      double precision s0,sinn(251),sji
      ls0m1=ls0-1
      DO i=1,n
         s0=0
         DO j1=1,ls0
            j=s0ind(j1)
            sji=si(j,i)
            s0=s0+sji
            sinn(j)=sji
         END DO
         s0=(s0+ls0m1)/ls0
         changed=0
         DO j1=1,lsi
            j=siind(j1)
            sji=si(j,i)
            if(sji.ge.s0) THEN
               sji=s0
               changed=1
            END IF
            sinn(j)=sji
         END DO
         sinn(nb1)=changed
         DO j=1,nb1
            sinew(j,i)=sinn(j)
         END DO
      END DO
      RETURN
      END

      subroutine mcorrlag(res,mask,n1,n2,n3,nv,sigma,mean,scorr,lag)

      implicit none
      integer n1,n2,n3,nv,lag(3)
      double precision scorr,res(nv,n1,n2,n3),sigma(n1,n2,n3),
     1       mean(n1,n2,n3)
      integer mask(n1,n2,n3)
      double precision vrm,zcorr,z,mi,mj
      integer i1,i2,i3,i4,l1,l2,l3,k,j1,j2,j3
      l1=lag(1)
      l2=lag(2)
      l3=lag(3)
      z=0.d0
      k=0
C  correlation in x
      do i1=1,n1-l1
         j1=i1+l1
         do i2=1,n2-l2
            j2=i2+l2
            do i3=1,n3-l3
               j3=i3+l3
               if (mask(i1,i2,i3)*mask(j1,j2,j3).eq.0) CYCLE
               vrm=sigma(i1,i2,i3)*sigma(j1,j2,j3)
               if(vrm.le.1e-10) CYCLE
               mi=mean(i1,i2,i3)
               mj=mean(j1,j2,j3)
               zcorr=(res(1,i1,i2,i3)-mi)*(res(1,j1,j2,j3)-mj)
               do i4=2,nv
                  zcorr=zcorr+(res(i4,i1,i2,i3)-mi)*
     1                         (res(i4,j1,j2,j3)-mj)
               enddo
               z=z+zcorr/vrm
               k=k+1
            enddo
         enddo
      enddo
      if( k.gt.0 ) then
         scorr=z/k/nv
      ELSE
         scorr=0.d0
      END IF
      return
      end

      subroutine msd(res,mask,n,nv,sigma,mean)
      implicit none
      integer n,nv
      double precision sigma(n),res(nv,n),mean(n)
      integer mask(n)
      integer i,iv
      double precision z,resi,zm,sigi
      DO i=1,n
         sigi=0.d0
         zm=0.d0
         if(mask(i).ne.0) THEN
            z=0.d0
            zm=0.d0
            DO iv=1,nv
               resi=res(iv,i)
               zm=zm+resi
               z=z+resi*resi
            END DO
            zm=zm/nv
            sigi=sqrt(z/nv-zm*zm)
         ENDIF
         mean(i)=zm
         sigma(i)=sigi
      END DO
      RETURN
      END

      subroutine mcorr(res,mask,n1,n2,n3,nv,sigma,mean,scorr,l1,l2,l3)

      implicit none
      integer n1,n2,n3,nv,l1,l2,l3,lag(3),n
      double precision scorr(l1,l2,l3),res(nv,n1,n2,n3),
     1       sigma(n1,n2,n3),mean(n1,n2,n3)
      integer mask(n1,n2,n3)
      integer i1,i2,i3
      double precision sci
      n=n1*n2*n3
      call msd(res,mask,n,nv,sigma,mean)
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call mcorrlag(res,mask,n1,n2,n3,nv,sigma,mean,sci,lag)
               scorr(i1,i2,i3)=sci
            END DO
         END DO
      END DO
      return
      end
      subroutine thcorr(w,n1,n2,n3,scorr,l1,l2,l3)

      implicit none
      integer n1,n2,n3,l1,l2,l3,lag(3)
      double precision scorr(l1,l2,l3),w(n1,n2,n3)
      integer i1,i2,i3
      double precision z,zcorr
      z=0.d0
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               z=z+w(i1,i2,i3)*w(i1,i2,i3)
            END DO
         END DO
      END DO
      Do i1=1,l1
         lag(1)=i1-1
         DO i2=1,l2
            lag(2)=i2-1
            DO i3=1,l3
               lag(3)=i3-1
               call thcorlag(w,n1,n2,n3,zcorr,lag)
               scorr(i1,i2,i3)=zcorr/z
               call rchkusr()
            END DO
         END DO
      END DO
      return
      end

      subroutine thcorlag(w,n1,n2,n3,scorr,lag)

      implicit none
      integer n1,n2,n3,lag(3)
      double precision scorr,w(n1,n2,n3)
      integer i1,i2,i3,c1,c2,c3,j1,j2,j3,l1,l2,l3
      double precision z
      c1=(n1-1)/2
      c2=(n2-1)/2
      c3=(n3-1)/2
      z=0.d0
      Do i1=-c1,c1
         j1=i1+c1+1
         l1=lag(1)-i1+c1+1
         if(l1.lt.1.or.l1.gt.n1) CYCLE
         DO i2=-c2,c2
            j2=i2+c2+1
            l2=lag(2)-i2+c2+1
            if(l2.lt.1.or.l2.gt.n2) CYCLE
            DO i3=-c3,c3
               j3=i3+c3+1
               l3=lag(3)-i3+c3+1
               if(l3.lt.1.or.l3.gt.n3) CYCLE
               z=z+w(j1,j2,j3)*w(l1,l2,l3)
            END DO
         END DO
      END DO
      scorr=z
      return
      end
      subroutine lconnect(segm,n1,n2,n3,i1,i2,i3,ind1,ind2,ind3,
     1                    mask)
C
C   assumes that we search for a connected region in segm==.TRUE.
C   that contains seed voxel (i1,i2,i3)
C   result: mask == .TRUE. if voxel is connected to seed
      implicit none
      integer n1,n2,n3,i1,i2,i3,ind1(*),ind2(*),ind3(*)
      integer mask(n1,n2,n3),segm(n1,n2,n3)
      logical final
      integer j1,j2,j3,k,l1,l2,l3,lind,lind0,ichecked
C     first find pixel close to (i1,i2) with segm(j1,j2)=0
      DO j1=1,n1
         DO j2=1,n2
            DO j3=1,n3
               mask(j1,j2,j3)=0
            END DO
         END DO
      END DO
      if(segm(i1,i2,i3).eq.0) THEN
         final=.FALSE.
         DO k=1,max(n1,n2,n3)
            DO l1=-k,k
               DO l2=-k,k
                  DO l3=-k,k
                     if(max(abs(l1),abs(l2),abs(l3)).ne.k) CYCLE
                     j1=i1+l1
                     if(j1.lt.1.or.j1.gt.n1) CYCLE
                     j2=i2+l2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     j3=i3+l3
                     if(j3.lt.1.or.j3.gt.n3) CYCLE
                     if(segm(j1,j2,j3).ne.0) THEN
                        final=.TRUE.
                        i1=j1
                        i2=j2
                        i3=j3
                     END IF
                     if(final) EXIT
                  END DO
                  if(final) EXIT
               END DO
               if(final) EXIT
            END DO
            if(final) EXIT
         END DO
      END IF
      mask(i1,i2,i3)=1
      ind1(1)=i1
      ind2(1)=i2
      ind3(1)=i3
      lind=1
      lind0=1
      ichecked=1
      final=.FALSE.
      DO while(.not.final)
         DO k=ichecked,lind0
            DO l1=-1,1
               DO l2=-1,1
                  DO l3=-1,1
                     if(l1.eq.0.and.l2.eq.0.and.i3.eq.0) CYCLE
                     j1=ind1(k)+l1
                     if(j1.lt.1.or.j1.gt.n1) CYCLE
                     j2=ind2(k)+l2
                     if(j2.lt.1.or.j2.gt.n2) CYCLE
                     j3=ind3(k)+l3
                     if(j3.lt.1.or.j3.gt.n3) CYCLE
              if(segm(j1,j2,j3).ne.0.and.mask(j1,j2,j3).eq.0) THEN
                        mask(j1,j2,j3)=1
                        lind=lind+1
                        ind1(lind)=j1
                        ind2(lind)=j2
                        ind3(lind)=j3
                     END IF
                  END DO
               END DO
            END DO
         END DO
         if(lind.eq.lind0) THEN
            final=.TRUE.
         ELSE
            ichecked=lind0
            lind0=lind
         END IF
      END DO
      RETURN
      END
      subroutine getmask(s0,n1,n2,n3,ns,level,msize,prop,s0m,mask)
      implicit none
      integer n1,n2,n3,ns,msize
      double precision  s0(n1,n2,n3,ns)
      double precision s0m(n1,n2,n3),prop,level
      integer mask(n1,n2,n3)
      integer i1,i2,i3,j,j1,j2,j3
      double precision z,anz,anz1
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               z=0.d0
               DO j=1,ns
                  z=z+s0(i1,i2,i3,j)
               END DO
               s0m(i1,i2,i3)=z/ns
            END DO
         END DO
      END DO
C
C   thats mean s0
C
      DO i1=1,n1
         DO i2=1,n2
            DO i3=1,n3
               anz=0
               anz1=0
               DO j1=max(1,i1-msize),min(n1,i1+msize)
                  DO j2=max(1,i2-msize),min(n2,i2+msize)
                     DO j3=max(1,i3-msize),min(n3,i3+msize)
                        if(s0m(j1,j2,j3).gt.level) anz1=anz1+1
                        anz=anz+1
                     END DO
                  END DO
               END DO
               mask(i1,i2,i3)=0
               if(anz1/anz.gt.prop) mask(i1,i2,i3)=1
            END DO
         END DO
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine means0(s0,n,ng0,level,ms0,mask)
C
C   calculate mean s0 value
C   generate mask
C   sweep s0 from si to generate  siq
C   calculate variance of siq
C
      integer n,ng0,mask(n),level
      double precision s0(ng0,n),ms0(n)
      integer i,k
      double precision z,thresh
      thresh = max(1,level*ng0)
      DO i=1,n
         z=0.d0
         DO k=1,ng0
            z=z+s0(k,i)
         END DO
         ms0(i) = z/ng0
         mask(i) = 0
         if(z.ge.thresh) mask(i) = 1
      END DO
      RETURN
      END

C
C __________________________________________________________________
C
      subroutine sweeps0(si,s0,n,ng0,ng1,level,siq,ms0,vsi,mask)
C
C   calculate mean s0 value
C   generate mask
C   sweep s0 from si to generate  siq
C   calculate variance of siq
C
      integer n,ng0,ng1,mask(n),level
      double precision si(ng1,n),s0(ng0,n)
      double precision siq(ng1,n),ms0(n),vsi(n)
      logical maskk
      integer i,k
      double precision s,z,z2,thresh,cv,s0mean,tvsi
      thresh = max(1,level*ng0)
      cv=ng1*(ng1-1)
      DO i=1,n
         z=0.d0
         DO k=1,ng0
            z=z+s0(k,i)
         END DO
         s0mean = z/ng0
         maskk = z.ge.thresh
         IF(maskk) THEN
            z=0.d0
            z2=0.d0
            DO k=1,ng1
               s=si(k,i)/s0mean
               if(s.gt.0.99d0) s=0.99d0
               z=z+s
               z2=z2+s*s
               siq(k,i)=s
            END DO
            tvsi=(ng1*z2-z)/cv
            if(tvsi.lt.1d-8) THEN
               maskk = .FALSE.
               tvsi=0.d0
            END IF
         ELSE
            tvsi=0.d0
            DO k=1,ng1
               siq(k,i)=1.d0
            END DO
         END IF
         ms0(i) = s0mean
         mask(i) = 0
         if(maskk) mask(i) = 1
         vsi(i) = tvsi
      END DO
      RETURN
      END
C
C __________________________________________________________________
C
      subroutine sweeps0p(si,s0,n,ng0,ng1,level,siq,ng2)
C
C   calculate mean s0 value
C   generate mask
C   sweep s0 from si to generate  siq
C   calculate variance of siq
C
      integer n,ng0,ng1,ng2,level
      double precision si(ng1,n),s0(ng0,n)
      double precision siq(ng2,n)
      logical maskk
      integer i,k
      double precision s,z,z2,thresh,cv,s0mean,tvsi,siqi(253)
      thresh = max(1,level*ng0)
      cv=ng1*(ng1-1)
      DO i=1,n
         z=0.d0
         DO k=1,ng0
            z=z+s0(k,i)
         END DO
         s0mean = z/ng0
         maskk = z.ge.thresh
         IF(maskk) THEN
            z=0.d0
            z2=0.d0
            DO k=1,ng1
               s=si(k,i)/s0mean
               s=min(s,0.99d0)
               z=z+s
               z2=z2+s*s
               siqi(k)=s
            END DO
            tvsi=(ng1*z2-z)/cv
            if(tvsi.lt.1d-8) THEN
               maskk = .FALSE.
               tvsi=0.d0
            END IF
         ELSE
            tvsi=0.d0
            DO k=1,ng1
               siqi(k)=1.d0
            END DO
         END IF
         siqi(ng1+1) = s0mean
         siqi(ng1+2) = tvsi
         if(maskk) THEN
            siqi(ng2) = 1
         ELSE
            siqi(ng2) = 0
         ENDIF
         DO k=1,ng2
            siq(k,i)=siqi(k)
         END DO
      END DO
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Hypergeometric 1F1 NIST HB 13.2.2, 13.2.39, 13.2(iv)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine hg1f1(a,b,z,n,fz)
      implicit none
      integer n
      double precision a,b,z(n),fz(n)
      integer i
      double precision x,y,d,eps,zi,ezi,ai,gofbai
      double precision gammaf
      external gammaf
      eps=1.d-15
      gofbai=gammaf(b)/gammaf(b-a)
      DO i=1,n
         d = 1.d0
         zi = z(i)
         IF(zi.lt.0) THEN
            ezi=exp(zi/2)
            ai=b-a
            if(zi.lt.-1400) THEN
               fz(i) = exp((-a)*log(-zi))*gofbai+5.64e-3+1.9e-3*b
C   add +5.64e-3+1.9e-3*b to keep the function monotone
               CYCLE
            END IF
         ELSE
            ezi=1.d0
            ai=a
         ENDIF
         x = ezi
         y = ezi
         DO WHILE (abs(y).gt.abs(x)*eps)
            y = -y*(ai+d-1.d0)/(b+d-1.d0)*zi/d
            x = x+y
            d = d+1.d0
         END DO
         fz(i) = ezi*x
      END DO
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     get rho from D (uses c(1,0,0,1,0,1) in case of inplausible D
C     use of estimated tensors (linearized model) as initial values
C     for nonlinear regression
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine d2rall(D,rho,nvox)
        implicit none
        integer nvox
        double precision D(6,nvox),rho(6,nvox)
        integer i,ierr
        double precision ew(3),ev(3,3),r1,r2,r3,r4,r5
        DO i=1,nvox
            call eigen3(D(1,i),ew,ev,ierr)
            if(ew(1).le.1.d-6) THEN
              rho(1,i) = 1.d-2
              rho(2,i) = 0.d0
              rho(3,i) = 0.d0
              rho(4,i) = 1.d-2
              rho(5,i) = 0.d0
              rho(6,i) = 1.d-2
            ELSE
              r1=sqrt(max(1d-12,D(1,i)))
              r2=D(2,i)/r1
              r3=D(3,i)/r1
              r4=sqrt(max(1d-12,D(4,i)-r2*r2))
              r5=(D(5,i)-r2*r3)/r4
              rho(6,i)=sqrt(max(1d-16,D(6,i)-r3*r3-r5*r5))
              rho(5,i)=r5
              rho(4,i)=r4
              rho(3,i)=r3
              rho(2,i)=r2
              rho(1,i)=r1
            END IF
        END DO
        RETURN
        END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     get D from rho
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine r2dall(rho,D,nvox)
      implicit none
      integer nvox,i
      double precision D(6,nvox),rho(6,nvox),eps
      double precision r1,r2,r3,r4,r5,r6
      eps=0.d-12
      DO i=1,nvox
          r1=rho(1,i)
          r2=rho(2,i)
          r3=rho(3,i)
          r4=rho(4,i)
          r5=rho(5,i)
          r6=rho(6,i)
          D(1,i)=r1*r1+eps
          D(2,i)=r1*r2
          D(3,i)=r1*r3
          D(4,i)=r2*r2+r4*r4+eps
          D(5,i)=r2*r3+r4*r5
          D(6,i)=r3*r3+r5*r5+r6*r6+eps
      END DO
C      call testreg2(D,rho,111)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     get D from rho
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rho2D(rho,D)
      implicit none
      double precision D(6),rho(6),eps
      eps=0.d-12
      D(1)=rho(1)*rho(1)+eps
      D(2)=rho(1)*rho(2)
      D(3)=rho(1)*rho(3)
      D(4)=rho(2)*rho(2)+rho(4)*rho(4)+eps
      D(5)=rho(2)*rho(3)+rho(4)*rho(5)
      D(6)=rho(3)*rho(3)+rho(5)*rho(5)+rho(6)*rho(6)+eps
C      call testreg2(D,rho,111)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     get D from rho without regularization
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rho2D0(rho,D)
      implicit none
      double precision D(6),rho(6),r1,r2,r3,r4,r5,r6
      r1=rho(1)
      r2=rho(2)
      r3=rho(3)
      r4=rho(4)
      r5=rho(5)
      r6=rho(6)
      D(1)=r1*r1
      D(2)=r1*r2
      D(3)=r1*r3
      D(4)=r2*r2+r4*r4
      D(5)=r2*r3+r4*r5
      D(6)=r3*r3+r5*r5+r6*r6
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     calculate fitted values of si
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sihat(th0i,Di,btb,F,nb)
      implicit none
      integer nb
      double precision th0i,Di(6),btb(6,nb),F(nb)
      integer j,k
      double precision z
      DO k=1,nb
        z=0.d0
        DO j=1,6
          z=z+btb(j,k)*Di(j)
        END DO
        F(k)=th0i*exp(-z)
      END DO
      RETURN
      END

      double precision function ddot3sq(phi,eta,g)
C  compute (g^T n(phi,eta))^2
      implicit none
      double precision phi,eta,g(3),sphi,z
      sphi=dsin(phi)
      z=g(1)*dcos(eta)*sphi+g(2)*dsin(eta)*sphi+g(3)*dcos(phi)
      ddot3sq=z*z
      RETURN
      END
      subroutine dgtddphi(phi,eta,g, d)
C  compute (g^T n(phi,eta))^2
C           d(g^T n(phi,eta))^2 / d phi
C           d(g^T n(phi,eta))^2 / d eta
C  as components of d
      implicit none
      double precision phi,eta,g(3),sphi,cphi,ceta,seta,z,ze,zp,d(3)
      sphi=dsin(phi)
      cphi=dcos(phi)
      seta=dsin(eta)
      ceta=dcos(eta)
      z=g(1)*ceta*sphi+g(2)*seta*sphi+g(3)*cphi
      zp=g(1)*ceta*cphi+g(2)*seta*cphi-g(3)*sphi
      ze=g(2)*ceta*sphi-g(1)*seta*sphi
      d(1)=z*z
      d(2)=2.d0*z*zp
      d(3)=2.d0*z*ze
      RETURN
      END
      subroutine besselq(x,n,fw)
      implicit none
      integer n
      double precision x(n),fw(n)
      external besseli
      double precision besseli
      integer i
      DO i=1,n
          x(i)=i*1.d-2
C  avoid transform to zero for numerical stability
          fw(i)=besseli(x(i),1.d0,2.d0)/besseli(x(i),0.d0,2.d0)
      END DO
      RETURN
      END
