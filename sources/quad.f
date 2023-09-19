       subroutine quad(Nrad0,Npoints)
       use integration_grid
       implicit double precision(a-h,o-z)
       double precision x(1000),y(1000),z(1000)

       PI=datan(1.0d0)*4.0d0 

        r=rr00

c         READ(*,*)Nrad
         CALL LEGZO(Nrad0,XR,WR)
         do i=1,nrad0
          wR(i)=2.0*r*wr(i)/(1.0-Xr(i))**2
          XR(i)=R*(Xr(i)+1.0)/(1.0-Xr(i))
c       if(XR(i).lt.0)then 
c        wR(i)=2.0*r*wr(i)/(1.0-Xr(i))**2
c        XR(i)=R*(Xr(i)+1.0)/(1.0-Xr(i))
c       else
c        a0=R
c        a1=2.0d0*Rmax
c        a2=(-Rmax-3.0d0*R)
c        a3=2.0d0*R
c        XR(i)=a0+a1*XR(i)+a2*XR(i)**2.d0+a3*XR(i)**3.d0
c        wR(i)=wr(i)*(a1+2.0d0*a2*XR(i)+3.0d0*a3*XR(i)**2)
c       end if
c        write(*,*)i,xr(i)
c       wR(i)=r*wr(i)/2.0
c       XR(i)=R*(Xr(i)+1.0)/2.0
       end do

c      WRITE(*,*)'Please enter angular points '
c      WRITE(*,*)'6,14,26,38,50,74,86,110,146,170,194,230,266,302'
c      WRITE(*,*)'350,434,590,770,974'
c      READ(*,*)Npoints
 
       if(npoints.eq.6)     CALL LD0006(X,Y,Z,W,N)
       if(npoints.eq.14)    CALL LD0014(X,Y,Z,W,N)
       if(npoints.eq.26)    CALL LD0026(X,Y,Z,W,N)
       if(npoints.eq.38)    CALL LD0038(X,Y,Z,W,N)
       if(npoints.eq.50)    CALL LD0050(X,Y,Z,W,N)
       if(npoints.eq.74)    CALL LD0074(X,Y,Z,W,N)
       if(npoints.eq.86)    CALL LD0086(X,Y,Z,W,N)
       if(npoints.eq.110)   CALL LD0110(X,Y,Z,W,N)
       if(npoints.eq.146)   CALL LD0146(X,Y,Z,W,N)
       if(npoints.eq.170)   CALL LD0170(X,Y,Z,W,N)
       if(npoints.eq.194)   CALL LD0194(X,Y,Z,W,N)
       if(npoints.eq.230)   CALL LD0230(X,Y,Z,W,N)
       if(npoints.eq.266)   CALL LD0266(X,Y,Z,W,N)
       if(npoints.eq.302)   CALL LD0302(X,Y,Z,W,N)
       if(npoints.eq.350)   CALL LD0350(X,Y,Z,W,N)
       if(npoints.eq.434)   CALL LD0434(X,Y,Z,W,N)
       if(npoints.eq.590)   CALL LD0590(X,Y,Z,W,N)
       if(npoints.eq.770)   CALL LD0770(X,Y,Z,W,N)
       if(npoints.eq.974)   CALL LD0974(X,Y,Z,W,N)
      
c angular mesh
       do iang=1,npoints
        xx=x(iang)
        yy=y(iang)
        zz=z(iang)
        th(iang)=dacos(zz)
        if(dsin(th(iang)).ne.0.0) then 
!        if(ABS(dsin(th(iang))).gt.1E-20) then 
         xs=xx/dsin(th(iang))
         if(xs.gt.1.0d0) xs=1.0d0
         if(xs.lt.-1.0d0) xs=-1.0d0
         ph(iang)=dacos(xs)
c        if(yy.lt.0.0) ph(iang)=ph(iang)+Pi
         if(yy.lt.0.0) ph(iang)=-ph(iang)
        else
         ph(iang)=0.0
        end if
       end do

c      vol=0.0
c      do i=1,nrad
c       do j=1,npoints
c        vol=vol+funct(xr(i),th(j),ph(j))*(xr(i)**2.0)*wr(i)*w(j)
c       end do
c      end do
c      vol=4.0*Pi*vol
c      write(*,*) '3D integ ', vol
c
       end

cfunction funct(R,th,ph)
c       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cPI=datan(1.0d0)*4.0d0 

cxx=exp(-R)
cfunct=xx*xx
cend

c###############################################################


        SUBROUTINE LEGZO(N,X,W)
C
C       =========================================================
C       Purpose : Compute the zeros of Legendre polynomial Pn(x)
C                 in the interval [-1,1], and the corresponding
C                 weighting coefficients for Gauss-Legendre
C                 integration
C       Input :   n    --- Order of the Legendre polynomial
C       Output:   X(n) --- Zeros of the Legendre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),W(N)
        N0=(N+1)/2
        DO 45 NR=1,N0
           Z=DCOS(3.1415926D0*(NR-0.25D0)/N)
10         Z0=Z
           P=1.0D0
           DO 15 I=1,NR-1
15            P=P*(Z-X(I))
           F0=1.0D0
           IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0.0D0
           F1=Z
           DO 20 K=2,N
              PF=(2.0D0-1.0D0/K)*Z*F1-(1.0D0-1.0D0/K)*F0
              PD=K*(F1-Z*PF)/(1.0D0-Z*Z)
              F0=F1
20            F1=PF
           IF (Z.EQ.0.0) GO TO 40
           FD=PF/P
           Q=0.0D0
           DO 35 I=1,NR
              WP=1.0D0
              DO 30 J=1,NR
                 IF (J.NE.I) WP=WP*(Z-X(J))
30            CONTINUE
35            Q=Q+WP
           GD=(PD-Q*FD)/P
           Z=Z-FD/GD
           IF (DABS(Z-Z0).GT.DABS(Z)*1.0D-15) GO TO 10
40         X(NR)=Z
           X(N+1-NR)=-Z
           W(NR)=2.0D0/((1.0D0-Z*Z)*PD*PD)
45         W(N+1-NR)=W(NR)
        RETURN
        END

