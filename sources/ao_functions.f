!! *********************************************************************** !!
!! BASIS FUNCTIONS AT A POINT, for analyses that work point by point      !!
!! instead of on the integration grid (DFT-DM1, OS-CENTROID):             !!
!!   gpoints   -- value and gradient of every basis function at a point   !!
!!   aofunct   -- value of one basis function at a point                  !!
!!   gxfunct, gyfunct, gzfunct -- x, y, z derivative of the basis         !!
!!                function selected in common /actual_bf/ (drho_xyz)      !!
!!   atsphere  -- radius of each atom's trust sphere (OS-CENTROID)        !!
!!   angle     -- cosine of the angle between two vectors                 !!
!! *********************************************************************** !!

      subroutine gpoints(x,y,z,gxchp,gychp,gzchp,aochp)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common/actual_bf/iact,jat,icenter
!! DFT-DM1 calls this from an OMP-parallelized loop (dft_dm1.f) -- each   !!
!! thread needs its own iact, not one shared across all of them.         !!
!$OMP THREADPRIVATE(/actual_bf/)
      common /nat/ nat,igr,ifg,idum(4)
      dimension gxchp(igr),gychp(igr),gzchp(igr),aochp(igr)

      do i=1,igr
       iact=i
       gxchp(iact)=gxfunct(x,y,z)
       gychp(iact)=gyfunct(x,y,z)
       gzchp(iact)=gzfunct(x,y,z)
       aochp(iact)=aofunct(iact,x,y,z)
      end do
      return
      end

!! gxfunct/gyfunct/gzfunct: thin wrappers around drho_xyz (dft_dm1.f)   !!

      function gxfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      gxfunct=drho_xyz(x,y,z,1)
      return
      end

      function aofunct(iact,x,y,z)
      use basis_set
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,idum(4)
      
      iactat=ihold(iact)
      f=0.d0
!! xx/yy/zz, not x/y/z: the arguments are gpoints' own point, which it   !!
!! reuses for every basis function                                      !!
      xx=x-coord(1,iactat)
      yy=y-coord(2,iactat)
      zz=z-coord(3,iactat)
      rr=dsqrt(xx**2+yy**2+zz**2)
      k=1
      do while(nprimbas(k,iact).ne.0)
       ipr=nprimbas(k,iact)
       nn=nlm(ipr,1)
       ll=nlm(ipr,2)
       mm=nlm(ipr,3)
       f=f+(xx**nn)*(yy**ll)*(zz**mm)*dexp(-expp(ipr)*(rr**2))*coefpb(ipr,iact)
       k=k+1
      enddo
      aofunct=f
      
      return
      end
      function gyfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      gyfunct=drho_xyz(x,y,z,2)
      return
      end

      function gzfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      gzfunct=drho_xyz(x,y,z,3)
      return
      end

      subroutine atsphere()
      use basis_set
      implicit real*8 (a-h,o-z)
      include 'parameter.h'
      common /quadrat/th(1000),ph(1000),w(1000),wr(500),Xr(500)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /nna/ xnonradi(maxnna),atsphradi(maxat),nna
      common /gridp/Nrad,Nang
      common/pha2/pha,phb,rr00

       xangle=cos(45.0d0*pi/180.0d0)
       tres=0.15d0 
       do i=1,nat
        atsphradi(i)=0.0d0
        do irad=nrad,1,-1
         if(xr(irad).gt.tres)then
          do iang=1,nang
           thx=th(iang)
           fix=ph(iang)
           rr=xr(irad)
           xabs0=rr*dsin(thx)*dcos(fix)
           yabs0=rr*dsin(thx)*dsin(fix)
           zabs0=rr*dcos(thx)
           xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
           yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
           zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
           x=xabs+coord(1,i)
           y=yabs+coord(2,i)
           z=zabs+coord(3,i)
           gx=dfunctxyz(1,x,y,z)
           gy=dfunctxyz(2,x,y,z)
           gz=dfunctxyz(3,x,y,z)
           if(angle(-xabs,-yabs,-zabs,gx,gy,gz).lt.xangle) go to 44
          end do  
         end if
        end do   
44     continue
       atsphradi(i)=xr(irad+1)
      end do 
      end 

       function angle(x,y,z,xx,yy,zz)
       implicit none
       double precision x,y,z,xx,yy,zz
       double precision angle

       angle=((x*xx+y*yy+z*zz)/(dsqrt((x)**2.+(y)**2.+(z)**2.)*dsqrt((xx)**2.+(yy)**2.+(zz)**2.)))

       end function angle
